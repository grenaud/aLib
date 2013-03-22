#include <iostream>
#include <vector>
#include <set>
#include <ctype.h>
#include <stdlib.h>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "utils.h"

using namespace std;
using namespace BamTools;

typedef struct{
    char bp;
    int offset;
} mdField;


typedef struct countMatMis{
    int match;
    int mismatch;
} countMatMis;

int asciiOffsetZero=48;
int maxQualScore=0;
int illuminaOffset=33;
char DUMMYCHAR='#';

inline int char2phred(char c){
    int toReturn=int(c)-illuminaOffset;
    if(toReturn>maxQualScore)
	maxQualScore=toReturn;
    return toReturn;
}

//This function converts the MD field into a vector of mdField structs
//we skip deletions in the read (ins in reference)
inline vector<mdField> mdString2Vector(const string & mdFieldToParse){
    vector<mdField> toReturn;
    int i=0;
    // int addToOffset=0;
    mdField toadd;
    

    toadd.offset=0;
    toadd.bp='N';

    while(int(mdFieldToParse.length()) != i){
	if(isdigit(mdFieldToParse[i])){
	    toadd.offset=toadd.offset*10+(int(mdFieldToParse[i])-asciiOffsetZero);
	}else{
	    //deletions in read (insertion in reference)
	    if(mdFieldToParse[i] == '^'){
		if(toadd.offset != 0){
		    toadd.bp=DUMMYCHAR;
		    toReturn.push_back(toadd);
		    toadd.offset=0;
		    toadd.bp='N';
		}

		i++;
		mdField toadd2;
		toadd2.offset=0;
		toadd2.bp='^';
		while(isalpha(mdFieldToParse[i])){
		    i++;
		    toadd2.offset++;
		}
		toReturn.push_back(toadd2);
		i--;
	    }else{
		toadd.bp=mdFieldToParse[i];
		toReturn.push_back(toadd);

		toadd.offset=0;
		toadd.bp='N';
	    }

	}
	i++;
    }
    return toReturn;
}


int main (int argc, char *argv[]) {
    bool pairedEnd=false;

    string usage=string(""+string(argv[0])+"   [aligned BAM file] [mask files] [cycles] "+
			"\nThis program computes the observed vs predicted tallies for the quality scores\n"+
			"arguments:\n"+
			"[mask files] : Comma delimited list of files containing the masked positions on the PhiX\n"+
			"[cycles]     : Number of cycles for the reads (use comma for paired-end e.g. 95,95)\n"+		       
			"\n");

    if(argc != 4 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    string bamfiletopen = string(argv[1]);
    set<int> maskedPos;    

    vector<string> fileMask=allTokens(string(argv[2]),',');

    for(unsigned int i=0;i<fileMask.size();i++){
	ifstream myFile; 
	string line;
	myFile.open(fileMask[i].c_str(), ios::in);

	if (myFile.is_open()){
	    while ( getline (myFile,line)){
		vector<string> lineTok=allTokens(line,'\t');
		maskedPos.insert(atoi(lineTok[1].c_str()));
	    }
	    myFile.close();
	}else{
	    cerr << "Unable to open file "<<fileMask[i]<<endl;
	    return 1;
	}
    }

    vector<string> cycles=allTokens(string(argv[3]),',');
    vector<countMatMis>  perCycleCountr122;
    vector< vector<countMatMis> > perCycleCountr1;
    vector< vector<countMatMis> > perCycleCountr1NM2;
    // vector< vector<count> > perCycleCountr2;
    // vector< vector<count> > perCycleCountr2NM2;


    if(cycles.size() == 1 ){
	pairedEnd=false;		
    }else{ //paired-end
	if(cycles.size() == 2){
	    pairedEnd=true;
	}else{
	    cerr << "Wrong number of cycles." << endl;
	    return 1;
	}
    }

    int cyclesFirst=0;
    int cyclesSecond=0;

    for(unsigned int j=0;j<cycles.size();j++){

	int numberOfCycles=atoi(cycles[j].c_str());
	    
	for(int k=0;k<numberOfCycles;k++){
	    vector<countMatMis> v1;
	    vector<countMatMis> v2;

	    for(int i=0;i<64;i++){
		countMatMis toadd;
		toadd.match=0;
		toadd.mismatch=0;
		v1.push_back(toadd);
		v2.push_back(toadd);
	    }

	    if(j == 0 ){
		cyclesFirst  = numberOfCycles;
	    }else{
		cyclesSecond = numberOfCycles;
	    }
	    perCycleCountr1.push_back(v1);
	    perCycleCountr1NM2.push_back(v2);	   
	}
    }


    BamReader reader;

    
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    vector<countMatMis> allCount;
    vector<countMatMis> editDist2Count;
    vector<countMatMis> countA;
    vector<countMatMis> countC;
    vector<countMatMis> countG;
    vector<countMatMis> countT;

    for(int i=0;i<64;i++){
	countMatMis toadd;
	toadd.match=0;
	toadd.mismatch=0;

	allCount.push_back(toadd);
	editDist2Count.push_back(toadd);
	countA.push_back(toadd);
	countC.push_back(toadd);
	countG.push_back(toadd);
	countT.push_back(toadd);
    }

    BamAlignment al;
    int editDist;
    string mdFieldString;
    string reconstructed;
    string reconstructedTemp;
    vector<CigarOp> cigarData;

    while ( reader.GetNextAlignment(al) ) {
	//initialize
	editDist=-1;
	mdFieldString="";
	reconstructed="";
	reconstructedTemp="";

	//skip unmapped
	if(!al.IsMapped())
	    continue;

	//get relevant data
	if(!al.GetTag("NM",editDist)){
	    cerr<<"Cannot get NM tag from "<<al.Name<<endl;
	    return 1;
	}
	if(!al.GetTag("MD",mdFieldString)){
	    cerr<<"Cannot get NM tag from "<<al.Name<<endl;
	    return 1;
	}
	
	cigarData=al.CigarData;
	for(unsigned int i=0;i<cigarData.size();i++){
	    reconstructedTemp+=string(cigarData[i].Length,cigarData[i].Type);
	}


	//get a vector representation of the MD field	

	vector<mdField> parsedMD=mdString2Vector(mdFieldString);
	// cout<<mdFieldString<<endl;
	// for(int i=0;i<parsedMD.size();i++){
	//     cout<<parsedMD[i].bp<<"\t"<<parsedMD[i].offset<<endl;
	// }
	//	continue;
	vector<int> positionsOnControl;
	int initialPositionControl=al.Position;

	//combine the CIGAR and MD into one single string
	unsigned int mdVectorIndex=0;

	for(unsigned int i=0;i<reconstructedTemp.size();i++){
	    if(reconstructedTemp[i] == 'M' ){ //only look at matches and indels	    
		
		if(mdVectorIndex<parsedMD.size()){ //still have mismatches

		    if(parsedMD[mdVectorIndex].offset == 0){ //we have reached a mismatch				

			if(parsedMD[mdVectorIndex].bp == DUMMYCHAR){ //no char to add, need to backtrack on the CIGAR
			    i--;
			}else{
			    reconstructed+=parsedMD[mdVectorIndex].bp;
			    positionsOnControl.push_back(initialPositionControl++);
			}
			mdVectorIndex++;
		    }else{ //wait until we reach a mismatch
			reconstructed+=reconstructedTemp[i];
			parsedMD[mdVectorIndex].offset--;
			positionsOnControl.push_back(initialPositionControl++);
		    }

		    //skipping all the positions with deletions on the read
		    while(parsedMD[mdVectorIndex].bp == '^'){ 
			initialPositionControl+=parsedMD[mdVectorIndex].offset;
			mdVectorIndex++;
		    }
		    
		}else{
		    reconstructed+=reconstructedTemp[i];
		    positionsOnControl.push_back(initialPositionControl++);
		}
	    }else{
		if(reconstructedTemp[i] == 'S' || reconstructedTemp[i] == 'I'){ //soft clipped bases and indels
		    reconstructed+=reconstructedTemp[i];
		    positionsOnControl.push_back(initialPositionControl);
		}
	    }
	}

	if(reconstructed.size() != al.QueryBases.size()){
	    cerr << "Could not recreate the sequence for read "<<al.Name << endl;
	    return 1;
	}

	if(positionsOnControl.size() != reconstructed.size()){
	    cerr << "Could not determine the positions for the read "<<al.Name << endl;
	    return 1;
	}

	// cout<<reconstructed<<endl;
	// for(int i=0;i<positionsOnControl.size();i++)
	//     cout<<reconstructed[i]<<"\t"<<positionsOnControl[i]<<endl;

	// int cyclesFirst;
	// int cyclesSecond;

	for(unsigned int i=0;i<reconstructed.size();i++){
	    //skip masked positions 
	    if(maskedPos.find(positionsOnControl[i]) != maskedPos.end())
		continue;
	    int cycle; //minus one for vectors in c++
	    if(al.IsSecondMate()){
		if(!pairedEnd){
		    cerr << "Cannot have paired reads like: "<<al.Name << " and be in unpaired mode"<<endl;
		    return 1;
		}
		if(al.IsReverseStrand()){
		    cycle=cyclesFirst+(cyclesSecond-i-1);
		}else{
		    cycle=cyclesFirst+i;
		}
	    }else{
		if(al.IsReverseStrand()){
		    cycle=cyclesFirst-i-1;
		}else{
		    cycle=i;
		}
	    }
	    
	    //skip insertions and softclip
	    if(reconstructed[i] == 'S' ||
	       reconstructed[i] == 'I' ){
		//reconstructed[i]='X';		
	    }else{
		//match
		if(reconstructed[i] == 'M'){
		    allCount[char2phred(al.Qualities[i])].match++;
		    perCycleCountr1[cycle][char2phred(al.Qualities[i])].match++;

		    if(editDist <=2 ){
			editDist2Count[char2phred(al.Qualities[i])].match++;
			perCycleCountr1NM2[cycle][char2phred(al.Qualities[i])].match++;
		    }

		    if(!al.IsReverseStrand()){
			if(al.QueryBases[i] == 'A')
			    countA[char2phred(al.Qualities[i])].match++;
			if(al.QueryBases[i] == 'C')
			    countC[char2phred(al.Qualities[i])].match++;
			if(al.QueryBases[i] == 'G')
			    countG[char2phred(al.Qualities[i])].match++;
			if(al.QueryBases[i] == 'T')
			    countT[char2phred(al.Qualities[i])].match++;
		    }else{
			if(al.QueryBases[i] == 'T')
			    countA[char2phred(al.Qualities[i])].match++;
			if(al.QueryBases[i] == 'G')
			    countC[char2phred(al.Qualities[i])].match++;
			if(al.QueryBases[i] == 'C')
			    countG[char2phred(al.Qualities[i])].match++;
			if(al.QueryBases[i] == 'A')
			    countT[char2phred(al.Qualities[i])].match++;
		    }
		}else{
		    //mismatch
		    if(reconstructed[i] == 'A' || reconstructed[i] == 'C' || reconstructed[i] == 'G' || reconstructed[i] == 'T' ){
			allCount[char2phred(al.Qualities[i])].mismatch++;
			perCycleCountr1[cycle][char2phred(al.Qualities[i])].mismatch++;

			if(editDist <=2 ){
			    editDist2Count[char2phred(al.Qualities[i])].mismatch++;
			    perCycleCountr1NM2[cycle][char2phred(al.Qualities[i])].mismatch++;
			}

			if(!al.IsReverseStrand()){
			    if(al.QueryBases[i] == 'A')
				countA[char2phred(al.Qualities[i])].mismatch++;
			    if(al.QueryBases[i] == 'C')
				countC[char2phred(al.Qualities[i])].mismatch++;
			    if(al.QueryBases[i] == 'G')
				countG[char2phred(al.Qualities[i])].mismatch++;
			    if(al.QueryBases[i] == 'T')
				countT[char2phred(al.Qualities[i])].mismatch++;
			}else{
			    if(al.QueryBases[i] == 'T')
				countA[char2phred(al.Qualities[i])].mismatch++;
			    if(al.QueryBases[i] == 'G')
				countC[char2phred(al.Qualities[i])].mismatch++;
			    if(al.QueryBases[i] == 'C')
				countG[char2phred(al.Qualities[i])].mismatch++;
			    if(al.QueryBases[i] == 'A')
				countT[char2phred(al.Qualities[i])].mismatch++;
			}

		    }else{
			cerr << "Could not recreate the sequence for read"<<al.Name << endl;
			return 1;
		    }
		}
	    }       
	}	    

	// cout<<reconstructed<<endl;

    }
    reader.Close();
    ofstream outfile;

    outfile.open(string(bamfiletopen+".baseobspred1").c_str());
    if ( !outfile.is_open() ) { cerr<<"Cannot write to "<<bamfiletopen<<".baseobspred1"<<endl; }
    for(int qual=0;qual<=maxQualScore;qual++){
	outfile<<qual<<"\t"<<editDist2Count[qual].match<<"\t"<<editDist2Count[qual].mismatch<<endl;
    }
    outfile.close();

    outfile.open(string(bamfiletopen+".baseobspred2").c_str());
    if ( !outfile.is_open() ) { cerr<<"Cannot write to "<<bamfiletopen<<".baseobspred1"<<endl; }
    for(int qual=0;qual<=maxQualScore;qual++){
    	outfile<<qual<<"\t"<<allCount[qual].match<<"\t"<<allCount[qual].mismatch<<endl;
    }
    outfile.close();

    outfile.open(string(bamfiletopen+".baseobspred3").c_str());
    if ( !outfile.is_open() ) { cerr<<"Cannot write to "<<bamfiletopen<<".baseobspred1"<<endl; }
    for(int qual=0;qual<=maxQualScore;qual++){
    	outfile<<qual<<"\t"<<countA[qual].match<<"\t"<<countA[qual].mismatch<<endl;
    }
    outfile.close();

    outfile.open(string(bamfiletopen+".baseobspred4").c_str());
    if ( !outfile.is_open() ) { cerr<<"Cannot write to "<<bamfiletopen<<".baseobspred1"<<endl; }
    for(int qual=0;qual<=maxQualScore;qual++){
    	outfile<<qual<<"\t"<<countC[qual].match<<"\t"<<countC[qual].mismatch<<endl;
    }
    outfile.close();

    outfile.open(string(bamfiletopen+".baseobspred5").c_str());
    if ( !outfile.is_open() ) { cerr<<"Cannot write to "<<bamfiletopen<<".baseobspred1"<<endl; }
    for(int qual=0;qual<=maxQualScore;qual++){
    	outfile<<qual<<"\t"<<countG[qual].match<<"\t"<<countG[qual].mismatch<<endl;
    }
    outfile.close();

    outfile.open(string(bamfiletopen+".baseobspred6").c_str());
    if ( !outfile.is_open() ) { cerr<<"Cannot write to "<<bamfiletopen<<".baseobspred1"<<endl; }
    for(int qual=0;qual<=maxQualScore;qual++){
    	outfile<<qual<<"\t"<<countT[qual].match<<"\t"<<countT[qual].mismatch<<endl;
    }
    outfile.close();


    //print per cycle




    
    int totalCycles=(cyclesFirst+cyclesSecond);
    outfile.open(string(bamfiletopen+".baseobspredcycle1").c_str());
    if ( !outfile.is_open() ) { cerr<<"Cannot write to "<<bamfiletopen<<".baseobspredcycle1"<<endl; }
    for(int qual=0;qual<=maxQualScore;qual++){
	outfile<<qual<<"\t";
	for(int k=0;k<totalCycles;k++){
	    if(k== (totalCycles-1))
		outfile<<perCycleCountr1[k][qual].match<<"\t"<<perCycleCountr1[k][qual].mismatch;
	    else
		outfile<<perCycleCountr1[k][qual].match<<"\t"<<perCycleCountr1[k][qual].mismatch<<"\t";
	}
	outfile<<endl;
    }
    outfile.close();


    outfile.open(string(bamfiletopen+".baseobspredcycle2").c_str());
    if ( !outfile.is_open() ) { cerr<<"Cannot write to "<<bamfiletopen<<".baseobspredcycle2"<<endl; }
    for(int qual=0;qual<=maxQualScore;qual++){
	outfile<<qual<<"\t";
	for(int k=0;k<totalCycles;k++){
	    if(k== (totalCycles-1))
		outfile<<perCycleCountr1NM2[k][qual].match<<"\t"<<perCycleCountr1NM2[k][qual].mismatch;
	    else
		outfile<<perCycleCountr1NM2[k][qual].match<<"\t"<<perCycleCountr1NM2[k][qual].mismatch<<"\t";
	}
	outfile<<endl;
    }
    outfile.close();





   
    return 0;
}

