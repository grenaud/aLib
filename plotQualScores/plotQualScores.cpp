#include <iostream>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

using namespace std;
using namespace BamTools;

short qualOffset=33;
// inline bool validBP(char c){
//     if(c == 'G')
// 	return true;
//     if(c == 'C')
// 	return true;
//     if(c == 'A')
// 	return true;
//     if(c == 'T')
// 	return true;
//     return false;
// }
short minQualScore=0;
short maxQualScore=50;

int main (int argc, char *argv[]) {
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<"plotQualScore input.bam"<<endl;
	return 1;
    }

    string bamfiletopen = string(argv[1]);
    BamReader reader;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    // if ( !reader.LocateIndex() ){
    // 	cerr << "warning: cannot locate index for file " << bamfiletopen<<endl;
    // 	//return 1;
    // }

    BamAlignment al;
    BamAlignment al2;
    
    bool unsurePEorSE=true;
    bool pe=true;
    int strLength=-1;
    int vecLengthToUse;

    map<short,unsigned long>  ** counterA;
    map<short,unsigned long>  ** counterC;
    map<short,unsigned long>  ** counterG;
    map<short,unsigned long>  ** counterT;
    
    int lengthIndex1=0;
    int lengthIndex2=0;
    string seqInd1;
    string seqInd2;
    string qualInd1;
    string qualInd2;
    int offsetInd2;

    while ( reader.GetNextAlignment(al) ) {
	if(unsurePEorSE){
	    strLength=al.QueryBases.length();
	    if(al.IsPaired()){
		pe=true;
		vecLengthToUse=2*strLength;		
	    }else{
		pe=false;
		vecLengthToUse=strLength;
	    }
	    string index1;
	    string index2;
	
	    if(al.HasTag("XI")){
		al.GetTag("XI",index1);
		vecLengthToUse+=index1.length();
		lengthIndex1=index1.length();
	    }

	    if(al.HasTag("XJ")){
		al.GetTag("XJ",index2);
		vecLengthToUse+=index2.length();
		lengthIndex2=index2.length();
	    }

	    counterA     = new map<short,unsigned long>  * [vecLengthToUse];
	    counterC     = new map<short,unsigned long>  * [vecLengthToUse];
	    counterG     = new map<short,unsigned long>  * [vecLengthToUse];
	    counterT     = new map<short,unsigned long>  * [vecLengthToUse];
	    for(int i=0;i<vecLengthToUse;i++){
		counterA[i]=new map<short,unsigned long>  ();
		counterC[i]=new map<short,unsigned long>  ();
		counterG[i]=new map<short,unsigned long>  ();
		counterT[i]=new map<short,unsigned long>  ();
		for(short k=minQualScore;k<=maxQualScore;k++){
		    (*counterA[i])[k]=0;
		    (*counterC[i])[k]=0;
		    (*counterG[i])[k]=0;
		    (*counterT[i])[k]=0;
		}
	    }
	    unsurePEorSE=false;
	}else{
	    if(pe  &&
	       !al.IsPaired()){
		cerr << "Cannot have unpaired reads in PE mode" << endl;
		return 1;
	    }

	    if(!pe  &&
	       al.IsPaired()){
		cerr << "Cannot have unpaired reads in SE mode" << endl;
		return 1;
	    }
	}
	
	if(al.QueryBases.length() !=  al.Qualities.length()){
	    cerr << "Cannot have different lengths for sequence and quality" << endl;
	    return 1;
	}
	if(al.QueryBases.length() !=  strLength){
	    cerr << "Cannot have different lengths for sequence and quality" << endl;
	    return 1;
	}

	if(pe){
	    if(al.IsFirstMate()){
		reader.GetNextAlignment(al2);
		if(al2.QueryBases.length() !=  al2.Qualities.length()){
		    cerr << "Cannot have different lengths for sequence and quality" << endl;
		    return 1;
		}

	    }else{
		cerr << "First read should be the first mate" << endl;
		return 1;	    
	    }
	}



	//cycle
	for(int i=0;i<al.QueryBases.length();i++){
	    short x=(short(al.Qualities[i])-qualOffset);
	    if(al.QueryBases[i] == 'A'){
	    	(*counterA[i])[x]++;
	    }
	    if(al.QueryBases[i] == 'C'){
	    	(*counterC[i])[x]++;
	    }
	    if(al.QueryBases[i] == 'G'){
	    	(*counterG[i])[x]++;
	    }
	    if(al.QueryBases[i] == 'T'){
	    	(*counterT[i])[x]++;
	    }
	}

	//The indices for al and al2 should hopefully be the same 
	if(lengthIndex1>0){
	    al.GetTag("XI",seqInd1);
	    al.GetTag("YI",qualInd1);
	    int j;

	    for(int i=0;i<lengthIndex1;i++){
		j=i+al.QueryBases.length();
		short x=(short(qualInd1[i])-qualOffset);
		if(seqInd1[i] == 'A'){
		    (*counterA[j])[x]++;
		}
		if(seqInd1[i] == 'C'){
		    (*counterC[j])[x]++;
		}
		if(seqInd1[i] == 'G'){
		    (*counterG[j])[x]++;
		}
		if(seqInd1[i] == 'T'){
		    (*counterT[j])[x]++;
		}
	    }
	}

	if(pe){
	    offsetInd2=al.QueryBases.length()+lengthIndex1+al2.QueryBases.length();
	    int j;
	    for(int i=0;i<al2.QueryBases.length();i++){
		j=i+al.QueryBases.length()+lengthIndex1;
		short x=(short(al2.Qualities[i])-qualOffset);
		if(al2.QueryBases[i] == 'A'){
		    (*counterA[j])[x]++;
		}
		if(al2.QueryBases[i] == 'C'){
		    (*counterC[j])[x]++;
		}
		if(al2.QueryBases[i] == 'G'){
		    (*counterG[j])[x]++;
		}
		if(al2.QueryBases[i] == 'T'){
		    (*counterT[j])[x]++;
		}
	    }
	}else{
	    offsetInd2=al.QueryBases.length()+lengthIndex1;
	}

	//The indices for al and al2 should hopefully be the same 
	if(lengthIndex2>0){
	    al.GetTag("XJ",seqInd2);
	    al.GetTag("YJ",qualInd2);
	    int j;

	    for(int i=0;i<lengthIndex2;i++){
		j=offsetInd2+i;
		short x=(short(qualInd2[i])-qualOffset);
		if(seqInd2[i] == 'A'){
		    (*counterA[j])[x]++;
		}
		if(seqInd2[i] == 'C'){
		    (*counterC[j])[x]++;
		}
		if(seqInd2[i] == 'G'){
		    (*counterG[j])[x]++;
		}
		if(seqInd2[i] == 'T'){
		    (*counterT[j])[x]++;
		}
	    }
	}

    }
    reader.Close();

    cout<<"cycle\t"<<"nuc\t";
    for(short k=minQualScore;k<maxQualScore;k++){
	cout<<k<<"\t";
    }
    cout<<maxQualScore<<endl;

    for(int i=0;i<vecLengthToUse;i++){
	cout<<(i+1)<<"\t";
	cout<<"A\t";
	for(short k=minQualScore;k<maxQualScore;k++){
	    cout<<(*counterA[i])[k]<<"\t";
	}
	cout<<(*counterA[i])[maxQualScore]<<endl;
	cout<<(i+1)<<"\t";
	cout<<"C\t";
	for(short k=minQualScore;k<maxQualScore;k++){
	    cout<<(*counterC[i])[k]<<"\t";
	}
	cout<<(*counterC[i])[maxQualScore]<<endl;
	cout<<(i+1)<<"\t";
	cout<<"G\t";
	for(short k=minQualScore;k<maxQualScore;k++){
	    cout<<(*counterG[i])[k]<<"\t";
	}
	cout<<(*counterG[i])[maxQualScore]<<endl;
	cout<<(i+1)<<"\t";
	cout<<"T\t";
	for(short k=minQualScore;k<maxQualScore;k++){
	    cout<<(*counterT[i])[k]<<"\t";
	}
	cout<<(*counterT[i])[maxQualScore]<<endl;	
    }



    
    return 0;
}

