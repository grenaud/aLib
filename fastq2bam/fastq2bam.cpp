/*
 * fastq2bam
 * Date: Feb-07-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>

#include "PutProgramInHeader.h"
#include "utils.h"
#include "FastQObj.h"
#include "FastQParser.h"

const uint32_t flagSingleReads =  4; // 00000100
const uint32_t flagFirstPair   = 77; // 01001101
const uint32_t flagSecondPair  =141; // 10001101

using namespace std;

int main (int argc, char *argv[]) {

    string bamoutfile;
    string index1    ="";
    string index2    ="";
    string readgroup ="";

    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<"This program converts fastq files into unaligned bam\nUsage: "<<string(argv[0])<<" [options] [fastq in r1] [fastq in r2]"<<endl;
        cout<<"Options:"<<endl;
        cout<<"\tMandatory:"<<endl;
        cout<<"\t\t-o [output bam]"<<endl;
        cout<<"\tOptional:"<<endl;
        cout<<"\t\t-i1 [index1]"<<endl;
        cout<<"\t\t-i2 [index2]"<<endl;
        cout<<"\t\t-r  [read group]"<<endl;

        return 1;
    }


    for(int i=1;i<(argc-2);i++){ //all but the last two arg
	if(strcmp(argv[i],"-o") == 0  ){
            bamoutfile=string(argv[i+1]);
            i++;
            continue;
        }

	if(strcmp(argv[i],"-i1") == 0  ){
            index1=string(argv[i+1]);
            i++;
            continue;
        }

	if(strcmp(argv[i],"-i2") == 0  ){
            index2=string(argv[i+1]);
            i++;
            continue;
        }

	if(strcmp(argv[i],"-r") == 0  ){
            readgroup=string(argv[i+1]);
            i++;
            continue;
        }

    }

    if(bamoutfile.empty()){
	cerr << "ERROR -o option is mandatory " << endl;
        return 1;
    }

    string fastqin1=argv[argc-2];
    string fastqin2=argv[argc-1];

    BamWriter writer;

    SamHeader header ;
    string pID          = string(argv[0]);
    string pName        = string(argv[0]);
    string pCommandLine = "";
    
    for(int i=0;i<(argc);i++){
        pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine);

    //adding @RG in read group
    if(!readgroup.empty()){
	
	SamReadGroupDictionary  srgd;
	SamReadGroup srg ( readgroup );	
	srgd.Add( srg );       
	header.ReadGroups=srgd;
    }

    RefVector references;

    if ( !writer.Open(bamoutfile,header,references) ) {
        cerr << "Could not open output BAM file "<<bamoutfile << endl;
        return 1;
    }
    
    FastQParser fqp1 (fastqin1);
    FastQParser fqp2 (fastqin2);
    unsigned int totalSeqs=0;
    while(fqp1.hasData()){
	if(!fqp2.hasData()){
	    cerr << "ERROR: Discrepency between fastq files " << endl;
	    return 1;
	}

	FastQObj * fo1=fqp1.getData();
	FastQObj * fo2=fqp2.getData();
	vector<string> def1=allTokens( *(fo1->getID()), ' '  );
	vector<string> def2=allTokens( *(fo2->getID()), ' ' );
	if(def1[0] != def2[0]){
	    cerr << "ERROR: Discrepency between fastq files, different names " << *(fo1->getID()) <<" and "<< *(fo2->getID()) <<endl;
	    return 1;
	}
	
	BamAlignment toWrite1;
	BamAlignment toWrite2;
	toWrite1.Name=def1[0];
	toWrite2.Name=def2[0];

	toWrite1.AlignmentFlag=flagFirstPair;
	toWrite2.AlignmentFlag=flagSecondPair;

	toWrite1.MapQuality=0;
	toWrite2.MapQuality=0;
		  
	toWrite1.QueryBases =  *(fo1->getSeq());
	toWrite1.Qualities  =  *(fo1->getQual());
		  
	toWrite2.QueryBases =  *(fo2->getSeq());
	toWrite2.Qualities  =  *(fo2->getQual());
	
	//add tags for indices and fake qualities for the indices
	if(!index1.empty()){
	    if(!toWrite1.AddTag("XI", "Z",index1) )                      {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite1.AddTag("YI", "Z",string(index1.length(),'!')) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite2.AddTag("XI", "Z",index1) )                      {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite2.AddTag("YI", "Z",string(index1.length(),'!')) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	}

	if(!index2.empty()){
	    if(!toWrite1.AddTag("XI", "Z",index2) )                      {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite1.AddTag("YI", "Z",string(index2.length(),'!')) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite2.AddTag("XI", "Z",index2) )                      {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite2.AddTag("YI", "Z",string(index2.length(),'!')) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	}

	if(!readgroup.empty()){
	    if(!toWrite1.AddTag("RG", "Z",readgroup) )                      {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	    if(!toWrite2.AddTag("RG", "Z",readgroup) )                      {cerr<<"Internal error, cannot add tag"<<endl; return 1; }
	}

	
	
	writer.SaveAlignment(toWrite1);
	writer.SaveAlignment(toWrite2);
		    
	totalSeqs++;
    }


    writer.Close();

    cerr<<"Wrote "<<totalSeqs<<" sequences successfully"<<endl;


    return 0;
}

