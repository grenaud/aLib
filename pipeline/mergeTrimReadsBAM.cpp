#include <iostream>
#include <string>
#include <cstring>
#include <sys/stat.h>

#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include "MergeTrimReads.h"
#include "PutProgramInHeader.h"

// #include "JSON.h"

#include "utils.h"

////////////////////////////////
// TO DO
//
////////////////////////////////

using namespace std;
//using namespace BamTools;
// using namespace __MergeTrimReads__;



// void initializeDefaultSequences(string configFile){
//     string line;
//     ifstream myFile;
//     string content="";
//     myFile.open(configFile.c_str(), ios::in);

//     if (myFile.is_open()){
// 	while ( getline (myFile,line)){
// 	    content+=line;
// 	}
// 	myFile.close();
//     }else{
// 	cerr << "Unable to open config file "<<configFile<<endl;
// 	exit(1);
//     }


//     JSONValue *value = JSON::Parse(content.c_str());
//     if (value == NULL){
// 	cerr<<"Failed to parse JSON file"<<endl;
// 	exit(1);
//     }

    
// }
string options_adapter_F_BAM="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG";
string options_adapter_S_BAM="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
string options_adapter_chimera_BAM="ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA";

inline bool isBamAlignEmpty(const BamAlignment & toTest){
    return ( toTest.Name.empty() &&
	     toTest.Length == 0 );    
}


int main (int argc, char *argv[]) {


    bool produceUnCompressedBAM=false;
    bool verbose=false;
    bool ancientDNA=false;
    bool keepOrig=false;

    string adapter_F=options_adapter_F_BAM;
    string adapter_S=options_adapter_S_BAM;
    string adapter_chimera=options_adapter_chimera_BAM;
    string key="";
    bool allowMissing=false;
    int trimCutoff=1;

    bool allowAligned=false;
    bool printLog=false;
    string logFileName;

    BamReader reader;
    BamWriter writer;

    string bamFile;
    string bamFileOUT="";

    string key1;
    string key2;
    
    const string usage=string(string(argv[0])+
			      "This program takes a BAM where each mate are consecutive and\ntrims and merges reads\n"+
			      +" [options] BAMfile"+"\n"+
			      //"\t"+"-p , --PIPE"+"\n\t\t"+"Read BAM from and write it to PIPE"+"\n"+
			      "\t"+"-o , --outfile" +"\t\t"+"Output (BAM format)."+"\n"+
			      "\t"+"-u            " +"\t\t"+"Produce uncompressed bam (good for pipe)"+"\n"+

			      //	"\t"+" , --outprefix" +"\n\t\t"+"Prefix for output files (default '"+outprefix+"')."+"\n"+
			      //"\t"+" , --SAM" +"\n\t\t"+"Output SAM not BAM."+"\n"+
			      "\t"+"--aligned" +"\t\t"+"Allow reads to be aligned (default "+boolStringify(allowAligned)+")"+"\n"+
			      "\t"+"-v , --verbose" +"\t\t"+"Turn all messages on (default "+boolStringify(verbose)+")"+"\n"+
			      "\t"+"--log [log file]" +"\t"+"Print a tally of merged reads to this log file (default only to stderr)"+"\n"+
			      
			      "\n\t"+"Paired End merging/Single Read trimming  options"+"\n"+
			      "\t\t"+"--ancientdna"+"\t\t\t\t"+"ancient DNA (default "+boolStringify(ancientDNA)+")"+"\n"+
			      "\t\t"+"            "+"\t\t\t\t"+" Allows for partial overlap"+"\n"+

			      "\t\t\t\t\t\t\tGood for merging ancient DNA reads into a single sequence\n\n"
			      "\t\t"+"--keepOrig"+"\t\t\t\t"+"Write original reads if they are trimmed or merged  (default "+boolStringify(keepOrig)+")"+"\n"+
			      "\t\t\t\t\t\t\tSuch reads will be marked as PCR duplicates\n\n"
			      "\t\t"+"-f , --adapterFirstRead" +"\t\t\t"+"Adapter that is observed after the forward read (def. Multiplex: "+options_adapter_F_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"-s , --adapterSecondRead" +"\t\t"+"Adapter that is observed after the reverse read (def. Multiplex: "+options_adapter_S_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"-c , --FirstReadChimeraFilter" +"\t\t"+"If the forward read looks like this sequence, the cluster is filtered out.\n\t\t\t\t\t\t\tProvide several sequences separated by comma (def. Multiplex: "+options_adapter_chimera_BAM.substr(0,30)+")"+"\n"+
			      "\t\t"+"-k , --key"+"\t\t\t\t"+"Key sequence with which each sequence starts. Comma separate for forward and reverse reads. (default '"+key+"')"+"\n"+
			      "\t\t"+"-i , --allowMissing"+"\t\t\t"+"Allow one base in one key to be missing or wrong. (default "+boolStringify(allowMissing)+")"+"\n"+
			      "\t\t"+"-t , --trimCutoff"+"\t\t\t"+"Lowest number of adapter bases to be observed for single Read trimming (default "+stringify(trimCutoff)+")");

    if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
    	cout<<"Usage:"<<endl;
    	cout<<""<<endl;
    	cout<<usage<<endl;
    	return 1;
    }

    for(int i=1;i<(argc-1);i++){ //all but the last arg

	if(strcmp(argv[i],"--log") == 0 ){
	    logFileName =string(argv[i+1]);
	    printLog=true;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-p") == 0 || strcmp(argv[i],"--PIPE") == 0 ){
	    cerr<<"This version no longer works with pipe, exiting"<<endl;
	    return 1;	    
	}

	if(strcmp(argv[i],"-u") == 0  ){
	    produceUnCompressedBAM=true;
	    continue;
	}

	if(strcmp(argv[i],"--aligned") == 0  ){
	    allowAligned=true;
	    continue;
	}



	if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	    bamFileOUT =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0 ){
	    verbose=true;
	    continue;
	}

	if(strcmp(argv[i],"--ancientdna") == 0 ){
	    ancientDNA=true;
	    continue;
	}

	if(strcmp(argv[i],"--keepOrig") == 0 ){
	    keepOrig=true;
	    continue;
	}

	if(strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--adapterFirstRead") == 0 ){
	    adapter_F =string(argv[i+1]);
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--adapterSecondRead") == 0 ){
	    adapter_S =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"--FirstReadChimeraFilter") == 0 ){
	    adapter_chimera =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-k") == 0 || strcmp(argv[i],"--keys") == 0 ){
	    key =string(argv[i+1]);
	    i++;
	    continue;
	}
	

	if(strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--allowMissing") == 0 ){
	    allowMissing=true;
	    continue;
	}

	if(strcmp(argv[i],"-t") == 0 || strcmp(argv[i],"--trimCutoff") == 0 ){
	    trimCutoff=atoi(argv[i+1]);
	    i++;
	    continue;
	}
	
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;	    
    }

    bamFile=argv[argc-1];
   //  initMerge();
//     set_adapter_sequences(adapter_F,
// 			  adapter_S,
// 			  adapter_chimera);
//     set_options(trimCutoff,allowMissing,mergeoverlap);
    if(key != ""){
	size_t found=key.find(",");
	if (found == string::npos){ //single end reads
	    key1=key;
	    key2="";
	} else{                     //paired-end
	    key1=key.substr(0,found);
	    key2=key.substr(found+1,key.length()-found+1);
	}
    }

    if( bamFileOUT == ""  ){
	cerr<<"The output must be a be specified, exiting"<<endl;
	return 1;
    }

    if ( !reader.Open(bamFile) ) {
    	cerr << "Could not open input BAM file  "<<bamFile << endl;
    	return 1;
    }
    SamHeader header = reader.GetHeader();

    

    string pID          = "mergeTrimReadsBAM";   
    string pName        = "mergeTrimReadsBAM";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));

    const RefVector references = reader.GetReferenceData();
    //we will not call bgzip with full compression, good for piping into another program to 
    //lessen the load on the CPU
    if(produceUnCompressedBAM) 
	writer.SetCompressionMode(BamWriter::Uncompressed);

    if ( !writer.Open(bamFileOUT,header,references) ) {
    	cerr << "Could not open output BAM file "<<bamFileOUT << endl;
    	return 1;
    }



    SamHeader sh=reader.GetHeader();
    //Up to the user to be sure that a sequence is followed by his mate
    // if(!sh.HasSortOrder() || 
    //    sh.SortOrder != "queryname"){
    // 	cerr << "Bamfile must be sorted by queryname" << endl;
    // 	return 1;
    // }
    
    MergeTrimReads mtr (adapter_F,adapter_S,adapter_chimera,
			key1,key2,
			trimCutoff,allowMissing,ancientDNA);


    BamAlignment al;
    BamAlignment al2;
    bool al2Null=true;
    
    while ( reader.GetNextAlignment(al) ) {

	
	if(al.IsMapped() || al.HasTag("NM") || al.HasTag("MD")  ){
	    if(!allowAligned){
		cerr << "Reads should not be aligned" << endl;
		return 1;
	    }else{
		//should we remove tags ?
	    }
	}


	if(al.IsPaired() && 
	   al2Null ){
	    al2=al;
	    al2Null=false;
	    continue;
	}else{
	    if(al.IsPaired() && 
	       !al2Null){

		bool  result =  mtr.processPair(al,al2);
		
		if( result ){//was merged
		    BamAlignment orig;
		    BamAlignment orig2;

		    if(keepOrig){
			orig2 = al2;
			orig  = al;
		    }

		    writer.SaveAlignment(al);

		    if(keepOrig){
			orig.SetIsDuplicate(true);
			orig2.SetIsDuplicate(true);
			writer.SaveAlignment(orig2);
			writer.SaveAlignment(orig);
		    }

		    //the second record is empty
		}else{
		    //keep the sequences as pairs

		    writer.SaveAlignment(al2);		    
		    writer.SaveAlignment(al);
		}

		//
		//  SINGLE END
		//
	    }else{ 
		BamAlignment orig;
		if(keepOrig){
		    orig =al;
		}
		mtr.processSingle(al);

		if(keepOrig){
		    //write duplicate
		    if(orig.QueryBases.length()  != al.QueryBases.length()){
			orig.SetIsDuplicate(true);
			writer.SaveAlignment(orig);
		    }
		}
		writer.SaveAlignment(al);



	    } //end single end
	    al2Null=true;
	}//second pair
		    

    } //while al
    reader.Close();
    writer.Close();

    cerr <<mtr.reportSingleLine()<<endl;

    if(printLog){
	ofstream fileLog;
	fileLog.open(logFileName.c_str());

	if (fileLog.is_open()){
	    fileLog <<mtr.reportMultipleLines() <<endl;

	}else{
	    cerr << "Unable to print to file "<<logFileName<<endl;
	}
	fileLog.close();
    }
   
    return 0;
}

