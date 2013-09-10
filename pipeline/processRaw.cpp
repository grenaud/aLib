#include <iostream>
#include <string>
#include <cstring>
#include <sys/stat.h>

#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>


#include "FilterBAMal.h"
#include "RGAssign.h"
#include "MergeTrimReads.h"

#include "PutProgramInHeader.h"

// #include "JSON.h"

#include "utils.h"

////////////////////////////////
// TO DO
//
////////////////////////////////

using namespace std;


static double rgScoreCutoff  = 80 ;             // two HQ bases can mismatch
static double fracConflict   = 20 ;             // top shall be 100x more likely
static double wrongness      = 30 ;             // flag wrong if the wrong pair is 1000x more likely
static int    mismatchesTrie = 2;
static int    maxErrorHits   = 20;



inline bool isBamAlignEmpty(const BamAlignment & toTest){
    return ( toTest.Name.empty() &&
	     toTest.Length == 0 );    
}

int main (int argc, char *argv[]) {

    BamReader reader;
    BamWriter writer;

    string bamFile;
    string bamFileOUT="";

    bool produceUnCompressedBAM=false;
    bool verbose=false;



    //MERGE VARS
    bool mergeoverlap=false;
    bool keepOrig=false;

    string adapter_F="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG";
    string adapter_S="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
    string adapter_chimera="ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA";

    string key="";
    bool allowMissing=false;
    int trimCutoff=1;

    bool allowAligned=false;
    bool printLog=false;
    string logFileName;

    string key1;
    string key2;


    
    //RG ASSIGN
    string index="";
    //     string outfile;
    bool   printSummary=false;
    string filenameSummary;

    bool   printError=false;
    string filenameError;

    ofstream ratioValuesOS;
    ofstream rgqualOS;

    bool ratioValuesFlag = false;
    bool rgqualFlag      = false;
    bool shiftByOne      = false;

    bool flag_ratioValues=false;
    bool flag_rgqual     =false;



    //FILTERING
    int minLength=-1;
    int maxLength=-1;
    //     string outfile;
    bool entropy      = false;
    bool frequency    = false;
    double compOrEntCutoff = 0.85;
    double cutoffLikelihood = 0.5;
    double cutoffAvgExpError = 0.01;
  
    ofstream likelihoodOS;
    bool     likelihoodFlag  = false;
    ofstream entropyOS;
    bool     entropyOSFlag  = false;
    ofstream frequencyOS;
    bool     frequencyOSFlag  = false;

    //     bool     verbose=false;
    //    bool usePercent=false;
    //    double bottomPercent=0.0;
    // bool definedCutoff=false;
    // bool definedExpCutoff=false;

    bool trimSeqs=false;
    //     bool produceUnCompressedBAM=false; 
    bool resetQC=false; 

    //report
    string reportFile="/dev/stderr";
    //     Report toprint;



    const string usage=string(string(argv[0])+
			      "This program takes a BAM file where each mate are consecutive and:\n\t- trims and merges reads\n\t-demultiplexes\n\t-filters low quality reads\n\n"+
			      "\nUsage:\n"+
			      " [options] BAMfile"+"\n"+
			      "\n\t"+"General I/O options"+"\n"+
			      //"\t"+"-p , --PIPE"+"\n\t\t"+"Read BAM from and write it to PIPE"+"\n"+
			      "\t"+"-o , --outfile" +"\t\t"+"Output (BAM format)."+"\n"+
			      "\t"+"-u            " +"\t\t"+"Produce uncompressed bam (good for pipe)"+"\n"+

			      //	"\t"+" , --outprefix" +"\n\t\t"+"Prefix for output files (default '"+outprefix+"')."+"\n"+
			      //"\t"+" , --SAM" +"\n\t\t"+"Output SAM not BAM."+"\n"+
			      "\t"+"--aligned" +"\t\t"+"Allow reads to be aligned (default "+boolStringify(allowAligned)+")"+"\n"+
			      "\t"+"-v , --verbose" +"\t\t"+"Turn all messages on (default "+boolStringify(verbose)+")"+"\n"+






			      
			      "\n\n\n\t"+"------------------------------------------------------------------------------------------------------------"+"\n"+

			      "\t"+"Paired End merging/Single Read trimming options"+"\n"+
			      "\t\t"+"--log\t\t\t[log file]" +"\t"+"Print a tally of merged reads to this log file (default only to stderr)"+"\n"+
			      "\t\t"+"--mergeoverlap"+"\t\t\t\t"+"Merge PE reads of molecules longer than read length that show a minimum overlap (default "+boolStringify(mergeoverlap)+")"+"\n"+
			      "\t\t\t\t\t\t\tGood for merging ancient DNA reads into a single sequence\n\n"
			      "\t\t"+"--keepOrig"+"\t\t\t\t"+"Write original reads if they are trimmed or merged  (default "+boolStringify(keepOrig)+")"+"\n"+
			      "\t\t\t\t\t\t\tSuch reads will be marked as PCR duplicates\n\n"
			      "\t\t"+"-f , --adapterFirstRead" +"\t\t\t"+"Adapter that is observed after the forward read (def. Multiplex: "+adapter_F.substr(0,30)+")"+"\n"+
			      "\t\t"+"-s , --adapterSecondRead" +"\t\t"+"Adapter that is observed after the reverse read (def. Multiplex: "+adapter_S.substr(0,30)+")"+"\n"+
			      "\t\t"+"-c , --FirstReadChimeraFilter" +"\t\t"+"If the forward read looks like this sequence, the cluster is filtered out.\n\t\t\t\t\t\t\tProvide several sequences separated by comma (def. Multiplex: "+adapter_chimera.substr(0,30)+")"+"\n"+
			      "\t\t"+"-k , --key"+"\t\t\t\t"+"Key sequence with which each sequence starts. Comma separate for forward and reverse reads. (default '"+key+"')"+"\n"+
			      "\t\t"+"-i , --allowMissing"+"\t\t\t"+"Allow one base in one key to be missing or wrong. (default "+boolStringify(allowMissing)+")"+"\n"+
			      "\t\t"+"-t , --trimCutoff"+"\t\t\t"+"Lowest number of adapter bases to be observed for single Read trimming (default "+stringify(trimCutoff)+")"+







			      "\n\n\n\t"+"------------------------------------------------------------------------------------------------------------"+"\n"+

			      "\t"+"Demultiplexing options"+"\n"+

			      "\t\tCutoffs options:"+"\n"
			      "\t\t\t"+"--rgqual"  +"\t[quality]"+"\t\t"+""+"Worst quality before flagging as unknown ["+stringify(rgScoreCutoff)+"]\n"+
			      "\t\t\t"+"--fracconf"+"\t[quality]"+"\t\t"+""+"Maximum quality difference considered a conflict ["+stringify(fracConflict)+"] \n"+
			      "\t\t\t"+"--wrongness"+"\t[quality]"+"\t\t"+""+"Mininum quality difference to flag as wrongly paired ["+stringify(wrongness)+"] \n"+
			      "\t\t\t"+"--mm"+"\t\t[mismatches]"+"\t\t"+""+"Maximum # of tolerated mismatches ["+stringify(mismatchesTrie)+"] \n"+

			      "\n\t\tRG assignment options:"+"\n"+
			      "\t\t\t"+"" +""+"--shift"+"\t"+"\t\t\t\t"+"Try shifting the index right by one at the cost of a mismatch"+"\n"+


			      "\n\t\tOutput options:"+"\n"+
			      
			      "\t\t"+"\tMandatory:"+"\n"+
			      "\t\t\t"+"-i"+","+"--index"+"\t[index]"+"\t\t\t"+"File describing index sequences used"+"\n"+
// 			      "\t\t\t"+"-o"+","+"--outfile"+"\t[outfile]"+"\t\t"+"Specify output file"+"\n"+
			      "\t\t"+"\tOptional:"+"\n"+
			      
			      "\t\t\t"+"--maxerr"+"\t[max err]"+"\t\t"+""+"Print  # wrongly of assigned RG in the error log (--error) ["+stringify(maxErrorHits)+"] \n"+
//                               "\t\t"+"-u" +"\t\t\t\t\t"           +"Produce uncompressed bam (good for pipe)"+"\n"+ 
			      "\t\t\t"+"-s"+","+"--summary"+"\t[summary file]"+"\t\t"+"Summarize the RG tally in this file"+"\n"+
			      "\t\t\t"+"-e"+","+"--error"  +"\t[error file]"+"\t\t"+"Summarize the indices that were not assigned to a RG"+"\n"+
			      "\t\t\t"+""+""+"--rgval"  +"\t\t[file]"+"\t\t\t\t"+"Write the rg qualities as a binary file"+"\n"+
			      "\t\t\t"+""+""+"--ratio"   +"\t\t[file]"+"\t\t\t"+"Write the likelihood ratios as a binary file"+"\n"+







			      "\n\n\n\t"+"------------------------------------------------------------------------------------------------------------"+"\n"+

			      "\t"+"Filtering options"+"\n"+


			      "\n\t\tInput options:"+"\n"+
			      "\t\t\t"+""+""+"-r"  +"\t\t"+"\t\t\t"+"Reset QC flags"+"\n"+
    
			      "\n\t\tOutput options:"+"\n"+
// 			      "\t\t\t"+""+""+"-u"  +"\t\t"+"\t\t\t"+"Produce uncompressed bam (good for pipe)"+"\n"+
			      "\t\t\t"+""+""+"--like"  +"\t\t[file]"+"\t\t\t"+"Write the sequence likelihoods as a binary file"+"\n"+
			      "\t\t\t"+""+""+"--ent"  +"\t\t[file]"+"\t\t\t"+"Write the sequence entropy as a binary file"+"\n"+
			      "\t\t\t"+""+""+"--freq"  +"\t\t[file]"+"\t\t\t"+"Write the sequence frequency complexity as a binary file"+"\n"+
			      "\t\t\t"+"-v"+","+"--verbose"  +"\t"+"\t\t\t"+"Print info on the stderr (Default: "+booleanAsString(verbose)+")\n"+
			      "\t\t\t"+""+" "+"--log"  +"\t"+"\t[log file]"+"\t\t"+"Print a report to this file (Default: stderr)\n"+

// 			      "\n\t\t"+"\tMandatory:"+"\n"+
// 			      "\t\t\t"+"-o"+" "+"--outfile"+"\t[outfile]"+"\t\t"+"Specify output file"+"\n"+

			      "\n\t\tFiltering options:"+"\n"+
			      "\t\t\t"+"-c" +","+"--cutoff"+"\t[cutoff]""\t\t"+"Sequence likelihood cutoff (Default: "+stringify(cutoffLikelihood)+")"+"\n"+
			      "\t\t\t"+"" +""+"--cutexp"+"\t[cutoff]""\t\t"+"Average of expectancy of base error cutoff (Default: "+stringify(cutoffAvgExpError)+")"+"\n"+

			      "\t\t\t"+"" +""+"--trim"+"\t\t\t\t\t"+"Try to trim from the 3' end (TO IMPLEMENT) (Default: "+booleanAsString(trimSeqs)+")"+"\n"+


			      "\t\t\t"+"" +""+"--min_length"+"\t[cutoff]"+"\t\t"+"Flag any sequence with strickly less than this min length as failed (Default: "+stringify(minLength)+")"+"\n"+
			      "\t\t\t"+"" +""+"--max_length"+"\t[cutoff]"+"\t\t"+"Flag any sequence with strickly more than this max length as failed (Default: "+stringify(maxLength)+")"+"\n"+
			      "\t\t\t"+"" +""+"--percent"+"\t[percentage]"+"\t\t"+"Flag any sequence the bottom % as failed, to use with only small datasets (Default: not used)"+"\n"+

			      "\n\t\tComplexity filter options:"+"\n"+
			      
			      "\t\t\t"+"-e"+","+"--entropy"     +"\t"+"\t\t\t"+"Apply sequence entropy filter (Default: "+booleanAsString(entropy)+")\n"+
			      "\t\t\t"+"-f"+","+"--frequency"   +"\t"+"\t\t\t"+"Apply base frequency filter  (Default: "+booleanAsString(frequency)+")\n"+
			      "\t\t\t"+""  +"" +"--comp_cutoff" +"\t"+"\t\t\t"+"Entropy value [0.0-2.0] or fraction [0.0-1.0] of most frequent base accepted (Default: "+stringify(compOrEntCutoff)+")\n"

   



);

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
	// General I/O options	
	if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	    bamFileOUT =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-u") == 0  ){
	    produceUnCompressedBAM=true;
	    continue;
	}

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


	if(strcmp(argv[i],"--aligned") == 0  ){
	    allowAligned=true;
	    continue;
	}




	if(strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0 ){
	    verbose=true;
	    continue;
	}		
	
	//"\t"+"Paired End merging/Single Read trimming options"+"\n"+


	if(strcmp(argv[i],"--mergeoverlap") == 0 ){
	    mergeoverlap=true;
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


	// "\t"+"Demultiplexing options"+"\n"+
	
	
	if(strcmp(argv[i],"--shift") == 0 ){
	    shiftByOne      = true;
	    continue;
	}

	if(strcmp(argv[i],"--rgval") == 0 ){
	    string temp =string(argv[i+1]);
	    rgqualOS.open(temp.c_str(), ios::out | ios::binary);
	    rgqualFlag      = true;
	    if (!rgqualOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    //setFileForRGQual(&rgqualOS);
	    flag_rgqual=true;
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"-u") == 0  ){ 
	    produceUnCompressedBAM=true; 
	    continue; 
	} 
	

	if(strcmp(argv[i],"--ratio") == 0 ){
	    string temp =string(argv[i+1]);
	    ratioValuesOS.open(temp.c_str(), ios::out | ios::binary);
	    ratioValuesFlag = true;

	    if (!ratioValuesOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    //setFileForRatio(&ratioValuesOS);
	    flag_ratioValues=true;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-e") == 0 || strcmp(argv[i],"--error") == 0 ){
	    printError=true;
	    filenameError =string(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--summary") == 0 ){
	    printSummary=true;
	    filenameSummary =string(argv[i+1]);
	    i++;
	    continue;
	}





	if(strcmp(argv[i],"--maxerr") == 0 ){
	    maxErrorHits =destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"--rgqual") == 0 ){
	    rgScoreCutoff =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"--fracconf") == 0 ){
	    fracConflict =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--wrongness") == 0 ){
	    wrongness =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--mm") == 0 ){
	    mismatchesTrie =destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--index") == 0 ){
	    index =string(argv[i+1]);
	    i++;
	    continue;
	}

	// if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	//     outfile=string(argv[i+1]);
	//     i++;
	//     continue;
	// }






	//"\t"+"Filtering options"+"\n"+

	if(strcmp(argv[i],"--freq") == 0 ){
	    string temp =string(argv[i+1]);
	    frequencyOS.open(temp.c_str(), ios::out | ios::binary);
	    frequencyOSFlag      = true;
	    if (!frequencyOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-u") == 0  ){ 
	    produceUnCompressedBAM=true; 
	    continue; 
	} 

	if(strcmp(argv[i],"-r") == 0  ){ 
	    resetQC=true; 
	    continue; 
	} 

	if(strcmp(argv[i],"--trim") == 0  ){ 
	    cerr<<"to implement"<<endl;
	    return 1;
	} 
	
	if(strcmp(argv[i],"--ent") == 0 ){
	    string temp =string(argv[i+1]);
	    entropyOS.open(temp.c_str(), ios::out | ios::binary);
	    entropyOSFlag      = true;
	    if (!entropyOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--like") == 0 ){
	    string temp =string(argv[i+1]);
	    likelihoodOS.open(temp.c_str(), ios::out | ios::binary);
	    likelihoodFlag      = true;
	    if (!likelihoodOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--log") == 0 ){
	    reportFile =string(argv[i+1]);
	    i++;
	    continue;
	}

	// if(  (strcmp(argv[i],"--percent") == 0)   ){
	//     bottomPercent =destringify<double>(argv[i+1]);
	//     usePercent=true;
	//     i++;
	//     continue;
	// }

	// if( (strcmp(argv[i],"-c") == 0) || (strcmp(argv[i],"--cutoff") == 0)   ){
	//     cutoffLikelihood =destringify<double>(argv[i+1]);
	//     definedCutoff=true;
	//     i++;
	//     continue;
	// }

	// if(  (strcmp(argv[i],"--cutexp") == 0)   ){
	//     cutoffAvgExpError =destringify<double>(argv[i+1]);
	//     definedExpCutoff=true;
	//     i++;
	//     continue;
	// }


	if( (strcmp(argv[i],"-e") == 0) || (strcmp(argv[i],"--entropy") == 0)   ){
	    entropy=true;
	    continue;
	}

	if( (strcmp(argv[i],"-f") == 0) || (strcmp(argv[i],"--frequency") == 0) ){
	    frequency=true;
	    continue;
	}



	if(strcmp(argv[i],"--min_length") == 0 ){
	    minLength =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--max_length") == 0 ){
	    maxLength =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	// if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	//     outfile=string(argv[i+1]);
	//     i++;
	//     continue;
	// }


	if( (strcmp(argv[i],"-v") == 0) || (strcmp(argv[i],"--verbose") == 0)   ){
	    verbose=true;
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

    

    string pID          = "processRaw";   
    string pName        = "processRaw";   
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
			trimCutoff,allowMissing,mergeoverlap);

    RGAssign rga =    RGAssign(  rgScoreCutoff  ,
				 fracConflict   ,
				 wrongness      ,

				 mismatchesTrie ,
				 maxErrorHits ,
				 shiftByOne,
				 
				 printSummary ,
				 printError ,
	       
				 flag_ratioValues,
				 flag_rgqual,

				 &rgqualOS,
				 &ratioValuesOS,


				 index
				 ) ;

    FilterBAMal fm (minLength,maxLength,cutoffAvgExpError,frequency,entropy,compOrEntCutoff,
		   likelihoodFlag,&likelihoodOS, entropyOSFlag, &entropyOS, frequencyOSFlag, &frequencyOS,verbose,resetQC);

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

		    rga.processSingleEndReads(al);//,writer,printError,unknownSeq,wrongSeq,conflictSeq);
		    fm.filterBAMAlign(&al);

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
		    rga.processPairedEndReads(al,al2); //,writer,printError,unknownSeq,wrongSeq,conflictSeq);
		    fm.filterBAMAlign(&al);
		    fm.filterBAMAlign(&al2);

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
		rga.processSingleEndReads(al);//,writer,printError,unknownSeq,wrongSeq,conflictSeq);
		fm.filterBAMAlign(&al);

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
   
    if(rgqualFlag)
	rgqualOS.close();

    if(ratioValuesFlag)
	ratioValuesOS.close();

    return 0;
}

