// vim:ts=8
#include <iostream>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <sstream>
#include <map>


#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include "RGAssign.h"
#include "PutProgramInHeader.h"


#include "utils.h"

using namespace std;
using namespace BamTools;


/****************************************/
/*                                      */
/*          STATIC VARIABLES            */
/*                                      */
/****************************************/



// ofstream * ratioValues;
// ofstream * rgqual;

static double rgScoreCutoff  = 80 ;             // two HQ bases can mismatch
static double fracConflict   = 20 ;             // top shall be 100x more likely
static double wrongness      = 30 ;             // flag wrong if the wrong pair is 1000x more likely
static int    mismatchesTrie = 2;
static int    maxErrorHits   = 20;



int main (int argc, char *argv[]) {

    BamReader reader;
    BamWriter writer;

    string bamFile;
    string bamFileOUT="";

    string index="";
    string outfile;
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

    bool produceUnCompressedBAM=false; 
    const string usage=string(string(argv[0])+
			      " [options] BAMfile"+"\n\n"+

			      "\tCutoffs options:"+"\n"
			      "\t\t"+"--rgqual"  +"\t[quality]"+"\t\t"+""+"Worst quality before flagging as unknown ["+stringify(rgScoreCutoff)+"]\n"+
			      "\t\t"+"--fracconf"+"\t[quality]"+"\t\t"+""+"Maximum quality difference considered a conflict ["+stringify(fracConflict)+"] \n"+
			      "\t\t"+"--wrongness"+"\t[quality]"+"\t\t"+""+"Mininum quality difference to flag as wrongly paired ["+stringify(wrongness)+"] \n"+
			      "\t\t"+"--mm"+"\t\t[mismatches]"+"\t\t"+""+"Maximum # of tolerated mismatches ["+stringify(mismatchesTrie)+"] \n"+

			      "\n\tRG assignment options:"+"\n"+
			      "\t\t"+"" +""+"--shift"+"\t"+"\t\t\t\t"+"Try shifting the index right by one at the cost of a mismatch"+"\n"+


			      "\n\tOutput options:"+"\n"+
			      
			      "\t"+"\tMandatory:"+"\n"+
			      "\t\t"+"-i"+","+"--index"+"\t[index]"+"\t\t\t"+"File describing index sequences used"+"\n"+
			      "\t\t"+"-o"+","+"--outfile"+"\t[outfile]"+"\t\t"+"Specify output file"+"\n"+
			      "\t"+"\tOptional:"+"\n"+
			      
			      "\t\t"+"--maxerr"+"\t[max err]"+"\t\t"+""+"Print  # wrongly of assigned RG in the error log (--error) ["+stringify(maxErrorHits)+"] \n"+
                              "\t\t"+"-u" +"\t\t\t\t\t"           +"Produce uncompressed bam (good for pipe)"+"\n"+ 
			      "\t\t"+"-s"+","+"--summary"+"\t[summary file]"+"\t\t"+"Summarize the RG tally in this file"+"\n"+
			      "\t\t"+"-e"+","+"--error"  +"\t[error file]"+"\t\t"+"Summarize the indices that were not assigned to a RG"+"\n"+
			      "\t\t"+""+""+"--rgval"  +"\t[file]"+"\t\t\t\t"+"Write the rg qualities as a binary file"+"\n"+
			      "\t\t"+""+""+"--ratio"   +"\t\t[file]"+"\t\t\t"+"Write the likelihood ratios as a binary file"+"\n"

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

    for(int i=1;i<(argc-1);i++){


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

	if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	    outfile=string(argv[i+1]);
	    i++;
	    continue;
	}


	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;             
    }

    if(outfile.size() == 0){
	cerr<<"The field -o is mandatory exiting"<<endl;
	return 1;             
    }

    if(index.size() == 0){
	cerr<<"The field -i is mandatory exiting"<<endl;
	return 1;             
    }

    

    bamFile=argv[argc-1];




    if ( !reader.Open(bamFile) ) {
    	cerr << "Could not open input BAM file  "<<bamFile << endl;
    	return 1;
    }

    SamHeader  myHeader=reader.GetHeader();
    SamProgram sp;
   
    string pID          = "assignRG";   
    string pName        = "assignRG";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }

    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));



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

    SamReadGroupDictionary  srgd;
    map<string,string>::const_iterator itRG;   
    for ( itRG  = rga.getRGS()->begin(); 
	  itRG != rga.getRGS()->end(); 
	  itRG++ ){
	SamReadGroup srg ( itRG->first );	
	srg.Description  = itRG->second; //description read in index file
	srgd.Add( srg );       	


	if(itRG->first == "conflict" || itRG->first == "unknown" || itRG->first == "wrong" ){
	    cerr<<"ERROR: The RG names cannot contain the words: \"conflict\" or \"unknown\" or \"wrong\""<<endl;
	    return 1;
	}
    }


    myHeader.ReadGroups=srgd;
    if(produceUnCompressedBAM)  
	writer.SetCompressionMode(BamWriter::Uncompressed); 

    if( !writer.Open(outfile,myHeader,reader.GetReferenceData() ) ) {
    	cerr << "Could not open output BAM file  "<<outfile << endl;
    	return 1;	
    }


    BamAlignment al;
    BamAlignment al2;

    // damn, this logic is convoluted...
    while( reader.GetNextAlignment(al) ) {
	while(1) {
            if( !reader.GetNextAlignment(al2) ) {
                // EOF, process the one leftover record
                rga.processSingleEndReads(al);//,writer,printError,unknownSeq,wrongSeq,conflictSeq);
		writer.SaveAlignment(al);
                break; 
            }
            // If it's paired, both should have the same index, and we
            // save some work.  Since the reads are probably not
            // ordered, check the names first
            if( al.IsPaired() && al.Name == al2.Name ) {
                rga.processPairedEndReads(al,al2); //,writer,printError,unknownSeq,wrongSeq,conflictSeq);
		writer.SaveAlignment(al);
		writer.SaveAlignment(al2);
                break ;
            } else {
                // no match, treat one(!) separately
                rga.processSingleEndReads(al );//,writer,printError,unknownSeq,wrongSeq,conflictSeq);
		writer.SaveAlignment(al);
                swap(al,al2) ;
            }
        }
    }

    reader.Close();
    writer.Close();




    //Print summary of RG assignment
    if(printSummary){
  
	ofstream fileSummary;
	fileSummary.open(filenameSummary.c_str());

	if (fileSummary.is_open()){

	    fileSummary<<rga.getSummaryString()<<endl;
	}else{
	    cerr << "Unable to print to file "<<filenameSummary<<endl;
	}
	fileSummary.close();
    }




    //Print over-represented sequences in conflict,unknown,wrong
    if(printError){

	ofstream fileError;
	fileError.open(filenameError.c_str());
	if (fileError.is_open())
        {
	    fileError<<rga.getErrorString()<<endl;
	}else{
	    cerr << "Unable to print to file "<<filenameError<<endl;
	}
	fileError.close();
    }



    //cleaning up
   
    if(rgqualFlag)
	rgqualOS.close();

    if(ratioValuesFlag)
	ratioValuesOS.close();

    return 0;
}

