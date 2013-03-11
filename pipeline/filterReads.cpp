/*
 * filterReads
 * Date: Nov-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>
#include <string>
#include <cstring>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "FilterBAMal.h"
#include "PutProgramInHeader.h"

#include "utils.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {


    int minLength=-1;
    int maxLength=-1;
    string outfile;
    bool entropy      = false;
    bool frequency    = false;
    double compOrEntCutoff = 0.85;
    double cutoffLikelihood = 0.5;
    double cutoffAvgExpError = 0.001;
  
    ofstream likelihoodOS;
    bool     likelihoodFlag  = false;
    ofstream entropyOS;
    bool     entropyOSFlag  = false;
    ofstream frequencyOS;
    bool     frequencyOSFlag  = false;
    bool     verbose=false;
    bool usePercent=false;
    double bottomPercent=0.0;
    bool definedCutoff=false;
    bool definedExpCutoff=false;

    bool trimSeqs=false;
    bool produceUnCompressedBAM=false; 

    const string usage=string("\t"+string(argv[0])+
			      " [options] BAMfile"+"\n\n"+
			      // "\tCutoffs options:"+"\n"
			      // "\t\t"+"--rgqual"+"\t[quality]"+"\t\t"+""+"Cutoffs for read group assignment quality (default:"+stringify(rgScoreCutoff)+")\n"+
			      // "\t\t"+"--fracconf"+"\t[fraction]"+"\t\t"+""+"Fraction of the second best RG probablity to the first to be\n\t\t\t\t\t\t\tlabeled conflict (default:"+stringify(fracConflict)+") \n"+
			      // "\t\t"+"--mm"+"\t\t[mismatches]"+"\t\t"+""+"Maximum # of tolerated mismatches (default:"+stringify(mismatchesTrie)+") \n"+


    
			      "\n\tOutput options:"+"\n"+
			      "\t\t"+""+""+"-u"  +"\t\t"+"\t\t\t"+"Produce uncompressed bam (good for pipe)"+"\n"+
			      "\t\t"+""+""+"--like"  +"\t\t[file]"+"\t\t\t"+"Write the sequence likelihoods as a binary file"+"\n"+
			      "\t\t"+""+""+"--ent"  +"\t\t[file]"+"\t\t\t"+"Write the sequence entropy as a binary file"+"\n"+
			      "\t\t"+""+""+"--freq"  +"\t\t[file]"+"\t\t\t"+"Write the sequence frequency complexity as a binary file"+"\n"+
			      "\t\t"+"-v"+""+"--verbose"  +"\t"+"\t\t\t"+"Print info on the stderr (Default: "+booleanAsString(verbose)+")\n"+
			      
			      "\n\t"+"\tMandatory:"+"\n"+
			      "\t\t"+"-o"+" "+"--outfile"+"\t[outfile]"+"\t\t"+"Specify output file"+"\n"+

			      "\n\tFiltering options:"+"\n"+
			      "\t\t"+"-c" +" "+"--cutoff"+"\t[cutoff]""\t\t"+"Sequence likelihood cutoff (Default: "+stringify(cutoffLikelihood)+")"+"\n"+
			      "\t\t"+"  " +" "+"--cutexp"+"\t[cutoff]""\t\t"+"Average of expectancy of base error cutoff (Default: "+stringify(cutoffAvgExpError)+")"+"\n"+

			      "\t\t"+"" +" "+"--trim"+"\t\t\t\t\t"+"Try to trim from the 3' end (TO IMPLEMENT) (Default: "+stringify(trimSeqs)+")"+"\n"+


			      "\t\t"+"" +""+"--min_length"+"\t[cutoff]"+"\t\t"+"Flag any sequence with strickly less than this min length as failed (Default: "+stringify(minLength)+""+"\n"+
			      "\t\t"+"" +""+"--max_length"+"\t[cutoff]"+"\t\t"+"Flag any sequence with strickly more than this max length as failed (Default: "+stringify(maxLength)+""+"\n"+
			      "\t\t"+"" +""+"--percent"+"\t[percentage]"+"\t\t"+"Flag any sequence the bottom % as failed, to use with only small datasets (Default: not used)"+"\n"+

			      "\n\tComplexity filter options:"+"\n"+
			      
			      "\t\t"+"-e"+","+"--entropy"     +"\t"+"\t\t\t"+"Apply sequence entropy filter (Default: "+booleanAsString(entropy)+")\n"+
			      "\t\t"+"-f"+","+"--frequency"   +"\t"+"\t\t\t"+"Apply base frequency filter  (Default: "+booleanAsString(frequency)+")\n"+
			      "\t\t"+""  +"" +"--comp_cutoff" +"\t"+"\t\t\t"+"Entropy value [0.0-2.0] or fraction [0.0-1.0] of most frequent base accepted (Default: "+stringify(compOrEntCutoff)+")\n"

			      //not used anymore
			      // "\t\t"+""+ "--clip"   +"\t"+"\t\t\t\t"+    "Clip number of bases (default 0, negative values = number of bases to trim from read end).",type='int',default=0)
			      // "\t\t"+""+ "--position"  +"\t"+"\t\t\t\t"+ "First base position to remove (1-based)",type='int',default=0)
			      // // "\n\tOutput options:"+"\n"+
			      // ""

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


	if(  (strcmp(argv[i],"--percent") == 0)   ){
	    bottomPercent =destringify<double>(argv[i+1]);
	    usePercent=true;
	    i++;
	    continue;
	}

	if( (strcmp(argv[i],"-c") == 0) || (strcmp(argv[i],"--cutoff") == 0)   ){
	    cutoffLikelihood =destringify<double>(argv[i+1]);
	    definedCutoff=true;
	    i++;
	    continue;
	}

	if(  (strcmp(argv[i],"--cutexp") == 0)   ){
	    cutoffAvgExpError =destringify<double>(argv[i+1]);
	    definedExpCutoff=true;
	    i++;
	    continue;
	}


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

	if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	    outfile=string(argv[i+1]);
	    i++;
	    continue;
	}


	if( (strcmp(argv[i],"-v") == 0) || (strcmp(argv[i],"--verbose") == 0)   ){
	    verbose=true;
	    continue;
	}


	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;             
    }

    if(usePercent && definedCutoff){
	cerr<<"Cannot defined a likelihood cutoff and a percentage"<<endl;
	return 1;             
    }


    // if(definedExpCutoff && (usePercent || definedCutoff) ){
    // 	cerr<<"Cannot defined a "<<endl;
    // 	return 1;             
    // }

    if(usePercent &&
       (bottomPercent<0 || bottomPercent>1 ) ){
	cerr<<"Cutoffs percentage must be between 0 and 1 for --percent"<<endl;
	return 1;    
    }

    if(outfile.size() == 0){
	cerr<<"The field -o is mandatory exiting"<<endl;
	return 1;             
    }

    if(entropy && frequency){
	cerr<<"Specify only one type of complexity filter"<<endl;
	return 1;             
    }

    if(frequency &&
       (compOrEntCutoff<0 || compOrEntCutoff>1 ) ){
	cerr<<"Cutoffs for frequency between 0 and 1 for --frequency"<<endl;
	return 1;             
    }

    if(entropy &&
       (compOrEntCutoff<0 || compOrEntCutoff>2 ) ){
	cerr<<"Cutoffs for entropy between 0 and 2 for --entropy"<<endl;
	return 1;             
    }

    if( (cutoffLikelihood<0 || cutoffLikelihood>1 ) ){
	cerr<<"The cutoff for likelihood must be between 0 and 1 for"<<endl;
	return 1;             
    }
       

    initializeLikelihoodScores();
    string bamfiletopen = string(argv[argc-1]);

    BamReader reader;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    if(usePercent){
	BamAlignment al2;
	vector<double> likelihoodsFound;
	while ( reader.GetNextAlignment(al2) ) {
	    likelihoodsFound.push_back(compLikelihoodSeq(&al2));
	}
	//	cout<<likelihoodsFound.size()<<endl;

	sort (likelihoodsFound.begin(), likelihoodsFound.end());

	cutoffLikelihood=likelihoodsFound[ int(bottomPercent*(likelihoodsFound.size()))  ];

	if ( !reader.Open(bamfiletopen) ) {
	    cerr << "Could not open input BAM files." << endl;
	    return 1;
	}

	//return 1;
    }

    setVariables(minLength,maxLength,cutoffLikelihood,frequency,entropy,compOrEntCutoff,
		 likelihoodFlag,&likelihoodOS, entropyOSFlag, &entropyOS, frequencyOSFlag, &frequencyOS,verbose);

    // cout<<"h"<<endl;
    // cout<<"outfile "<<outfile<<endl;
    // cout<<"#"<<myHeader.ToString()<<endl;    
    // SamProgram sp;
    // sp.ID          = "filterReads";   
    // sp.Name        = "filterReads";   
    // sp.CommandLine = "";
    // cout<<"b"<<endl;
    // for(int i=0;i<(argc-1);i++){
    // 	sp.CommandLine += (string(argv[i])+" ");
    // }    
    // //need to add version
    // //sp.Version        ;
    // myHeader.Programs.Add(sp);
    // cout<<"#"<<endl;
    // cout<<"#"<<myHeader.ToString()<<endl;
    SamHeader  myHeader=reader.GetHeader();

    BamWriter writer;
    //SamHeader  myHeader=;
    // cout<<"outfile "<<outfile<<endl;
    // cout<<"#"<<reader.GetHeader().ToString()<<"#"<<endl;
    // cout<<reader.GetReferenceData().size()<<endl;
    string pID          = "filterReads";   
    string pName        = "filterReads";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine);


    if(produceUnCompressedBAM)  
	writer.SetCompressionMode(BamWriter::Uncompressed); 

    if( !writer.Open(outfile,
		     myHeader,
		     reader.GetReferenceData() 
		     ) ) {
    	cerr << "Could not open output BAM file  "<<outfile << endl;
    	return 1;	
    }

    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	filterBAMAlign(&al);
	writer.SaveAlignment(al);
    }

    reader.Close();
    writer.Close();

    if(likelihoodFlag)
	likelihoodOS.close();


    return 0;
}

