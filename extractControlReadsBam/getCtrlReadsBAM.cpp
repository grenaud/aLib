#include <iostream>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "utils.h"

using namespace std;
using namespace BamTools;

int main (int argc, char *argv[]) {

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cerr<<"This program reads a series of BAM files (or a single one)"<<endl;
	cerr<<"and produces an BAM file as output with only the reads "<<endl;
	cerr<<"whose XI field is the sequence specified in [CONTROL INDEX]"<<endl;
	cerr<<"This program is used normally to extract phiX reads for QC"<<endl;
	cerr<<"Usage:"<<endl;
	cerr<<""<<endl;
	//	cerr<<"getendposition [READ GROUP] [CONTROL INDEX] s_sequence.bam out.bam"<<endl;
	cerr<<argv[0]<<" [CONTROL INDEX(es)] out.bam in1.bam in2.bam ..."<<endl;
	cerr<<"if multiple indices are used, use commas to separate them"<<endl;
	return 1;
    }

    //    string readgroup    = string(argv[1]);
    string ctrlindex      = string(argv[1]);
    string outputFilename = string(argv[2]);
    BamWriter writer;
    vector<string> ctrlindexVec = allTokens(ctrlindex,',');

    for(int argument=3;argument<argc;argument++){
	string bamfiletopen = string(argv[argument]);
	BamReader reader;
	if ( !reader.Open(bamfiletopen) ) {
	    cerr << "Could not open input BAM files." << endl;
	    return 1;
	}


	vector<RefData>  testRefData=reader.GetReferenceData();
	const SamHeader header = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();

	if(argument==3){ //first file
	    if ( !writer.Open(outputFilename, header, references) ) {
		cerr << "Could not open output BAM file" << endl;
		return 1;
	    }
	}

	BamAlignment al;
	while ( reader.GetNextAlignment(al) ) {
	    string rgTag;
	    string firstIndex;

	    ///al.GetTag("RG",rgTag);
	    if(!al.GetTag("XI",firstIndex)){
		cerr << "Could not get XI tag for read" <<al.Name << " in file " << bamfiletopen <<endl;
		return 1;
	    }

	    bool matchedIndex=false;
	    //probably not optimal if there are a lot of possibilities but for less than 10 we should be ok
	    for(unsigned int i=0;i<ctrlindexVec.size();i++){

		if(firstIndex    == ctrlindexVec[i] ){
		    matchedIndex =  true;
		    break;
		}

	    }

	    if( matchedIndex ){
		writer.SaveAlignment(al);	   
	    }

	}
	reader.Close();
    }
    writer.Close();
   
    return 0;
}

