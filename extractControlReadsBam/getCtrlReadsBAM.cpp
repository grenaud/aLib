#include <iostream>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

using namespace std;
using namespace BamTools;

int main (int argc, char *argv[]) {
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"This program reads a series of BAM files (or a single one)"<<endl;
	cout<<"and produces an BAM file as output with only the reads "<<endl;
	cout<<"whose XI field is the sequence specified in [CONTROL INDEX]"<<endl;
	cout<<"This program is used normally to extract phiX reads for QC"<<endl;
	cout<<"Usage:"<<endl;
	cout<<""<<endl;
	//	cout<<"getendposition [READ GROUP] [CONTROL INDEX] s_sequence.bam out.bam"<<endl;
	cout<<argv[0]<<" [CONTROL INDEX] out.bam in1.bam in2.bam ..."<<endl;
	return 1;
    }

    //    string readgroup    = string(argv[1]);
    string ctrlindex      = string(argv[1]);
    string outputFilename = string(argv[2]);
    BamWriter writer;

    for(int argument=3;argument<argc;argument++){
	string bamfiletopen = string(argv[argument]);
	BamReader reader;
	if ( !reader.Open(bamfiletopen) ) {
	    cerr << "Could not open input BAM files." << endl;
	    return 1;
	}

	// if ( !reader.LocateIndex() ){
	// 	cerr << "warning: cannot locate index for file " << bamfiletopen<<endl;
	// 	//return 1;
	// }


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

	    if(//rgTag      == readgroup &&
	       firstIndex == ctrlindex ){
		writer.SaveAlignment(al);
		// cout<<rgTag<<endl;
		// cout<<firstIndex<<endl;
	   
	    }
	}
	reader.Close();
    }
    writer.Close();
   
    return 0;
}

