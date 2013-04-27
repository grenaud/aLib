/*
 * countClustersBAM
 * Date: Aug-02-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

using namespace std;
using namespace BamTools;

int main (int argc, char *argv[]) {
    bool ignorePCRDup=false;
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<"countClustersBAM [options] file1.bam [file2.bam file3.bam ..]\nOptions:\n\t-i ignore those reads marked as PCR duplicates"<<endl;
	return 1;
    }

    for(int argument=1;argument<argc;argument++){
	if(string(argv[argument]) == "-i"){
	    ignorePCRDup=true;
	}
    }

    unsigned long counter=0;

    for(int argument=1;argument<argc;argument++){
	string bamfiletopen = string(argv[argument]);

	BamReader reader;

	if ( !reader.Open(bamfiletopen) ) {
	    cerr << "Could not open input BAM files." << endl;
	    return 1;
	}



	BamAlignment al;
	//    while ( reader.GetNextAlignment(al) ) {
	while ( reader.GetNextAlignmentCore(al) ) {
	    //skip PCR duplicates if the flag was set
	    if(ignorePCRDup && al.IsDuplicate())
		continue;

	    if(al.IsPaired()){
		if(al.IsFirstMate()){
		    counter++;
		}	   
	    }else{
		counter++;
	    }
	}
	reader.Close();
    }

    cout<<counter<<endl;

    return 0;
}

