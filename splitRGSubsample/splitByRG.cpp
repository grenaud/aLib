/*
 * splitByRG
 * Date: Feb-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "utils.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

     if( (argc!= 4 && argc !=5 && argc !=6) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cerr<<"Usage:splitByRG [in bam] [rg Tally] [out prefix] (optional target)"<<endl<<"this program will subsample a BAM file per read group for a certain target\nFor example splitByRG in.bam tally.txt out will create\nout.rg1.bam\nout.rg2.bam\n"<<endl;
    	return 1;
    }


     string bamfiletopen      = string(argv[1]);
     string rgTally           = string(argv[2]);
     string bamDirOutPrefix   = string(argv[3]);
     
     int target            =  200000;
     int maxTarget         = 1000000;

     if(argc==5){
	 target    = destringify<int> ( string(argv[4]) );	 
     }

     if(argc==6){
	 target    = destringify<int> ( string(argv[4]) );	 
	 maxTarget = destringify<int> ( string(argv[5]) );	 
     }


     cerr<<"minimum fragments:\t"<<target<<endl;
     cerr<<"target  fragments:\t"<<maxTarget<<endl;

     string line;
     ifstream myFileTally;
     map<string,double> rg2Fraction;

     myFileTally.open(rgTally.c_str(), ios::in);
     cerr<<"Retained groups:\n"<<endl;
     cerr<<"RG\t#mapped\tfraction retained"<<endl;
     cerr<<"-----------------------------------"<<endl;

     if (myFileTally.is_open()){
	 while ( getline (myFileTally,line)){
	     vector<string> tokens = allTokens(line,'\t');
	     if(tokens.size() > 6)
		 if( tokens[1] == "pass" && 
		    (tokens[0] != "\"\""    && 
		     tokens[0] != "control" && 
		     tokens[0] != "TOTAL") ){
		     //cout<<tokens[0]<<"\t"<<tokens[5]<<endl;
		     int count = destringify<int>(tokens[5]);

		     if(count>target){

			 if(count>=maxTarget){
			     rg2Fraction[  tokens[0] ] = double(maxTarget)/double(count);
			     cout<<tokens[0]<<"\t"<<count<<"\t"<<double(maxTarget)/double(count)<<endl;
			 }else{
			     cout<<tokens[0]<<"\t"<<count<<"\t"<<1.0<<endl;
			     rg2Fraction[  tokens[0] ] = 1.0;
			 }
		     }
		 }
	 }
	 myFileTally.close();
     }else{
	 cerr << "Unable to open file "<<rgTally<<endl;
	 return 1;
     }



     map<string,BamWriter *> rg2BamWriter;
     
     // if(!isDirectory(bamDirOut)){
     // 	 cerr<<"ERROR: the out directory does not exist"<<endl;
     // 	return 1;
     // }

     BamReader reader;

     if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
     }
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    vector<RefData>  refData=reader.GetReferenceData();

    SamReadGroupDictionary 	srgd=header.ReadGroups;
    for(SamReadGroupConstIterator srgci=srgd.ConstBegin();
	srgci<srgd.ConstEnd();
	srgci++){
	//cout<<*srgci<<endl;
	const SamReadGroup rg = (*srgci);
	//cout<<rg.ID<<endl;
	if( rg2Fraction.find(rg.ID) != rg2Fraction.end() ){
	    rg2BamWriter[rg.ID] = new  BamWriter();
	    rg2BamWriter[rg.ID]->Open(bamDirOutPrefix+"."+rg.ID+".bam",header,references); 
	}
	//cout<<bamDirOutPrefix+"."+rg.ID+".bam"<<endl;
    }
    //    return 1;

    //    BamWriter unmapped;

    // cout<<header.ToString()<<endl;
    // return 1;

    // if ( !unmapped.Open(bamDirOutPrefix+".unmapped.bam",header,references) ) {
    // 	cerr << "Could not open output BAM file "<< bamDirOutPrefix+".unmapped.bam" << endl;
    // 	return 1;
    // }

    //    cout<<"reading"<<endl;

    BamAlignment al;
    unsigned int total=0;
    while ( reader.GetNextAlignment(al) ) {


	if(al.HasTag("RG") &&
	   al.IsMapped() ){
	    string rgTag;
	    al.GetTag("RG",rgTag);
	    //cout<<rgTag<<endl;
	    if(rg2BamWriter.find(rgTag) == rg2BamWriter.end()){ //new: ignore completely
	
		
	    }else{
		if( randomProb() <= rg2Fraction[  rgTag ] ){
		    rg2BamWriter[rgTag]->SaveAlignment(al);	 
		    //cout<<"wrote "<<rgTag<<endl;
		}   else{
		    //cout<<"skipped "<<rgTag<<endl;
		}	   
	    }
	}// else{
	//     string rgTag="unknown";	    
	//     //cout<<rgTag<<endl;
	//     if(rg2BamWriter.find(rgTag) == rg2BamWriter.end()){ //new
	// 	cerr<<"Found new RG "<<rgTag<<endl;
	// 	rg2BamWriter[rgTag] = new  BamWriter();
	//  	if ( !rg2BamWriter[rgTag]->Open(bamDirOutPrefix+"."+rgTag+".bam",header,references) ) {
	//  	    cerr     << "Could not open output BAM file "<< bamDirOutPrefix<<"."<<rgTag<<".bam" << endl;
	//  	    return 1;
	//  	}
	// 	rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	//     }else{
	// 	rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	//     }

	//     // cerr << "Cannot get RG tag for " << al.Name<<endl;
	//     // return 1;
	// }

	total++;
    } //while al

    reader.Close();
    // writer.Close();
    
    // unmapped.Close();

    map<string,BamWriter *>::iterator rg2BamWriterIt;
    for (rg2BamWriterIt =rg2BamWriter.begin(); 
	 rg2BamWriterIt!=rg2BamWriter.end(); 
	 rg2BamWriterIt++){
	rg2BamWriterIt->second->Close();
    }
    cerr<<"Wrote succesfully "<<total<<" reads"<<endl;


    return 0;
}

