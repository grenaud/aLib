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

#include "PutProgramInHeader.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

     if( (argc!= 3) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cerr<<"Usage:splitByRG [in bam] [out prefix]"<<endl<<"this program creates one bam file per RG in the with the outprefix\nFor example splitByRG in.bam out will create\nout.rg1.bam\nout.rg2.bam\n"<<endl;
    	return 1;
    }


     string bamfiletopen = string(argv[1]);
     // if(!strEndsWith(bamfiletopen,".bam")){

     // }
     string bamDirOutPrefix    = string(argv[2]);
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

    SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    vector<RefData>  refData=reader.GetReferenceData();
    string pID          = "splitByRG";   
    string pName        = "splitByRG";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
        pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));


    SamReadGroupDictionary 	srgd=header.ReadGroups;
    for(SamReadGroupConstIterator srgci=srgd.ConstBegin();
	srgci<srgd.ConstEnd();
	srgci++){
	//cout<<*srgci<<endl;
	const SamReadGroup rg = (*srgci);
	//cout<<rg.ID<<endl;
	rg2BamWriter[rg.ID] = new  BamWriter();
	rg2BamWriter[rg.ID]->Open(bamDirOutPrefix+"."+rg.ID+".bam",header,references); 
    }



    BamAlignment al;
    unsigned int total=0;
    while ( reader.GetNextAlignment(al) ) {

	// al.SetIsFailedQC(false);
	// writer.SaveAlignment(al);
	// if(al.IsMapped () ){
	//     if(rg2BamWriter.find(refData[al.RefID].RefName) == rg2BamWriter.end()){ //new
	// 	rg2BamWriter[refData[al.RefID].RefName] = new  BamWriter();
	// 	if ( !rg2BamWriter[refData[al.RefID].RefName]->Open(bamDirOutPrefix+"."+refData[al.RefID].RefName+".bam",header,references) ) {
	// 	    cerr     << "Could not open output BAM file "<< bamDirOutPrefix<<"."<<refData[al.RefID].RefName<<".bam" << endl;
	// 	    return 1;
	// 	}
	
	//     }else{
	// 	rg2BamWriter[refData[al.RefID].RefName]->SaveAlignment(al);
	//     }
	// }else{
	//     unmapped.SaveAlignment(al);
	// }
	if(al.HasTag("RG")){
	    string rgTag;
	    al.GetTag("RG",rgTag);
	    //cout<<rgTag<<endl;
	    if(rg2BamWriter.find(rgTag) == rg2BamWriter.end()){ //new
		cerr<<"Found new RG "<<rgTag<<endl;
		rg2BamWriter[rgTag] = new  BamWriter();
	 	if ( !rg2BamWriter[rgTag]->Open(bamDirOutPrefix+"."+rgTag+".bam",header,references) ) {
	 	    cerr     << "Could not open output BAM file "<< bamDirOutPrefix<<"."<<rgTag<<".bam" << endl;
	 	    return 1;
	 	}
		rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	    }else{
		rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	    }
	}else{
	    string rgTag="unknown";	    
	    //cout<<rgTag<<endl;
	    if(rg2BamWriter.find(rgTag) == rg2BamWriter.end()){ //new
		cerr<<"Found new RG "<<rgTag<<endl;
		rg2BamWriter[rgTag] = new  BamWriter();
	 	if ( !rg2BamWriter[rgTag]->Open(bamDirOutPrefix+"."+rgTag+".bam",header,references) ) {
	 	    cerr     << "Could not open output BAM file "<< bamDirOutPrefix<<"."<<rgTag<<".bam" << endl;
	 	    return 1;
	 	}
		rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	    }else{
		rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	    }

	    // cerr << "Cannot get RG tag for " << al.Name<<endl;
	    // return 1;
	}

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

