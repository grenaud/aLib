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

     if( (argc!= 2) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cerr<<"Usage:insertSize [in bam]"<<endl<<"\n"<<endl;
    	return 1;
    }


     string bamfiletopen = string(argv[1]);
     // if(!strEndsWith(bamfiletopen,".bam")){

     // }
     //string bamDirOutPrefix    = string(argv[2]);
     map< string, map<int,unsigned int> > rg2Counter;
     

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
	const SamReadGroup rg = (*srgci);
	map<int,unsigned int> tempMap;
	rg2Counter[rg.ID] = tempMap;
    }




    BamAlignment al;
    unsigned int total=0;
    while ( reader.GetNextAlignment(al) ) {
	int insertSize=-1;
	//TODO is mapped ???
	if(al.IsPaired()   ){
	    if(al.IsMapped() && al.IsProperPair()   ){
		if( al.IsFirstMate()  ){
		    if(al.InsertSize>1)
			insertSize=al.InsertSize;
		}
	    }

	}else{
	    insertSize=al.Length;
	}
	
	if(insertSize==-1)
	    continue;

	if(al.HasTag("RG")){
	    string rgTag;
	    al.GetTag("RG",rgTag);
	    //cout<<rgTag<<endl;
	    if(rg2Counter.find(rgTag) == rg2Counter.end()){ //new
		cerr<<"Found new RG "<<rgTag<<endl;
		map<int,unsigned int> tempMap;
		rg2Counter[rgTag]   = tempMap;
		//rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	    }else{
		//rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	    }
	    rg2Counter[rgTag][insertSize]++;

	}else{
	    // string rgTag="unknown";	    
	    // //cout<<rgTag<<endl;
	    // if(rg2BamWriter.find(rgTag) == rg2BamWriter.end()){ //new
	    // 	cerr<<"Found new RG "<<rgTag<<endl;
	    // 	rg2BamWriter[rgTag] = new  BamWriter();
	    // 	if ( !rg2BamWriter[rgTag]->Open(bamDirOutPrefix+"."+rgTag+".bam",header,references) ) {
	    // 	    cerr     << "Could not open output BAM file "<< bamDirOutPrefix<<"."<<rgTag<<".bam" << endl;
	    // 	    return 1;
	    // 	}
	    // 	rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	    // }else{
	    // 	rg2BamWriter[rgTag]->SaveAlignment(al);	    	   
	    // }
	    // cerr << "Cannot get RG tag for " << al.Name<<endl;
	    // return 1;
	}

	total++;
    } //while al

    reader.Close();
    // writer.Close();
    
    // unmapped.Close();
    //rg2CounterIt
    
    vector<string> keysrg = allKeysMap(rg2Counter);
    for(unsigned int i=0;i<keysrg.size();i++){
	vector<int> insertSizes = allKeysMap( rg2Counter[ keysrg[i] ] );
	sort (insertSizes.begin(), insertSizes.end());
	cout<<"#"<<keysrg[i]<<"\t"<<0<<endl;
	for(unsigned int j=0;j<insertSizes.size();j++){
	    cout<<insertSizes[j]<<"\t"<<rg2Counter[ keysrg[i] ][ insertSizes[j] ]<<endl;
	}
	
    }

    // map<string, map<int,unsigned int> >::iterator rg2CounterIt;
    // for (rg2CounterIt =rg2Counter.begin(); 
    // 	 rg2CounterIt!=rg2Counter.end(); 
    // 	 rg2CounterIt++){
    // 	map<int,unsigned int>::iterator counterIt;
    // 	for (counterIt =rg2Counter.begin(); 
    // 	     counterIt!=rg2Counter.end(); 
    // 	     counterIt++){
	    
    // 	}
    // 	//rg2CounterIt->second->Close();
    // }
    
    cerr<<"Looked at "<<total<<" reads"<<endl;


    return 0;
}

