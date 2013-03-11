#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream> 

using namespace std;

int main (int argc, char *argv[]) {
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<"./getCtrlReadsFASTQ [FIRST CYCLE FIRST INDEX] [CONTROL INDEX] [FIRST CYCLE SECOND INDEX] [OUT PREFIX] s_sequence.txt"<<endl;
	cout<<"ex:"<<endl;
	cout<<"./getCtrlReadsFASTQ 76 TTGCCGC 159 s_sequence.txt"<<endl;
	return 1;
    }

    int    firstCycleIndex  = atoi(argv[1]);
    string ctrlindex        = string(argv[2]);
    int    ctrlindexLen     = ctrlindex.length();
    int    secondCycleIndex = atoi(argv[3]);
    string outprefix        = string(argv[4]);
    string fastqFileToOpen  = string(argv[5]);
    ifstream fastqin;
    string line;
    string id;
    string sequence;
    bool readIsCtrl=false;
    long lineCounter=0;
    ofstream freads;
    ofstream rreads;

    freads.open(string(outprefix+"r1.txt").c_str());
    rreads.open(string(outprefix+"r2.txt").c_str());
    if (!freads.is_open()){
	cerr << "Problem writing forward reads for "<<outprefix<<endl;
	return 1;       
    }
    if (!rreads.is_open()){
	cerr << "Problem writing forward reads for "<<outprefix<<endl;
	return 1;       
    }
    // cout<<"ok2"<<endl;

    // cout<<"firstCycleIndex "<< firstCycleIndex <<endl;
    // cout<<"ctrlindex       "<< ctrlindex <<endl;
    // cout<<"ctrlindexLen    "<< ctrlindexLen <<endl;
    // cout<<"secondCycleIndex "<< secondCycleIndex <<endl;
    // cout<<"fastqFileToOpen "<< fastqFileToOpen <<endl;

    fastqin.open(fastqFileToOpen.c_str(), ios::in);    // open the streams
    if (fastqin.is_open()) {
	while ( getline (fastqin,line) ){

	    switch(lineCounter){
	    case 0:
		if(line[0] != '@'){
		    cerr << "Wrong defline "<<line<<endl;
		    return 1;       
		}
		id=line;
		break;
	    case 1:
		sequence=line;
		break;
	    case 3:
		if(sequence.substr(firstCycleIndex-1,ctrlindexLen) == ctrlindex){		    
		    // cout<<sequence<<endl<<line<<endl;
		    // cout<<"F"<<endl;
		    freads<<id<<endl
			  <<sequence.substr(0,firstCycleIndex-1)<<endl
			  <<"+"<<endl
			  <<line.substr(0,firstCycleIndex-1)<<endl;
		    //cout<<"R"<<endl;
		    rreads<<id<<endl
			  <<sequence.substr(firstCycleIndex+ctrlindexLen-1,secondCycleIndex-firstCycleIndex-ctrlindexLen)<<endl
			  <<"+"<<endl
			  <<line.substr(firstCycleIndex+ctrlindexLen-1,secondCycleIndex-firstCycleIndex-ctrlindexLen)<<endl;



		}
		lineCounter=-1;
		break;
	    default:
		break;
	    }		    
	    lineCounter++;    
	    
  
	}

	fastqin.close();
	freads.close();
	rreads.close();

    }else {
	cerr << "Unable to open file "<<fastqFileToOpen<<endl;
	return 1;       
    }
   
    return 0;
}

