#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <sstream>

using namespace std;

typedef struct { 
    char baseC;
    int  qualC;
} base;



inline base bin2base(unsigned char c ){
    base toreturn;
    toreturn.baseC="ACGT"[c % 4];
    c /= 4;
    toreturn.qualC=(c % 64);
    return toreturn;
}





int main (int argc, char *argv[]) {
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<"estimateControlReadsBam [FIRST CYCLE OF INDEX] [CONTROL INDEX] [LANE#] [TILE#] [DIRECTORY WITH BASECALLS]"<<endl;
	cout<<endl;
	cout<<"estimateControlReadsBam 77 TTGCCGC 1 20 /mnt/solexa/120217_SOLEXA-GA02_00048_PEdi_AW_CB_MM/Data/Intensities/BaseCalls/"<<endl;
	return 1;
    }

    int firstCycleIDX       = atoi(  argv[1]);
    string ctrlSeqINDEX     = string(argv[2]);
    int ctrlSeqINDEXlength  =ctrlSeqINDEX.length();
    string lane                = string(argv[3]);
    string tile                = string(argv[4]);
    string directory        = string(argv[5]);
    stringstream ss;
    fstream mybclfile [firstCycleIDX+ctrlSeqINDEXlength-1];
    unsigned int numberOfClustersFirst=0 ;
    unsigned int numberOfCtrlClusters =0 ;


    //open files
    for(int cycle=firstCycleIDX;cycle<=(firstCycleIDX+ctrlSeqINDEXlength-1);cycle++){
	ss.str(std::string());
	ss<<cycle;
	string bclFile=directory+"/L00"+lane+"/C"+ ss.str()+".1/s_"+lane+"_"+tile+".bcl";

	mybclfile[cycle-firstCycleIDX].open(bclFile.c_str(),ios::in|ios::binary);
	if (!mybclfile[cycle-firstCycleIDX]) {
	    cout<<"Unable to read file "<<bclFile<<endl;
	}
	unsigned int numberOfClusters ;
	mybclfile[cycle-firstCycleIDX].read((char*)&numberOfClusters, sizeof (int));

	if(numberOfClustersFirst == 0 ){
	    numberOfClustersFirst=numberOfClusters;
	}else{
	    if(numberOfClustersFirst!=numberOfClusters){
		cerr<<"Number of clusters in "<<bclFile<<" is different from the other cycles, exiting\n";
		return 1;
	    }
	}
    }


    for(unsigned int i=0;i<numberOfClustersFirst;i++){
	string sequence="";
	string quality ="";

	for(int idxcycle=0;idxcycle<ctrlSeqINDEXlength;idxcycle++){
	    char toread;
	    //cout<<idxcycle<<endl;
	    mybclfile[idxcycle].read(&toread, sizeof (char));
	    base returned=bin2base(toread);
	    sequence   +=returned.baseC;
	    quality    +=char(33+returned.qualC);
	}
	// cout<<sequence<<endl;
	// cout<<quality<<endl;
	// cout<<i<<endl;
	// cout<<ctrlSeqINDEX<<endl;
	if(sequence == ctrlSeqINDEX){
	    numberOfCtrlClusters++;
	}
    }

    cout<<lane<<"\t"<<tile<<"\t"<<numberOfCtrlClusters<<"\t"<<numberOfClustersFirst<<"\t"<<100.0*(double(numberOfCtrlClusters)/double(numberOfClustersFirst))<<"%"<<endl;
    /*
    //    string filebcl=argv[1];
    fstream mybclfile (argv[1],ios::in|ios::binary);
    //ifstream myFile ("data.bin", ios::in | ios::binary);

    // cout<<argv[1]<<endl;
    if (!mybclfile) {
	cout<<"Unable to read file "<<argv[1]<<endl;
    }
    unsigned int numberOfClusters ;
    mybclfile.read((char*)&numberOfClusters, sizeof (int));
    cout<<"numberOfClusters "<<numberOfClusters<<endl;

    char toread;
    int counter=0;
    while(mybclfile.read(&toread, sizeof (char))){
	base returned=bin2base(toread);
	cout<<"base "<< returned.baseC<<endl;
	cout<<"qual "<< returned.qualC<<endl;
	counter++;
    }
    // cout<<"counter "<<counter<<endl;
    */
    return 0;
}

