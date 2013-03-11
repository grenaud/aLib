#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <iomanip>

using namespace std;

//taken from http://www.umanitoba.ca/afs/Plant_Science/psgendb/local/pkg/CASAVA_v1.8.2/src/c++/include/alignment/BclReader.hh
static const int blockSize = 25;
static const int imageWidth = 2048;
static const int blocksPerlLine = (imageWidth + blockSize - 1) / blockSize;
static const int imageHeight = 20000;

int main (int argc, char *argv[]) {

    unsigned long totalClusters=0;
    fstream myclocsfile (argv[1],ios::in|ios::binary);
    if (!myclocsfile) {
	cout<<"Unable to read file "<<argv[1]<<endl;
    }
    char toread;
    myclocsfile.read(&toread,sizeof(char)); //version
    unsigned int numberOfTotalBins ;
    myclocsfile.read((char*)&numberOfTotalBins, sizeof (int)); //number of bins
    unsigned int binIndex=0 ;
    double xOffset = 0;
    double yOffset = 0;


    while(binIndex<numberOfTotalBins){ //for each bin

	unsigned char numberOfrecordsInbin ;
	myclocsfile.read((char*)&numberOfrecordsInbin, sizeof (char));

	totalClusters+=int(numberOfrecordsInbin);
	if(numberOfrecordsInbin == 0){
	    binIndex++;
	    continue;
	}

	xOffset = float(blockSize) * (binIndex % blocksPerlLine);
	yOffset = float(blockSize) * (binIndex / blocksPerlLine);

	for(int i=0;i<int(numberOfrecordsInbin);i++){ //for each block in the bin
	    unsigned char x ;
	    myclocsfile.read((char*)&x, sizeof (char));
	    unsigned char y ;
	    myclocsfile.read((char*)&y, sizeof (char));	    
	}

	binIndex++;
	
    }
    cout<<totalClusters<<endl;

    return 0;
}

