/*
 * FilterBAMal
 * Date: Nov-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */



#include "FilterBAMal.h"


int minLength=-1;
int maxLength=-1;
double likeQual[64];
double expError[64];

bool frequency;
bool entropy;
bool compOrEntCutoff;
double cutoffLikelihood = 0.5;

ofstream * likelihoodOS;
bool     likelihoodFlag  = false;
ofstream * entropyOS;
bool     entropyOSFlag   = false;
ofstream * frequencyOS;
bool     frequencyOSFlag   = false;
bool     verbose=false;



void initializeLikelihoodScores(){
    for(int i=0;i<64;i++){
	if(i == 0)
	    likeQual[i]    = -1.0; // need to fix this to minimize precision loss
	else
	    likeQual[i]    = log10(1.0- pow(10.0,i/-10.0) ); // need to fix this to minimize precision loss	
	//cout<<"likeQual["<<i<<"] "<<likeQual[i]<<endl;
	expError[i]    = pow(10.0,i/-10.0); 
	// cout<<"expError["<<i<<"] "<<expError[i]<<endl;
    }
}

void setVariables(int _minLength,int _maxLength,double _cutoffLikelihood,bool _frequency,bool _entropy,bool _compOrEntCutoff, 
		  bool _likelihoodFlag,ofstream * _likelihoodOS,bool _entropyOSFlag,ofstream  * _entropyOS,bool _frequencyOSFlag,ofstream  * _frequencyOS,bool _verbose){
   
    
    minLength        = _minLength;
    maxLength        = _maxLength;    

    frequency        = _frequency;
    entropy          = _entropy;
    compOrEntCutoff  = _compOrEntCutoff;
    likelihoodFlag   = _likelihoodFlag;
    cutoffLikelihood = _cutoffLikelihood;
    likelihoodOS     = _likelihoodOS;

    entropyOS        = _entropyOS;
    entropyOSFlag    = _entropyOSFlag;
    frequencyOS      = _frequencyOS;
    frequencyOSFlag  = _frequencyOSFlag;
    verbose          = _verbose;

}


/* This subroutine finds the most common nucleotide and checks if
 * the fraction of the count of this nucleotide and the total
 * resolved nucleotide is less than a certain cutoff otherwise it flags the BamAlignment
 * as failed.
 */
void complexFrequency(BamAlignment * al){
    int countBp [4];
    int total=0;
    int maxBp=-1;
    char c;

    for(int j=0;j<4;j++)
	countBp[j]=0;

    for(int i=0;i<int(al->QueryBases.size());i++){
	c=al->QueryBases[i];
	if(c == 'N'){ continue;     } //do not care about 'N'
	if(c == 'A'){ countBp[0]++; }
	if(c == 'C'){ countBp[1]++; }
	if(c == 'G'){ countBp[2]++; }
	if(c == 'T'){ countBp[3]++; }
	total++;
    }


    for(int j=0;j<4;j++)
	if(countBp[j] > maxBp)
	    maxBp =countBp[j];

    if(total == 0) {//just Ns
	al->SetIsFailedQC(true);
    }else{
	
	double finalFreq=(double(maxBp)/double(total));
	if(frequencyOSFlag)
	    frequencyOS->write( (char *)&finalFreq, sizeof(finalFreq));

	if(finalFreq  > compOrEntCutoff)
	    if(frequency)
		al->SetIsFailedQC(true);
    }
}



/* This subroutine computes the entropy as defined by:
 *  H = sum p_i * log(p_i) for i each individual character
 *  it flags the BamAlignment if it below a certain cutoff
 * 
 */
void complexEntropy(BamAlignment * al){
    int countBp [4];
    int total=0;
    char c;

    for(int j=0;j<4;j++)
	countBp[j]=0;

    for(int i=0;i<int(al->QueryBases.size());i++){
	c=al->QueryBases[i];
	if(c == 'N'){ continue;     } //do not care about 'N'
	if(c == 'A'){ countBp[0]++; }
	if(c == 'C'){ countBp[1]++; }
	if(c == 'G'){ countBp[2]++; }
	if(c == 'T'){ countBp[3]++; }
	total++;
    }


    if(total == 0) {//just Ns
	al->SetIsFailedQC(true);
    }else{
	double entropyCalc=0.0;
	for(int j=0;j<4;j++){
	    if(countBp[j] != 0 ){
		entropy +=  -1.0*( (double(countBp[j])/double(total)) * (log( (double(countBp[j])/double(total)) )/log(2.0)) ) ;
	    }
	}

	if(entropyOSFlag)
	    entropyOS->write( (char *)&entropyCalc, sizeof(entropyCalc) );

	if(entropyCalc < compOrEntCutoff)
	    if(entropy)
		al->SetIsFailedQC(true);	
    }

}


double compLikelihoodSeq(BamAlignment * al){
    int qualOffset=33;

    vector<int> quals1;
    for(int i=0;i<int(al->Qualities.length());i++){
	if(al->QueryBases[i] != 'N')
	    quals1.push_back( max( (int(char( al->Qualities[i] ))-qualOffset),2 ))  ; //since qual scores less than 2 do not make sense
    }

    double totalLike=0;
    for(int i=0;i<int(quals1.size());i++){
	totalLike += likeQual[ quals1[i] ];
    }
    double likeSeq= pow(10.0,totalLike);
    return likeSeq;
}


double computeExpectationError(BamAlignment * al){
    int qualOffset=33;

    vector<int> quals1;
    for(int i=0;i<int(al->Qualities.length());i++){
	if(al->QueryBases[i] != 'N')
	    quals1.push_back( max( (int(char( al->Qualities[i] ))-qualOffset),2 ))  ; //since qual scores less than 2 do not make sense
    }

    double totalExp=0;
    for(int i=0;i<int(quals1.size());i++){
	totalExp += expError[ quals1[i] ];
    }
    // cout<<totalExp<<endl;
    double likeSeq= totalExp; ///double( al->QueryBases.length() );
    return likeSeq;
}


void filterBAMAlign(BamAlignment * al){

    // if(al->IsMapped()) {
    // 	cerr<<"Read cannot be mapped"<<endl;
    // 	exit(1);
    // }
    //reset the flag
    al->SetIsFailedQC(false);

    if( (minLength != -1) &&
	(int(al->QueryBases.size()) < minLength)
	){
	al->SetIsFailedQC(true);
    }

    if( (maxLength != -1) &&
	(int(al->QueryBases.size()) > maxLength)
	){
	al->SetIsFailedQC(true);
    }

    // int qualOffset=33;

    // vector<int> quals1;
    // //cout<<al->Qualities<<endl;
    // for(int i=0;i<int(al->Qualities.length());i++){
    // 	//cout<<(int(char( al->Qualities[i] ))-qualOffset)<<endl;
    // 	quals1.push_back( max( (int(char( al->Qualities[i] ))-qualOffset),2 ))  ; //since qual scores less than 2 do not make sense
    // }

    // double totalLike=0;
    // for(int i=0;i<int(quals1.size());i++){
    // 	totalLike += likeQual[ quals1[i] ];
    // }

    double likeSeq=compLikelihoodSeq(al);

    
    //cout<<pow(10.0,totalLike)<<endl;
    //if(pow(10.0,totalLike) <0.2)
    //cout<<pow(10.0,totalLike)<<"\t"<<al->Qualities<<endl;
    if(verbose){             
	cerr<<al->Name<<"\tlike:\t"<<likeSeq<<"\tcutoff:\t"<<cutoffLikelihood<<"\tfailed:"<<boolStringify(likeSeq < cutoffLikelihood)<<endl;
    }


    double avgExp=computeExpectationError(al);
    // if(!al->IsPaired() ){
    // 	cout<<al->QueryBases.size()<<"\t"<<avgExp<<"\t"<<log10(avgExp)<<endl;
    // }


    if(avgExp > 1.0)
	al->SetIsFailedQC(true);
    
    //to uncomment
    // if(likeSeq < cutoffLikelihood)
    // 	al->SetIsFailedQC(true);

    if(likelihoodFlag){
	likelihoodOS->write( (char *)&likeSeq, sizeof(likeSeq));
    }

    if(entropy || entropyOSFlag  )
	complexEntropy(al);

    if(frequency || frequencyOSFlag)
	complexFrequency(al);

}
