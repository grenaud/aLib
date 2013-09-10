/*
 * FilterBAMal
 * Date: Nov-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */



#include "FilterBAMal.h"


//! Initialize likelihood values
/*!

  Initialize likelihood values, called by constructor

*/
void FilterBAMal::initializeLikelihoodScores(){
    for(int i=0;i<64;i++){
	if(i == 0)
	    likeQual[i]    = -1.0; // need to fix this to minimize precision loss
	else
	    likeQual[i]    = log10(1.0- pow(10.0,i/-10.0) ); // need to fix this to minimize precision loss	
	//cout<<"likeQual["<<i<<"] "<<likeQual[i]<<endl;
	expError[i]    = pow(10.0,i/-10.0); 
	//cout<<"expError["<<i<<"] "<<expError[i]<<endl;
    }
}



//! Constructor
/*!

  It calls initializeLikelihoodScores()  to initialize likelihood values and copies values cutoffs
  \param _minLength Minimum sequence length required not to flag as fail
  \param _maxLength Maximum sequence length required not to flag as fail
  \param _cutoffAvgExpError Threshold for the expected number of mismatches
  \param _frequency If we use sequence frequency/complexity as cutoff
  \param _entropy   If we use sequence entropy as cutoff
  \param _compOrEntCutoff Cutoff for either frequency/complexity or entropy
  \param _likelihoodFlag Whether we write the expected number of mismatches to a file
  \param _likelihoodOS ofstream pointer to file where the expected number of mismatches will be written
  \param _entropyOSFlag Whether we write the entropy to a file
  \param _entropyOS ofstream pointer to file where the entropy will be written
  \param _frequencyOSFlag Whether we write the frequency/complexity to a file
  \param _frequencyOS ofstream pointer to file where the frequency/complexity will be written
  \param _verbose Whether we will write a message for each read
  \param _resetQC If set to true, we will flag every sequence as not failed from the beginning
*/
FilterBAMal::FilterBAMal(int _minLength,int _maxLength,double _cutoffAvgExpError,bool _frequency,bool _entropy,bool _compOrEntCutoff, bool _likelihoodFlag,ofstream * _likelihoodOS,bool _entropyOSFlag,ofstream  * _entropyOS,bool _frequencyOSFlag,ofstream  * _frequencyOS,bool _verbose,bool _resetQC):

    minLength        ( _minLength),
    maxLength        ( _maxLength),    

    frequency        ( _frequency),
    entropy          ( _entropy),
    compOrEntCutoff  ( _compOrEntCutoff),
    likelihoodFlag   ( _likelihoodFlag),
    //cutoffLikelihood ( _cutoffLikelihood),
    cutoffAvgExpError( _cutoffAvgExpError),

    likelihoodOS     ( _likelihoodOS),

    entropyOS        ( _entropyOS),
    entropyOSFlag    ( _entropyOSFlag),
    frequencyOS      ( _frequencyOS),
    frequencyOSFlag  ( _frequencyOSFlag),
    verbose          ( _verbose),

    resetQC          ( _resetQC)
    //    repToPrint       ( _repToPrint)

{
   
    initializeLikelihoodScores();

    repToPrint.totalSeq  = 0;
    repToPrint.qcLength  = 0;

    repToPrint.qcFailClx = 0;
    repToPrint.qcFailEnt = 0;
    repToPrint.qcFailExp = 0;
    repToPrint.qcBefore  = 0;
}



//! Destructor
/*!

*/
FilterBAMal::~FilterBAMal(){

}



//! Destructor
/* This subroutine finds the most common nucleotide and checks if
   the fraction of the count of this nucleotide and the total
   resolved nucleotide is less than a certain cutoff otherwise it flags the BamAlignment
   as failed.
   \param al input pointer to BamAlignment object that needs to be written
 */
void FilterBAMal::complexFrequency(BamAlignment * al){
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
	repToPrint.qcFailClx++;
    }else{
	
	double finalFreq=(double(maxBp)/double(total));
	if(frequencyOSFlag)
	    frequencyOS->write( (char *)&finalFreq, sizeof(finalFreq));

	if(finalFreq  > compOrEntCutoff)
	    if(frequency){
		al->SetIsFailedQC(true);
		repToPrint.qcFailClx++;
	    }
    }
}




//! To compute sequence entropy
/* 
   This subroutine computes the entropy as defined by:
   H = sum p_i * log(p_i) for i each individual character
   it flags the BamAlignment if it below a certain cutoff.
   
   \param al input pointer to BamAlignment object that needs to be written
 */
void FilterBAMal::complexEntropy(BamAlignment * al){
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
	repToPrint.qcFailEnt++;
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
	    if(entropy){
		al->SetIsFailedQC(true);	
		repToPrint.qcFailEnt++;
	    }
    }

}



// //! To compute sequence entropy
// /* 
//    This subroutine computes the entropy as defined by:
//    H = sum p_i * log(p_i) for i each individual character
//    it flags the BamAlignment if it below a certain cutoff.
   
//    \param al input pointer to BamAlignment object that needs to be written
//  */
// double FilterBAMal::compLikelihoodSeq(BamAlignment * al){
//     int qualOffset=33;

//     vector<int> quals1;
//     for(int i=0;i<int(al->Qualities.length());i++){
// 	//if(al->QueryBases[i] != 'N')
// 	quals1.push_back( max( (int(char( al->Qualities[i] ))-qualOffset),2 ))  ; //since qual scores less than 2 do not make sense
//     }

//     double totalLike=0;
//     for(int i=0;i<int(quals1.size());i++){
// 	totalLike += likeQual[ quals1[i] ];
//     }
//     double likeSeq= pow(10.0,totalLike);
//     return likeSeq;
// }


//! Compute the expectation of error
/* This subroutine sums the expectation of error for a given sequence/quality and returns it.
   Since each quality score is the expectation of error for a base, the expectation for
   an entire sequence is given by the sum.
   
   \param al input pointer to BamAlignment object that needs to be written
   \return Expectation of error
 */
double FilterBAMal::computeExpectationError(BamAlignment * al){
    int qualOffset=33;

    vector<int> quals1;
    for(int i=0;i<int(al->Qualities.length());i++){
	//if(al->QueryBases[i] != 'N')
	quals1.push_back( max( (int(char( al->Qualities[i] ))-qualOffset),2 ))  ; //since qual scores less than 2 do not make sense	
    }

    double totalExp=0;
    for(int i=0;i<int(quals1.size());i++){
	totalExp += expError[ quals1[i] ];
    }
    //    cout<<totalExp<<endl;
    //double likeSeq= totalExp;///double( al->QueryBases.length() );
    double likeSeq= totalExp/double( al->QueryBases.length() );
    return likeSeq;
}


//! subroutine to flag as QC fail BamAlignment
/* This subroutine will use the different filters
   defined in the constructor.
   
   \param al input pointer to BamAlignment object that needs to be written
 */
void FilterBAMal::filterBAMAlign(BamAlignment * al){

    // if(al->IsMapped()) {
    // 	cerr<<"Read cannot be mapped"<<endl;
    // 	exit(1);
    // }
    repToPrint.totalSeq++;

    if(al->IsFailedQC()){
	repToPrint.qcBefore++;
    }

    //reset the flag
    if(resetQC){
	al->SetIsFailedQC(false);
    }

    

    if( (minLength != -1) &&
	(int(al->QueryBases.size()) < minLength)
	){
	al->SetIsFailedQC(true);
	repToPrint.qcLength++;
    }

    if( (maxLength != -1) &&
	(int(al->QueryBases.size()) > maxLength)
	){
	al->SetIsFailedQC(true);
	repToPrint.qcLength++;
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

//     double likeSeq=compLikelihoodSeq(al);

    
    //cout<<pow(10.0,totalLike)<<endl;
    //if(pow(10.0,totalLike) <0.2)
    //cout<<pow(10.0,totalLike)<<"\t"<<al->Qualities<<endl;


    double avgExp=computeExpectationError(al);
    // if(!al->IsPaired() ){
    //cout<<al->QueryBases.size()<<"\t"<<avgExp<<"\t"<<-10*log10(avgExp)<<endl;
    // }

    if(verbose){             
	cerr<<al->Name<<"\t"<<al->QueryBases<<"\t"<<al->Qualities<<"\t#"<<avgExp<<"#\t"<<-10*log10(avgExp)<<endl;
    }



    if( (avgExp) > cutoffAvgExpError ){
	//cout<<al->Qualities<<endl;
	al->SetIsFailedQC(true);
	repToPrint.qcFailExp++;
    }else{
	//cout<<al->Name<<"\t"<<al->Qualities<<endl;
    }
    
    //to uncomment
    // if(likeSeq < cutoffLikelihood)
    // 	al->SetIsFailedQC(true);

    if(likelihoodFlag){
	likelihoodOS->write( (char *)&avgExp, sizeof(avgExp));
    }

    if(entropy || entropyOSFlag  )
	complexEntropy(al);
    

    if(frequency || frequencyOSFlag)
	complexFrequency(al);



}


//! subroutine to print a log
/* This subroutine produces a string detailing what was done by the filterBAMAlign()
   subroutine.
   
   \return The log as a string
 */
string FilterBAMal::printLog(){
    stringstream st;
    st <<"Total reads                                               : "<<repToPrint.totalSeq<<endl;
    st <<"Reads that were previously flagged as QC failed           : "<<repToPrint.qcBefore<<endl;
    st <<"Reads that are QC failed due to length                    : "<<repToPrint.qcLength<<endl;
    st <<"Reads that are QC failed due to complexity                : "<<repToPrint.qcFailClx<<endl;
    st <<"Reads that are QC failed due to entropy                   : "<<repToPrint.qcFailEnt<<endl;
    st <<"Reads that are QC failed due to expectation of mismatches : "<<repToPrint.qcFailExp<<endl;
    return st.str();
}
