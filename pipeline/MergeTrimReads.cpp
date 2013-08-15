#include "MergeTrimReads.h"


// #define DEBUG2
// //#define DEBUGSR

//#define DEBUGADAPT
// #define DEBUGOVERLAP
// // //#define DEBUGTOTALOV
// // #define DEBUGPARTIALOV
// // #define CONSBASEPROB

#define DEBUGPR
#define DEBUGSCORE

double max_prob_N = 0.25;
// double likelihoodChimera   = -15.0;
double maxLikelihoodRatio = 0.95;

// size_t min_length = 5;
// double likelihoodAdapterSR = -0.4;
/* static const double likelihoodAdapterPR = -1.0; */
// double likelihoodAdapterPR = -120.0;
// const double max_prob_N = 0.25;

// const size_t min_length = 5;
// const int    qualOffset = 33;

const size_t min_length =5;
const int    qualOffset =33;

const char* const  chimInit[]= {
    "ACACTCTTTCCCTACACGTCTGAACTCCAG",
    "ACACTCTTTCCCACACGTCTGAACTCCAGT",
    "ACACTCTTTCCCTACACACGTCTGAACTCC",
    "CTCTTTCCCTACACGTCTGAACTCCAGTCA",
    "GAAGAGCACACGTCTGAACTCCAGTCACII",
    "GAGCACACGTCTGAACTCCAGTCACIIIII",
    "GATCGGAAGAGCACACGTCTGAACTCCAGT",
    "AGATCGGAAGAGCACACGTCTGAACTCCAG",
    "AGAGCACACGTCTGAACTCCAGTCACIIII",
    "ACACGTCTGAACTCCAGTCACIIIIIIIAT",
    "GTGCACACGTCTGAACTCCAGTCACIIIII",
    "AGCACACGTCTGAACTCCAGTCACIIIIII",
    "CGTATGCCGTCTTCTGCTTGAAAAAAAAAA"};

string options_adapter_F="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG";
string options_adapter_S="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";


vector<string> adapter_chimeras (chimInit,chimInit+13);




// bool   initialized       = false;
// double cutoff_merge_trim = 0.80;
size_t maxadapter_comp  = 30;
size_t min_overlap_seqs = 10;
// double cutoff_merge_seqs_early = 0.95;
// double cutoff_merge_seqs = 0.90;

// //  Key variables ///
 bool handle_key           = false;
string keys0="";
string keys1="";
int len_key1=0;
int len_key2=0;


size_t  options_trimCutoff   = 1;
bool    options_mergeoverlap = false;
bool    options_allowMissing = false;

double likeMatch[64];
double likeMismatch[64];

double probForQual[64];
double likeRandomMatch;    // 1/4
double likeRandomMisMatch; // 3/4

static inline char revComp(char c){
    if(c ==    'A')
	return 'T';

    if(c ==    'C')
	return 'G';

    if(c ==    'G')
	return 'C';

    if(c ==    'T')
	return 'A';

    if(c ==    'N')
	return 'N';
    
    cerr<<"Wrong base in revComp function = "<<c<<endl;
    exit(1);
}

static inline string revcompl(const string seq){
    string toReturn="";
    int seqL=seq.length();
    int i=0;
    while(i<seqL){
	toReturn+=revComp(seq[seqL-i-1]);
	i++;
    }
    return toReturn;
}

// void    setLikelihoodScores(double likelihoodChimera_,
// 			    double likelihoodAdapterSR_,
// 			    double likelihoodAdapterPR_){

//     likelihoodChimera   = likelihoodChimera_;
//     likelihoodAdapterSR = likelihoodAdapterSR_;
//     likelihoodAdapterPR = likelihoodAdapterPR_;
    
// }

void initMerge(){

    
    likeRandomMatch       = log( 1.0/4.0 )/log(10);
    likeRandomMisMatch    = log( 3.0/4.0 )/log(10);

#ifdef DEBUG2
    cout<<"likeRandomMisMatch "<<likeRandomMisMatch<<endl;
    cout<<"likeRandomMatch    "<<likeRandomMatch<<endl;
#endif

    //probability scores for qscores less than 2 make no sense
    //since random DNA has by default a probability of error of 0.75
    //anyway, we set the min for qual scores at 2
    for(int i=0;i<2;i++){

	
	likeMatch[i]        = log1p(    -pow(10.0,2.0/-10.0) )/log(10);
	    
        likeMismatch[i]     = log  (     pow(10.0,2.0/-10.0)/3.0 )/log(10);// 2.0/(-10.0*1.0);  

	probForQual[i]      = max(double(1.0)-pow(double(10.0),double(2.0)/double(-10.0)),
				  max_prob_N);

    }

    for(int i=2;i<64;i++){
        // if(i == 0)
        //     likeMatch[i]    = -3.0; // this is vrong, hope it's never accessed
        // else
	likeMatch[i]        = log1p(    -pow(10.0,i/-10.0) )/log(10);
	    
        likeMismatch[i]     = log  (     pow(10.0,i/-10.0)/3.0  )/log(10); //i/(-10.0*1.0);  

	probForQual[i] = max(double(1.0)-pow(double(10.0),double(i)/double(-10.0)),
			     max_prob_N);

#ifdef DEBUG2
        cout<<"qual = "<<i<<endl;
        cout<<likeMatch[i]<<endl;
        cout<<likeMismatch[i]<<endl;
#endif

    }


}

static inline string convert_logprob_quality(vector<int> logScores){
    string toReturn="";
    for(unsigned int i=0;i<logScores.size();i++){
	toReturn+=char(max(qualOffset,min(126,logScores[i]+qualOffset)));
    }
    return toReturn;
}


// static inline vector<double>  convert_logprob_prob(vector<int> logScores){
//     vector<double> toReturn;
//     for(unsigned int i=0;i<logScores.size();i++){
// 	toReturn.push_back( max(double(1.0)-pow(double(10.0),logScores[i]/double(-10.0)),
// 				max_prob_N) );
//     }
//     return toReturn;
// }

static inline double randomGen(){
    static bool initialized = false;

    if(initialized == false){
	struct timeval time;
	gettimeofday(&time,NULL);
	srand(time.tv_usec);
	initialized=true;
	return  (double(rand())/double(RAND_MAX));       
    }else{
	return  (double(rand())/double(RAND_MAX));       
    }
}

static inline baseQual cons_base_prob(baseQual  base1,baseQual base2){
    const string dnaAlphabet  = "ACGT";
    
    vector<char> bases;
    if(base1.base == base2.base){
	bases.push_back(base1.base);
    }else{
	bases.push_back(base1.base);
	bases.push_back(base2.base);
    }

    double aprob1 = log10(1.0-base1.prob)-log10(3.0); 
    double lprob1 = log10(base1.prob);
    double lprob2 = log10(base2.prob);
    double aprob2 = log10(1.0-base2.prob)-log10(3.0);

#ifdef CONSBASEPROB
    cerr<<base1.base<<"\t"<<base1.prob<<"\t"<<base2.base<<"\t"<<base2.prob<<"\t"<<endl;
#endif

    double total_prob = 0.0;
    for(int i=0;i<4;i++){
	double help = 0.0;
	if(base1.base == dnaAlphabet[i])
	    help+=lprob1;
	else 
	    help+=aprob1;
	if(base2.base == dnaAlphabet[i])
	    help+=lprob2;
	else
	    help+=aprob2;
	total_prob+=pow(10.0,help);
    }
    total_prob = log10(total_prob);

    baseQual toReturn;
    toReturn.base='N' ;
    toReturn.qual= -1;

    for(unsigned int i=0;i<bases.size();i++){
	double thelp=0.0;
	for(int k=0;k<4;k++){
	    if(bases[i] != dnaAlphabet[k]){
		double help=0.0;
		if(base1.base == dnaAlphabet[k])
		    help+=lprob1;
		else 
		    help+=aprob1;

		if(base2.base == dnaAlphabet[k])
		    help+=lprob2;
		else
		    help+=aprob2;
		thelp+=pow(double(10.0),help);
	    }
	}
	thelp=log10(thelp);
	int hqual = int( floor( min(60.0,-10.0*(thelp-total_prob)) + 0.5) );
       	if( (hqual > toReturn.qual) || ((hqual == toReturn.qual) && (randomGen() >=0.5)) ){
	    toReturn.base = bases[i];
	    toReturn.qual = hqual;
	}
    }

    return toReturn;
}


void set_adapter_sequences(const string& forward, const string& reverse, const string& chimera, int max_comp){
    options_adapter_F.assign( forward.length(), 0 );
    options_adapter_S.assign( reverse.length(), 0 );
    transform(forward.begin(), forward.end(), options_adapter_F.begin(), ::toupper);
    transform(reverse.begin(), reverse.end(), options_adapter_S.begin(), ::toupper);

    string tempChimera( chimera.length(), 0 );
    transform(chimera.begin(), chimera.end(), tempChimera.begin(), ::toupper);

    adapter_chimeras.clear();
    string token=returnFirstToken(&tempChimera,",");
    bool foundAdapterF=false;
    while(token.length()!=0){
	if(options_adapter_F == token)
	    foundAdapterF=true;
	adapter_chimeras.push_back(token);
	token=returnFirstToken(&tempChimera,",");	    
    }
    if(!foundAdapterF)
	adapter_chimeras.push_back(options_adapter_F);
    vector<string>::iterator it;
    it=adapter_chimeras.begin() ; 
    while(it < adapter_chimeras.end()  ){

	it->erase( remove( it->begin(), it->end(), ' ' ), it->end() ); //remove white spaces

	if( ( it->length() > 0) && (it->length()  <= maxadapter_comp)){
	    while ( it->length() <= maxadapter_comp ){
		*it += "I";
	    }
	    it++;
	}else{
	    if( it->empty() ){
		adapter_chimeras.erase(it);
	    }else{
		it++;
	    }
	}
    }
  
    
    while ( options_adapter_F.length() <= maxadapter_comp){
	options_adapter_F += "I";
    }

    while ( options_adapter_S.length() <= maxadapter_comp){
	options_adapter_S += "I";
    }


}
   

string returnFirstToken(string * toparse,string delim){
    size_t found;
    string toreturn;
    found=toparse->find(delim);
    if (found!=string::npos){
	toreturn=toparse->substr(0,found);
	toparse->erase(0,found+delim.length());
	return toreturn;
    }
    toreturn=(*toparse);
    toparse->erase(0,toparse->length());    
    return toreturn;
}

void set_keys(const string& key1, const string& key2){

    keys0=key1;
    if(key2==""){
	keys1=key1;
    }else{
	keys1=key2;
    }

    len_key1=keys0.length();
    len_key2=keys1.length();

    transform(keys0.begin(), keys0.end(),keys0.begin(), ::toupper);
    transform(keys1.begin(), keys1.end(),keys1.begin(), ::toupper);

    if( (len_key1 > 0) || (len_key2 > 0))
	handle_key=true;
    else
	handle_key=false;
	
}

void set_options(int trimcutoff,bool allowMissing,bool mergeoverlap){
    options_trimCutoff   = trimcutoff;
    options_allowMissing = allowMissing;
    options_mergeoverlap = mergeoverlap;
}


static inline double detectChimera(const string      & read,
				   const vector<int> & qual,
				   const string      & chimeraString,
				   unsigned int        offsetChimera=0){

    double likelihoodMatch=0.0;
    unsigned int maxidx=min(read.length(),chimeraString.length()-offsetChimera);

    if(maxidx <= 0)
	return -DBL_MAX;

    for(unsigned i=0;i<maxidx;i++){
	// cerr<<"i= "<<i<<endl;
	if(          read[i]                        == chimeraString[i+offsetChimera] || 
		     chimeraString[i+offsetChimera] == 'I' ){
	    likelihoodMatch  +=    likeMatch[ qual[i] ];
	}else{
	    likelihoodMatch  += likeMismatch[ qual[i] ];
	}
	// cerr<<"i= "<<likelihoodMatch<<endl;
    }
    //add likelihood of remaining bases
    //    cout<<double(read.length() - maxidx )<<endl;
    likelihoodMatch  += double(read.length() - maxidx ) * likeRandomMatch;
    
    return likelihoodMatch;
}





static inline double measureOverlap(const string      & read1,
				    const vector<int> & qual1,
				    const string      & read2,
				    const vector<int> & qual2,
				    const int         & maxLengthForPair,
				    unsigned int      offsetRead=0,				    
				    double *          iterations =0 ,
				    int  *            matches=0){
    
#ifdef DEBUGOVERLAP
    string comparedread1;
    string comparedread2;
#endif

    int i1=0; //start index read 1
    int i2=maxLengthForPair-offsetRead;
    double likelihoodMatch=0.0;

    for(int i=0;i<int(offsetRead);i++){
	
	if(     i1<int(read1.size()) && 
		i2<int(read2.size()) && 
		i2>=0){
		
#ifdef DEBUGOVERLAP
	    // cerr<<"overlap "<<" read1["<<(i1)<<"] "<<read1[i1]<<"\tread2["<<i2<< "] "<<read2[i2]<<endl;
	    comparedread1+=read1[i1];
	    comparedread2+=read2[i2];
#endif

	    if(read1[i1] == read2[i2] ){	    
		(*matches)++;
		likelihoodMatch  +=    likeMatch[    min(qual1[i1],qual2[i2])  ];	    
	    }else{
		likelihoodMatch  +=    likeMismatch[ min(qual1[i1],qual2[i2])  ];	    
	    }
	    likelihoodMatch  += likeRandomMatch; 

	    
	    i1++;
	    i2++;
	    continue;
	}

	
	if( i1<int(read1.size()) ){

#ifdef DEBUGOVERLAP
	    //cerr<<"overlap "<<" read1["<<(i1)<<"] "<<read1[i1]<<"\tread2["<<i2<< "] "<<read2[i2]<<endl;
	    comparedread1+=read1[i1];
	    comparedread2+="-";
#endif

	    likelihoodMatch  += likeRandomMatch; 

	
	    i1++;
	    i2++;
	    continue;
	}

	if( i2<int(read2.size()) && 
	    i2>=0 ){
		
#ifdef DEBUGOVERLAP
	    //cerr<<"overlap "<<" read1["<<(i1)<<"] "<<read1[i1]<<"\tread2["<<i2<< "] "<<read2[i2]<<endl;
	    comparedread1+="-";
	    comparedread2+=read2[i2];
#endif
	    likelihoodMatch  += likeRandomMatch; 
	    
	    i1++;
	    i2++;
	    continue;
	}

	cerr<<"Wrong state"<<endl;
	exit(1);
	    
    }	
	

#ifdef DEBUGOVERLAP
    cerr<<"comparedread1 "<<comparedread1<<endl;
    cerr<<"comparedread2 "<<comparedread2<<endl;
    cerr<<"result        "<<likelihoodMatch<<endl;
#endif

    return likelihoodMatch;
}


// static inline double measureOverlap(const string      & read1,
// 				    const vector<int> & qual1,
// 				    const string      & read2,
// 				    const vector<int> & qual2,
// 				    const int startRead1,				   
// 				    const int startRead2,				   
// 				    int	maxLength,
// 				    double * iterations=0,
// 				    int  *  matches=0){
//     if(maxLength < 0)
// 	return -DBL_MAX;
    
//     double likelihoodMatch=0.0;

//     //index for both reads:
//     int i1=startRead1;
//     int i2=startRead2;

// #ifdef DEBUGOVERLAP
//     string comparedread1;
//     string comparedread2;
// #endif

//     //iterate over r1
//     for(int i=0;i<maxLength;i++){


// #ifdef DEBUGOVERLAP
// 	//cerr<<"overlap "<<" read1["<<(i1)<<"] "<<read1[i1]<<"\tread2["<<i2<< "] "<<read2[i2]<<endl;
// 	comparedread1+=read1[i1];
// 	comparedread2+=read2[i2];
// #endif

// 	//first base
// 	if(read1[i1] == read2[i2] ){	    
// 	    (*matches)++;
// 	    likelihoodMatch  +=    likeMatch[    min(qual1[i1],qual2[i2])  ];	    
// 	}else{
// 	    likelihoodMatch  +=    likeMismatch[ min(qual1[i1],qual2[i2])  ];	    
// 	}

// 	//second base
// 	likelihoodMatch  += likeRandomMatch; 


// 	(*iterations)++;

// 	i1++;
// 	i2++;
//     }

// #ifdef DEBUGOVERLAP
//     cerr<<"comparedread1 "<<comparedread1<<endl;
//     cerr<<"comparedread2 "<<comparedread2<<endl;
//     cerr<<"result        "<<likelihoodMatch<<endl;
// #endif

//     return likelihoodMatch;

// }


    
/*! Computes the likelihood of seeing the adapter in read 
    at index offsetRead  
*/  
static inline double detectAdapter(const string      & read,
				   const vector<int> & qual,
				   const string      & adapterString,
				   unsigned int      offsetRead=0,
				   double *          iterations =0 ,
				   int  *            matches=0){
    
    double likelihoodMatch=0.0;
    
    unsigned maxIterations=min ( min( (read.length()-offsetRead),
			       adapterString.length())  , maxadapter_comp  );


#ifdef DEBUGADAPT
    string comparedread;
#endif
    unsigned  i=0;
    unsigned  indexRead=offsetRead;

    while(i<maxIterations){
	indexRead = i+offsetRead;

#ifdef DEBUGADAPT
	comparedread+=read[indexRead];
#endif

	
	if(         read[indexRead]           == adapterString[i] || 
		    adapterString[indexRead]  == 'I' ){ //match
	    (*matches)++;
	    likelihoodMatch  +=    likeMatch[ qual[indexRead] ];
	}else{ //mismatch
	    likelihoodMatch  += likeMismatch[ qual[indexRead] ];
	}
	(*iterations)++;
	i++;

    }

    int extraBases = max(0,int(read.length())-int(indexRead)-1);

#ifdef DEBUGADAPT
    cerr<<"detectAdapter1 "<<read<<"\n              "<<comparedread <<"\n              "<<adapterString<<"\tlm:"<<likelihoodMatch<<"\trl:"<<read.length()<<"\tor:"<<offsetRead<<"\tadd:"<<"\tmi:"<<maxIterations<<"\teb:"<<extraBases<<endl;
#endif

   

    likelihoodMatch += double( extraBases ) * likeRandomMatch;


#ifdef DEBUGADAPT
    cerr<<"detectAdapter2 "<< comparedread <<"\n              "<<adapterString<<"\tlm:"<<likelihoodMatch<<"\trl:"<<read.length()<<"\tor:"<<offsetRead<<endl;
    //cerr<<double(read.length()-indexRead)<<endl;
#endif

    return likelihoodMatch;
}


static inline int edits(const string & seq1,const string & seq2){
    int lmin = min(seq1.length(),seq2.length());
    int lmax = max(seq1.length(),seq2.length());
    int dist = lmax-lmin;

    for(int pos=0;pos<lmin;pos++){
	if (seq1[pos] != seq2[pos]) 
	    dist+=1;	
    }

    return dist;
}

merged process_SR(string  read1, 
		  string  qual1){

    merged toReturn;
    // int qualOffset=33;


    //check for the key (if needed)
    if( handle_key ){
	if ((read1.substr(0,len_key1) == keys0)  || //perfect match
	    (options_allowMissing && (edits(read1.substr(0,len_key1),keys0) == 1))  //1mm in first key	     
	    ){ //1mm in second key
	    read1 = read1.substr(len_key1,read1.size()-len_key1);
	    qual1 = qual1.substr(len_key1,read1.size()-len_key1);
	}else{
	    if(options_allowMissing && 
	       (read1.substr(0,len_key1-1) == keys0.substr(1,len_key1-1))  ){
		read1 = read1.substr(len_key1-1,read1.size()-(len_key1-1));
		qual1 = qual1.substr(len_key1-1,read1.size()-(len_key1-1));
	    }else{
		toReturn.code    ='K';
		toReturn.sequence="";
		toReturn.quality ="";		    
		return toReturn;	       
	    }
	}
    }
    //end check for the key (if needed)


    if( (read1.length() != qual1.length()) ){
	cerr<<"MergeTrimReads: The reads and qualities must have equal lengths"<<endl;
	exit(1);
    }

    vector<int> qualv1;
    for(unsigned int i=0;i<qual1.length();i++){
	qualv1.push_back( max( (int(char( qual1[i] ))-qualOffset),2) );
    }
    

    //start detecting chimera //
    double logLikelihoodTotal     = double(read1.length())*likeRandomMatch;
    int logLikelihoodTotalIdx     = -1;

    // double lowChimeraLike=-DBL_MAX;
    //finding best match
    for(unsigned int indexChimera=0;indexChimera<adapter_chimeras.size();indexChimera++){
	//lowChimeraLike = max(lowChimeraLike,detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 0 ) ); 
	
	if( detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 0 )  > logLikelihoodTotal ){
	    toReturn.code    ='D';
	    toReturn.sequence="";
	    toReturn.quality ="";		    
	    return toReturn;
	}

	//lowChimeraLike = max(lowChimeraLike,detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 1 ) ); //try an off by 1 match

	if( detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 1 )  > logLikelihoodTotal ){
	    toReturn.code    ='D';
	    toReturn.sequence="";
	    toReturn.quality ="";
	    return toReturn;
	}

	// cout<<"res "<<adapter_chimeras[indexChimera]<<"\t"<<detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera],1 )<<endl;
    }

    // end detecting chimera //



    //start detecting adapter //
    // double       lowAdapterLike   =-DBL_MAX;
    // int       indexAdapterBest = read1.size();
    
    // double logLikelihoodTotal     = double(read1.length())*likeRandomMatch ;
    // int logLikelihoodTotalMatches =0;

    //second best hit
    double sndlogLikelihoodTotal      = -DBL_MAX;
    int sndlogLikelihoodTotalIdx      = -1;
    // int sndlogLikelihoodTotalMatches  =  0;



    for(unsigned int indexAdapter=0;
	indexAdapter<(read1.length()-options_trimCutoff);
	indexAdapter++){
	double iterations=0;
	double logLike=detectAdapter( read1 , qualv1 , options_adapter_F,indexAdapter,&iterations );



	if(logLikelihoodTotal < logLike){
    	    sndlogLikelihoodTotal        = logLikelihoodTotal;
    	    sndlogLikelihoodTotalIdx     = logLikelihoodTotalIdx;
    	    // sndlogLikelihoodTotalMatches = logLikelihoodTotalMatches;
	    
    	    logLikelihoodTotal           = logLike;
    	    logLikelihoodTotalIdx        = indexAdapter;
    	    //logLikelihoodTotalMatches    = matches;

    	}else{ 
    	    if(sndlogLikelihoodTotal < logLike){ //not more likely than first but more likely than second
    		sndlogLikelihoodTotal        = logLike;
    		sndlogLikelihoodTotalIdx     = indexAdapter;
		//    		sndlogLikelihoodTotalMatches = matches;
    	    }
    	}

#ifdef DEBUGSR
	 cerr<<indexAdapter<<"\t"<<tm<<endl;
#endif
    }
#ifdef DEBUGSR
     cerr<<logLikelihoodTotal<<"\t"<<endl;
#endif
     // cerr<<lowAdapterLike<<endl;

    // if (lowAdapterLike > likelihoodAdapterSR  ) {
    if( logLikelihoodTotalIdx != -1    &&                       //the status quo is not the most likely
	(logLikelihoodTotal/sndlogLikelihoodTotal)      < maxLikelihoodRatio ){
     
	
        read1 = read1.substr(0,logLikelihoodTotalIdx);
        qual1 = qual1.substr(0,logLikelihoodTotalIdx);



	if( read1.length() < min_length){
	    toReturn.code    ='D';
	    toReturn.sequence="";	   
	    toReturn.quality ="";
	    return toReturn;	
	}else{
	    toReturn.code    =' ';
	    toReturn.sequence=read1;	   
	    toReturn.quality =qual1;	
	    return toReturn;
	}


    }


    toReturn.code    =' ';
    toReturn.sequence="";	   
    toReturn.quality ="";
    return toReturn;
  
}

merged process_PE( string  read1,  string  qual1,
		   string  read2,  string  qual2){


    merged toReturn;

    if( handle_key && read1.length() > 0){

	if (((read1.substr(0,len_key1) == keys0) && (read2.substr(0,len_key2) == keys1)) || //perfect match
	    (options_allowMissing && (edits(read1.substr(0,len_key1),keys0) == 1) && (read2.substr(0,len_key2) == keys1)) || //1mm in first key
	    (options_allowMissing && (read1.substr(0,len_key1) == keys0) && (edits(read2.substr(0,len_key2),keys1) == 1))
	    ){ //1mm in second key
	    read1 = read1.substr(len_key1,read1.size()-len_key1);
	    qual1 = qual1.substr(len_key1,read1.size()-len_key1);
	    read2 = read2.substr(len_key2,read2.size()-len_key2);
	    qual2 = qual2.substr(len_key2,read2.size()-len_key2);
	}else{
	    if(options_allowMissing && 
	       ( (len_key1>0?(read1.substr(0,len_key1-1) == keys0.substr(1,len_key1-1)):true) && (read2.substr(0,len_key2) == keys1)) ){
		read1 = read1.substr(len_key1-1,read1.size()-(len_key1-1));
		qual1 = qual1.substr(len_key1-1,read1.size()-(len_key1-1));
		read2 = read2.substr(len_key2,read2.size()-len_key2);
		qual2 = qual2.substr(len_key2,read2.size()-len_key2);
	    }else{
		if(options_allowMissing && 
		   (read1.substr(0,len_key1) == keys0) && (len_key2>0?(read2.substr(0,len_key2-1) == keys1.substr(1,len_key2-1)):true)){
		    read1 = read1.substr(len_key1,read1.size()-len_key1);
		    qual1 = qual1.substr(len_key1,read1.size()-len_key1);
		    read2 = read2.substr(len_key2-1,read2.size()-(len_key2-1));
		    qual2 = qual2.substr(len_key2-1,read2.size()-(len_key2-1));		    
		}else{
		    toReturn.code    ='K';
		    toReturn.sequence="";
		    toReturn.quality ="";		    
		    return toReturn;
		}

	    }
	}
    }

    if( (read1.size() != qual1.size()) ||
	(read2.size() != qual2.size()) ){
	cerr<<"MergeTrimReads: The reads and qualities must have equal lengths"<<endl;
	exit(1);
    }

    vector<int> qualv1;
    for(unsigned int i=0;i<qual1.length();i++){
	qualv1.push_back( max( (int(char( qual1[i] ))-qualOffset),2) );
    }

    vector<int> qualv2;
    for(unsigned int i=0;i<qual2.length();i++){
	qualv2.push_back( max( (int(char( qual2[i] ))-qualOffset),2) );
    }
    
    int lengthRead1 = int(read1.size());
    int lengthRead2 = int(read2.size());

    double logLikeSecondRead      = double(lengthRead2)*likeRandomMatch;
    double logLikelihoodTotal     = double(lengthRead1)*likeRandomMatch   + logLikeSecondRead;


    int logLikelihoodTotalIdx     = -1;
    int logLikelihoodTotalMatches =0;

    //start detecting chimera //
    double lowChimeraLike=-DBL_MAX;
    //finding best match
    for(unsigned int indexChimera=0;indexChimera<adapter_chimeras.size();indexChimera++){
	//cerr<<"indexChimera1 "<<lowChimeraLike<<endl;
	//lowChimeraLike = max(lowChimeraLike,detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 0 ) ); 
	//lowChimeraLike = max(lowChimeraLike,detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 1 ) ); //try an off by 1 match
	double logLikeChimera = (detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 0 ) + logLikeSecondRead );
	if( logLikeChimera  > logLikelihoodTotal  ){
	    //cout<<"logLikelihoodTotal "<<logLikelihoodTotal<<"\t"<<logLikeChimera<<"\t"<<read1<<"\t"<<adapter_chimeras[indexChimera]<<endl;
	    toReturn.code    ='D';
	    toReturn.sequence="";
	    toReturn.quality ="";		    
	    return toReturn;
	}



	logLikeChimera = (detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 1 ) + logLikeSecondRead );
	if( logLikeChimera > logLikelihoodTotal  ){  //try an off by 1 match
	    toReturn.code    ='D';
	    toReturn.sequence="";
	    toReturn.quality ="";		    
	    return toReturn;
	}


	
	//cerr<<"indexChimera2 "<<lowChimeraLike<<endl;
	// cout<<"res "<<adapter_chimeras[indexChimera]<<"\t"<<detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera],1 )<<endl;
    }

    //cerr<<"TEST "<<lowChimeraLike<<"\t"<<likelihoodChimera<<endl;
   
    // end detecting chimera //


    //computing rev compl for read 2
    string  read2_rev=  revcompl(read2);
    string  qual2_rev=           qual2;

    reverse(qual2_rev.begin(), 
	    qual2_rev.end());
    
    vector<int> qualv2_rev (qualv2);

    reverse(qualv2_rev.begin(), 
	    qualv2_rev.end());

   




    //For a given putative index, we compute the likelihood of seeing:
    //
    // 
    //         L(adapter in read2)+  L(likelihood overlap) + L(adapter in read1)
    //                                            |--------> (adapter)
    //read1                         ---------------------------------
    //read2      ---------------------------------
    //                     <--------| (adapter 2)
    // let :
    // logLike1 = L(adapter in read1)
    // logLike2 = L(adapter in read2)
    // logLike3 = L(adapter in overlapping region)

    // logLike1 = L(seeing adapter bases in read1) + L(remaining bases)
    // logLike2 = L(seeing adapter bases in read2) + L(remaining bases)
    // logLike3 = L(bases matching read1/read2 overlap) + L(bases read1/read2 overlap) + L(any remaining bases
    
    // The default case is simply the likelihood of both :
    // 
  


    //second best hit
    double sndlogLikelihoodTotal      = -DBL_MAX;
    int sndlogLikelihoodTotalIdx      = -1;
    int sndlogLikelihoodTotalMatches  =  0;

#ifdef DEBUGPR

    cerr<<"fst: "<<read1<<endl<<"raw: "<<read2<<endl<<"rev: "<<read2_rev<<endl<<endl;
#endif
  

    
#ifdef DEBUGSCORE
    cerr<<"mergeo "<<"\t"<<logLikelihoodTotal<<"\t"<<logLikelihoodTotalIdx<<"\t"<<logLikelihoodTotalMatches<<endl;
    //cout<<logLikelihoodTotal<<"\t"<<logLikelihoodTotalMatches<<endl;
#endif
    

    int maxLengthForPair=max(lengthRead1,lengthRead2);
    //int lengthDiffR1_R2 =  (lengthRead1-lengthRead2);

    for(int indexAdapter=0; //let index adapters be the index of the potential P7/P5 adapter
	indexAdapter<(2*maxLengthForPair-min_overlap_seqs);
	indexAdapter++){
	double logLike=0.0;

 	double iterations=0;
 	int    matches   =0;

	if(indexAdapter > maxLengthForPair) //partial overlap
	    if(!options_mergeoverlap) //no point in continuing 
		break;
	    
	//adapters
 	if(indexAdapter<lengthRead1)
 	    logLike  +=  detectAdapter(    read1 ,     qualv1     , options_adapter_F,indexAdapter,&iterations, &matches);

 	if(indexAdapter<lengthRead2)
 	    logLike  +=  detectAdapter(    read2 ,     qualv2     , options_adapter_S,indexAdapter,&iterations, &matches );

	
	//overlap
	logLike  +=	 measureOverlap(   read1 ,
					   qualv1     ,
					   read2_rev ,
					   qualv2_rev,  
					   maxLengthForPair,					      
					   indexAdapter,
					   &iterations,
					   &matches);


#ifdef DEBUGPR
    	cerr<<"idx: "<<indexAdapter<<"\ttl:"<<logLike<<"\tlr:"<<"\tit:"<<iterations<<"\tma:"<<matches<<endl;
    	//cerr<<logLike1<<"\t"<<logLike2<<"\t"<<logLike3<<endl;
#endif
	
    	if(logLikelihoodTotal < logLike){
    	    sndlogLikelihoodTotal        = logLikelihoodTotal;
    	    sndlogLikelihoodTotalIdx     = logLikelihoodTotalIdx;
    	    sndlogLikelihoodTotalMatches = logLikelihoodTotalMatches;
	    
    	    logLikelihoodTotal           = logLike;
    	    logLikelihoodTotalIdx        = indexAdapter;
    	    logLikelihoodTotalMatches    = matches;
    	}else{ 
    	    if(sndlogLikelihoodTotal < logLike){ //not more likely than first but more likely than second
    		sndlogLikelihoodTotal        = logLike;
    		sndlogLikelihoodTotalIdx     = indexAdapter;
    		sndlogLikelihoodTotalMatches = matches;
    	    }
    	}
        

    }



    
    
    
#ifdef DEBUGSCORE
    cerr<<"mergeo "<<"\t"<<logLikelihoodTotal<<"\t"<<logLikelihoodTotalIdx<<"\t"<<logLikelihoodTotalMatches<<"\trt:"<<(sndlogLikelihoodTotal/logLikelihoodTotal)<<"\t"<<sndlogLikelihoodTotal<<"\t"<<(logLikelihoodTotal/sndlogLikelihoodTotal)<<endl;
    //cout<<logLikelihoodTotal<<"\t"<<logLikelihoodTotalMatches<<endl;
#endif
    

    


    baseQual b1;
    baseQual b2;
    //test for merging:
    if( logLikelihoodTotalIdx != -1    &&                       //the status quo is not the most likely
	(logLikelihoodTotal/sndlogLikelihoodTotal)      < maxLikelihoodRatio &&
	logLikelihoodTotalMatches >= min_overlap_seqs  ){       //sufficient # of mismatches for partial overlap (artifically to min_overlap_seqs for complete overlap)



	int i1=0; //start index read 1
	int i2=maxLengthForPair-logLikelihoodTotalIdx;

	string        newSeq  = "";
	vector<int>   newQual ;

	for(int i=0;i<int(logLikelihoodTotalIdx);i++){



	    if(     i1<int(read1.size()) && 
		    i2<int(read2.size()) && 
		    i2>=0){
		
		b1.base = read1[i1];
		b1.prob = probForQual[ qualv1[i1] ];
		b1.qual = -1;

		b2.base = read2_rev[i2];
		b2.prob = probForQual[ qualv2_rev[i2] ]; 
		b2.qual = -1;

		baseQual RT    = cons_base_prob(b1,b2);
		newSeq         +=  RT.base ;
		newQual.push_back( RT.qual );

		i1++;
		i2++;
		continue;
	    }

	
	    if( i1<int(read1.size()) ){

		newSeq              +=	 read1[i1] ;
		newQual.push_back(      qualv1[i1] );

		i1++;
		i2++;
		continue;
	    }

	    if( i2<int(read2.size()) && 
		i2>=0 ){
		

		newSeq              +=	read2_rev[i2] ;
		newQual.push_back(     qualv2_rev[i2] );
	    
		i1++;
		i2++;
		continue;
	    }

	    cerr<<"Wrong state"<<endl;
	    exit(1);
	    
	}



	toReturn.code    =' ';
	toReturn.sequence=newSeq;	   
	toReturn.quality =convert_logprob_quality(newQual);	

	// cout<<toReturn.sequence<<endl;
	return toReturn;

	
    }else{

	toReturn.code    =' ';
	toReturn.sequence="";	   
	toReturn.quality ="";
	return toReturn;    
    }
    


}




