#include "MergeTrimReads.h"


// #define DEBUG2
// //#define DEBUGSR
// #define DEBUGPR
// // #define DEBUGADAPT
// #define DEBUGOVERLAP
// //#define DEBUGTOTALOV
// #define DEBUGPARTIALOV
// #define CONSBASEPROB






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

void initMerge(){

    
    likeRandomMatch       = log1p( -1.0/4.0 )/log(10);
    likeRandomMisMatch    = log1p( -3.0/4.0 )/log(10);
    
    //probability scores for qscores less than 2 make no sense
    //since random DNA has by default a probability of error of 0.75
    //anyway, we set the min for qual scores at 2
    for(int i=0;i<2;i++){

	
	likeMatch[i]        = log1p( -pow(10.0,2.0/-10.0) )/log(10);
	    
        likeMismatch[i]     = 2.0/-10.0;  

	probForQual[i]      = max(double(1.0)-pow(double(10.0),double(2.0)/double(-10.0)),
				  max_prob_N);

    }

    for(int i=2;i<64;i++){
        // if(i == 0)
        //     likeMatch[i]    = -3.0; // this is vrong, hope it's never accessed
        // else
	likeMatch[i]        = log1p( -pow(10.0,i/-10.0) )/log(10);
	    
        likeMismatch[i]     = i/-10.0;  

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
   

static string returnFirstToken(string * toparse,string delim){
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
    unsigned maxidx=min(read.length(),chimeraString.length()-offsetChimera);

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
    
    return likelihoodMatch;
}



static inline double measureOverlap(const string      & read1,
				    const vector<int> & qual1,
				    const string      & read2,
				    const vector<int> & qual2,
				    const int startRead1,				   
				    const int startRead2,				   
				    int	maxLength,
				    double * iterations=0,
				    int  *  matches=0){
    if(maxLength < 0)
	return -DBL_MAX;
    
    // unsigned int maxidx=endRead ;
    double likelihoodMatch=0.0;
    int i1=startRead1;
    int i2=startRead2;

#ifdef DEBUGOVERLAP
    //cerr<<"st1:"<<startRead1<<"\t"<<"st2:"<<startRead2<<"\tml:"<<maxLength<<endl;
    string comparedread1;
    string comparedread2;
#endif

    //iterate over r1
    for(int i=0;i<maxLength;i++){


#ifdef DEBUGOVERLAP
	//cerr<<"overlap "<<" read1["<<(i1)<<"] "<<read1[i1]<<"\tread2["<<i2<< "] "<<read2[i2]<<endl;
	comparedread1+=read1[i1];
	comparedread2+=read2[i2];
#endif
	
	if(read1[i1] == read2[i2] ){	    
	    (*matches)++;
	    likelihoodMatch  +=    likeMatch[    min(qual1[i1],qual2[i2])  ];
	}else{
	    likelihoodMatch  +=    likeMismatch[ min(qual1[i1],qual2[i2])  ];
	}
	// iterations++;
	(*iterations)++;
	// i++;
	i1++;
	i2++;
    }

#ifdef DEBUGOVERLAP
    cerr<<"comparedread1 "<<comparedread1<<endl;
    cerr<<"comparedread2 "<<comparedread2<<endl;
    cerr<<"result        "<<likelihoodMatch<<endl;
#endif

    return likelihoodMatch;

}
    

static inline double detectAdapter(const string      & read,
				   const vector<int> & qual,
				   const string      & adapterString,
				   unsigned int offsetRead=0,
				   double * iterations =0 ,
				   int  *  matches=0){

    double likelihoodMatch=0.0;
    unsigned maxidx=min ( min(read.length()-offsetRead,
			      adapterString.length())  , maxadapter_comp  );
    // cerr<<read<<endl;
    // cerr<<adapterString<<endl;

    unsigned int i1;

#ifdef DEBUGADAPT
    string comparedread;
#endif

    for(unsigned i=0;i<maxidx;i++){
	i1 = i+offsetRead;

#ifdef DEBUGADAPT
	//cerr<<"da "<<read[i1] <<"\t"<<adapterString[i]<<"\t"<<likelihoodMatch<<endl;
	comparedread+=read[i1];
#endif
	if(         read[i1]  == adapterString[i] || 
		    adapterString[i1]  == 'I' ){ //match
	    (*matches)++;
	    likelihoodMatch  +=    likeMatch[ qual[i1] ];
	}else{ //mismatch
	    likelihoodMatch  += likeMismatch[ qual[i1] ];
	}
	(*iterations)++;
    }
    
#ifdef DEBUGADAPT
	cerr<<"detectAdapter "<< comparedread <<"\n              "<<adapterString<<"\t"<<likelihoodMatch<<endl;
#endif

    return likelihoodMatch;
}

// static inline double detectAdapterRev(const string      & read,
// 		     const vector<int> & qual,
// 		     const string      & adapterString,
// 		     unsigned int offsetRead){

//     double likelihoodMatch=0.0;
//     unsigned maxidx=min ( min(read.length()-offsetRead,
// 			      adapterString.length())  , maxadapter_comp  );
//     // cerr<<read<<endl;
//     // cerr<<adapterString<<endl;
//     unsigned int i2;

//     double iterations=0;
//     for(unsigned i=0;i<maxidx;i++){
// 	i2 = read.size()-i-offsetRead-1;

// #ifdef DEBUGADAPT
// 	cerr<<"dar "<<i2<<"\t"<<read[i2] <<"\t"<<adapterString[i]<<"\t"<<likelihoodMatch<<endl;
// #endif
// 	if(         read[i2]  == adapterString[i] || 
// 		    adapterString[i2]  == 'I' ){ //match
// 	    likelihoodMatch  +=    likeMatch[ qual[i2] ];
// 	}else{ //mismatch
// 	    likelihoodMatch  += likeMismatch[ qual[i2] ];
// 	}
// 	iterations++;
//     }
    
//     return likelihoodMatch/iterations;
// }

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
    double lowChimeraLike=-DBL_MAX;
    //finding best match
    for(unsigned int indexChimera=0;indexChimera<adapter_chimeras.size();indexChimera++){
	lowChimeraLike = max(lowChimeraLike,detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 0 ) ); 
	lowChimeraLike = max(lowChimeraLike,detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 1 ) ); //try an off by 1 match
	// cout<<"res "<<adapter_chimeras[indexChimera]<<"\t"<<detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera],1 )<<endl;
    }

    if(lowChimeraLike > likelihoodChimera){
	toReturn.code    ='D';
	toReturn.sequence="";
	toReturn.quality ="";		    
	return toReturn;
    }
    // end detecting chimera //



    //start detecting adapter //
    double       lowAdapterLike   =-DBL_MAX;
    unsigned int indexAdapterBest = read1.size();

    for(unsigned int indexAdapter=0;
	indexAdapter<(read1.length()-options_trimCutoff);
	indexAdapter++){
	double iterations=0;
	double tm=detectAdapter( read1 , qualv1 , options_adapter_F,indexAdapter,&iterations );
	tm=tm/iterations;
	if( lowAdapterLike   <  tm){
	    lowAdapterLike   =  tm;
	    indexAdapterBest =  indexAdapter;
	}
#ifdef DEBUGSR
	 cerr<<indexAdapter<<"\t"<<tm<<endl;
#endif
    }
#ifdef DEBUGSR
     cerr<<lowAdapterLike<<"\t"<<indexAdapterBest<<endl;
#endif
     // cerr<<lowAdapterLike<<endl;

    if (lowAdapterLike > likelihoodAdapterSR  ) {
       
	
        read1 = read1.substr(0,indexAdapterBest);
        qual1 = qual1.substr(0,indexAdapterBest);



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


	// cerr<<read1<<"\t"<<indexAdapterBest<<endl;
	// cerr<<lowAdapterLike<<endl;
    }


    //exit(1);
    //cerr<<
    //end detecting adapter //

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
    

    //start detecting chimera //
    double lowChimeraLike=-DBL_MAX;
    //finding best match
    for(unsigned int indexChimera=0;indexChimera<adapter_chimeras.size();indexChimera++){
	//cerr<<"indexChimera1 "<<lowChimeraLike<<endl;
	lowChimeraLike = max(lowChimeraLike,detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 0 ) ); 
	lowChimeraLike = max(lowChimeraLike,detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera], 1 ) ); //try an off by 1 match
	//cerr<<"indexChimera2 "<<lowChimeraLike<<endl;
	// cout<<"res "<<adapter_chimeras[indexChimera]<<"\t"<<detectChimera( read1 , qualv1 , adapter_chimeras[indexChimera],1 )<<endl;
    }

    //cerr<<"TEST "<<lowChimeraLike<<"\t"<<likelihoodChimera<<endl;

    if(lowChimeraLike > likelihoodChimera){
	toReturn.code    ='D';
	toReturn.sequence="";
	toReturn.quality ="";		    
	return toReturn;
    }
    // end detecting chimera //


    //computing rev compl for read 2
    string  read2_rev=  revcompl(read2);
    string  qual2_rev=           qual2;

    reverse(qual2_rev.begin(), 
	    qual2_rev.end());
    
    vector<int> qualv2_rev (qualv2);

    reverse(qualv2_rev.begin(), 
	    qualv2_rev.end());

    double logLikeTotalOverlap    = -DBL_MAX;
    double logLikeTotalOverlapIdx = 0;
    int    logLikeTotalOverlapMatches=0;


    //second best hit
    double sndlogLikeTotalOverlap    = -DBL_MAX;
    double sndlogLikeTotalOverlapIdx = 0;
    int    sndlogLikeTotalOverlapMatches=0;

    // double logLikePartialOverlap  =-DBL_MAX;
    // double logLikePartialOverlapIdx=0;







    //start  detecting partial overlap//

    // double       lowAdapterLike   =-DBL_MAX;


#ifdef DEBUGPR

    cerr<<"fst: "<<read1<<endl<<"raw: "<<read2<<endl<<"rev: "<<read2_rev<<endl<<endl;
#endif
    int size1 = int(read1.size());
    int size2 = int(read2.size());

    //int minLengthForPair=min(size1,size2);
    int maxLengthForPair=max(size1,size2);

    int lengthDiffR1_R2 =  (size1-size2);
    // int lengthDiffR2_R1 =  (size2-size1);



    for(int indexAdapter=0; //let index adapters be the index of the potential P7/P5 adapter
	indexAdapter<(2*maxLengthForPair-min_overlap_seqs);
	indexAdapter++){


#ifdef DEBUGPR
	cerr<<"idx: "<<indexAdapter<<endl;
#endif
    

	double iterations=0;
	int    matches   =0;

	double logLike1=0.0; //p7 likelihood
	double logLike2=0.0; //p5 likelihood
	double logLike3=0.0; //overlap likelihood

	if(indexAdapter<size1)
	    logLike1  =  detectAdapter(    read1 ,     qualv1     , options_adapter_F,indexAdapter,&iterations, &matches);
	if(indexAdapter<size2)
	    logLike2  =  detectAdapter(    read2 ,     qualv2     , options_adapter_S,indexAdapter,&iterations, &matches );

	
	if(indexAdapter > maxLengthForPair){
	    if(!options_mergeoverlap) //no point in continuing 
		break;
	 

	    //read1      ---------------------------------
	    //read2                  ---------------------------------

	    int startr1=indexAdapter - maxLengthForPair+max(0,lengthDiffR1_R2);
	    int startr2=max( (maxLengthForPair-indexAdapter-max(0,lengthDiffR1_R2)),0);
	    int endr1  =int(read1.size());
	    int lengthIt =endr1 -startr1;



	    logLike3      =  measureOverlap(   read1 ,     qualv1     ,  read2_rev , qualv2_rev,  
					       startr1,					      
					       startr2,
					       lengthIt,
					       &iterations,
					       &matches);

	} else{


	    //read1                         ---------------------------------
	    //read2      ---------------------------------
	    //cerr<<"A " <<lengthDiffR1_R2<<"\t"<<(indexAdapter-maxLengthForPair)<<endl;
	    int startr1=max(0,indexAdapter- size2);
	    int startr2=max( (maxLengthForPair-indexAdapter-max(0,lengthDiffR1_R2)),0);
	    int endr1  =min(int(indexAdapter),int(read1.size()));
	    int lengthIt =endr1 -startr1;
	    
	    //cerr<<"123 "<<startr1<<"\t"<<startr2<<"\t"<<lengthIt<<endl;
	    logLike3      =  measureOverlap(   read1 ,     qualv1     ,  read2_rev , qualv2_rev,  
					       startr1,
					       //is the max length minus the index of the adapter
					       //plus a potential offset for the case where r2 is shorter
					       startr2,
					       //Length of match
					       //the matching has to stop on r1 either at the size of the read or index of the adapter
					       //the matching starts at max(0,lengthDiffR1_R2), 
					       //hence the # of required comparisons is:
					       lengthIt,
					       &iterations,
					       &matches );

	    //if( indexAdapter <= (min_overlap_seqs+1))
	    //matches = min_overlap_seqs;//here, very small overlaps can have an arbitrary amount of matches
	}


	// double totalL=(logLike1+logLike2+logLike3)/iterations;
	//cout<<likeRandomMatch<<"\t"<<likeRandomMisMatch<<endl;
	double likelihoodRandom  = (   (iterations)*likeRandomMisMatch );

	double likelihoodMerge   = (logLike1+logLike2+logLike3);
	double totalL=
	    -2.0*( likelihoodRandom )
	    + 
	    2.0*(likelihoodMerge);
	
	//observing a high likelihood for short matches is easy, however, observing it over 
	//longer matches is more unlikely. Hence, we need to scale for the matches
	
	//http://en.wikipedia.org/wiki/Likelihood-ratio_test
	//null model  : this is noise
	//alternative : we observe overlap/primers

#ifdef DEBUGPR
	cerr<<"idx: "<<indexAdapter<<"\t"<<totalL<<"\t"<<( likelihoodRandom )<<"\t"<<( likelihoodMerge )<<"\t"<<iterations<<"\t"<<matches<<endl;
#endif
	
	if(logLikeTotalOverlap < totalL){
	    sndlogLikeTotalOverlap        = logLikeTotalOverlap;
	    sndlogLikeTotalOverlapIdx     = logLikeTotalOverlapIdx;
	    sndlogLikeTotalOverlapMatches = logLikeTotalOverlapMatches;

	    logLikeTotalOverlap        = totalL;
	    logLikeTotalOverlapIdx     = indexAdapter;
	    logLikeTotalOverlapMatches = matches;
	}else{ 
	    if(sndlogLikeTotalOverlap < totalL){ //not more likely than first but more likely than second
		sndlogLikeTotalOverlap        = totalL;
		sndlogLikeTotalOverlapIdx     = indexAdapter;
		sndlogLikeTotalOverlapMatches = matches;
	    }
	}
    }





    
    
    
#ifdef DEBUGPR   
    cerr<<"mergeo "<<"\t"<<logLikeTotalOverlap<<"\t"<<logLikeTotalOverlapIdx<<"\t"<<logLikeTotalOverlapMatches<<endl;
#endif
    
    cout<<logLikeTotalOverlap<<endl;
    //exit(1);
    //cout<<pow(10.0,logLikeTotalOverlap)<<endl;

    

    if( logLikeTotalOverlap         > likelihoodAdapterPR     && //sufficient likelihood
        logLikeTotalOverlapMatches >= min_overlap_seqs  ){       //sufficient # of mismatches for partial overlap (artifically to min_overlap_seqs for complete overlap)

	// cout<<-10.0*log10(logLikeTotalOverlap/sndlogLikeTotalOverlap)<<endl;
	// if( (logLikeTotalOverlap/sndlogLikeTotalOverlap)> 0.8){
	//     cout<<read1<<endl;
	//     cout<<read2<<endl;
	//     cout<<read2_rev<<endl;
	//     cout<<logLikeTotalOverlapIdx<<endl;
	//     cout<<sndlogLikeTotalOverlapIdx<<endl;

	//     exit(1);
	// }
	int indexAdapter=	logLikeTotalOverlapIdx;

	int startr1;
	int startr2;
	int endr1  ;
	int lengthIt;

	if(indexAdapter > maxLengthForPair){//partial overlap

	     startr1=indexAdapter - maxLengthForPair+max(0,lengthDiffR1_R2);
	     startr2=max( (maxLengthForPair-indexAdapter-max(0,lengthDiffR1_R2)),0);
	     endr1  =int(read1.size());
	     lengthIt =endr1 -startr1;
	     
	}else{
	    
	    startr1=max(0,indexAdapter- size2);
	    startr2=max( (maxLengthForPair-indexAdapter-max(0,lengthDiffR1_R2)),0);
	    endr1  =min(int(indexAdapter),int(read1.size()));
	    lengthIt =endr1 -startr1;
	    
	}
	    
	//logLikeTotalOverlapIdx
	string newSeq         (indexAdapter,'N');
	vector<int>   newQual (indexAdapter, 0 );
	

	//vector<int> (qualv1.begin(),qualv1.begin()+qualv1.size() - logLikePartialOverlapIdx );
	//string(read1.substr(0, read1.size()-logLikePartialOverlapIdx )+read2_rev); //new seq is part of read 1 and all of read 2
	//newQual.insert(newQual.end(),qualv2_rev.begin(),qualv2_rev.end());
	
	//coyping before the overlap
	int pos=0;
	for(pos=0;pos<startr1;pos++){
	    newSeq[pos]     =  read1[pos];
	    newQual[pos]  = qualv1[pos];	    
	}


	int posr1=startr1;
	int posr2=startr2;

	baseQual b1;
	baseQual b2;
	for(int idx=0;idx<lengthIt;idx++){

	    b1.base =  read1[posr1];
	    b1.prob = probForQual[ qualv1[posr1] ];
	    b1.qual = -1;

	    b2.base =  read2_rev[posr2];
	    b2.prob = probForQual[ qualv2_rev[posr2] ]; 
	    b2.qual = -1;

	    baseQual RT    = cons_base_prob(b1,b2);
	    newSeq[pos]    = RT.base; //since newseq starts with read1
	    newQual[pos++] = RT.qual;

	    posr1++;
	    posr2++;		
	}


	//coyping after the overlap
	for(int pos2=posr2;pos2<size2;pos2++){		
	    newSeq[pos]     =  read2_rev[pos2];
	    newQual[pos++]  = qualv2_rev[pos2];
	}


	// cerr<<indexAdapter<<"\t"<<"#"<<newSeq<<"#"<<endl;
	// exit(1);

	toReturn.code    =' ';
	toReturn.sequence=newSeq;	   
	toReturn.quality =convert_logprob_quality(newQual);	

	// cerr<<indexAdapter<<"\t"<<"#\n"<<newSeq<<"\n"<<toReturn.quality<<endl;
	// exit(1);

	return toReturn;

	
    }else{
	//exit(1);   
	toReturn.code    =' ';
	toReturn.sequence="";	   
	toReturn.quality ="";
	return toReturn;    
    }
    








//     if(max(logLikePartialOverlap,logLikeTotalOverlap) > likelihoodAdapterPR){
// 	//computing new sequence
// 	string newSeq;
// 	vector<int> newQual;

// #ifdef DEBUGTOTALOV
// 	cerr<< logLikePartialOverlap << "\t"<<logLikeTotalOverlap <<endl;
// #endif

// 	if(logLikePartialOverlap < logLikeTotalOverlap){ //total overlap
// 	    //TODO: case where read 1/2 is shorter
// #ifdef DEBUGTOTALOV
// 	    cerr<<"total overlap "<<logLikeTotalOverlapIdx<<endl;
// #endif

// 	    newSeq  = string(read1.substr(0, logLikeTotalOverlapIdx ) );

// #ifdef DEBUGTOTALOV
// 	    cerr<<"newSeq "<<newSeq<<endl;
// #endif


// 	    vector<int>   newQual = vector<int> (qualv1.begin(),qualv1.begin()+ logLikeTotalOverlapIdx );

// #ifdef DEBUGTOTALOV
// 	    cerr<<"newQual "<<newQual.size()<<endl;
// #endif

// 	    baseQual b1;
// 	    baseQual b2;
	    
// 	    for(int pos=0;pos<logLikeTotalOverlapIdx;pos++){
// 		int posr1=pos;
// 		int posr2=read1.size()-logLikeTotalOverlapIdx+pos;

// #ifdef DEBUGTOTALOV
// 	    cerr<<"pos "<<pos<<" posr1 "<<posr1<<" posr2 "<<posr2<<endl;
// 	    cerr<<"r1 "<<read1[posr1]<<endl;
// 	    cerr<<"r2 "<<read2_rev[posr2]<<endl;

// #endif

// 		b1.base =  read1[posr1];
// 		b1.prob = probForQual[ qualv1[posr1] ];
// 		b1.qual = -1;

// 		b2.base =  read2_rev[posr2];
// 		b2.prob = probForQual[ qualv2_rev[posr2] ]; 
// 		b2.qual = -1;

// 		baseQual RT  = cons_base_prob(b1,b2);
// 		newSeq[pos]  = RT.base; //since newseq starts with read1
// 		newQual[pos] = RT.qual;
// 	    }
// #ifdef DEBUGTOTALOV
// 	    cerr<<"newSeq "<<newSeq<<endl;
// #endif
	    
// 	    toReturn.code    =' ';
// 	    toReturn.sequence=newSeq;	   
// 	    toReturn.quality =convert_logprob_quality(newQual);	
// 	    return toReturn;

// 	}else{//partial overlap

// 	    //logLikeTotalOverlapIdx
// 	    newSeq  = string(read1.substr(0, read1.size()-logLikePartialOverlapIdx )+read2_rev); //new seq is part of read 1 and all of read 2
// 	    vector<int>   newQual = vector<int> (qualv1.begin(),qualv1.begin()+qualv1.size() - logLikePartialOverlapIdx );
// 	    newQual.insert(newQual.end(),qualv2_rev.begin(),qualv2_rev.end());

// 	    baseQual b1;
// 	    baseQual b2;

// #ifdef DEBUGPARTIALOV	    
// 	    cerr<<"PARTIAL RAW  "<<logLikePartialOverlapIdx<<endl<<read1<<endl<<qual1<<endl<<read2_rev<<endl<<qual2_rev<<endl<<newSeq<<"\t"<<endl;
// #endif
// 	    //TODO: case where read 1/2 is shorter

// 	    //calling consensus on the overlapping chunk
// 	    for(int pos=0;pos<logLikePartialOverlapIdx;pos++){
// 		int posr1=read1.size()-logLikePartialOverlapIdx+pos;
// 		int posr2=pos;

// #ifdef DEBUGPARTIALOV
// 		cerr<<"pos "<<pos<<" read1["<<posr1<<"] = "<<read1[posr1]<<" ("<<qualv1[posr1]<<") read2_rev["<<posr2<<"] = "<<read2_rev[posr2]<<" ("<<qualv2_rev[posr2]<<")"<<endl;
// #endif
// 		b1.base =  read1[posr1];
// 		b1.prob = probForQual[ qualv1[posr1] ];
// 		b1.qual = -1;

// 		b2.base =  read2_rev[posr2];
// 		b2.prob = probForQual[ qualv2_rev[posr2] ]; 
// 		b2.qual = -1;

// 		baseQual RT    = cons_base_prob(b1,b2);
// 		newSeq[posr1]  = RT.base; //since newseq starts with read1
// 		newQual[posr1] = RT.qual;

// #ifdef DEBUGPARTIALOV
// 		cerr<<newSeq[posr1]<<"\t"<<newQual[posr1]<<endl;
// #endif
// 	    }
	    
// #ifdef DEBUGPARTIALOV	    
// 	    cerr<<"PARTIAL TREAT  "<<newSeq<<endl<<vectorToString(newQual)<<endl<<convert_logprob_quality(newQual)<<endl;
// #endif


// 	    toReturn.code    =' ';
// 	    toReturn.sequence=newSeq;	   
// 	    toReturn.quality =convert_logprob_quality(newQual);	
// 	    return toReturn;
// 	}

//     }else{

// 	toReturn.code    =' ';
// 	toReturn.sequence="";	   
// 	toReturn.quality ="";
// 	return toReturn;    
//     }
}


// int main(){

    
//     return 0;
// }
//     if(0){
// 	cout<<edits("ACGTACGT","ACGTCCGT")<<endl;
// 	string qualString=string("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJK");
// 	vector<int> logprob=convert_quality_logprob(qualString);
// 	for(int i=0;i<logprob.size();i++)
// 	    cout<<"logprob["<<i<<"]"<<logprob[i]<<endl;
// 	cout<<convert_logprob_quality(logprob)<<endl;
// 	vector<double> prob=convert_logprob_prob(logprob);
// 	for(int i=0;i<prob.size();i++)
// 	    cout<<"prob["<<i<<"]"<<prob[i]<<endl;
// 	cout<<"TTTTTTTTTAAAAAAAAAGGGGGGGGGGGCCCCCCCCCC"<<endl;
// 	cout<<revcompl("TTTTTTTTTAAAAAAAAAGGGGGGGGGGGCCCCCCCCCC")<<endl;
// 	// while(1){
// 	// 	cout<<randomGen()<<endl;
// 	// }
// 	baseQual b1;
// 	baseQual b2;
// 	baseQual res;
    
// 	b1.base='A' ;
// 	b1.prob= 0.999;    
// 	b2.base='A' ;
// 	b2.prob= 0.99;
// 	res = cons_base_prob(b1,b2);
// 	cout<<res.base<<"\t"<<res.qual<<endl;

// 	b1.base='A' ;
// 	b1.prob= 0.99;    
// 	b2.base='C' ;
// 	b2.prob= 0.999;
// 	res = cons_base_prob(b1,b2);
// 	cout<<res.base<<"\t"<<res.qual<<endl;

// 	b1.base='A' ;
// 	b1.prob= 0.99;    
// 	b2.base='N' ;
// 	b2.prob= 0.25;
// 	res = cons_base_prob(b1,b2);
// 	cout<<res.base<<"\t"<<res.qual<<endl;

// 	b1.base='T' ;
// 	b1.prob= 0.9;    
// 	b2.base='N' ;
// 	b2.prob= 0.9;
// 	res = cons_base_prob(b1,b2);
// 	cout<<res.base<<"\t"<<res.qual<<endl;

// 	b1.base='T' ;
// 	b1.prob= 0.9;    
// 	b2.base='C' ;
// 	b2.prob= 0.9;
// 	res = cons_base_prob(b1,b2);
// 	cout<<res.base<<"\t"<<res.qual<<endl;

// 	vector<double>  prob1  (5,0.9);
// 	vector<double>  prob2  ;
// 	cout<<"identity "<<quality_ident("ACGTT",prob1,"ACGTC",prob2)<<endl;
// 	cout<<"identity "<<quality_ident("ACGTT",prob1,"ACGTI",prob2)<<endl;

// 	string read1 ="NCCGNCGCCGGGTACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
// 	string read2 ="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
// 	double qual1 [] = {0.25, 0.98004737685031118, 0.94988127663727273, 0.99205671765275716, 0.25, 0.98999999999999999, 0.99205671765275716, 0.9748811356849042, 0.99748811356849043, 0.98999999999999999, 0.9748811356849042, 0.98741074588205835, 0.99601892829446503, 0.98999999999999999, 0.98415106807538888, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
// 	double qual2 [] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
// 	int pqual1[]= {0, 17, 13, 21, 0, 20, 21, 16, 26, 20, 16, 19, 24, 20, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
// 	int pqual2[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
// 	seqQual sq1;

// 	sq1.sequence     =read1;
// 	for(int i=0;i<65;i++){
// 	    sq1.probabilities.push_back( qual1[i]);
// 	    sq1.logProbs     .push_back(pqual1[i]);
// 	}


// 	seqQual sq2;

// 	sq2.sequence     =read2;
// 	for(int i=0;i<65;i++){
// 	    sq2.probabilities.push_back( qual2[i]);
// 	    sq2.logProbs     .push_back(pqual2[i]);
// 	}
// 	seqQual cm = check_merge(sq1,sq2);
// 	cout<<cm.sequence<<endl;
// 	cout<<cm.quality<<endl;

// 	// AAGTAGAAAGGAGGAAAAACTGAAGAGATTACTTACAAATAACTTTGATCTCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTAAGTA
// 	// AGAGATCAAAGTTATTTGTAAGTAATCTCTTCAGTTTTTCCTCCTTTCTACTTGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAATTAGTGTAGA
// 	// 8=??=?=A@??=AACA:BF?EBECAEB@E@B>CAAAAC?C@DA:D?A@A?H@DDBFD>B<?EC?>>DA???C<AC<A>F<<B@D>D;BAB?>@;>
// 	// ><==<?B<>?@<B<??BA@@AEBC@B;?;CA@@C=A>?@<>?<=>>A=E>>BC=<<B=@>;@>?@C>>>=?>;?>>;=9;?>?>>BA??8>9>>8
// 	// ('', 'AAGTAGAAAGGAGGAAAAACTGAAGAGATTACTTACAAATAACTTTGATCTCT', ']]]]]]]]]]]]]]]]\\]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]')
//     }
//     //    set_adapter_sequences(options_adapter_F, options_adapter_S,"ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA",maxadapter_comp);
//     set_adapter_sequences("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG",
// 			  "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",
// 			  "ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA",maxadapter_comp);

    
//     set_options(options_trimCutoff,options_allowMissing,true);
//     set_keys("","");
//     merged res;
//     string r1="TGTGCAGAGTTTATGCAAATGACACTGGCCTTTAGTAATGTTTCTAAAAATGGAACTCTGCAAGATTTAATATCACTGTGGTTGCCCGAAGAGAGTCATG";	
//     string q1=";<2@<<<;3CB<><<<??A:9=AD=@G@9??7@?>D??F>=:D=<@;DAFB8=8=?B9@=A><<C>C@>C=?;=>===?=@89<76<>;D@=>;958=49";
//     string r2="ATGACTCTCTTCGGGCAACCACAGTGATATTAAATCTTGCAGAGTTCCATTTTTAGAAACATTACTAAAGGCCAGTGTCATTTGCATAAACTCTGCACAC";
//     string q2="564;>@9>>?:;>A<<@E=;EAD9E@D;CF>ACC@:@A@=E?@EBAC?B??:F<=8?@?@B=9>@BC@CA<<8<@@8A8>;49;:>8@<?9:6973865/";


     //res =process_PE(r1,q1,r2,q2);
    //cout<<res.code<<endl;
    //cout<<res.sequence<<"\t"<<res.quality<<endl;

//     // string r1="AAGTAGAAAGGAGGAAAAACTGAAGAGATTACTTACAAATAACTTTGATCTCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTAAGTA";
//     // string q1="8=\?\?=?=A@\?\?=AACA:BF?EBECAEB@E@B>CAAAAC?C@DA:D?A@A?H@DDBFD>B<?EC?>>DA\?\??C<AC<A>F<<B@D>D;BAB?>@;>";
//     // string r2="AGAGATCAAAGTTATTTGTAAGTAATCTCTTCAGTTTTTCCTCCTTTCTACTTGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAATTAGTGTAGA";
//     // string q2="><==<?B<>?@<B<??BA@@AEBC@B;?;CA@@C=A>?@<>?<=>>A=E>>BC=<<B=@>;@>?@C>>>=?>;?>>;=9;?>?>>BA??8>9>>8";

//     // string r1;
//     // string q1;
//     // string r2;
//     // string q2;
//     // merged res;  

//     // r1="TTTGCAATTACCATTTGTGAGGTGTTTAAATAAATGTACATGTACGAGCTTCCAGCTTCACAAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACAA";
//     // q1="<8::<<?=ED?=>AEBA8>C>=D?E?>EBAAABA=;?=BDB>>D>?AB9>C>;=:>A?7A<F8<>D><=<8F=B=D<<:<9?;87=<B8?:8;9667350";
//     // r2="TTTGTGAAGCTGGAAGCTCGTACATGTACATTTATTTAAACACCTCACAAATGGTAATTGCAAAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGA";
//     // q2="587;:??I9DC>6A>CDD:@DC=A;AF?@DAA@@B?:I@A;?9;<?F<ABE>A<:@@@?F?>DCB=<;8>E<?=>=7;B=:>=>988?:<=<<=;8976:";

//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;

//     // r1="GAAATAAGTGCTTCATTAGTAATAAAGGCTCTAGGTTATGTTACCATGTAACCTTTTGTTCAAACATGAGATCGGAAGAGCACACGTCTGAACTCCAGTC";
//     // q1="879=:>>:9<@BG>9;?==BEB?CDEF@=?=?>DEA<CB=:?E9?@8C9@=A8>A<@><>9?<C9<<E@A?;6>A=:=9=9<?=<>>8@==498737878";
//     // r2="CATGTTTGAACAAAAGGTTACATGGTAACATAATCTAGAGCCTTTATTACTAATGAAGCACTTATTTCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG";
//     // q2=":583:B7;@H>CCA@B>=A>>BDA<==H>EB;D>AC=>@??:?>??A?=9<BBC>AED<=89;A=>:??8>>8?B=99A@9:><8=8;@:;<:<9776:9";


//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;
//     // r1="GCATGTGGTCGGGGCAGCAGCCTGGGCCAGGGATGGCCGTACCGCGCCTGGCTCCACCACAGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACAAC";
//     // q1="::<==@<9=@=8=?AC@?A?A@?7A?<@BA>9G?B??@8::9AB=;D=<@>@A9==?<A>>D?;=BD8;>??>@A8>;74>;9;C>D:==85=85:,2:3";
//     // r2="TCTGTGGTGGAGCCAGGCGCGGTACGGCCATCCCTGGCCCAGGCTGCTGCCCCGACCACATGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT";
//     // q2=":87=89=?<<CADBD<D>??=@DDADD:@BA<BBC?:<=>E?;;@9<>@AA@>B>8:B7;BC<@<B<DC9@A9?=8>=98:BC878:9<=:8:<<9<853";


//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;
//     // r1="TTTTTACTTGTGAAAGTACAGCAAATTTAGTTCCTCTCATGAATGGACAGAAAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACAACCGGTATCT";
//     // q1=";9>6@?=BD6A9??C??@??;<B>:6B:BAB@7A<>9?@@:ADD7>;@=98BD?>A<D@==<<=@B;9A<B;=>8A=B===<@>EC>@7>A:6C<60804";
//     // r2="ATTTTCTGTCCATTCATGAGAGGAACTAAATTTGCTGTACTTTCACAAGTAAAAAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG";
//     // q2="6=:<8>>:@;BF<@@>B?CD@:DA@AEA9B?FD;=?@@<;@>B?BA@BE?CD@@<A=AF>><=CFA:;A;<:C@5A;B?=@>?<9=><=99?;96;=6<8";


//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;
//     // r1="GCCCTTTGACCTCTCTCTCCTAACTATGGTCACGCGGCCAAAATGTTCAGCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACAACCGGTATCTCGT";
//     // q1="6999=B>:@;@C?=AD?AA?@E@;B@<@AD?D?@=??@CA@;DD;I@=<>=>BEC<>A>??<>H;:>A89<:;?@>=<;<@B=<=?C>9>K:B<45<51:";
//     // r2="TGCTGAACATTTTGGCCGCGTGACCATAGTTAGGAGAGAGAGGTCAAAGGGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCG";
//     // q2="669?9<6>C?HBBA>=B?=>?>F7?@@IB?=?>=ADB@?>DAC<;A<DFA<;D==D9>BA;?G@<B=<=B?>>@;=@A?<>@>=8:5A8742@989::/5";


//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;

//     // r1="GCAAATCTCTGAACAGCATAACCTATGCTGTTGAGAACATGTGGTATGTCAGTCAGTCCTAACTTACTGACCGTTCTGGATAACTATGAAATTTCAATAG";
//     // q1=":7<:@A?A=?CC?>E@9=@B?@@=E?=<A?;A;G?<A9DA<>E9DE<@<=@A:>@:C;9=A<;?:A<9=@<<?@><<<<;>>9><A3:=??6;6598235";
//     // r2="GCTGTTGAAATTTCATAGTTATCCAGAACGGTCAGTAAGTTAGGACTGACTGACATACCACATGTTCTCAACAGCATAGGTTATGCTGTTCAGAGATTTG";
//     // q2="6765AAB=>FB<@@CA@>DG@=@=H??J?=D>@EGAC?9C=ACDG>A<A?;B?6?CD?>9;8>=;B<?9B=8C9;>@@8;:>;;79;;<<4:;<9>5878";

//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;
//     // r1="AATGTCAACCGATAGGATATTCTTTAATGATATCATGGCAGCATGGGCTATCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACAACCGGTATCTC";
//     // q1=";;99B=3>?=@D=C>F;?A@A@G@>BAD<A>=@@@>;:@A<@BA>?7>:>C8<CG;:@7<@C=BA@=C<=>==6<>?:6;:7B6?6<==>;9;;@8<52/";
//     // r2="ACGATAGCCCATGCTGCCATGATATCATTAAAGAATATCCTATCGGTTGACATTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGT";
//     // q2="957-7<;=B;C<D@?A<BDB<C@CA>DABEFA=G<CB<A@B=@<@GBF>;?9DAA9@>><>A>>BG<>=<8=;=<=:8<E><C=7@>8:7><::6944;2";


//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;

//     // r1="GAAGCTATCTGCACATGCAAGTAGCCACAAAGTTGTAGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCGGTATCTCGTATGCCGTCTTCTG";
//     // q1=";=68?9<C:D:6@B?AC??@=7D@9?E@C?><=ACDG9>A><B>@EC=C>?<A6?:<@A?@C@<C=?@@C=>:7A5>B:;C=;;99<==C5;<>6:6530";
//     // r2="ACTACAACTTTGTGGCTACTTGCATGTGCAGATAGCTTCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTTTC";
//     // q2="9359=<>?AA=@F:@@C@>@AA=D@?AE@C;AC@A>A>=B:DC:A?>E<=@7@@8;=CAA=;@@@D:A><CB::=@:<;;>9?:;:>76;==21>>;:<3";

//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;
//     // r1="AAAAATAAATAGTTGCTGTTTATAAGCCACTCAGGGCATAGTATTTTATTATAGCAGCTTGAGCAGACTAAGGCAATGGAGCATCAGTATCAGTCTACTA";
//     // q1="99:<?<7B=??>D@?@@D@BBGEFEB?:A@;:>;B9>>@@>C>>?CA:BC=?C>=8??8BB?=@==?<=BC=C><CB>:>D8==;8=@;;;;=65=627/";
//     // r2="TTCTGTGATTAACATATTACATCAGATTAGTAGACTGATACTGATGCTCCATTGCCTTAGTCTGCTCAGGCTGCTATAATAAAATACTATGCCCTGAGTG";
//     // q2=":=59???;B@B<<BB=BEDA<=A?<FD>A<@<<D@DCBAG>@@@=9@?>9CE@D;<BA>;@9:@><:=@B>;?9F=;@>:D998898;?>96:8<8=@82";

//     // res =process_PE(r1,q1,r2,q2);
//     // //cout<<res.code<<endl;
//     // cout<<res.sequence<<"\t"<<res.quality<<endl;

//     // r1="TGGTGCCGGAGCACACCAGCCTCAGCAACAGAGACCCCATCTCAAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACAACCGGTATCTCGTATGCCGT";
//     // q1="8<7==9;ABD><D>?<=AD?;A?GB@>=@B?BCA>@>=>?@7>CD:BB;<A=A>=EBC;=A?>A?A;=<<>:;@@@<?:=A;?<=?A?@;:B9;9;7842";
//     // r2="TTGAGATGGGGTCTCTGTTGCTGAGGCTGGTGTGCTCCGGCACCAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATC";
//     // q2="8856;>;C>>:F:>@@B<EDB@;A8?@B>>B=8@>@=>A>BD8AEA>@B?7C?@8E;=9C;<=@9EA>DF@;;=:B>?@<<<:39>@9;<@8;69;<4<3";

//     return 0;
// }
