#include "MergeTrimReads.h"

bool   initialized       = false;
double cutoff_merge_trim = 0.80;
size_t maxadapter_comp=30;
size_t min_overlap_seqs=10;
double cutoff_merge_seqs_early = 0.95;
double cutoff_merge_seqs = 0.90;

//  Key variables ///
bool handle_key           = false;
bool options_allowMissing = false;
string keys0="";
string keys1="";
int len_key1=0;
int len_key2=0;
size_t options_trimCutoff =1;
bool options_mergeoverlap = false;


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
    
    cerr<<"Wrong base in revComp function "<<endl;
    exit(1);
}

static inline int edits(const string seq1,const string seq2){
    int lmin = min(seq1.length(),seq2.length());
    int lmax = max(seq1.length(),seq2.length());
    int dist = lmax-lmin;
    int  pos=0;
    while(pos<lmin){
	if (seq1[pos] != seq2[pos]) 
	    dist+=1;
	pos++;
    }
    return dist;
}

//! Function to convert quality string to log scores (int) on the phred scale (0-40)
/*!
  \param qualstring : string containing the quality scores
  \return Vector of ints containing the log scores
*/
static inline vector<int>  convert_quality_logprob(const string qualstring){
    vector<int> toReturn;
    unsigned int i=0;
    while(i<qualstring.length()){
	// cout<<(int(qualstring[i])-offset)<<endl;
	toReturn.push_back(int(qualstring[i])-offset);
	i++;
    }
    return toReturn;
}

//! Function to convert log scores to quality string
/*!
  \param qualstring : string containing the quality scores
  \return Vector of ints containing the log scores
*/
static inline string convert_logprob_quality(vector<int> logScores){
    string toReturn="";
    unsigned int i=0;
    while(i<logScores.size()){	
	toReturn+=char(max(33,min(126,logScores[i]+offset)));
	i++;
    }
    return toReturn;
}

//! Function to convert log scores to probablilities for identity calculation
/*!
  \param qualstring : string containing the quality scores
  \return Vector of ints containing the log scores
*/
static inline vector<double>  convert_logprob_prob(vector<int> logScores){
    vector<double> toReturn;
    unsigned int i=0;
    while(i<logScores.size()){
	//cout<<max(double(1.0)-pow(double(10.0),logScores[i]/double(-10.0)),max_prob_N)<<endl;
	toReturn.push_back(max(double(1.0)-pow(double(10.0),logScores[i]/double(-10.0)),max_prob_N));
	    //char(max(33,min(126,logScores[i]+offset)));
	i++;
    }
    return toReturn;

}

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

//! Function to retrieve the reverse complement
/*!
  \param qualstring : string containing the quality scores
  \return Vector of ints containing the log scores
*/
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


// specialized for qualities for second sequence
static double quality_ident(
        string::const_iterator seqA, size_t lenA, vector<double>::const_iterator qualA,
        string::const_iterator seqB, size_t lenB, vector<double>::const_iterator qualB)
{
    double ident = 0.0;
    
    //switch seq1 seq2 around?
    bool flip = lenB < lenA ;
    string::const_iterator& seq1 = flip?seqB:seqA ;
    string::const_iterator& seq2 = flip?seqA:seqB ;
    int len1 = flip?lenB:lenA;
    vector<double>::const_iterator& qual1 = flip?qualB:qualA ;
    vector<double>::const_iterator& qual2 = flip?qualA:qualB ;

    for( int pos=0; pos != len1 ; ++pos ) {
        if( seq1[pos] != seq2[pos] )
            ident += min(qual1[pos],qual2[pos]);
    }

    return len1 ? 1 - ident/len1 : 0;
}

// specialized for no qualities for second sequence
static double quality_ident(
        string::const_iterator seq1, size_t len1, vector<double>::const_iterator qual1,
        string::const_iterator seq2, size_t len2, size_t maxcomp)
{
    double ident = 0.0;
    maxcomp = min(maxcomp, len1 ) ;
    
    if( len2 >= maxcomp ) {
        for( size_t pos=0; pos != maxcomp; ++pos ) {
            if( (seq1[pos] != seq2[pos]) && (seq2[pos] != 'I'))
                ident += qual1[pos];
        }
    }else{
	cerr<<"MergeTrimReads:Quality_ident: Call with invalid arguments"<<endl;
	exit(1);	
    }

    return maxcomp ? 1 - ident/maxcomp : 0;
}

static seqQual check_merge(const seqQual& read1, 
			   const seqQual& read2){
    seqQual toReturn;
        
    unsigned int lread1 = read1.sequence.length();
    unsigned int lread2 = read2.sequence.length();

    if(lread1 != read1.probabilities.size() ||
       lread1 != read1.logProbs.size() ){
	cerr<<"MergeTrimReads: check_merge read1 cannot have variables of different lengths"<<endl;
	exit(1);	
    }

    if(lread2 != read2.probabilities.size() ||
       lread2 != read2.logProbs.size() ){
	cerr<<"MergeTrimReads: check_merge read2 cannot have variables of different lengths"<<endl;
	exit(1);	
    }
       

    if( (lread1 > 0) && (lread2 > 0) ){
	double oident = quality_ident(read2.sequence.begin(),read2.sequence.length(),read2.probabilities.begin()
                                     ,read1.sequence.begin(),read1.sequence.length(),read1.probabilities.begin());
	if(oident > cutoff_merge_trim){
	    toReturn.sequence=read1.sequence;
	    toReturn.logProbs=vector<int>(read1.logProbs);	    

	    unsigned int  pos=lread1;
	    while(pos<lread2){
		toReturn.sequence        +=   read2.sequence[pos];
		toReturn.logProbs.push_back(  read2.logProbs[pos]);
		pos++;
	    }
	    pos=0;

	    while(pos<lread2){
		if( (toReturn.sequence[pos] == 'I') || 
		    ((toReturn.sequence[pos] == 'N') && (read2.sequence[pos] != 'N') && (read2.sequence[pos] != 'I'))
		    ){
		    toReturn.sequence[pos]      = read2.sequence[pos];
		    toReturn.logProbs[pos]      = read2.logProbs[pos];
		}else{
		    if( (pos < lread1) && (read1.sequence[pos] != 'N') && (read2.sequence[pos] != 'N') && (read1.sequence[pos] != 'I') && (read2.sequence[pos] != 'I')){
			baseQual b1;
			baseQual b2;
			b1.base = read1.sequence[pos];
			b1.prob = read1.probabilities[pos];
			b1.qual = -1;

			b2.base = read2.sequence[pos];
			b2.prob = read2.probabilities[pos];
			b2.qual = -1;
			baseQual b = cons_base_prob(b1,b2);
			toReturn.sequence[pos] = b.base;
			toReturn.logProbs[pos] = b.qual;
		    }
		}
		pos++;
	    }

	}
    }


    toReturn.quality=convert_logprob_quality(toReturn.logProbs);
    return toReturn;
}



merged process_PE(string read1,string qual1,string read2,string qual2){
    merged toReturn;
    int lread1 = read1.length();
    int lread2 = read2.length();

    if( (read1.length() != qual1.length()) ||
	(read2.length() != qual2.length()) ){
	cerr<<"The reads and qualities must have equal lengths"<<endl;
	exit(1);
    }
	
    //check key, if key sequence specified
    //and matches remove key from both ends
    if( handle_key && read1.length() > 0){

	if (
	    ((read1.substr(0,len_key1) == keys0) && (read2.substr(0,len_key2) == keys1)) || //perfect match
	    (options_allowMissing && (edits(read1.substr(0,len_key1),keys0) == 1) && (read2.substr(0,len_key2) == keys1)) || //1mm in first key
	    (options_allowMissing && (read1.substr(0,len_key1) == keys0) && (edits(read2.substr(0,len_key2),keys1) == 1))
	    ){ //1mm in second key
	    read1 = read1.substr(len_key1,lread1-len_key1);
	    qual1 = qual1.substr(len_key1,lread1-len_key1);
	    read2 = read2.substr(len_key2,lread2-len_key2);
	    qual2 = qual2.substr(len_key2,lread2-len_key2);
	}else{
	    if(options_allowMissing && 
	       ( (len_key1>0?(read1.substr(0,len_key1-1) == keys0.substr(1,len_key1-1)):true) && (read2.substr(0,len_key2) == keys1)) ){
		read1 = read1.substr(len_key1-1,lread1-(len_key1-1));
		qual1 = qual1.substr(len_key1-1,lread1-(len_key1-1));
		read2 = read2.substr(len_key2,lread2-len_key2);
		qual2 = qual2.substr(len_key2,lread2-len_key2);
	    }else{
		if(options_allowMissing && 
		   (read1.substr(0,len_key1) == keys0) && (len_key2>0?(read2.substr(0,len_key2-1) == keys1.substr(1,len_key2-1)):true)){
		    read1 = read1.substr(len_key1,lread1-len_key1);
		    qual1 = qual1.substr(len_key1,lread1-len_key1);
		    read2 = read2.substr(len_key2-1,lread2-(len_key2-1));
		    qual2 = qual2.substr(len_key2-1,lread2-(len_key2-1));		    
		}else{
		    toReturn.code    ='K';
		    toReturn.sequence="";
		    toReturn.quality ="";		    
		    return toReturn;
		}

	    }
	}
    }

    vector<int>    lqual1 = convert_quality_logprob(qual1);
    vector<double> pqual1 = convert_logprob_prob(lqual1);

    //check adapter chimeras
    if( adapter_chimeras.size() > 0){
	unsigned int indexChimera=0;
	while(indexChimera<adapter_chimeras.size()){

	    if( quality_ident(read1.begin(),read1.length(),pqual1.begin()
                             ,adapter_chimeras[indexChimera].begin(),adapter_chimeras[indexChimera].length()
                             ,maxadapter_comp) > cutoff_merge_trim){
		read1 = "";
		read2 = "";
		break;
	    }else{
		if (quality_ident(read1.begin(),read1.length(),pqual1.begin()
                                 ,adapter_chimeras[indexChimera].begin()+1,adapter_chimeras[indexChimera].length()-1
                                 ,maxadapter_comp) > cutoff_merge_trim){ // consider losing first base
		    read1 = "";
		    read2 = "";
		    break;
		}
	    }
	    indexChimera++;
	}
	if(read1 == "" || read2 == "" ){
	    toReturn.code    ='D';
	    toReturn.sequence="";
	    toReturn.quality ="";		    
	    return toReturn;
	}
    }

    if( (read1.length() != qual1.length()) ||
	(read2.length() != qual2.length()) ){
	cerr<<"Internal error: the reads and qualities must have equal lengths"<<endl;
	exit(1);
    }


    lread1 = read1.length();
    lread2 = read2.length();
    // convert quality strings to probabiliities
    vector<int>    lqual2 = convert_quality_logprob(qual2);
    vector<double> pqual2 = convert_logprob_prob(lqual2);

    
    string         rread2 = revcompl(read2);
    vector<double> rqual2 (pqual2);
    reverse(rqual2.begin(),rqual2.end());

    vector<int> rlqual2 (lqual2);
    reverse(rlqual2.begin(),rlqual2.end());
    int mlength = min(lread1,lread2);
    int clength = max(lread1,lread2);

    bool have_merged = false;
    seqQual r1;
    seqQual r2;
    seqQual resultCM;
    int idxmax;

    int start=mlength;
    while(start>=0){

	double cval1 = quality_ident(read1.begin()+start, lread1-start, pqual1.begin()+start,
				     options_adapter_F.begin(), options_adapter_F.length(), maxadapter_comp);
	double cval2 = quality_ident(read2.begin()+start, lread2-start, pqual2.begin()+start,
				     options_adapter_S.begin(), options_adapter_S.length(), maxadapter_comp);
	r1.sequence        = read1.substr(0,start);
	r1.probabilities   = vector<double> (pqual1.begin(),pqual1.begin()+start);
	r1.logProbs        = vector<int>    (lqual1.begin(),lqual1.begin()+start);
	idxmax=max(0,lread2-start);
	r2.sequence        = rread2.substr(idxmax,lread2-idxmax);
	r2.probabilities   = vector<double> (rqual2.begin()+idxmax,rqual2.end());
	r2.logProbs        = vector<int>    (rlqual2.begin()+idxmax,rlqual2.end());


	resultCM=check_merge(r1,r2);
	
	if( (resultCM.sequence != "") && (max(cval1,cval2) > cutoff_merge_trim)){
	    have_merged = true; //utterly useless 
	    toReturn.code    = ' ';
	    toReturn.sequence=resultCM.sequence;
	    toReturn.quality =resultCM.quality;

	    return toReturn;
	}
	start--;
    }


    if(!have_merged){
	int start=mlength+1;
	while(start<(clength+1)){ // sequence left for asymmetric read length
	    double cval1;
	    double cval2;
	    string s1touse="";
	    string s2touse="";

	    if( lread1>start ){
		s1touse=read1.substr(start,lread1-start);
		cval1 = quality_ident(read1.begin()+start, lread1-start, pqual1.begin()+start,
				      options_adapter_F.begin(), options_adapter_F.length(), maxadapter_comp);
	    }else
		cval1=0.0;

	    if( lread2>start ){
		s2touse=read2.substr(start,lread2-start);
		cval2 = quality_ident(read2.begin()+start, lread2-start, pqual2.begin()+start,
				      options_adapter_S.begin(), options_adapter_S.length(), maxadapter_comp);
	    }else
		cval2=0.0;

	    
	    seqQual help;
	    help.sequence="";
	    int ipos=0;
	    while(ipos<max(0,start-lread2)){
		help.sequence     +="I";
		help.probabilities.push_back(0.0);
		help.logProbs     .push_back(0  );
		ipos++;
	    }
	    r1.sequence        = read1.substr(0,start);
	    r1.probabilities   = vector<double> (pqual1.begin(),pqual1.begin()+min(start,int(pqual1.size())));
	    r1.logProbs        = vector<int>    (lqual1.begin(),lqual1.begin()+min(start,int(lqual1.size())));

	    idxmax=max(0,lread2-start);
	    help.sequence    +=rread2.substr(idxmax,lread2-idxmax);
	    help.probabilities.insert(help.probabilities.end(),rqual2.begin()+idxmax,rqual2.end());
	    help.logProbs     .insert(help.logProbs.end()     ,rlqual2.begin()+idxmax,rlqual2.end());
	  

	    resultCM=check_merge(r1,help);

	    if( resultCM.sequence != "" && 
		( (max(cval1,cval2)>cutoff_merge_trim) || 
		  ( (s1touse          == s2touse) &&  
		    (s1touse.length() == 0) ) )){
		    have_merged = true; 
		    if( resultCM.sequence.length() < min_length){
			toReturn.code     ='D';
			toReturn.sequence ="";
			toReturn.quality  ="";			
		    }else{
			toReturn.code     =' ';
			toReturn.sequence = resultCM.sequence;
			toReturn.quality  = resultCM.quality ;						
		    }
		    return toReturn;
	    }
	    start++;
	}
    }


    //insert might be longer than read length. try to merge sequences with overlap
    if(!have_merged && options_mergeoverlap){
	double max_value1 =-1.0;
	int max_pos1   =-1;
	size_t start=0;
	while(start<(lread1-min_overlap_seqs)){
	
	    if (lread1-start <= lread2){
		double cval = quality_ident(read1.begin()+start, lread1-start, pqual1.begin()+start,
					    rread2.begin(), rread2.length(), rqual2.begin());
		if( cval > cutoff_merge_seqs_early){
		    max_value1=cval;
		    max_pos1  =start;
		    break;
		}else{
		    if( (max_value1 == double(-1.0)) || (cval > max_value1)){
			max_value1=cval;
			max_pos1=start;
		    }
		}
	    }
	    start++;
	}


	if (max_value1 > cutoff_merge_seqs){
	    string        new_lseq  = string(read1.substr(0,max_pos1)+rread2);
	    vector<int>   new_lqual = vector<int> (lqual1.begin(),lqual1.begin()+max_pos1);
	    new_lqual.insert(new_lqual.end(),rlqual2.begin(),rlqual2.end());
	    baseQual b1;
	    baseQual b2;

	    int pos=max_pos1;
	    while(pos<mlength){
		b1.base=new_lseq [pos];
		b1.prob=rqual2[pos-max_pos1];
		b1.qual=-1;
		b2.base=read1[pos];
		b2.prob=pqual1[pos];
		b2.qual=-1;
		baseQual RT = cons_base_prob(b1,b2);
		new_lseq[pos]=RT.base;
		new_lqual[pos]=RT.qual;
		pos++;
	    }
	    
	    toReturn.code    =' ';
	    toReturn.sequence=new_lseq;	   
	    toReturn.quality =convert_logprob_quality(new_lqual);	
	    return toReturn;
	}
    } //end if(!have_merged && options_mergeoverlap){


    toReturn.code    =' ';
    toReturn.sequence="";	   
    toReturn.quality ="";	
    return toReturn;
}





merged process_SR(string read1,string qual1){
    merged toReturn;
    vector<double> emptyQual2;

    size_t lread1 = read1.length();
    if( (read1.length() != qual1.length()) ){
	cerr<<"MergeTrimReads: The reads and qualities must have equal lengths"<<endl;
	exit(1);
    }

    if( handle_key ){
	if ((read1.substr(0,len_key1) == keys0)  || //perfect match
	     (options_allowMissing && (edits(read1.substr(0,len_key1),keys0) == 1))  //1mm in first key	     
	    ){ //1mm in second key
	    read1 = read1.substr(len_key1,lread1-len_key1);
	    qual1 = qual1.substr(len_key1,lread1-len_key1);
	}else{
	    if(options_allowMissing && 
	       (read1.substr(0,len_key1-1) == keys0.substr(1,len_key1-1))  ){
		read1 = read1.substr(len_key1-1,lread1-(len_key1-1));
		qual1 = qual1.substr(len_key1-1,lread1-(len_key1-1));
	    }else{
		toReturn.code    ='K';
		toReturn.sequence="";
		toReturn.quality ="";		    
		return toReturn;	       
	    }
	}
    }

    if( (read1.length() != qual1.length()) ){
	cerr<<"MergeTrimReads : The reads and qualities must have equal lengths"<<endl;
	exit(1);
    }

    lread1 = read1.length();
    // convert quality strings to probabiliities
    vector<int>    lqual1 = convert_quality_logprob(qual1);
    vector<double> pqual1 = convert_logprob_prob(lqual1);
    
    //check adapter chimeras
    if( adapter_chimeras.size() > 0){
	vector<int>    lqual1 = convert_quality_logprob(qual1);
	vector<double> pqual1 = convert_logprob_prob(lqual1);


	unsigned int indexChimera=0;
	while(indexChimera<adapter_chimeras.size()){
	   
	    if( quality_ident(read1.begin(),read1.length(),pqual1.begin(),
                              adapter_chimeras[indexChimera].begin(), adapter_chimeras[indexChimera].length(),
                              maxadapter_comp) > cutoff_merge_trim){
		read1 = "";
		break;
	    }else{
		if (quality_ident(read1.begin(),read1.length(),pqual1.begin(),
                                  adapter_chimeras[indexChimera].begin()+1, adapter_chimeras[indexChimera].length()-1, 
                                  maxadapter_comp-1) > cutoff_merge_trim){ // consider losing first base
		    read1 = "";
		    break;
		}
	    }
	    indexChimera++;
	}
	if( read1 == ""  ){
	    toReturn.code    ='D';
	    toReturn.sequence="";
	    toReturn.quality ="";		    
	    return toReturn;
	}
    }

    int adapter_pos = -2;
    double max_value =-1;
    int max_pos   =-1;
    size_t start=0;
    while(start<(lread1)){
	double cval = quality_ident(read1.begin()+start, lread1-start, pqual1.begin()+start,
				    options_adapter_F.begin(), options_adapter_F.length(), maxadapter_comp);

	if( cval > cutoff_merge_seqs_early){
	    max_value=cval;
	    max_pos  =start;
	    break;
	}else{
	    if(((cval > cutoff_merge_trim) && (cval > max_value) && (lread1-start > min_length)) || 
	       ((lread1-start <= min_length) && (max_value == -1))
	       ){
		max_value =cval;
		max_pos   =start;
	    }
	}	
	start++;
    }

    if( (max_value > cutoff_merge_trim) && 
	((lread1-max_pos) >= options_trimCutoff) ){
	adapter_pos = max_pos;
    
        read1 = read1.substr(0,adapter_pos);
	if( read1.length() < min_length){
	    toReturn.code    ='D';
	    toReturn.sequence="";	   
	    toReturn.quality ="";
	    return toReturn;	
	}else{
	    qual1 = qual1.substr(0,adapter_pos);
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





void set_options(int trimcutoff,bool allowMissing,bool mergeoverlap){
    options_trimCutoff = trimcutoff;
    options_allowMissing = allowMissing;
    options_mergeoverlap = mergeoverlap;
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




// int main(){
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


//     res =process_PE(r1,q1,r2,q2);
//     //cout<<res.code<<endl;
//     cout<<res.sequence<<"\t"<<res.quality<<endl;

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
