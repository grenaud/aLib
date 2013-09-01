#ifndef MergeTrimReads_h
#define MergeTrimReads_h

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <cfloat>
#include <math.h>
#include <sys/time.h>

#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>

#include "utils.h"

using namespace std;
using namespace BamTools;



typedef struct{
    char   code;
    string sequence;
    string quality;
} merged;

typedef struct {
    char   base;
    int    qual;
    double prob;
} baseQual;

class MergeTrimReads{
 private:
    //VARIABLES
    double maxLikelihoodRatio  ;
    
    double likelihoodChimera  ;
    double likelihoodAdapterSR ;
    
    double likelihoodAdapterPR ;
    /* static const double likelihoodAdapterPR = -1.0; */
    bool initialized ;

    const size_t min_length   ;
    const int    qualOffset    ;

    /* static const double max_prob_N = 0.25; */
    /* extern double cutoff_merge_trim; */
    size_t maxadapter_comp;

    size_t min_overlap_seqs;


    /* //  Key variables /// */
    bool handle_key;
    bool options_allowMissing;
    string keys0;
    string keys1;
    int len_key1;
    int len_key2;
    size_t options_trimCutoff;
    bool options_mergeoverlap;
    double max_prob_N ;
    /* extern size_t min_length ; */

    //Chimera options and adapter
    char*  chimInit[];/* = { */
     
    vector<string> adapter_chimeras ;
    string options_adapter_F;
    string options_adapter_S;
    


    // //  Key variables ///
    /*     bool handle_key; */
    /*     string keys0; */
    /*     string keys1; */
    /*     int len_key1 */
    /*     int len_key2; */


    //likelihood variables
    double likeMatch[64];
    double likeMismatch[64];
    
    double probForQual[64];
    double likeRandomMatch;    // 1/4
    double likeRandomMisMatch; // 3/4

    //vector<string> adapter_chimeras;

    //FUNCTIONS
    string returnFirstToken(string * toparse,string delim);
    char revComp(char c);
    string revcompl(const string seq);
    inline string convert_logprob_quality(vector<int> logScores);
    inline double randomGen();
    inline baseQual cons_base_prob(baseQual  base1,baseQual base2);


    void    setLikelihoodScores(double likelihoodChimera_,
				double likelihoodAdapterSR_,
				double likelihoodAdapterPR_);

    void set_options(int trimcutoff=1,bool allowMissing=false,bool mergeoverlap=false);
    void set_adapter_sequences(const string& forward, const string& reverse, const string& chimera);
    void set_keys(const string& key1="", 
		  const string& key2="");

    void initMerge();
    merged process_PE(string  read1,string  qual1,string read2,string qual2);
    merged process_SR(string  read1, string qual1);
    double detectChimera(const string      & read,
			 const vector<int> & qual,
			 const string      & chimeraString,
			 unsigned int        offsetChimera=0);
    double measureOverlap(const string      & read1,
			  const vector<int> & qual1,
			  const string      & read2,
			  const vector<int> & qual2,
			  const int         & maxLengthForPair,
			  unsigned int      offsetRead=0,				    
			  //double *          iterations =0 ,
			  int  *            matches=0);
    double detectAdapter(const string      & read,
			 const vector<int> & qual,
			 const string      & adapterString,
			 unsigned int        offsetRead=0,
			 int              *  matches=0);

    int edits(const string & seq1,const string & seq2);
    void sanityCheckLength(const string & seq,const string & qual);
    bool checkKeySingleEnd(string & read1,string & qual1,merged & toReturn);
    bool checkKeyPairedEnd(string & read1,string & qual1,
			   string & read2,string & qual2,
			   merged & toReturn);
    bool checkChimera(const string & read1,
		      const vector<int> & qualv1,
		      merged & toReturn, 
		      const double & logLikelihoodTotal);
    void string2NumericalQualScores(const string & qual,vector<int> & qualv);
    void computeBestLikelihoodSingle(const string      & read1,
				     const vector<int> & qualv1,
				     double & logLikelihoodTotal,
				     int &    logLikelihoodTotalIdx,
				     double & sndlogLikelihoodTotal,
				     int &    sndlogLikelihoodTotalIdx);
    void computeBestLikelihoodPairedEnd(const string &      read1,
					const vector<int> & qualv1,
					    
					const string &      read2,
					const vector<int> & qualv2,
					    
					const string &      read2_rev,
					const vector<int> & qualv2_rev,
					    
					const int & lengthRead1,
					const int & lengthRead2,
					const int & maxLengthForPair,

					double & logLikelihoodTotal,
					int    & logLikelihoodTotalIdx,
					int    & logLikelihoodTotalMatches,
					    
					double & sndlogLikelihoodTotal,
					int    & sndlogLikelihoodTotalIdx,
					int    & sndlogLikelihoodTotalMatches);


    void computeConsensusPairedEnd( const string & read1,
				    const vector<int> &   qualv1,
					
				    const string & read2_rev,
				    const vector<int> & qualv2_rev,
					
							      
				    const double & logLikelihoodTotal,
				    const int    & logLikelihoodTotalIdx,
				    const int    & logLikelihoodTotalMatches,
					
				    const double & sndlogLikelihoodTotal,
				    const int    & sndlogLikelihoodTotalIdx,
				    const int    & sndlogLikelihoodTotalMatches,
					
				    const int & maxLengthForPair,
				    merged & toReturn);


    string sortUniqueChar(string v);
    bool set_extra_flag( BamAlignment &al, int32_t f );

    const string MERGEDBAMFLAG ;
    const int32_t TRIMMEDFLAG       ;
    const int32_t MERGEDFLAG        ;
    const int32_t TRIMMEDMERGEDFLAG ;
 
    int count_all ;
    int count_fkey ;
    int count_merged ;
    int count_merged_overlap ;
    int count_trimmed ;
    int count_nothing ;
    int count_chimera ;
    vector<string> checkedTags;

 public:
    MergeTrimReads (const string& forward, const string& reverse, const string& chimera,
		    const string& key1="", const string& key2="",
		    int trimcutoff=1,bool allowMissing=false,bool mergeoverlap=false);

    MergeTrimReads(const MergeTrimReads & other);
    ~MergeTrimReads();
    MergeTrimReads & operator= (const MergeTrimReads & other);
    

    pair<BamAlignment,BamAlignment> processPair(const BamAlignment & al,const BamAlignment & al2);
    BamAlignment                    processSingle(const BamAlignment & al);

    string reportSingleLine();
    string reportMultipleLines();

};
#endif
