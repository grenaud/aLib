#ifndef RGAssign_h
#define RGAssign_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <sys/time.h>
//#include <ctype.h>
#include <string>
#include <map>
#include <set>


#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>

#include "JSON.h"

#include "PrefixTree.h"

using namespace std;
using namespace BamTools;

typedef struct {
    map<string,string>  namesMap;

    vector<string>  names;
    vector<string>  indices1;
    vector<string>  indices2;    
    bool isDoubleIndex;
    int mlindex1;
    int mlindex2;   
} indexData;

typedef struct{
    string predictedGroup;    
    bool conflict;
    bool wrong;
    bool unknown;

    double logRatioTopToSecond;            // ~~ conflict
    double logLikelihoodScore;             // ~~ unknown
    double topWrongToTopCorrect;        // ~~ wrong
    int numberOfMismatches;          //total # of mismatches
} rgAssignment;

struct tallyForRG{    
    unsigned int assigned;
    unsigned int conflict;
    unsigned int unknown;
    unsigned int wrong;    
};

class RGAssign{
 private:
    string toUpperCase(string toCheck);
    bool isValidDNA(string tocheck);
    void setFileForRatio(ofstream * streamFile);
    void setFileForRGQual(ofstream * streamFile);

    indexData intern_readIndex(string filename);
    map<string,string>  readIndexFile(string filename,int mismatchesTrie,bool _shiftByOne);
    rgAssignment assignReadGroup(string  & index1,string & index1q,string & index2,string & index2q,double rgScoreCutoff,double fracConflict,int mismatchesTrie);
    void deallocate();

    // XXX really, a global?!
    map<string,tallyForRG> namesMap; //a map name of RG to count of how many observed

    PrefixTree<string> * trieKnownString;
    string dashes;
    map<string,string> rg;

    PrefixTree<int> * indTrie1;
    PrefixTree<int> * indTrie2;
    indexData values;               // SRSLY?!


    ofstream * rgqual;
    ofstream * ratioValues;
    
    // array containing the log10 of likelihoods corresponding to quality
    // scores (these are negative numbers)
    double likeMatch[64];
    double likeMismatch[64];
    

    map<string,string> rgs;
    
    map<string,int> unknownSeq;
    map<string,int> wrongSeq;
    map<string,int> conflictSeq;

    double rgScoreCutoff  ;
    double fracConflict   ;
    double wrongness      ;
    int    mismatchesTrie ;
    int    maxErrorHits   ;

    bool shiftByOne;

    bool   printSummary;
    bool   printError;

    bool flag_ratioValues;
    bool flag_rgqual;


    double oplus( double x, double y );
    void checkRGname(string tocheck);
    inline double computeLike(const string & indexRef,const string & indexRead,const vector<int> * quals);
    inline int computeMM(const string & indexRef,const string & indexRead);
    string getCWD();
    inline string get_string_field( BamAlignment &al, const char* name );
    void updateRecord( BamAlignment &al, const rgAssignment &rg );
    inline bool containsNoNs(const string & sN); 
    void initializeKnownIndices(PrefixTree<string> * trieKnownString,string configFile);
    void check_thresholds( rgAssignment &rg ) ;


    void printUnfoundToFile(vector< pair<string,int> > * unfound,stringstream & fileError);

 public:

    RGAssign( double rgScoreCutoff_  ,
	      double fracConflict_   ,
	      double wrongness_      ,
	      int    mismatchesTrie_ ,
	      int    maxErrorHits_ ,
	      bool shiftByOne_   ,

	      bool   printSummary_ ,
	      bool   printError_ ,


	      bool flag_ratioValues_,
	      bool flag_rgqual_,

		    
	  
	      ofstream * rgqual_,
	      ofstream * ratioValues_,
  

	      string 	  index
	      ) ;

    RGAssign(const RGAssign & other);
    ~RGAssign();
    RGAssign & operator= (const RGAssign & other);

    void processPairedEndReads( BamAlignment &al, BamAlignment &al2);
    void processSingleEndReads( BamAlignment &al);


    const map<string,string> *  getRGS();
    string getSummaryString();
    string getErrorString();

};
/* static double cutoffRatioLogLike=0.5; */


#endif
