#ifndef MergeTrimReads_h
#define MergeTrimReads_h

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <sys/time.h>


using namespace std;



static const size_t min_length = 5;
static const int offset = 33;
static const double max_prob_N = 0.25;

extern double cutoff_merge_trim;
extern size_t maxadapter_comp;
extern size_t min_overlap_seqs;
extern double cutoff_merge_seqs_early;
extern double cutoff_merge_seqs;

//  Key variables ///
extern bool handle_key;
extern bool options_allowMissing;
extern string keys0;
extern string keys1;
extern int len_key1;
extern int len_key2;
extern size_t options_trimCutoff;
extern bool options_mergeoverlap;

//Chimera options and adapter
static const char* const chimInit[]= {
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

static vector<string> adapter_chimeras (chimInit,chimInit+13);
static string options_adapter_F="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG";
static string options_adapter_S="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
// static string returnFirstToken(string * toparse,string delim);

typedef struct {
    char   base;
    int    qual;
    double prob;
} baseQual;

typedef struct {
    string         sequence;
    string         quality;
    vector<double> probabilities;
    vector<int>    logProbs;
} seqQual;

typedef struct{
    char   code;
    string sequence;
    string quality;
} merged;


void set_options(int trimcutoff=1,bool allowMissing=false,bool mergeoverlap=false);
void set_adapter_sequences(const string& forward="", const string& reverse="", const string& chimera="",int max_comp=30);
void set_keys(const string& key1, const string& key2="");

merged process_PE(string read1,string qual1,string read2,string qual2);
merged process_SR(string read1,string qual1);

/* def convert_quality_logprob(qualstring): */
/* def revcompl(seq): */
/* def cons_base_prob(base1,base2,prob1,prob2): */
/* static void process_PE(const string  read1,const string qual1,const string read2,const string qual2); */
/* static void process_SR(const string  read1,const string qual1); */

/* class MergeTrimReads{ */
/* private: */

/* public: */

/* }; */
#endif
