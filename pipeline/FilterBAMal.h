/*
 * FilterBAMal
 * Date: Nov-16-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef FilterBAMal_h
#define FilterBAMal_h

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include <math.h>
#include <stdlib.h>

#include "utils.h"

using namespace std;
using namespace BamTools;

typedef struct Report {
    unsigned int totalSeq;
    unsigned int qcBefore;
    unsigned int qcFailExp;
    unsigned int qcFailEnt;
    unsigned int qcFailClx;
    unsigned int qcLength;
} Report;

void initializeLikelihoodScores();
void setVariables(int _minLength,int _maxLength,double _cutoffLikelihood,double _cutoffAvgExpError,bool _frequency,bool _entropy,bool _compOrEntCutoff, 
		  bool _likelihoodFlag,ofstream * _likelihoodOS,bool _entropyOSFlag,ofstream  * _entropyOS,bool _frequencyOSFlag,ofstream  * _frequencyOS,bool verbose,bool _resetQC,Report * repToPrint);
void filterBAMAlign(BamAlignment * al);
double compLikelihoodSeq(BamAlignment * al);


#endif
