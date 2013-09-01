#include "RGAssign.h"
// vim:ts=8

// #define DEBUG
// #define DEBUG2
//#define DEBUG3


// void setFileForRatio(ofstream * streamFile){
//     flag_ratioValues=true;
//     ratioValues=streamFile;
// }

// void setFileForRGQual(ofstream * streamFile){
//     flag_rgqual=true;
//     rgqual=streamFile;
// }



// compare by likelihood (second component).  since likelihoods are
// negative, the meaning is flipped
struct comparePair 
{
    bool operator() (pair<int,double> i,pair<int,double> j) { 
        return (i.second>j.second); 
    }
} ;


struct compareNameRG {
    bool operator() (pair<string,int> i,pair<string,int> j) {
        return (i.second>j.second);
    }
};

struct compareNameTally {
    bool operator() (pair<string,tallyForRG> i,pair<string,tallyForRG> j) { 
        return ( (i.second.assigned+i.second.unknown+i.second.conflict+i.second.wrong)
		 >
		 (j.second.assigned+j.second.unknown+j.second.conflict+j.second.wrong) );
    }
};


// Returns log10( pow(10,x)+pow(10,y) ), but does so without causing
// overflow or loss of precision.
double RGAssign::oplus( double x, double y )
{
    return x > y 
        ? x + log1p( pow( 10, y-x ) ) / log(10)
        : y + log1p( pow( 10, x-y ) ) / log(10) ;
}



string RGAssign::toUpperCase(string toCheck){
    string toReturn=toCheck;
    transform(toReturn.begin(), toReturn.end(), toReturn.begin(), ::toupper);
    return toReturn;
}

bool RGAssign::isValidDNA(string tocheck){
    for(unsigned int  i=0;i<tocheck.size();i++){
	if(tocheck[i] != 'A' &&
	   tocheck[i] != 'C' &&
	   tocheck[i] != 'G' &&
	   tocheck[i] != 'T' )
	    return false;	  
    }
    return true;
}

void RGAssign::checkRGname(string tocheck){
    if(tocheck == "unknown") {  cerr<<"Error: Cannot use reserved name \"unknown\" as read group name"<<endl;  exit(1);}
    if(tocheck == "conflict"){  cerr<<"Error: Cannot use reserved name \"conflict\" as read group name"<<endl; exit(1);}
    if(tocheck == "wrong")   {  cerr<<"Error: Cannot use reserved name \"wrong\" as read group name"<<endl;    exit(1);}
    //fine
}

indexData RGAssign::intern_readIndex(string filename){
    string line;
    ifstream myFile;

    indexData toReturn;
    toReturn.mlindex1=0;
    toReturn.mlindex2=0;

    bool isFirstLine  =true;


    //initialize the values for the likelihood of matches or mismatches 
    for(int i=0;i<64;i++){
	if(i == 0)
	    likeMatch[i]    = -3.0; // this is vrong, hope it's never accessed
	else
	    likeMatch[i]    = log1p( -pow(10.0,i/-10.0) )/log(10);
	
	likeMismatch[i]     = i/-10.0;	
#ifdef DEBUG2
	cout<<"qual = "<<i<<endl;
	cout<<likeMatch[i]<<endl;
	cout<<likeMismatch[i]<<endl;
#endif
    }

    //reading the files

    myFile.open(filename.c_str(), ios::in);
    if (myFile.is_open()){

	while ( getline (myFile,line)){
	    line+=' ';

	    if(isFirstLine){
		if(line[0] == '#'){
		    unsigned int i=0;
		    int numberOfFields=0;
		    bool inWS=true;
		    while(i<line.length()){			
			if( isspace(line[i])){			    
			    inWS=true;
			}else{
			    if(inWS){
				numberOfFields++;
			    }
			    inWS=false;			    
			}
			i++;
		    }
		    
		    if(numberOfFields==2){ 
			toReturn.isDoubleIndex=false; 
		    }else{
			if(numberOfFields==3 || numberOfFields==5){
			    toReturn.isDoubleIndex=true; 
			}else{
			    cerr << "Must have 2, 3 or 5 fields"<<endl;
			    exit(1);
			}			
		    }

		}else{
		    cerr << "First line must begin with #"<<endl;
		    exit(1);
		}
		isFirstLine=false;
	    }else{
		int i=0;
		int fieldIndex=0;
		bool inWS=false;
		int lastOneNW=0;
		string foundName;
		while(i<int(line.length())){		
		    
		    if( isspace(line[i]) && i==0){
			cerr<<line<<endl;
			cerr << "First character cannot be a space"<<endl;
			exit(1);
		    }			    
		    if( isspace(line[i]) ){			    
			if(!inWS){ //found a field

			    //first field, first index
			    if(fieldIndex==0){
				toReturn.indices1.push_back(toUpperCase(line.substr(lastOneNW,i-lastOneNW)));

				if(toReturn.mlindex1 < (i-lastOneNW)){
				    toReturn.mlindex1 =(i-lastOneNW);
				}

			    }else{
				//second field, either name of single ind or second index
				if(fieldIndex==1){
				    if(toReturn.isDoubleIndex){
					toReturn.indices2.push_back(toUpperCase(line.substr(lastOneNW,i-lastOneNW)));
					if(toReturn.mlindex2 < (i-lastOneNW)){
					    toReturn.mlindex2 =(i-lastOneNW);
					}
				    }else{
					foundName=line.substr(lastOneNW,i-lastOneNW);
					//duplicated names ?					
					if(toReturn.namesMap.find(  foundName  ) !=  toReturn.namesMap.end()){
					    cerr<<"Warning: The sequence name is duplicated "<<foundName<<endl;
					    //exit(1);
					}else{
					    toReturn.namesMap[ foundName ] = ""; 
					}

					toReturn.names.push_back( foundName );

				    }
                                }else if(fieldIndex==2){
                                    //sequence name when two indices
                                    if(toReturn.isDoubleIndex){
                                        //duplicated names
                                        foundName=line.substr(lastOneNW,i-lastOneNW);
					
                                        if(toReturn.namesMap.find(  foundName  ) !=  toReturn.namesMap.end()){
                                            cerr<<"Warning: The sequence name is duplicated "<<foundName<<endl;
                                            //exit(1);
                                        }else{
                                            toReturn.namesMap[ foundName ] = ""; 
                                        }

                                        toReturn.names.push_back( foundName );
                                    }else{
                                        //it's a comment for single index
                                        toReturn.namesMap[ foundName ] +=  line.substr(lastOneNW,i-lastOneNW);
                                        // cerr<<"Single index file cannot have 3 fields"<<endl;
                                        // exit(1);
                                    }
                                }else{
                                    //it's a comment again
                                    toReturn.namesMap[ foundName ] +=  line.substr(lastOneNW,i-lastOneNW);
				}

			    }
			    fieldIndex++;
			}
			inWS=true;		    
			
		    }else{
			if(inWS)
			    lastOneNW=i;
			inWS=false;			    
		    }
		    i++;		
		} //ending while(i<line.length()){		
	    }  // ending else firstline

	}  // ending while myFile.good() ){
	myFile.close();
    }else{ 
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }

    //checking for size
    // cout<<toReturn.indices1.size()<<endl;
    // cout<<toReturn.indices2.size()<<endl;
    // cout<<toReturn.names.size()<<endl;
    if(toReturn.isDoubleIndex)
	if((toReturn.indices1.size() != toReturn.indices2.size()) ){
	    cerr << "Size of the fields inconsistent "<<filename<<endl;
	    exit(1);
	}


    if(toReturn.indices1.size() != toReturn.names.size() ){
	cerr << "Size of the fields inconsistent "<<filename<<endl;
	exit(1);
    }


    //checking for valid dna    
    for(unsigned int i=0;i<toReturn.indices1.size();i++){
	if(!isValidDNA(toReturn.indices1[i])){
	    cerr << "Index " << toReturn.indices1[i] <<" is not a valid DNA sequence"<<endl;
	    exit(1);
	}
	if(toReturn.isDoubleIndex)
	    if(!isValidDNA(toReturn.indices2[i])){
		cerr << "Index " << toReturn.indices2[i] <<" is not a valid DNA sequence"<<endl;
		exit(1);
	    }
    }
    return toReturn;
}


map<string,string>  RGAssign::readIndexFile(string filename,int mismatchesTrie,bool _shiftByOne){
    values = intern_readIndex(filename);
    shiftByOne=_shiftByOne;
    
    indTrie1 = new PrefixTree<int>  ();    
    if(values.isDoubleIndex){
	indTrie2 = new PrefixTree<int>  ();
    }
    // PrefixTree<int> pairsOfIndex;
    for(int i=0;i<int(values.names.size());i++){

	indTrie1->insertIntoTree( values.indices1[i].c_str() , i);
	if(values.isDoubleIndex){
	    indTrie2->insertIntoTree( values.indices2[i].c_str() , i);
	    //pairsOfIndex->insertIntoTree( values.indices1[i].c_str()+values.indices2[i].c_str() , i);	    
	}	
    }
    
    
    
    //detect conflicts ?
    cerr<<"Conflicts for index1:"<<endl;
    for(int i=0;i<int(values.names.size());i++){
	vector< int > matchesind ;
	indTrie1->searchMismatch(values.indices1[i].c_str(),&matchesind,mismatchesTrie);
	string toPrint="";

	for(unsigned int j=0;j<matchesind.size();j++){
	    if(matchesind[j] != i)
		toPrint+=values.names[ matchesind[j]  ]+" ";
	}

	if(toPrint.length() > 0){
	    if(values.names[i] != toPrint)
		cerr<<values.indices1[i]<<" from #"<<values.names[i]<<"# causes a conflict with #"<<toPrint<<"#"<<endl;
	}
    }

    if(values.isDoubleIndex){
	cerr<<"Conflicts for index2:"<<endl;
	for(int i=0;i<int(values.names.size());i++){

	    vector< int > matchesind;
	    indTrie2->searchMismatch(values.indices2[i].c_str(),&matchesind,mismatchesTrie);

	    string toPrint="";

	    for(unsigned int j=0;j<matchesind.size();j++){
		if(matchesind[j] != i)
		    toPrint+=values.names[ matchesind[j]  ]+" ";
	    }

	    if(toPrint.length() > 0){
		cerr<<values.indices1[i]<<" from "<<values.names[i]<<" causes a conflict with "<<toPrint<<endl;
	    }
	}

	cerr<<"Conflicts for pairs:"<<endl;
	for(int i=0;i<int(values.names.size());i++){
	    vector< int > matchesind1;
	    vector< int > matchesind2;

	    indTrie1->searchMismatch(values.indices1[i].c_str(),&matchesind1,mismatchesTrie);	    
	    indTrie2->searchMismatch(values.indices2[i].c_str(),&matchesind2,mismatchesTrie);

	    string toPrint="";

	    for(unsigned int j=0;j<matchesind1.size();j++){
		//found in second as well
		if( find(matchesind2.begin(),
			 matchesind2.end(),
			 matchesind1[j]) != matchesind2.end() 
		    &&
		    matchesind1[j] != i)
		    toPrint+=values.names[ matchesind1[j]  ]+" ";
	    }

	    if(toPrint.length() > 0){
		cerr<<values.indices1[i]+"#"+values.indices2[i]<<" from "<<values.names[i]<<" causes a conflict with "<<toPrint<<endl;
	    }
	}
	
    }

    return values.namesMap ;
}

// likelihood for an index read given the correct index sequence
// this is effectively the sum of qualities of mismatching bases; a
// negative number is our world.
inline double RGAssign::computeLike(const string & indexRef,const string & indexRead,const vector<int> * quals){
#ifdef DEBUG2
    cerr<<"computeLike() "<<indexRef<<"\t"<<indexRead<<endl;
#endif

    double toReturn=0.0;
    for(unsigned int i=0;i<min(indexRef.length(),indexRead.length());i++){

	
	if( indexRef[i] == indexRead[i] ){    
	    toReturn+=likeMatch[    (*quals)[i] ]; 
#ifdef DEBUG2
	    cerr<<"i="<<i<<"\tmatch="<<indexRef[i]<<"\t"<<indexRead[i]<<"\t"<<(*quals)[i]<<"\t"<<likeMatch[    (*quals)[i] ]<<endl;
#endif

	}else{
	    toReturn+=likeMismatch[ (*quals)[i] ]; 
#ifdef DEBUG2
	    cerr<<"i="<<i<<"\tmismatch="<<indexRef[i]<<"\t"<<indexRead[i]<<"\t"<<(*quals)[i]<<"\t"<<likeMismatch[    (*quals)[i] ]<<endl;
#endif

	}

#ifdef DEBUG
	cerr<<"i="<<i<<"\t"<<indexRef[i]<<"\t"<<indexRead[i]<<"\t"<<toReturn<<endl;
#endif

	
    }

    return toReturn;
}



//computes mismatches between index ref and index from the read
inline int RGAssign::computeMM(const string & indexRef,const string & indexRead){
    int toReturn=0;
    for(unsigned int i=0;i<min(indexRef.length(),indexRead.length());i++){
	
	if( indexRef[i] != indexRead[i] )   
	    toReturn++;
	    
    }

    return toReturn;
}

rgAssignment RGAssign::assignReadGroup(string &index1, 
				       string &index1q, 
				       string &index2,
				       string &index2q, 
				       double rgScoreCutoff,
				       double fracConflict,
				       int mismatchesTrie){
    //BEGIN DEBUG
    // cout<<"DEBUG"<<endl;
    // vector< int > * test1=new vector<int>();
    // indTrie1->searchMismatch("NTGNNNN",test1,2);
    // //    vector< matches<int> > * test2=indTrie2->searchForWordMismatch("CCGGTAC",0);
    // cout<<test1->size()<<endl;
    // //    cout<<test2->size()<<endl;

    // for(unsigned int j=0;j<test1->size();j++){
    // 	// list< int  >::const_iterator iter;		    
    // 	// for (iter  = (*test1)[j].listOfDeflinesOfMatches->begin(); 
    // 	//      iter != (*test1)[j].listOfDeflinesOfMatches->end(); 
    // 	//      iter++){
    // 	//     cout<<"test "<<*iter<<"\t"<<values.indices1[*iter]<<"\t"<<values.names[*iter]<<endl;
    // 	// }
    // 	cout<<"test "<< (*test1)[j] <<"\t"<<values.indices1[  (*test1)[j] ]<<"\t"<<values.names[ (*test1)[j] ]<<endl;
    // }
    // delete test1;
    //exit(1);
    //END DEBUG

    rgAssignment toReturn;
    int qualOffset=33;

    vector<int> quals1;
    vector<int> quals2;

    //Find matches using the prefix trie

    vector< int > matchesind1, matchesind2;
    indTrie1->searchMismatch(      index1.c_str(),&matchesind1,mismatchesTrie);

    
    if(shiftByOne){       
	indTrie1->searchMismatch(     ("N"+index1.substr(0, index1.size() -1)     ).c_str(),&matchesind1,mismatchesTrie);
	indTrie1->searchMismatch(         (index1.substr(1, index1.size() -1)+"N" ).c_str(),&matchesind1,mismatchesTrie);
    }

    for(unsigned int i=0;i<index1q.length();i++){
	char temp_3 = char(index1q[i]);   
	int  temp_2 = (int(temp_3)-qualOffset);
	int  temp_  =  max(temp_2 ,2);
	quals1.push_back( temp_  )   ; //since qual scores less than 2 do not make sense
    }

    if(!index2.empty()){
	indTrie2->searchMismatch(  index2.c_str(),&matchesind2,mismatchesTrie);

	if(shiftByOne){       
	    indTrie2->searchMismatch( ("N"+index2.substr(0, index2.size() -1)     ).c_str(),&matchesind2,mismatchesTrie);
	    indTrie2->searchMismatch(     (index2.substr(1, index2.size() -1)+"N" ).c_str(),&matchesind2,mismatchesTrie);
	}

	for(unsigned int i=0;i<index2q.length();i++){
	    char temp_3 = char(index2q[i]);
	    int  temp_2 = (int(temp_3)-qualOffset);
	    int  temp_  = max( temp_2,2);

	    quals2.push_back( temp_ )   ;	//since qual scores less than 2 do not make sense
	}
    }


    set<int> foundIndices; //set of indices found

    double like1;
    double like2;

    vector< pair<int,double>  > sortedLikelihoodAll; //likelihood for both indices
    vector< pair<int,double>  > sortedLikelihood1;   //likelihood for index 1
    vector< pair<int,double>  > sortedLikelihood2;   //likelihood for index 1

    //putting all the index # into the set
    for(unsigned int j=0;j<matchesind1.size();j++)
        foundIndices.insert ( matchesind1[j]  );   

    for(unsigned int j=0;j<matchesind2.size();j++)
        foundIndices.insert ( matchesind2[j]  );   


    //for every # in the set, compute likelihood    
    for (set<int>::iterator sit=foundIndices.begin(); sit!=foundIndices.end(); sit++){

	like1      =  computeLike(index1,values.indices1[ *sit ],&quals1);
	
	if(shiftByOne){ //check if shifting by one improves the likelihood
	    like1      =  max(like1,
			      max(
				  computeLike( ("N"+index1.substr(0, index1.size() -1)     ), values.indices1[ *sit ],&quals1),
				  computeLike( (    index1.substr(1, index1.size() -1)+"N" ), values.indices1[ *sit ],&quals1)
				  )			      
			      );
	}

	sortedLikelihood1.push_back(     make_pair (*sit,like1) );

	if(!index2.empty()){

	    like2  =  computeLike(index2,values.indices2[ *sit ],&quals2);

	    if(shiftByOne){ //check if shifting by one improves the likelihood
		like2      =  max(like2,
				  max(
				      computeLike( ("N"+index2.substr(0, index2.size() -1)     ), values.indices2[ *sit ], &quals2),
				      computeLike( (    index2.substr(1, index2.size() -1)+"N" ), values.indices2[ *sit ], &quals2)
				      )
				  );
	    }

	    sortedLikelihood2.push_back( make_pair (*sit,like2) );
	}else{
	    like2=0.0;
	}

	sortedLikelihoodAll.push_back( make_pair (*sit,like1+like2) ) ;      
    }

    //if nothing was found, unknown
    if(foundIndices.empty() ){
	toReturn.predictedGroup.clear();
	return toReturn;
    }
    //At this point, we will return a RG 

    //sorting by likelihood
    sort (sortedLikelihood1.begin(),   sortedLikelihood1.end(),   comparePair()); 
    sort (sortedLikelihood2.begin(),   sortedLikelihood2.end(),   comparePair()); 
    sort (sortedLikelihoodAll.begin(), sortedLikelihoodAll.end(), comparePair()); 

#ifdef DEBUG
    cerr<<endl;
    cerr<<"first:"<<endl;
    for(unsigned int j=0;j<sortedLikelihood1.size();j++){
    	cerr<< values.names[ sortedLikelihood1[j].first ]<<"\t"<< sortedLikelihood1[j].second<<endl;
    }
    cerr<<"second:"<<endl;
    for(unsigned int j=0;j<sortedLikelihood2.size();j++){
    	cerr<< values.names[ sortedLikelihood2[j].first ]<<"\t"<< sortedLikelihood2[j].second<<endl;
    }
    cerr<<"all:"<<endl;
    for(unsigned int j=0;j<sortedLikelihoodAll.size();j++){
    	cerr<< values.names[ sortedLikelihoodAll[j].first ]<<"\t"<< sortedLikelihoodAll[j].second<<"\t"<<pow(10.0,sortedLikelihoodAll[j].second)<<endl;
    }
    cerr<<endl;
    exit(1);
#endif

    // DETECT WRONGS
    // Look for a wrong index.  Suppose the top indices do not match up,
    // then those give the likelihood for being wrong, which we compare
    // to that of the top correct pair.
    //
    // Suppose they do, then we have to find the highest scoring wrong
    // pair.  That's either the top first index with the runner-up
    // second index, or vice versa.  Could actually be both to some
    // extent, so we add them.

    if( (sortedLikelihood1.size()>1) && 
	(sortedLikelihood2.size()>1)   ){
        if(sortedLikelihood1[0].first != sortedLikelihood2[0].first) {
            // mismatch, we found the wrong pair
            toReturn.topWrongToTopCorrect = sortedLikelihoodAll[0].second
                                          - sortedLikelihood1[0].second
                                          - sortedLikelihood2[0].second ;
        }else {
            // we compare one correct pair to two potentially wrong
            // ones; add 0.3 to make it fair
            toReturn.topWrongToTopCorrect = sortedLikelihoodAll[0].second + 0.3 -
                oplus( sortedLikelihood1[0].second + sortedLikelihood2[1].second 
                     , sortedLikelihood1[1].second + sortedLikelihood2[0].second ) ;
        }
    }else{
	toReturn.topWrongToTopCorrect = 0x7fffffff ; // +infinity for practical purposes
    }

    // DETECT CONFLICTS
    // Checking likelihood of inferior hits; if the ratio is too low,
    // it's a conflict.  Checking the second best only is already a
    // useable approximation, adding all is more appropriate (and a bit
    // more expensive).
    double probRG    = sortedLikelihoodAll[0].second;
    toReturn.predictedGroup = values.names[ sortedLikelihoodAll[0].first ];
    toReturn.logLikelihoodScore = probRG;

    if(sortedLikelihoodAll.size() > 1){
	double probRG2nd    = sortedLikelihoodAll[1].second;
        for( size_t i = 2 ; i != sortedLikelihoodAll.size() ; ++i )
            probRG2nd = oplus( probRG2nd, sortedLikelihoodAll[i].second ) ; //oplus= log10( pow(10,x)+pow(10,y) )

	toReturn.logRatioTopToSecond = probRG2nd - oplus(probRG,probRG2nd) ;
	double temporaryD=-10.0*toReturn.logRatioTopToSecond;
	if(flag_ratioValues ) // && toReturn.logRatioTopToSecond > -5) //to avoid very small values
	    ratioValues->write( (char *)&(temporaryD), sizeof(toReturn.logRatioTopToSecond));
    }else{
	toReturn.logRatioTopToSecond = 1;
    }

    if(flag_rgqual){
	double temporaryD=-10.0*toReturn.logLikelihoodScore;
	rgqual->write( (char *)&(temporaryD), sizeof(toReturn.logLikelihoodScore));
    }


#ifdef DEBUG3    

    toReturn.numberOfMismatches      = computeMM(index1,values.indices1[ sortedLikelihoodAll[0].first ]);
    if(!index2.empty())
	toReturn.numberOfMismatches += computeMM(index2,values.indices2[ sortedLikelihoodAll[0].first ]);
    
    //cerr<<endl;
    // if(toReturn.numberOfMismatches!=0){
    // 	cerr<<toReturn.logLikelihoodScore<<endl;
    // cerr<<toReturn.topWrongToTopCorrect<<endl;
    // cerr<<toReturn.logRatioTopToSecond<<endl;
    //cerr<<"MM"<<toReturn.numberOfMismatches<<endl;
    //}
    //cerr<<toReturn.numberOfMismatches<<"\t"<<toReturn.logLikelihoodScore<<endl;
#endif
    
    return toReturn;
}





string RGAssign::getCWD(){
   char temp[1000];
   return ( getcwd(temp, 1000) ? string( temp ) : string("") );
}


inline string RGAssign::get_string_field( BamAlignment &al, const char* name ) 
{
    if(al.HasTag(name)) {
        char ttype;
        if( !al.GetTagType(name,ttype) ) {
            cerr << "Unable to get tag ("<<name<<") type for read  " << al.Name<<endl;
            exit(1);
        }
        if( ttype=='Z' || ttype=='H' ) {
            string tagInfo;
            if(!al.GetTag(name,tagInfo)){
                cerr << "Unable to edit " << name << " tag" << endl;
                exit(1);     
            }
            return tagInfo;
        } else if( ttype=='A' ) {
            int val;
            if( !al.GetTag(name,val) ){
                cerr << "Unable to edit " << name << " tag" << endl;
                exit(1);     
            }
            return string(1,val) ;
        }
    }
    return string() ;
}

static void getIndices( const BamAlignment &al,string & index1,string & index1Q,string & index2,string & index2Q){
    if(!al.GetTag("XI",index1) ){ 	
	cerr << "Cannot retrieve XI field  "<<al.Name << endl;
	exit(1); 
    }
    if(!al.GetTag("YI",index1Q)){ 	
	cerr << "Cannot retrieve YI field  "<<al.Name << endl;
	exit(1); 
    }

    if(!al.GetTag("XJ",index2)) {
        index2 ="";
        index2Q="";
    }
    else if(!al.GetTag("YJ",index2Q)) {
	    cerr << "Cannot retrieve YJ field  "<<al.Name << endl;
	    exit(1); 
    }
}


void RGAssign::updateRecord( BamAlignment &al, const rgAssignment &rg )
{
    // get old ZQ field, remove "ICW"
    string zq = get_string_field(al, "ZQ");
    string::iterator p = zq.begin(), q = zq.begin(), e = zq.end() ;
    while( p != e ) {
        if( *p != 'I' && *p != 'C' && *p != 'W' ) 
        {
            *q = *p ;
            ++q ;
        }
        ++p ;
    }
    zq.erase( q, e ) ;

    string predictedGroup=rg.predictedGroup;
    //will overwrite the RG
    bool assigned=true;
    if( rg.predictedGroup.empty() ) {
	predictedGroup="unknown";
    }
    bool incrInTally=false;
    if( rg.conflict ){ zq += 'C' ; assigned=false; if(!incrInTally){ namesMap[predictedGroup].conflict++; incrInTally=true;}  }
    if( rg.wrong    ){ zq += 'W' ; assigned=false; if(!incrInTally){ namesMap[predictedGroup].wrong++;    incrInTally=true;}  }
    if( rg.unknown  ){ zq += 'I' ; assigned=false; if(!incrInTally){ namesMap[predictedGroup].unknown++;  incrInTally=true;}  }
    if(assigned){
	namesMap[predictedGroup].assigned++;
    }


    if( rg.predictedGroup.empty() ) {
        al.RemoveTag("RG");
        al.RemoveTag("Z0");
        al.RemoveTag("Z1");
        al.RemoveTag("Z2");
        zq += 'I';
	//namesMap[ "unknown" ] ++;
	//        al.EditTag("RG","Z","unknown");
    } else    {
	//namesMap[ predictedGroup ] ++;
        al.EditTag("RG","Z",rg.predictedGroup);

        al.EditTag("Z0","i",(int)round(-10 * rg.logLikelihoodScore));

        if( rg.logRatioTopToSecond <= 0 )
            al.EditTag("Z1","i",(int)round(-10 * rg.logRatioTopToSecond));
        else
            al.RemoveTag("Z1") ;

        if( (rg.topWrongToTopCorrect) <= 0) 
            al.EditTag("Z2","i",(int)round(-10 * rg.topWrongToTopCorrect));
        else
            al.RemoveTag("Z2") ;

    }

    // store new ZQ field and set FailedQC flag if it isn't empty
    al.SetIsFailedQC( !zq.empty() ) ;
    if( zq.empty() ) 
	al.RemoveTag("ZQ") ;
    else 
	al.EditTag( "ZQ", "Z", zq ) ;
}

inline bool RGAssign::containsNoNs(const string & sN){
    return (sN.find("N") == string::npos);
}

void RGAssign::initializeKnownIndices(PrefixTree<string> * trieKnownString,string configFile){
    //p7 300
    string line;
    ifstream myFile;
    string content="";
    myFile.open(configFile.c_str(), ios::in);

    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    content+=line;
	}
	myFile.close();
    }else{
	cerr << "Unable to open config file "<<configFile<<endl;
	exit(1);
    }


    JSONValue *value = JSON::Parse(content.c_str());
    if (value == NULL){
	cerr<<"Failed to parse JSON file"<<endl;
	exit(1);
    }

    JSONObject root;
    root = value->AsObject();
    if(root.find(L"indices") == root.end()){
	cerr<<"Failed to parse JSON file, needs a indices field"<<endl;
	exit(1);
    }
    

    JSONValue *jsonIndices    = root.at(L"indices");
    JSONObject jsonIndicesObj = jsonIndices->AsObject();


    //p7 indices
    if(jsonIndicesObj.find(L"p7indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p7indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp7    = jsonIndicesObj.at(L"p7indices");
    JSONObject jsonIndicesp7Obj = jsonIndicesp7->AsObject();
    JSONArray arrayp7 = jsonIndicesp7Obj[L"p7index"]->AsArray();
    for (unsigned int i = 0; i < arrayp7.size(); i++){
	JSONObject temp=	arrayp7[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p7indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p7indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p7#"    +tempID);
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p7REVC#"+tempID);
    }


    //p5 indices
    if(jsonIndicesObj.find(L"p5indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p5indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp5    = jsonIndicesObj.at(L"p5indices");
    JSONObject jsonIndicesp5Obj = jsonIndicesp5->AsObject();
    JSONArray arrayp5 = jsonIndicesp5Obj[L"p5index"]->AsArray();
    for (unsigned int i = 0; i < arrayp5.size(); i++){
	JSONObject temp=	arrayp5[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p5indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p5indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());

	trieKnownString->insertIntoTree( tempSeq.c_str() , "p5#"    +tempID);
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p5REVC#"+tempID);
    }




    //Other sequences after the primming site
    string IS4="AGATCTC";
    trieKnownString->insertIntoTree( IS4.c_str() ,                    "IS4");
    trieKnownString->insertIntoTree( reverseComplement(IS4).c_str() , "REVC#IS4");


}


void RGAssign::printUnfoundToFile(vector< pair<string,int> > * unfound,stringstream & fileError){

    for(int i=0;i<min(int(unfound->size()),maxErrorHits);i++){	       
	//Searching in known strings
	vector<string> temp = allTokens((*unfound)[i].first,'#');
	vector<string> temp2;
	for(unsigned int j=0;j<temp.size();j++){
	    vector<string> * temp3=new vector<string>();
	    vector<string> * temp4=new vector<string>();

	    trieKnownString->searchMismatch(temp[j].c_str(),temp3,0);
	    trieKnownString->searchMismatch( ( "N"+temp[j].substr(0, temp[j].size() -1)     ).c_str(),temp4,1);
	    trieKnownString->searchMismatch( (     temp[j].substr(1, temp[j].size() -1)+"N" ).c_str(),temp4,1);

	    //adding a tag before the shifted ones
	    for(unsigned int k=0;k<temp4->size();k++)
		(*temp4)[k]="SHFT#"+(*temp4)[k];
	    temp3->insert( temp3->end(), temp4->begin(), temp4->end() );

	    if(temp3->size() == 0 && temp4->size() == 0){
		temp2.push_back( "?");
	    }else{
		temp2.push_back( vectorToString(*temp3,","));  
	    }


	    delete temp3;
	    delete temp4;

	    
	}
	
	fileError<< 
	    vectorToString( temp,"\t" )<<"\t"<<
	    (*unfound)[i].second<<"\t"<<
	    vectorToString( temp2,"\t" )<< endl;
    }
}

void RGAssign::check_thresholds( rgAssignment &rg ) {
    rg.unknown  = -10 * rg.logLikelihoodScore > rgScoreCutoff ;
    rg.conflict = -10 * rg.logRatioTopToSecond < fracConflict && rg.logRatioTopToSecond < 0 ;
    rg.wrong    = -10 * rg.topWrongToTopCorrect > wrongness ;
}

void RGAssign::processSingleEndReads( BamAlignment &al){ //, BamWriter &writer, bool printError, map<string,int> &unknownSeq, map<string,int> &wrongSeq, map<string,int> &conflictSeq)
    //{
    string index1;
    string index1Q;
    string index2;
    string index2Q;

    getIndices(al,index1,index1Q,index2,index2Q);

    rgAssignment rgReturn=assignReadGroup(index1,index1Q,index2,index2Q,rgScoreCutoff,fracConflict,mismatchesTrie);
    check_thresholds( rgReturn ) ;

    updateRecord(al,rgReturn);
    //writer.SaveAlignment(al);

    //record unresolved indices
    if(printError){
	string keyIndex;
	if(index2.empty()){
	    keyIndex=index1;
	}else{
	    keyIndex=index1+"#"+index2;
	}

	if( rgReturn.conflict ) conflictSeq[ keyIndex ] ++;
	if( rgReturn.unknown  ) unknownSeq [ keyIndex ] ++;
	if( rgReturn.wrong    ) wrongSeq   [ keyIndex ] ++;
    }

}

void RGAssign::processPairedEndReads( BamAlignment &al, BamAlignment &al2){//, BamWriter &writer, bool printError, map<string,int> &unknownSeq, map<string,int> &wrongSeq, map<string,int> &conflictSeq)

    string index1;
    string index1Q;
    string index2;
    string index2Q;

    string sindex1;
    string sindex1Q;
    string sindex2;
    string sindex2Q;

    //retrieve indices
    getIndices(al,index1,index1Q,index2,index2Q);
    //check to see if the other indices are the same just for fun
    getIndices(al2,sindex1,sindex1Q,sindex2,sindex2Q);

    if(index1 !=sindex1 ){cerr<<"Seq#1 has a different index 1 than seq #2, exiting "        <<al.Name<<" vs "<<al2.Name<< endl; exit(1);}
    if(index1Q!=sindex1Q){cerr<<"Seq#1 has a different index 1 quality than seq #2, exiting "<<al.Name<<" vs "<<al2.Name<< endl; exit(1);}
    if(index2 !=sindex2 ){cerr<<"Seq#1 has a different index 2 than seq #2, exiting "        <<al.Name<<" vs "<<al2.Name<< endl; exit(1);}
    if(index2Q!=sindex2Q){cerr<<"Seq#1 has a different index 2 quality than seq #2, exiting "<<al.Name<<" vs "<<al2.Name<< endl; exit(1);}

    rgAssignment rgReturn = assignReadGroup(index1,index1Q,index2,index2Q,rgScoreCutoff,fracConflict,mismatchesTrie);
    check_thresholds( rgReturn ) ;

    updateRecord(al, rgReturn);
    updateRecord(al2,rgReturn);
//     writer.SaveAlignment(al);
//     writer.SaveAlignment(al2);

    // record unresolved indices
    if(printError){
        string keyIndex;
        if(index2.empty()){
            keyIndex=index1;
        }else{
            keyIndex=index1+"#"+index2;
        }

        if( rgReturn.conflict ) conflictSeq[ keyIndex ] += 2;
        if( rgReturn.unknown  ) unknownSeq [ keyIndex ] += 2;
        if( rgReturn.wrong    ) wrongSeq   [ keyIndex ] += 2;
    }
}


const map<string,string> *  RGAssign::getRGS(){
    return &rgs;
}

RGAssign::RGAssign( double rgScoreCutoff_  ,
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
		    ) :

    rgScoreCutoff(rgScoreCutoff_),
    fracConflict(fracConflict_),
    wrongness(wrongness_),
    mismatchesTrie(mismatchesTrie_),
    maxErrorHits(maxErrorHits_),
    shiftByOne(shiftByOne_),
    printSummary(printSummary_),
    printError(printError_),
 
    flag_ratioValues(flag_ratioValues_),
    flag_rgqual(flag_rgqual_)

{
    rgs = readIndexFile(index,mismatchesTrie,shiftByOne);
    map<string,string>::const_iterator itRG;   

    for ( itRG  = rgs.begin(); 
	  itRG != rgs.end(); 
	  itRG++ ){
	namesMap[ itRG->first ].assigned =0;
	namesMap[ itRG->first ].unknown  =0;
	namesMap[ itRG->first ].conflict =0;
	namesMap[ itRG->first ].wrong    =0;
    }

    namesMap[ "unknown" ].assigned=0;
    namesMap[ "unknown" ].unknown=0;
    namesMap[ "unknown" ].conflict=0;
    namesMap[ "unknown" ].wrong=0;

    if(flag_ratioValues){
	ratioValues = ratioValues_;
    }

    if(flag_rgqual){
	rgqual = rgqual_;
    }



    dashes = "--------------------------------";
    if(printError){
	trieKnownString = new PrefixTree<string>();
	initializeKnownIndices(trieKnownString,getCWD()+"/../webForm/config.json");
	// //debug
	// vector<string> * temp3=new vector<string>();
	// vector<string> * temp4=new vector<string>();
	// string t="AAGGTCT";
	// trieKnownString->searchMismatch(t.c_str(),temp3,0);
	// trieKnownString->searchMismatch( ( "N"+t.substr(0, t.size() -1)     ).c_str(),temp4,1);
	// trieKnownString->searchMismatch( (     t.substr(1, t.size() -1)+"N" ).c_str(),temp4,1);
	// cout<<"t3 "<<vectorToString(*temp3)<<endl;
	// cout<<"t4 "<<vectorToString(*temp4)<<endl;

	// return 1;
	// //end debug
    }
}


// void deallocate() {
//     // cout<<"deallocate"<<endl;
//     delete indTrie1;
//     delete indTrie2;
// }

string RGAssign::getSummaryString(){
    map<string,tallyForRG>::iterator it;   
    unsigned int totalRG=0;	
    unsigned int totalAssignRG=0;	

    vector< pair<string,tallyForRG> > toprintVec;
    for ( it=namesMap.begin() ; it != namesMap.end(); it++ ){
	toprintVec.push_back(  make_pair( it->first , it->second ) );
	totalRG+=it->second.assigned+it->second.unknown+it->second.conflict+it->second.wrong;
    }

    sort (toprintVec.begin(),   toprintVec.end(),   compareNameTally() );

    stringstream fileSummary;
    fileSummary << "RG\ttotal\ttotal%\tassigned\tassigned%\tunknown\tunknown%\tconflict\tconflict%\twrong\twrong%"<<endl;
    fileSummary<<dashes<<endl;
    for(unsigned int i=0;i<toprintVec.size();i++){		
	unsigned int totalForRQ=toprintVec[i].second.assigned+toprintVec[i].second.unknown+toprintVec[i].second.conflict+toprintVec[i].second.wrong;

	fileSummary << toprintVec[i].first << "\t" << totalForRQ << "\t"
		    << 100.0*double(totalForRQ)/double(totalRG) << "%\t" ;

	fileSummary  << toprintVec[i].second.assigned << "\t"
		     << 100.0*double(toprintVec[i].second.assigned)/double(totalForRQ) << "%\t" ;

	fileSummary  << toprintVec[i].second.unknown << "\t"
		     << 100.0*double(toprintVec[i].second.unknown)/double(totalForRQ) << "%\t" ;

	fileSummary <<  toprintVec[i].second.conflict << "\t"
		    << 100.0*double(toprintVec[i].second.conflict)/double(totalForRQ) << "%\t" ;

	fileSummary <<  toprintVec[i].second.wrong << "\t"
		    << 100.0*double(toprintVec[i].second.wrong)/double(totalForRQ) << "%\n" ;
		
	if(toprintVec[i].first != "unknown" )
	    totalAssignRG+=toprintVec[i].second.assigned;
    }

    fileSummary<<dashes<<endl;
    fileSummary<<"ASSIGNED:\t"<< totalAssignRG<<"\t"<<100.0*double(totalAssignRG)/double(totalRG)<<"%"<<endl;
    fileSummary<<"PROBLEMS:\t"<< (totalRG-totalAssignRG)<<"\t"<<100.0*double(totalRG-totalAssignRG)/double(totalRG)<<"%"<<endl;

    fileSummary<<"TOTAL:\t"<<totalRG<<"\t100.0%"<<endl;

    return fileSummary.str();

}

string RGAssign::getErrorString(){
    vector< pair<string,int> > conflictToPrint( conflictSeq.begin(), conflictSeq.end() ) ;
    vector< pair<string,int> > unknownToPrint(  unknownSeq.begin(),  unknownSeq.end() ) ;
    vector< pair<string,int> > wrongToPrint(    wrongSeq.begin(),    wrongSeq.end() ) ;
    
    sort (conflictToPrint.begin(),   conflictToPrint.end(),   compareNameRG() ); 
    sort (unknownToPrint.begin(),    unknownToPrint.end(),    compareNameRG() ); 
    sort (wrongToPrint.begin(),      wrongToPrint.end(),      compareNameRG() ); 
    stringstream fileError;
    fileError<<      dashes<<endl<<"Conflict:"<<endl<<dashes<<endl;
    printUnfoundToFile(&conflictToPrint,fileError);
    
    fileError<<endl<<dashes<<endl<<"Unknown:" <<endl<<dashes<<endl;
    printUnfoundToFile(&unknownToPrint,fileError);
    
    fileError<<endl<<dashes<<endl<<"Wrong:"   <<endl<<dashes<<endl;
    printUnfoundToFile(&wrongToPrint,fileError);

    return fileError.str();
}

RGAssign::~RGAssign() {
    // cout<<"deallocate"<<endl;
    delete indTrie1;
    delete indTrie2;
    if(printError){
	delete trieKnownString;
    }

}
