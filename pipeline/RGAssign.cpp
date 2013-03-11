#include "RGAssign.h"


// #define DEBUG
// #define DEBUG2

PrefixTree<int> * indTrie1;
PrefixTree<int> * indTrie2;
indexData values;

ofstream * ratioValues;
ofstream * rgqual;
bool flag_ratioValues=false;
bool flag_rgqual=false;



//array containing the likelihood
double likeMatch[64];
double likeMismatch[64];

bool shiftByOne=false;

void setFileForRatio(ofstream * streamFile){
    flag_ratioValues=true;
    ratioValues=streamFile;
}

void setFileForRGQual(ofstream * streamFile){
    flag_rgqual=true;
    rgqual=streamFile;
}



bool comparePair (pair<int,double> i,pair<int,double> j) { 
    return (i.second>j.second); 
}

string toUpperCase(string toCheck){
    string toReturn=toCheck;
    transform(toReturn.begin(), toReturn.end(), toReturn.begin(), ::toupper);
    return toReturn;
}

bool isValidDNA(string tocheck){
    for(unsigned int  i=0;i<tocheck.size();i++){
	if(tocheck[i] != 'A' &&
	   tocheck[i] != 'C' &&
	   tocheck[i] != 'G' &&
	   tocheck[i] != 'T' )
	    return false;	  
    }
    return true;
}


void checkRGname(string tocheck){
    if(tocheck == "unknown") {  cerr<<"Error: Cannot use reserved name \"unknown\" as read group name"<<endl;  exit(1);}
    if(tocheck == "conflict"){  cerr<<"Error: Cannot use reserved name \"conflict\" as read group name"<<endl; exit(1);}
    if(tocheck == "wrong")   {  cerr<<"Error: Cannot use reserved name \"wrong\" as read group name"<<endl;    exit(1);}
    //fine
}

indexData intern_readIndex(string filename){
    string line;
    ifstream myFile;

    indexData toReturn;
    toReturn.mlindex1=0;
    toReturn.mlindex2=0;

    bool isFirstLine  =true;


    //initialize the values for the likelihood of matches or mismatches 
    for(int i=0;i<64;i++){
	if(i == 0)
	    likeMatch[i]    = -3.0; // need to fix this to minimize precision loss
	else
	    likeMatch[i]    = log10(1.0- pow(10.0,i/-10.0) ); // need to fix this to minimize precision loss
	
	likeMismatch[i]     = log10(pow(10.0,i/-10.0) );	
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
			if(numberOfFields==3){
			    toReturn.isDoubleIndex=true; 
			}else{
			    cerr << "Must have 2 or 3 fields"<<endl;
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
				}else{
				    if(fieldIndex==2){
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

					// cerr<<"Wrong field index"<<endl;
					// exit(1);
				    }
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

void deallocate(bool isDouble){
    // cout<<"deallocate"<<endl;
    delete indTrie1;
    if(isDouble)
	delete indTrie2;

}

pair<bool,map<string,string> > readIndexFile(string filename,int mismatchesTrie,bool _shiftByOne){
    values = intern_readIndex(filename);
    shiftByOne=_shiftByOne;
    
    indTrie1 = new PrefixTree<int>  ();    
    if(values.isDoubleIndex){
	indTrie2 = new PrefixTree<int>  ();
    }

    for(int i=0;i<int(values.names.size());i++){

	indTrie1->insertIntoTree( values.indices1[i].c_str() , i);
	if(values.isDoubleIndex)
	    indTrie2->insertIntoTree( values.indices2[i].c_str() , i);
	
    }

   

    //detect conflicts ?
    cerr<<"Conflicts for index1:"<<endl;
    for(int i=0;i<int(values.names.size());i++){
	//vector< matches<int> > * matchesind=indTrie1->searchMismatch(values.indices1[i].c_str(),mismatchesTrie);
	vector< int > * matchesind=new vector< int >();
	indTrie1->searchMismatch(values.indices1[i].c_str(),matchesind,mismatchesTrie);
	string toPrint="";

	for(unsigned int j=0;j<matchesind->size();j++){
	    // list< int  >::const_iterator iter;		    
	    // for (iter  = (*matchesind)[j].listOfDeflinesOfMatches->begin(); 
	    // 	 iter != (*matchesind)[j].listOfDeflinesOfMatches->end(); 
	    // 	 iter++){
	    // 	if(*iter != i)
	    // 	    toPrint+=values.names[i]+" ";
	    // }	    
	    if(matchesind->at(j) != i)
		toPrint+=values.names[ matchesind->at(j)  ]+" ";
	}

	if(toPrint.length() > 0){
	    if(values.names[i] != toPrint)
		cerr<<values.indices1[i]<<" from #"<<values.names[i]<<"# causes a conflict with #"<<toPrint<<"#"<<endl;
	}
	
	delete(matchesind);
    }

    if(values.isDoubleIndex){
	cerr<<"Conflicts for index2:"<<endl;
	for(int i=0;i<int(values.names.size());i++){

	    //vector< matches<int> > * matchesind=indTrie2->searchForWordMismatch(values.indices2[i].c_str(),mismatchesTrie);
	    vector< int > * matchesind=new vector< int >();
	    indTrie2->searchMismatch(values.indices2[i].c_str(),matchesind,mismatchesTrie);

	    string toPrint="";

	    for(unsigned int j=0;j<matchesind->size();j++){
		// list< int  >::const_iterator iter;		    
		// for (iter  = (*matchesind)[j].listOfDeflinesOfMatches->begin(); 
		//      iter != (*matchesind)[j].listOfDeflinesOfMatches->end(); 
		//      iter++){
		//     if(*iter != i)
		// 	toPrint+=values.names[i]+" ";
		// }
		if(matchesind->at(j) != i)
		    toPrint+=values.names[ matchesind->at(j)  ]+" ";
	    }

	    if(toPrint.length() > 0){

		cerr<<values.indices1[i]<<" from "<<values.names[i]<<" causes a conflict with "<<toPrint<<endl;
	    }
	    delete(matchesind);
	}

    }


    return make_pair(values.isDoubleIndex, values.namesMap) ;
}

inline double computeLike(const string & indexRef,const string & indexRead,const vector<int> * quals){
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

rgAssignment assignReadGroup(bool doubleInd,string & index1,string & index1q,string & index2,string & index2q,double rgScoreCutoff,double fracConflict,int mismatchesTrie){

  
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
    toReturn.predictedGroup  = "";
    toReturn.conflict        = false;
    toReturn.wrong           = false;
    toReturn.unknown         = false;
    toReturn.likelihoodScore = 0.0;

    int qualOffset=33;

    vector<int> quals1;
    vector<int> quals2;

    //Find matches using the prefix trie

    vector< int > * matchesind1=new vector< int > ();
    indTrie1->searchMismatch(      index1.c_str(),matchesind1,mismatchesTrie);

    
    if(shiftByOne){       
	indTrie1->searchMismatch(     ("N"+index1.substr(0, index1.size() -1)     ).c_str(),matchesind1,mismatchesTrie);
	indTrie1->searchMismatch(         (index1.substr(1, index1.size() -1)+"N" ).c_str(),matchesind1,mismatchesTrie);
    }

    for(unsigned int i=0;i<index1q.length();i++)
	quals1.push_back(     max( (int(char(index1q[i]))-qualOffset),2)  )   ; //since qual scores less than 2 do not make sense

    vector< int > * matchesind2;
    if(doubleInd){
	matchesind2=new vector< int > ();	
	indTrie2->searchMismatch(  index2.c_str(),matchesind2,mismatchesTrie);

	if(shiftByOne){       
	    indTrie2->searchMismatch( ("N"+index2.substr(0, index2.size() -1)     ).c_str(),matchesind2,mismatchesTrie);
	    indTrie2->searchMismatch(     (index2.substr(1, index2.size() -1)+"N" ).c_str(),matchesind2,mismatchesTrie);
	}

	for(unsigned int i=0;i<index2q.length();i++)
	    quals2.push_back( max( (int(char(index2q[i]))-qualOffset),2)  )   ;	//since qual scores less than 2 do not make sense
    }


    set<int> foundIndices; //set of indices found

    double like1;
    double like2;

    vector< pair<int,double>  > sortedLikelihoodAll; //likelihood for both indices
    vector< pair<int,double>  > sortedLikelihood1;   //likelihood for index 1
    vector< pair<int,double>  > sortedLikelihood2;   //likelihood for index 1



    //putting all the index # into the set
    for(unsigned int j=0;j<matchesind1->size();j++){
	 for(unsigned int j=0;j<matchesind1->size();j++)
	     foundIndices.insert ( matchesind1->at(j)  );   
    }
    delete(matchesind1);


    if(doubleInd){
	 for(unsigned int j=0;j<matchesind2->size();j++)
	     foundIndices.insert ( matchesind2->at(j)  );   
	delete(matchesind2);
    }


    //for every # in the set, compute likelihood    
    for (set<int>::iterator sit=foundIndices.begin(); sit!=foundIndices.end(); sit++){

	like1      =  computeLike(index1,values.indices1[ *sit ],&quals1);

	// cerr<<endl<<"1 "<<endl<<index1<<endl;
	// cerr<<values.names[ *sit ]<<endl;
	// cerr<<values.indices1[ *sit ]<<endl;

	// cerr<<like1<<endl;
	// cerr<<computeLike( ("N"+index1.substr(0, index1.size() -1)     ), values.indices1[ *sit ],quals1)<<endl;
	// cerr<<computeLike( (    index1.substr(1, index1.size() -1)+"N" ), values.indices1[ *sit ],quals1)<<endl;


	if(shiftByOne){ //check if shifting by one improves the likelihood

	    like1      =  max(like1,
			      max(
				  computeLike( ("N"+index1.substr(0, index1.size() -1)     ), values.indices1[ *sit ],&quals1),
				  computeLike( (    index1.substr(1, index1.size() -1)+"N" ), values.indices1[ *sit ],&quals1)
				  )			      
			      );
	}

	sortedLikelihood1.push_back(     make_pair (*sit,like1) );

	if(doubleInd){

	    like2  =  computeLike(index2,values.indices2[ *sit ],&quals2);

	    // cerr<<endl<<"2 "<<endl<<index2<<endl;
	    // cerr<<("N"+index2.substr(0, index2.size() -1)     )<<endl;
	    // cerr<<(    index2.substr(1, index2.size() -1)+"N" )<<endl;
	    // cerr<<values.names[ *sit ]<<endl;
	    // cerr<<values.indices2[ *sit ]<<endl;
	    // cerr<<like2<<endl;
	    // cerr<<computeLike( ("N"+index2.substr(0, index2.size() -1)     ), values.indices2[ *sit ], quals2)<<endl;
	    // cerr<<computeLike( (    index2.substr(1, index2.size() -1)+"N" ), values.indices2[ *sit ], quals2)<<endl;


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
	toReturn.unknown         = true;
	return toReturn;
    }



    //sorting by likelihood
    sort (sortedLikelihood1.begin(),   sortedLikelihood1.end(),   comparePair); 
    sort (sortedLikelihood2.begin(),   sortedLikelihood2.end(),   comparePair); 
    sort (sortedLikelihoodAll.begin(), sortedLikelihoodAll.end(), comparePair); 


    // // //DEBUG
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
#endif

    //exit(1);
    // // //END DEBUG

    //DETECT WRONGS

    //Checking if first index is first for first and second index
    //have the assignment as the first choice
    //if not, 'wrong'
    if(!sortedLikelihood1.empty() &&
       !sortedLikelihood2.empty() ){
	double bestLike1 = sortedLikelihood1[0].second; //best likelihood for index 1

	//we need to do this if there are multiple best hits to avoid false "wrong"
	bool foundInFirstIndex1=false; //best likelihood for index 1
	for(unsigned int j=0;j<sortedLikelihood1.size();j++){

	    if(sortedLikelihood1[j].second != bestLike1)
		break;
	    if( (sortedLikelihoodAll[0].first ==  sortedLikelihood1[j].first ) )
		foundInFirstIndex1=true;
	}

	double bestLike2 = sortedLikelihood2[0].second; //best likelihood for index 2
	bool foundInFirstIndex2=false;
	for(unsigned int j=0;j<sortedLikelihood2.size();j++){

	    if(sortedLikelihood2[j].second != bestLike2)
		break;
	    if( (sortedLikelihoodAll[0].first ==  sortedLikelihood2[j].first ) )
		foundInFirstIndex2=true;
	}

	if( !foundInFirstIndex1 || !foundInFirstIndex2 ){ //meaning that the best index is not found for either the first or second
	    toReturn.wrong=true;
	    return toReturn;
	}    
    }
    



    //DETECT CONFLICTS
    //Checking likelihood of second best hit
    //if the ratio is too low, it's a conflict
    double probRG    = pow(10.0,sortedLikelihoodAll[0].second);

    if(sortedLikelihoodAll.size() > 1){
	double probRG2nd    = pow(10.0,sortedLikelihoodAll[1].second);
	toReturn.ratioTopToSecond=(probRG2nd/probRG);
	//cout<<"ratio "<<	toReturn.ratioTopToSecond<<endl;
	if(flag_ratioValues && toReturn.ratioTopToSecond > 0.00001)//to avoid very small values
	    ratioValues->write( (char *)&toReturn.ratioTopToSecond, sizeof(toReturn.ratioTopToSecond));
	if( (probRG2nd/probRG) > fracConflict){ //if the second has a probability score within 80% of the second one, flag as conflict
	    toReturn.conflict=true;
	    return toReturn;
	}
    }else{
	toReturn.ratioTopToSecond=0.0;
    }




    //DETECT UNKNOWNS
    double probRGAll = 0.0;
    for(unsigned int j=0;j<sortedLikelihoodAll.size();j++){
	probRGAll+=pow(10.0,sortedLikelihoodAll[j].second);
    }


    double scoreBayes1=-10.0*log10(1.0- (probRG*( probRG/probRGAll) ) );
#ifdef DEBUG
    cerr<<"probRG "<<probRG<<endl;
    cerr<<"probRGAll "<<probRGAll<<endl;
    cerr<<"rg qual "<<scoreBayes1<<endl;
#endif


    toReturn.likelihoodScore = scoreBayes1;
    if(flag_rgqual)
	rgqual->write( (char *)&toReturn.likelihoodScore, sizeof(toReturn.likelihoodScore));
    
    //cout<<"scoreBayes1 "<<scoreBayes1<<endl;
    if(scoreBayes1 < rgScoreCutoff){
	toReturn.unknown = true;
	return toReturn;
    }


    toReturn.predictedGroup = values.names[ sortedLikelihoodAll[0].first ];

    // delete(matchesind1);
    // if(doubleInd){
    // 	delete(matchesind2);
    // }



    return toReturn;







}
