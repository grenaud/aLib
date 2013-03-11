#include <iostream>
#include <string>
#include <cstring>
#include <sys/stat.h>

#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include "MergeTrimReads.h"
#include "PutProgramInHeader.h"

////////////////////////////////
// TO DO
//
////////////////////////////////

using namespace std;
using namespace BamTools;
// using namespace __MergeTrimReads__;


string options_adapter_F_BAM="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG";
string options_adapter_S_BAM="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
string options_adapter_chimera_BAM="ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA";
int maxadapterComp_BAM = 30;

const string MERGEDBAMFLAG = "FF";
const int32_t TRIMMEDFLAG       = 1;
const int32_t MERGEDFLAG        = 2;
const int32_t TRIMMEDMERGEDFLAG = 3;

string boolStringify(const bool b){
    return b ? "true" : "false";
}

string intStringify(int i) {
    stringstream s;
    s << i;
    return s.str();
}

inline string sortUniqueChar(string v){
    // vector<char> v( input.begin(), input.end() );
    string::iterator it;
    sort (v.begin(), v.end());
    char   previousChar='!';
    string toReturn    ="";
    for (it=v.begin(); it!=v.end();it++ ){
	if(previousChar == '!' ||
	   previousChar != *it){
	    previousChar=*it;
	    toReturn+=previousChar;
	}	
    }
    return toReturn;
}

bool isaDirectory(const string& directory){
    struct stat sb;
    return !stat(directory.c_str(), &sb) && S_ISDIR(sb.st_mode); 
}


int main (int argc, char *argv[]) {

    bool produceUnCompressedBAM=false;
    bool verbose=false;
    bool mergeoverlap=false;
    string adapter_F=options_adapter_F_BAM;
    string adapter_S=options_adapter_S_BAM;
    string adapter_chimera=options_adapter_chimera_BAM;
    string key="";
    bool allowMissing=false;
    int trimCutoff=1;
    
    bool printLog=false;
    string logFileName;

    BamReader reader;
    BamWriter writer;

    string bamFile;
    string bamFileOUT="";
    vector<string> checkedTags;
    checkedTags.push_back("RG");
    checkedTags.push_back("XI");
    checkedTags.push_back("YI");
    checkedTags.push_back("XJ");
    checkedTags.push_back("YJ");

    const string usage=string(string(argv[0])+
			      "This program takes a BAM where each mate are consecutive and\ntrims and merges reads\n"+
			      +" [options] BAMfile"+"\n"+
			      //"\t"+"-p , --PIPE"+"\n\t\t"+"Read BAM from and write it to PIPE"+"\n"+
			      "\t"+"-o , --outfile" +"\t\t"+"Output (BAM format)."+"\n"+
			      "\t"+"-u            " +"\t\t"+"Produce uncompressed bam (good for pipe)"+"\n"+

			      //	"\t"+" , --outprefix" +"\n\t\t"+"Prefix for output files (default '"+outprefix+"')."+"\n"+
			      //"\t"+" , --SAM" +"\n\t\t"+"Output SAM not BAM."+"\n"+
			      "\t"+"-v , --verbose" +"\t\t"+"Turn all messages on (default "+boolStringify(verbose)+")"+"\n"+
			      "\t"+"--log [log file]" +"\t"+"Print a tally of merged reads to this log file (default only to stderr)"+"\n"+
			      
			      "\n\t"+"Paired End merging/Single Read trimming  options"+"\n"+
			      "\t\t"+"--mergeoverlap"+"\t\t\t\t"+"Merge PE reads of molecules longer than read length that show a minimum overlap (default "+boolStringify(mergeoverlap)+")"+"\n"+
			      "\t\t\t\t\t\tGood for merging ancient DNA reads into a single sequence\n\n"
			      "\t\t"+"-f , --adapterFirstRead" +"\t\t\t"+"Adapter that is observed after the forward read (def. Multiplex: "+options_adapter_F_BAM .substr(0,maxadapterComp_BAM)+")"+"\n"+
			      "\t\t"+"-s , --adapterSecondRead" +"\t\t"+"Adapter that is observed after the reverse read (def. Multiplex: "+options_adapter_S_BAM.substr(0,maxadapterComp_BAM)+")"+"\n"+
			      "\t\t"+"-c , --FirstReadChimeraFilter" +"\t\t"+"If the forward read looks like this sequence, the cluster is filtered out.\n\t\t\t\t\t\t\tProvide several sequences separated by comma.(def. Multiplex: "+options_adapter_chimera_BAM.substr(0,maxadapterComp_BAM)+")"+"\n"+
			      "\t\t"+"-k , --key"+"\t\t\t\t"+"Key sequence with which each sequence starts. Comma separate for forward and reverse reads. (default '"+key+"')"+"\n"+
			      "\t\t"+"-i , --allowMissing"+"\t\t\t"+"Allow one base in one key to be missing or wrong. (default "+boolStringify(allowMissing)+")"+"\n"+
			      "\t\t"+"-t , --trimCutoff"+"\t\t\t"+"Lowest number of adapter bases to be observed for single Read trimming (default "+intStringify(trimCutoff)+")");

    if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
    	cout<<"Usage:"<<endl;
    	cout<<""<<endl;
    	cout<<usage<<endl;
    	return 1;
    }

    for(int i=1;i<(argc-1);i++){ //all but the last arg

	if(strcmp(argv[i],"--log") == 0 ){
	    logFileName =string(argv[i+1]);
	    printLog=true;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-p") == 0 || strcmp(argv[i],"--PIPE") == 0 ){
	    cerr<<"This version no longer works with pipe, exiting"<<endl;
	    return 1;	    
	}

	if(strcmp(argv[i],"-u") == 0  ){
	    produceUnCompressedBAM=true;
	    continue;
	}

	if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	    bamFileOUT =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"--verbose") == 0 ){
	    verbose=true;
	    continue;
	}

	if(strcmp(argv[i],"--mergeoverlap") == 0 ){
	    mergeoverlap=true;
	    continue;
	}

	if(strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--adapterFirstRead") == 0 ){
	    adapter_F =string(argv[i+1]);
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--adapterSecondRead") == 0 ){
	    adapter_S =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"--FirstReadChimeraFilter") == 0 ){
	    adapter_chimera =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-k") == 0 || strcmp(argv[i],"--keys") == 0 ){
	    key =string(argv[i+1]);
	    i++;
	    continue;
	}
	

	if(strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--allowMissing") == 0 ){
	    allowMissing=true;
	    continue;
	}

	if(strcmp(argv[i],"-t") == 0 || strcmp(argv[i],"--trimCutoff") == 0 ){
	    trimCutoff=atoi(argv[i+1]);
	    i++;
	    continue;
	}
	
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;	    
    }

    bamFile=argv[argc-1];

    set_adapter_sequences(adapter_F,
			  adapter_S,
			  adapter_chimera,
			  maxadapterComp_BAM);
    set_options(trimCutoff,allowMissing,mergeoverlap);
    
    if(key != ""){
	size_t found=key.find(",");
	if (found == string::npos){ //single end reads
	    set_keys(key);
	} else{                     //paired-end
	    set_keys(key.substr(0,found),
		     key.substr(found+1,key.length()-found+1));
	}
    }

    if( bamFileOUT == ""  ){
	cerr<<"The output must be a be specified, exiting"<<endl;
	return 1;
    }

    if ( !reader.Open(bamFile) ) {
    	cerr << "Could not open input BAM file  "<<bamFile << endl;
    	return 1;
    }
    SamHeader header = reader.GetHeader();

    string pID          = "mergeTrimReadsBAM";   
    string pName        = "mergeTrimReadsBAM";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine);

    const RefVector references = reader.GetReferenceData();
    //we will not call bgzip with full compression, good for piping into another program to 
    //lessen the load on the CPU
    if(produceUnCompressedBAM) 
	writer.SetCompressionMode(BamWriter::Uncompressed);

    if ( !writer.Open(bamFileOUT,header,references) ) {
    	cerr << "Could not open output BAM file "<<bamFileOUT << endl;
    	return 1;
    }



    SamHeader sh=reader.GetHeader();
    //Up to the user to be sure that a sequence is followed by his mate
    // if(!sh.HasSortOrder() || 
    //    sh.SortOrder != "queryname"){
    // 	cerr << "Bamfile must be sorted by queryname" << endl;
    // 	return 1;
    // }

    int count_all = 0;
    int count_fkey = 0;
    int count_merged = 0;
    int count_merged_overlap = 0;
    int count_trimmed = 0;
    int count_nothing = 0;
    int count_chimera = 0;

    string read1;
    string read2;
    string qual1;
    string qual2;

    BamAlignment al;
    BamAlignment al2;
    bool al2Null=true;
    
    while ( reader.GetNextAlignment(al) ) {

	if(al.IsMapped() || al.HasTag("NM") || al.HasTag("MD")  ){
	    cerr << "Reads should not be aligned" << endl;
	    return 1;
	}


	if(al.IsPaired() && 
	   al2Null ){
	    al2=al;
	    al2Null=false;
	    continue;
	}else{
	    if(al.IsPaired() && 
	       !al2Null){
		if(al.Name != al2.Name ){
		    cerr << "Seq#1 has a different id than seq #2, exiting " << endl;
		    return 1;
		} 
		count_all ++;
		if(al.IsFirstMate()  &&
		   al2.IsSecondMate() ){
		    read1 =string(al.QueryBases);
		    read2 =string(al2.QueryBases);
		    qual1 =string(al.Qualities);
		    qual2 =string(al2.Qualities);		    
		}else{
		    if(al2.IsFirstMate()  &&
		       al.IsSecondMate() ){
			read1 =string(al2.QueryBases);
			qual1 =string(al2.Qualities);
			read2 =string(al.QueryBases);
			qual2 =string(al.Qualities);		    
		    }else{
			cerr << "Seq#1 must be the first mate for seq #2, exiting " << endl;
			return 1;
		    }
		}

		if(qual1 == "*"){
		    qual1=string(read1.length(),'0');
		}
		if(qual2 == "*"){
		    qual2=string(read1.length(),'0');
		}


		merged result=process_PE(read1,qual1,
					 read2,qual2);

		if(result.code != ' '){ 
		    string prevZQ1="";
		    string prevZQ2="";

		    al.SetIsFailedQC(true);
		    al.GetTag("ZQ",prevZQ1);		    
		    prevZQ1+=result.code;
		    if(al.HasTag("ZQ") ){ //this is done because bamtools was not intelligent enough to understand that "ZQ:A" becomes "ZQ:Z" when you add a char, oh well.. 
			al.RemoveTag("ZQ");
			if(al.HasTag("ZQ") ){
			    cerr << "Failed to remove tag for "<< al.Name<<endl;
			    return 1;
			}
		    }

		    if(prevZQ1 != "")
			if(!al.AddTag("ZQ","Z",sortUniqueChar(prevZQ1))){
			    cerr << "Error while editing tags new tag11:"<<prevZQ1 <<"#"<< endl;
			    return 1;
			}

		   

		    al2.SetIsFailedQC(true);
		    al2.GetTag("ZQ",prevZQ2);
		    prevZQ2+=result.code;
		    if(al2.HasTag("ZQ") ){ 
			al2.RemoveTag("ZQ");
			if(al2.HasTag("ZQ") ){ 
			    cerr << "Failed to remove tag for "<< al2.Name<< endl;
			    return 1;
			}
		    }
		    if(prevZQ2 != "")
		    if(!al2.AddTag("ZQ","Z",sortUniqueChar(prevZQ2))){
			cerr << "Error while editing tags new tag21:" << prevZQ2<<"#"<<endl;
			return 1;
		    }

		  
		    
		    if( result.code == 'K'){
			count_fkey ++;
		    }else{
			if( result.code  == 'D'){
			    count_chimera++;
			}
		    }
		}

		if(result.sequence != ""){ //new sequence
		  BamAlignment toWrite (al);//build from the previous one
		  string towriteZQ="";
		  al.GetTag("ZQ",towriteZQ);  //get from the first one
		  // if(!toWrite.RemoveTag("ZQ")){ 
		  //     cerr << "Failed to remove tag for new "<< toWrite.Name<< endl;
		  //     return 1;
		  // }
		  if( result.code == ' '){
		      if( result.sequence.length() > max(read1.length(),read2.length())){
			  count_merged_overlap ++;			  
			  if(!al.AddTag(MERGEDBAMFLAG,"i",MERGEDFLAG)){
			      cerr << "Unable to add tag" << endl;
			      return 1;
			  }
		      }else{
			  count_merged++;
			  if(!al.AddTag(MERGEDBAMFLAG,"i",TRIMMEDMERGEDFLAG)){
			      cerr << "Unable to add tag" << endl;
			      return 1;
			  }
		      }
		  }
		  
		  toWrite.AlignmentFlag=4;
		  toWrite.MapQuality=0;
		  
		  toWrite.QueryBases = result.sequence;
		  toWrite.Qualities  = result.quality;
		  toWrite.SetIsMapped(false);
		

		  toWrite.Position    =-1;
		  toWrite.MatePosition=-1;
		  toWrite.SetIsFailedQC( al.IsFailedQC() && al2.IsFailedQC() );
		  toWrite.TagData=al.TagData; //copy tag info
		  //toWrite.RemoveTag("ZQ");
		  if(toWrite.HasTag("ZQ") ){ 
		      toWrite.RemoveTag("ZQ");
		      if(toWrite.HasTag("ZQ") ){ 
			  cerr << "Failed to remove tag for new "<< toWrite.Name<< endl;
			  return 1;
		      }
		  }

		  if( toWrite.QueryBases.length()  < min_length){
		      toWrite.SetIsFailedQC( true );		   
		      // if(!al.EditTag("ZQ","Z",string("L"))){
		      // 	  cerr << "Error while editing tags" << endl;
		      // 	  return 1;
		      // }
		      towriteZQ+="L";
		  }

		  /////////////////////////////
		  //       Fixing tags       //
		  /////////////////////////////
		  string dummy1;
		  string dummy2;
		  

		  //paranoid check to make sure our tags are identical
		  for(size_t idx=0;idx<checkedTags.size();idx++)
		      if(al.GetTag(checkedTags[idx],dummy1)){
			  if(al2.GetTag(checkedTags[idx],dummy2)){
			      if(dummy1 != dummy2){
				  cerr << "Value for "<<checkedTags[idx]<<" cannot differ between mates " << endl;
				  return 1;
			      }
			      //fine otherwise
			  }else{
			      cerr << "One read has been assigned a  "<<checkedTags[idx]<<" tag but not the other " << endl;
			      return 1;

			  }
		      }

		  //The new read has the same ZQ tag as al. at this point
		  if(al.HasTag("ZQ") || al2.HasTag("ZQ") ){
		      if( al2.HasTag("ZQ") ){
			  if(!al.GetTag("ZQ",dummy1)) {
			      cerr << "Failed to get ZQ field from read 1" << endl;
			      return 1;
			  }
			  if(!al.GetTag("ZQ",dummy2)) {
			      cerr << "Failed to get ZQ field from read 2" << endl;
			      return 1;
			  }
			  //we then need to add the ZQ from the second read
			  if(dummy1 != dummy2){
			      towriteZQ+=dummy2;
			  }	       	    
		      }
		  }


		  if(towriteZQ != ""){
		      if(!toWrite.AddTag("ZQ","Z",sortUniqueChar(towriteZQ))){
			  cerr << "Error while editing tags new tag20:"<<towriteZQ <<"#"<< endl;
			  return 1;
		      }
		  }		    
		  writer.SaveAlignment(toWrite);
		    
		}else{
		    if( result.code == ' ')
			count_nothing++;
		    //keep the sequences as pairs
		    writer.SaveAlignment(al2);
		    writer.SaveAlignment(al);
		}





		//
		//  SINGLE END
		//
	    }else{ 
		count_all ++;

		read1 =string(al.QueryBases);
		qual1 =string(al.Qualities);
		if(qual1 == "*"){
		    qual1=string(read1.length(),'0');
		}
	

		merged result=process_SR(read1,qual1);
		
		if(result.code != ' '){ 
		    string prevZQ1="";

		    al.SetIsFailedQC(true);
		    al.GetTag("ZQ",prevZQ1);
		    prevZQ1+=result.code;
		    if(al.HasTag("ZQ") ){ 
			al.RemoveTag("ZQ");
			if(al.HasTag("ZQ") ){ 
			    cerr << "Failed to remove tag for "<< al.Name<< endl;
			    return 1;
			}
		    }

		    if(prevZQ1 != ""){
			if(!al.EditTag("ZQ","Z",sortUniqueChar(prevZQ1))){
			    cerr << "Error while editing tags new tag11:"<<prevZQ1 <<"#"<< endl;
			    return 1;
			}
		    }	
	    
		    if( result.code == 'K'){
			count_fkey ++;
		    }else{
			if( result.code  == 'D'){
			    count_chimera++;
			}
		    }
		}

		if(result.sequence != ""){ //new sequence

		    BamAlignment toWrite (al);//build from the previous al

		    toWrite.MapQuality=0;


		    if(!al.AddTag(MERGEDBAMFLAG,"i",TRIMMEDFLAG)){
		      cerr << "Unable to add tag" << endl;
		      return 1;
		    }

		    
		    toWrite.QueryBases = result.sequence;
		    toWrite.Qualities  = result.quality;
		    toWrite.SetIsMapped(false);
		    toWrite.TagData=al.TagData; //copy tag info

		    toWrite.Position    =-1;
		    toWrite.MatePosition=-1;		    
		    if( result.code == ' ')
			count_trimmed++;
		    writer.SaveAlignment(toWrite);

		}else{
		    if( result.code == ' ')
			count_nothing++;

		    writer.SaveAlignment(al);
		}



	    } //end single end
	    al2Null=true;
	}//second pair
	
	if(verbose && (count_all%10000 == 0) && (count_all > 0)){
	    cerr<<"<==> Total "<<count_all<<"; Merged (trimming) "<<count_merged<<"; Merged (overlap) "<<count_merged_overlap<<"; Kept PE/SR "<<count_nothing<<"; Trimmed SR "<<count_trimmed<<"; Adapter dimers/chimeras "<<count_chimera<<"; Failed Key "<<count_fkey<<endl;	
	}
	    

    } //while al
    reader.Close();
    writer.Close();

    cerr <<"Total "<< count_all<<"; Merged (trimming) "<<count_merged <<"; Merged (overlap) "<<count_merged_overlap <<"; Kept PE/SR "<< count_nothing<<"; Trimmed SR "<<count_trimmed <<"; Adapter dimers/chimeras "<<count_chimera <<"; Failed Key "<<count_fkey <<endl;

    if(printLog){
	ofstream fileLog;
	fileLog.open(logFileName.c_str());

	if (fileLog.is_open()){
	    fileLog <<"Total reads :"            <<count_all<<           "\t"<<100.0*double(count_all)           /double(count_all)<<"%"<<endl;
	    fileLog <<"Merged (trimming) "       <<count_merged<<        "\t"<<100.0*double(count_merged)        /double(count_all)<<"%"<<endl;
	    fileLog <<"Merged (overlap) "        <<count_merged_overlap<<"\t"<<100.0*double(count_merged_overlap)/double(count_all)<<"%"<<endl;
	    fileLog <<"Kept PE/SR "              <<count_nothing<<       "\t"<<100.0*double(count_nothing)       /double(count_all)<<"%"<<endl;
	    fileLog <<"Trimmed SR "              <<count_trimmed<<       "\t"<<100.0*double(count_trimmed)       /double(count_all)<<"%"<<endl;
	    fileLog <<"Adapter dimers/chimeras " <<count_chimera<<       "\t"<<100.0*double(count_chimera)       /double(count_all)<<"%"<<endl;
	    fileLog <<"Failed Key "              <<count_fkey<<          "\t"<<100.0*double(count_fkey)          /double(count_all)<<"%"<<endl;

	}else{
	    cerr << "Unable to print to file "<<logFileName<<endl;
	}
	fileLog.close();
    }
   
    return 0;
}

