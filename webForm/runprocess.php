<html><head>
<title>Illumina Form</title>
<body>
<?php

require "email.php";

ini_set('display_errors', 'On');
error_reporting(E_ALL);


if( ! file_exists(getcwd()."/config.xml")  ){
    echo "Configuration file not found ".getcwd()."/config.xml";
    exit(1);
 }

$xmlconf = simplexml_load_file( getcwd()."/config.xml" );

$emailAddrToSend=$xmlconf->emailAddrToSend;
$genomedirectory=$xmlconf->genomedirectory;
$outputdirectory=$xmlconf->outputdirectory;



// $emailAddrToSend=array("gabriel_renaud@eva.mpg.de");
// $genomedirectory="/mnt/solexa/Genomes/";
// $outputdirectory="/mnt/solexa/tmp/";
#$outputdirectory="/tmp/";


//convert this to XML and put it into config.xml ...
$protocol2chimera = array(  'IllMultiplex'  => array("Illumina Multiplex library",'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT','ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA'),
			    'IllTruSeq'     => array("Illumina TruSeq RNA/DNA library",'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT','ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA'),
			    'IllPE'         => array("Illumina Paired End library",'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG','AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT','AGATCTAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGT,GAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCG,AAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTG'),
			    'IllSR'         => array("Illumina Single Read library",'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG','',''),
			    'IllsRNA'       => array("Illumina smallRNA/DpnI library",'TCGTATGCCGTCTTCTGCTTG','',''),
			    'IllsRNA2'      => array("Illumina smallRNA II",'ATCTCGTATGCCGTCTTCTGCTTG','','TCTCGTATGCCGTCTTCTGCTTG,TCGTATGCCGTCTTCTGCTTG,GTATGCCGTCTTCTGCTTG'),
			    'IllsRNA2Index' => array("Illumina smallRNA II (5p Index)",'ATCTCGTATGCCGTCTTCTGCTTG','','TCTCGTATGCCGTCTTCTGCTTG,TCGTATGCCGTCTTCTGCTTG,GTATGCCGTCTTCTGCTTG'),
			    'IllSAGE'       => array("Illumina SAGE (NlaIII) library",'TCGTATGCCGTCTTCTGCTTG','','GTATGCCGTCTTCTGCTTG,CCGTCTTCTGCTTGTCGTATGCCGTCTTCTGCTTG,CGTATGCCGTCTTCTGCTTG,CCGTCTTCTGCTTGAAAAAAAAA,TCGGACTGTAGAACTCT'),
			    #'454Std': ("454 Standard library (converted)",'CTGAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGGATCTCGTATGCCGTCTTCTGCTTG','CTGAGACAGGGAGGGAACAGATGGGACACGCAGGGATGAGATGGGTGTAGATCTCGGTGGTCGCCGTATCATT','ACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGGATCTCGTATGCCGTCTTCTGCTTG'),
			    #'454Nea': ("454 Neandertal library (converted)",'GTCAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGGATCTCGTATGCCGTCTTCTGCTTG','GTCAGACAGGGAGGGAACAGATGGGACACGCAGGGATGAGATGGGTGTAGATCTCGGTGGTCGCCGTATCATT','ACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGGATCTCGTATGCCGTCTTCTGCTTG'),
			    #'454StdIndex': ("454 Standard library (Index converted)",'CTGAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGGIIIIIIIATCTCGTATGCCGTCTTCTGCTTG','CTGAGACAGGGAGGGAACAGATGGGACACGCAGGGATGAGATGGGTGTAGATCTCGGTGGTCGCCGTATCATT','ACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGGIIIIIIIATCTCGTATGCCGTCTTCTGCTTG'),
			    #'454NeaIndex': ("454 Neandertal library (Index converted)",'GTCAGACACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGGIIIIIIIATCTCGTATGCCGTCTTCTGCTTG','GTCAGACAGGGAGGGAACAGATGGGACACGCAGGGATGAGATGGGTGTAGATCTCGGTGGTCGCCGTATCATT','ACGCAACAGGGGATAGGCAAGGCACACAGGGGATAGGIIIIIIIATCTCGTATGCCGTCTTCTGCTTG'),			 
			    'IllMultiplexCR1' =>array("Illumina CR1 Multiplex (Neandertal group)",'GAGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG','GAGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',''),
			    'IllMultiplexCR2' =>array("Illumina CR2 Multiplex (Neandertal group)",'AGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG','AGACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',''),
			    'IllMultiplex_ssDNA' =>array("Illumina ssDNA Multiplex (ssDNA ligation protocol with 5nt deletion)",'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG','GGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT','') );

$p7_block300 = array("TCGCAGG",
		     "CTCTGCA",
		     "CCTAGGT",
		     "GGATCAA",
		     "GCAAGAT",
		     "ATGGAGA",
		     "CTCGATG",
		     "GCTCGAA",
		     "ACCAACT",
		     "CCGGTAC",
		     "AACTCCG",
		     "TTGAAGT",
		     "ACTATCA",
		     "TTGGATC",
		     "CGACCTG",
		     "TAATGCG",
		     "AGGTACC",
		     "TGCGTCC",
		     "GAATCTC",
		     "CATGCTC",
		     "ACGCAAC",
		     "GCATTGG",
		     "GATCTCG",
		     "CAATATG",
		     "TGACGTC",
		     "GATGCCA",
		     "CAATTAC",
		     "AGATAGG",
		     "CCGATTG",
		     "ATGCCGC",
		     "CAGTACT",
		     "AATAGTA",
		     "CATCCGG",
		     "TCATGGT",
		     "AGAACCG",
		     "TGGAATA",
		     "CAGGAGG",
		     "AATACCT",
		     "CGAATGC",
		     "TTCGCAA",
		     "AATTCAA",
		     "CGCGCAG",
		     "AAGGTCT",
		     "ACTGGAC",
		     "AGCAGGT",
		     "GTACCGG",
		     "GGTCAAG",
		     "AATGATG",
		     "AGTCAGA",
		     "AACTAGA",
		     "CTATGGC",
		     "CGACGGT",
		     "AACCAAG",
		     "CGGCGTA",
		     "GCAGTCC",
		     "CTCGCGC",
		     "CTGCGAC",
		     "ACGTATG",
		     "ATACTGA",
		     "TACTTAG",
		     "AAGCTAA",
		     "GACGGCG",
		     "AGAAGAC",
		     "GTCCGGC",
		     "TCAGCTT",
		     "AGAGCGC",
		     "GCCTACG",
		     "TAATCAT",
		     "AACCTGC",
		     "GACGATT",
		     "TAGGCCG",
		     "GGCATAG",
		     "TTCAACC",
		     "TTAACTC",
		     "TAGTCTA",
		     "TGCATGA",
		     "AATAAGC",
		     "AGCCTTG",
		     "CCAACCT",
		     "GCAGAAG",
		     "AGAATTA",
		     "CAGCATC",
		     "TTCTAGG",
		     "CCTCTAG",
		     "CCGGATA",
		     "GCCGCCT",
		     "AACGACC",
		     "CCAGCGG",
		     "TAGTTCC",
		     "TGGCAAT",
		     "CGTATAT",
		     "GCTAATC",
		     "GACTTCT",
		     "GTACTAT",
		     "CGAGATC",
		     "CGCAGCC");





$p7_block = array("ACAGTG",
		  "GATCAG",
		  "ATCACG",
		  "CGATGT",
		  "CTTGTA",
		  "GGCTAC",
		  "TGACCA",
		  "AAAGCA",
		  "AAATGC",
		  "AAGCGA",
		  "AAGGAC",
		  "AATAGG",
		  "ACCCAG",
		  "ACTCTC",
		  "AGAAGA",
		  "AGCATC",
		  "AGGCCG",
		  "ATACGG",
		  "ATCCTA",
		  "ATCTAT",
		  "ATGAGC",
		  "CATTTT",
		  "CCGCAA",
		  "CTCAGA",
		  "GAATAA",
		  "GCCGCG",
		  "GCTCCA",
		  "GGCACA",
		  "GGCCTG",
		  "TCGGCA",
		  "TCTACC",
		  "TGCCAT",
		  "TGCTGG",
		  "AGGTTT",
		  "AGTCAA",
		  "AGTTCC",
		  "ATGTCA",
		  "CCGTCC",
		  "GTAGAG",
		  "GTGAAA",
		  "GTGGCC",
		  "GTTTCG",
		  "CGTACG",
		  "GAGTGG",
		  "GGTAGC",
		  "ACTTGA",
		  "CAGATC",
		  "GCCAAT",
		  "TAGCTT",
		  "TTAGGC",
		  "AACCGCC",
		  "AACGAAC",
		  "AACGCCT",
		  "AACGGTA",
		  "AACTAGT",
		  "AACTGAG",
		  "AAGAATT",
		  "AAGATAG",
		  "AAGCTCT",
		  "AAGTCTG",
		  "AATAACC",
		  "AATCCGT",
		  "ACCGATT",
		  "ACCGTAG",
		  "ACCTCAT",
		  "ACCTTGC",
		  "ACGACCT",
		  "ACGATTC",
		  "ACGCGGC",
		  "ACGGAGG",
		  "ACGTAAC",
		  "ACTACTG",
		  "ACTCGTT",
		  "ACTGCGC",
		  "AGACCTC",
		  "AGACTAG",
		  "AGAGACC",
		  "AGAGCGT",
		  "AGATATG",
		  "AGATTCT",
		  "AGCAAGC",
		  "AGCAGTT",
		  "AGCGCTG",
		  "AGTATAC",
		  "ATAAGTC",
		  "ATAATGG",
		  "ATACTCC",
		  "ATAGAAG",
		  "ATCTCCG",
		  "ATGCAGT",                                           
		  "ATGGTAT",
		  "ATTATCT",
		  "ATTCGAC",
		  "ATTGCTA",
		  "CAACCGG",
		  "CAACTAA",
		  "AATCTTC",
		  "ACCAACG",
		  "AGATGGC",
		  "CCAGGTT",
		  "CCGTTAG",
		  "CGCCTCT",
		  "CTTGCGG",
		  "GGCGGAG",
		  "TGGACGT",
		  "AACCATG",
		  "CAGGAAG",
		  "CATACCT",
		  "CCAATCC",
		  "CCGGCGT",
		  "CGCATAG",
		  "CGTAATC",
		  "CGTTGGT",
		  "CTATACG",
		  "GACCTAC",
		  "GATATTG",
		  "AAGACGC",
		  "GCAGTAT",
		  "GGTCCGC",
		  "GTCGACT",                                           
		  "GTTAGAT",
		  "TAACTCG",
		  "TGCTTCC",
		  "TGGCGCT",
		  "AATGGCG",
		  "ACCAGAC",
		  "ACGCCAG",
		  "ACTAAGT",
		  "AGAACCG",
		  "ATCGTTC",
		  "CAACGTC");


$p5_block = array("TCGCAGG",
		  "CTCTGCA",
		  "CCTAGGT",
		  "GGATCAA",
		  "GCAAGAT",
		  "ATGGAGA",
		  "CTCGATG",
		  "GCTCGAA",
		  "ACCAACT",
		  "CCGGTAC",
		  "AACTCCG",
		  "TTGAAGT",
		  "ACTATCA",
		  "TTGGATC",
		  "CGACCTG",
		  "TAATGCG",
		  "AGGTACC",
		  "TGCGTCC",
		  "GAATCTC",
		  "CATGCTC",
		  "ACGCAAC",
		  "GCATTGG",
		  "GATCTCG",
		  "CAATATG",
		  "TGACGTC",
		  "GATGCCA",
		  "CAATTAC",
		  "AGATAGG",
		  "CCGATTG",
		  "ATGCCGC",
		  "CAGTACT",
		  "AATAGTA",
		  "CATCCGG",
		  "TCATGGT",
		  "AGAACCG",
		  "TGGAATA",
		  "CAGGAGG",
		  "AATACCT",
		  "CGAATGC",
		  "TTCGCAA",
		  "AATTCAA",
		  "CGCGCAG",
		  "AAGGTCT",
		  "ACTGGAC",
		  "AGCAGGT",
		  "GTACCGG",
		  "GGTCAAG",
		  "AATGATG",
		  "AGTCAGA",
		  "AACTAGA",
		  "CTATGGC",
		  "CGACGGT",
		  "AACCAAG",
		  "CGGCGTA",
		  "GCAGTCC",
		  "CTCGCGC",
		  "CTGCGAC",
		  "ACGTATG",
		  "ATACTGA",
		  "TACTTAG",
		  "AAGCTAA",
		  "GACGGCG",
		  "AGAAGAC",
		  "GTCCGGC",
		  "TCAGCTT",
		  "AGAGCGC",
		  "GCCTACG",
		  "TAATCAT",
		  "AACCTGC",
		  "GACGATT",
		  "TAGGCCG",
		  "GGCATAG",
		  "TTCAACC",
		  "TTAACTC",
		  "TAGTCTA",
		  "TGCATGA",
		  "AATAAGC",
		  "AGCCTTG",
		  "CCAACCT",
		  "GCAGAAG",
		  "AGAATTA",
		  "CAGCATC",
		  "TTCTAGG",
		  "CCTCTAG",
		  "CCGGATA",
		  "GCCGCCT",
		  "AACGACC",
		  "CCAGCGG",
		  "TAGTTCC",
		  "TGGCAAT",
		  "CGTATAT",
		  "GCTAATC",
		  "GACTTCT",
		  "GTACTAT",
		  "CGAGATC",
		  "CGCAGCC" );

$is4seq="AGATCTC";

function getGenomes(){
    global $genomedirectory;
    $myDirectory = opendir($genomedirectory);
    $filelist=array();
    
    while($entryName = readdir($myDirectory)) {
	#$item = $directoryToCheck."/".$entryName;
	if($entryName != "." and $entryName != ".."){
	    $filelist[] = $entryName; 
	}
    }   

    closedir($myDirectory);
    #$stringFound=$filelist[0]['name'];
    #$stringFound=substr($stringFound,1,strlen($stringFound)-3);
    return $filelist;
}

#variables
#echo "-------------";
#var_dump($_POST);
#echo "-------------";

if ( isset( $_POST["step"] ) and $_POST["step"] >= 1 and $_POST["step"]<= 9 ) {
    call_user_func( "processStep" . (int)$_POST["step"] );
} else {
    $runid=array_keys($_POST);
    $runid=$runid[0];
    displayStep1($runid);
}


function processStep1() {
    displayStep1();
}

function processStep2() {

    #if ( isset( $_POST["submitButton"] ) and $_POST["submitButton"] =="< Back" ) {
    displayStep2();
    #} else {
	#displayThanks();
    #}
}

function processStep3() {
    displayStep3();
}

function processStep4() {
    displayStep4();
}

function processStep5() {
    displayStep5();
}

function processStep6() {
    displayStep6();
}


function processStep7() {
    displayStep7();
}


function processStep8() {
    displayStep8();
}



function processStep9() {
    displayStep9();
}




function displayStep1($runid) {

    $illuminareaddir="/mnt/solexa/";

    echo "<h1>Analysis for ".$runid."</h1>";
    #gathering run information

    #echo $illuminareaddir;
    $rundirectory=$illuminareaddir."/".$runid;
    $runxml=$rundirectory."/RunInfo.xml";
    $runinformation=array("runid"        => $runid,
			  "cyclesread1"  => 0,
			  "cyclesread2"  => 0,
			  "cyclesindx1"  => 0,
			  "cyclesindx2"  => 0,
			  "LaneCount"    => 0, 
			  "SurfaceCount" => 0, 
			  "SwathCount"   => 0, 
			  "TileCount"    => 0);
    
    //////////////////////////////////////////
    //GETTING SOME INFO FROM THE RunInfo.xml
    //////////////////////////////////////////
    if(file_exists($runxml)){
	
	$xmlp = simplexml_load_file($runxml);
	//getting number of cycles
	foreach($xmlp->Run->Reads->Read as $readelem){
	    
	    if($readelem["IsIndexedRead"]  == "N"){
		
		if($runinformation["cyclesread1"] == 0){//first read
		    $runinformation["cyclesread1"]=(int)$readelem["NumCycles"];
		}else{//second
		    $runinformation["cyclesread2"]=(int)$readelem["NumCycles"];
		}
		
	    }else{
		if($readelem["IsIndexedRead"]  == "Y"){
		    
		    if($runinformation["cyclesindx1"] == 0){//first indx
			$runinformation["cyclesindx1"]=(int)$readelem["NumCycles"];
		    }else{//second
			$runinformation["cyclesindx2"]=(int)$readelem["NumCycles"];
		    }
		    
		}else{
		    echo "ERROR: cannot determine type of read in RunInfo.xml";
		    exit;	 
		}
	    }
	}
	
	//getting flowcell layout
	$runinformation["LaneCount"]    = (int)$xmlp->Run->FlowcellLayout["LaneCount"];
	$runinformation["SurfaceCount"] = (int)$xmlp->Run->FlowcellLayout["SurfaceCount"];
	$runinformation["SwathCount"]   = (int)$xmlp->Run->FlowcellLayout["SwathCount"];
	$runinformation["TileCount"]    = (int)$xmlp->Run->FlowcellLayout["TileCount"];
	
	

    }else{
	echo "Warning: The run has no RunInfo.xml file: ".$runxml;
    }

    #var_dump($runinformation);


    echo "<h3>Step 1: Run Information</h3>";
    echo "<form action=\"runprocess.php\" method=\"post\">";
    echo "<input type=\"hidden\" name=\"runid\" value=\"".$runid."\" />\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"2\" />";
    

    echo "<label for=\"email\">Your email</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"75\" name=\"email\"><br />";

    echo "<BR>Number of cycles:<BR>\n";
    echo "<label for=\"cyclesread1\">Cycles for read#1</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"cyclesread1\" value=\"".$runinformation["cyclesread1"]."\"><br />";

    echo "<label for=\"cyclesindx1\">Cycles for index#1</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"cyclesindx1\" value=\"".$runinformation["cyclesindx1"]."\"><br />";

    echo "<label for=\"cyclesread2\">Cycles for read#2</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"cyclesread2\" value=\"".$runinformation["cyclesread2"]."\"><br />";

    echo "<label for=\"cyclesindx2\">Cycles for index#2</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"cyclesindx2\" value=\"".$runinformation["cyclesindx2"]."\"><br />";

    echo "<BR>Flowcell Layout:<BR>\n";
    echo "<label for=\"LaneCount\">Lane count</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"LaneCount\" value=\"".$runinformation["LaneCount"]."\" readonly><br />";

    echo "<label for=\"SurfaceCount\">Lane count</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"SurfaceCount\" value=\"".$runinformation["SurfaceCount"]."\" readonly><br />";

    echo "<label for=\"SwathCount\">Lane count</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"SwathCount\" value=\"".$runinformation["SwathCount"]."\" readonly><br />";

    echo "<label for=\"TileCount\">Lane count</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"TileCount\" value=\"".$runinformation["TileCount"]."\" readonly><br />";
    
    #echo "<label for=\"TileCount\">Lane count</label>:\n";

    if($runinformation["LaneCount"] != 1){
	echo "<br>Which lanes  ? <br>\n";

	for ($i=1; $i<=$runinformation["LaneCount"]; $i++){
	    echo $i.":<input type=\"checkbox\" value=\"".$i."\" name=\"lanes".$i."\"><br/>\n";
	}
    }else{
	echo "<br>Which lanes  ? <br>\n";
	echo "1:<input type=\"checkbox\" value=\"1\" name=\"lanes1\" checked ><br/>\n";
    }

    echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";
    echo "</form>\n";

}

function displayStep2() {
    #var_dump($_POST);

    echo "<form action=\"runprocess.php\" method=\"post\">\n";

    //////////////////////////////////
    //BEGIN checking step 1 variables
    ////////////////////////////////
    //lanes

    $arrayfield=explode("_",$_POST["runid"]);    
    $seqtype="unknown";
    #echo "field 1#".$arrayfield[1]."#";
    $sequencers=array("M00518"     =>   "miseq",
		      "SN7001204"  =>   "hiseq");


    #echo  array_keys($sequencers);
    if(isset($sequencers[ $arrayfield[1] ])) {
	$seqtype=$sequencers[ $arrayfield[1] ];
    }else{

	#exit;
    }

    $lanestoanalyze=array();

    if(isset($_POST["lanes1"]) ){ array_push($lanestoanalyze,1);    }
    if(isset($_POST["lanes2"]) ){ array_push($lanestoanalyze,2);    }
    if(isset($_POST["lanes3"]) ){ array_push($lanestoanalyze,3);    }
    if(isset($_POST["lanes4"]) ){ array_push($lanestoanalyze,4);    }
    if(isset($_POST["lanes5"]) ){ array_push($lanestoanalyze,5);    }
    if(isset($_POST["lanes6"]) ){ array_push($lanestoanalyze,6);    }
    if(isset($_POST["lanes7"]) ){ array_push($lanestoanalyze,7);    }
    if(isset($_POST["lanes8"]) ){ array_push($lanestoanalyze,8);    }


    if(count($lanestoanalyze) == 0){
	echo "No lanes selected";
	exit;
    }

    if($_POST["cyclesread1"] == 0 &&
       $_POST["cyclesread1"] == 0 &&
       $_POST["cyclesindx1"] == 0 &&
       $_POST["cyclesindx2"] == 0 ){
	echo "Cycles are all null";
	exit;
    }
       
    if(!isset($_POST["email"]) ||
       !strstr($_POST["email"],"@") ){
	echo "Please enter a valid email";
	exit;
    }

    $runinformation=array("runid"        => $_POST["runid"],
			  "cyclesread1"  => $_POST["cyclesread1"],
			  "cyclesread2"  => $_POST["cyclesread2"],
			  "cyclesindx1"  => $_POST["cyclesindx1"],
			  "cyclesindx2"  => $_POST["cyclesindx2"],
			  "LaneCount"    => $_POST["LaneCount"],
			  "SurfaceCount" => $_POST["SurfaceCount"],
			  "SwathCount"   => $_POST["SwathCount"],
			  "TileCount"    => $_POST["TileCount"],
			  "email"        => $_POST["email"],
			  "sequencer"    => $seqtype,
			  "lanes"        => $lanestoanalyze);
    #echo "\nvar";
    #echo var_dump($runinformation);
    #echo "\nend var\n";
    //////////////////////////////////
    //END checking step 1 variables
    ////////////////////////////////




    #$arrayfield=explode("_",$entryName);
    #$seqtype="unknown";
    #if(isset($sequencers[ $arrayfield[1] ])) {
    #$seqtype=$sequencers[ $arrayfield[1] ];
    #}

    echo "<form action=\"runprocess.php\" method=\"post\">\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"3\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    #echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";
    echo "<h3>Step 2: Basecalling</h3>";

    if($seqtype == "hiseq" || $seqtype == "ga" ){
	echo "Bustard:<input type=\"checkbox\" value=\"bustard\" name=\"bustard\"\"><br/>\n";
	echo "freeIbis:<input type=\"checkbox\" value=\"freeibis\" name=\"freeibis\"\" checked><br/>\n";
    }else{
	echo "Bustard:<input type=\"checkbox\" value=\"bustard\" name=\"bustard\"\" checked><br/>\n";
	echo "freeIbis:<input type=\"checkbox\" value=\"freeibis\" name=\"freeibis\"\"><br/>\n";
    }
    echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";

    echo "</form>\n";

	
}

function displayStep3() {
    global $protocol2chimera;
    
    //////////////////////////////////
    //BEGIN checking step 2 variables
    ////////////////////////////////

    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));
    $runinformation["freeibis"] = isset($_POST["freeibis"]);
    $runinformation["bustard"]  = isset($_POST["bustard"]);
    


    //////////////////////////////////
    //END checking step 2 variables
    ////////////////////////////////

    echo "<h3>Step 3: Merging/trimming adapters</h3>";

    echo "<form action=\"runprocess.php\" method=\"post\">\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"4\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    #echo "#";
    #echo var_dump($runinformation);
    #echo var_dump($protocol2chimera);
    echo "Select the protocol used :<br/>\n";
    echo "<select name=\"protocol\">\n";
    foreach(array_keys($protocol2chimera) as $protocol){
	#echo "PROT ".var_dump($protocol)."<BR>\n";
	echo "<option value=".$protocol.">". $protocol2chimera[$protocol][0] ."</option>\n";
    }
    echo "</select>\n";
    echo "<BR><BR><img src=\"images/diagramOverlapseq.gif\" alt=\"merge diagram\"  height=\"500\" width=\"600\">\n";
    echo "<BR><BR>Merge paired reads if<BR>\n";
    echo "<input type=\"radio\" name=\"mergeoverlap\" value=\"False\" checked>after adapter trimming, they overlap completely<br>\n";
    echo "<input type=\"radio\" name=\"mergeoverlap\" value=\"True\" >The above plus if they show partial overlap (recommended for ancient DNA)<br>\n";

    echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";
    echo "</form>\n";
    #echo "#";   
    
    
}

function displayStep4() {
    global $protocol2chimera;

    ///////////////////////////////////
    //BEGIN checking step 3 variables//
    ///////////////////////////////////
    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));

    #var_dump($_POST);
    //protocol
    $runinformation["adapter1"]= $protocol2chimera[ $_POST["protocol"] ][1];
    $runinformation["adapter2"]= $protocol2chimera[ $_POST["protocol"] ][2];
    $runinformation["chimeras"]= $protocol2chimera[ $_POST["protocol"] ][3];
    $runinformation["protocol"]= $_POST["protocol"];
    $runinformation["mergeoverlap"]= ($_POST["mergeoverlap"]=="True");

    //merging
    
    //////////////////////////////////
    //END checking step 3 variables //
    //////////////////////////////////
    $keypossible=False;
    $key="";
    if($runinformation["protocol"] == "IllMultiplexCR1"){
	$key="ACTC";
	$keypossible=True;
    }elseif($runinformation["protocol"] == "IllMultiplexCR2"){
	$key="GTCT";
	$keypossible=True;
    }



    if($keypossible){
	echo "<h3>Step 4: Key and Quality flagging</h3>\n";
	echo "<form action=\"runprocess.php\" method=\"post\">\n";   
	echo "The library protocol allows for reading a key in the beginning of a read, specify the key:<BR>\n";
	if($runinformation["cyclesread2"] == 0 ){ //single-end
	    echo "<input type=\"text\" name=\"key1\" value=\"".$key."\" size=\"5\"> <BR>\n";
	    echo "<input type=\"hidden\" name=\"key2\" value=\"\" />\n";
	}else{
	    echo "<input type=\"text\" name=\"key1\" value=\"\"         size=\"5\"> <BR>\n";
	    echo "<input type=\"text\" name=\"key2\" value=\"".$key."\" size=\"5\"> <BR>\n";
	}
	echo "<BR><BR>";
    }else{
	echo "<h3>Step 4: Quality flagging</h3>";
	echo "<form action=\"runprocess.php\" method=\"post\">\n";
	echo "<input type=\"hidden\" name=\"key1\" value=\"\" />\n";
	echo "<input type=\"hidden\" name=\"key2\" value=\"\" />\n";
    }





    echo "<input type=\"hidden\" name=\"step\" value=\"5\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    echo "Basic filtering:<BR>\n";
    echo "Filter reads with a high likelihood of error:<input type=\"checkbox\" value=\"True\" name=\"filterseqlike\" checked ><br/>\n";
    #echo "Instead of filtering, trim the low quality bases at the 3' end instead<input type=\"checkbox\" value=\"1\" name=\"lanes1\" checked ><br/>\n";

    echo "Likelihood cutoff : <input type=\"text\" name=\"seqlikecutoff\" value=\"0.5\" size=\"5\"><BR><BR>\n";
    echo "Additional filtering:<BR>\n";
    echo "<input type=\"radio\" name=\"addfilters\" value=\"False\" checked>Do not use additional filters<br>\n";
    echo "<input type=\"radio\" name=\"addfilters\" value=\"entropy\">Apply sequence entropy [0.0-2.0]  filter at:  <input type=\"text\" name=\"entropycutoff\" value=\"0.85\" size=\"5\"> <BR>\n";
    echo "<input type=\"radio\" name=\"addfilters\" value=\"frequency\">Apply base frequency [0.0-1.0] filter at: <input type=\"text\" name=\"frequencycutoff\" value=\"0.1\" size=\"5\"> <BR>\n";    

    echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";

    echo "</form>\n";

}


function displayStep5() {
    #var_dump($_POST);
    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));

    ///////////////////////////////////
    //BEGIN checking step 4 variables//
    ///////////////////////////////////
    if((int)$_POST["seqlikecutoff"]<0 || (int)$_POST["seqlikecutoff"]>1){
	echo "Sequence likelihood cannot be less than 0 or more than 1";
	exit;
    }

    if($_POST["addfilters"] == "entropy"){
	if((int)$_POST["entropycutoff"]<0 || (int)$_POST["entropycutoff"]>2){
	    echo "Sequence entropy cutoff cannot be less than 0 or more than 2";
	    exit;
	}
    }

    if($_POST["addfilters"] == "frequency"){
	if((int)$_POST["frequencycutoff"]<0 || (int)$_POST["frequencycutoff"]>1){
	    echo "Sequence base frequency cutoff cannot be less than 0 or more than 1";
	    exit;
	}
    }

    $runinformation["key1"]    = $_POST["key1"];
    $runinformation["key2"]    = $_POST["key2"];

    $runinformation["filterseqlike"]    = ($_POST["filterseqlike"] == "True");
    $runinformation["seqlikecutoff"]   = $_POST["seqlikecutoff"];
    $runinformation["filterentropy"]   = ($_POST["addfilters"] == "entropy");
    $runinformation["filterfrequency"] = ($_POST["addfilters"] == "frequency");		    
    $runinformation["entropycutoff"]   = $_POST["entropycutoff"];
    $runinformation["frequencycutoff"] = $_POST["frequencycutoff"];
		    

    //////////////////////////////////
    //END checking step 4 variables //
    //////////////////////////////////


    echo "<h3>Step 5: Read group assignment</h3>";
    echo "<form action=\"runprocess.php\" method=\"post\">\n";

    echo "<input type=\"hidden\" name=\"step\" value=\"6\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
		    
    ?>
	Put your indices here : (see format below)<br />
	<br />
	<br />

	<textarea rows="20" cols="100" name="indextext" wrap="physical" placeholder="Paste your indices here"></textarea>
	<br>
	<input type="submit" name="submitButton" id="nextButton" value="Next &gt;" />
        <br>
	<input type="reset" value="Clear fields" />
	<br>	<br>

	format (tab separated): <br />
	double index: <br />	
	<pre>
	#index p7 p5
	RG1  23	4
	RG2  24	3
	RG3  25	2
	...
	</pre>
	single index: <br />	
	<pre>
        #index p7
	RG1  45
	RG2  46
	...
	</pre>



	<?php

    echo "</form>\n";

}



function displayStep6() {
    global $p7_block300;
    global $p7_block;
    global $p5_block;


   # echo "test6";
    

    ///////////////////////////////////
    //BEGIN checking step 5 variables//
    ///////////////////////////////////
    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));
    $indextext  = $_POST["indextext"];
    $stringToPrint="";

    $indexOfLines=0;
    $foundControl=0;

    if(!$indextext ){
	echo "ERROR: please enter your indices in the text field";
	exit;	 
    }

    
    #CHECKING EACH LINE
    foreach(explode("\n",$indextext) as $line){
	
	$line=trim($line);
	#skip empty lines
	if(strlen($line) == 0){
	    continue 1;
	}

	#$arrayfield=explode("\t",$line);
	$arrayfield=preg_split('/\s+/', $line);



#	if($indextype       == "single"){
	if($runinformation["cyclesread2"] == 0 ){
	    if(count($arrayfield) != 2){
		echo "ERROR: For single index, lines must have 2 fields check line \"".$line."\"";
		exit;
	    }
	}else{
	    if(count($arrayfield) != 3){
		echo "ERROR: For double index, lines must have 3 fields check line \"".$line."\"";
		exit;
	    }
	}


	if($indexOfLines == 0){
	    #CHECKING HEADER
	    if($arrayfield[0] != "#index"){
		echo "ERROR: the first field of the header must be #index (case sensitive)";
		exit;
	    }

	    if($arrayfield[1] != "p7"){
		echo "ERROR: the second field of the header must be p7 (case sensitive)";
		exit;
	    }

	    #if($indextype  == "double"){
	    if($runinformation["cyclesread2"] != 0 ){

		if($arrayfield[2] != "p5"){
		    echo "ERROR: the first field of the header must be p5 (case sensitive)";
		    exit;
		}
	    }

	}else{
	    #CHECKING REMAINING FIELDS
	    #if($indextype       == "single"){

	    if(strstr($arrayfield[0],"\"") ||
	       strstr($arrayfield[0],"'")  || 
	       strstr($arrayfield[0],"\\")  || 
	       strstr($arrayfield[0],"/")   ){	       
		echo "ERROR: The first field cannot have quotes or (back)slashes characters for line ".$line;
		exit;
	    }

	    if($runinformation["cyclesread2"] == 0 ){

		if( ($arrayfield[1] > count($p7_block)) &&
		    ( ($arrayfield[1]<301) || ($arrayfield[1]>(300+count($p7_block300)) ) )  ){
		    echo "ERROR: index for p7 is not within the expected bound for ".$line . "";
		    exit;	    
		}

		if( ($arrayfield[1]>=301) && ($arrayfield[1]<=(300+count($p7_block300)) )   ){
		    $stringToPrint.=$p7_block300[ $arrayfield[1] - 301] ."\t".$arrayfield[0]."\n"; #-1 for zero base array
		}else{
		    $stringToPrint.=$p7_block[ $arrayfield[1] - 1] ."\t".$arrayfield[0]."\n";
		}
	    }else{


		if( ($arrayfield[1] > count($p7_block)) &&
		    ( ($arrayfield[1]<301) || ($arrayfield[1]>(300+count($p7_block300)) ) )  ){
		    echo "ERROR: index for p7 is not within the expected bound for ".$line . "";
		    exit;	    
		}

		if($arrayfield[2] > count($p5_block)){
		    echo "ERROR: index for p5 is greater than size of the array in line ".$line;
		    exit;	    
		}

		if( ($arrayfield[1]>=301) && ($arrayfield[1]<=(300+count($p7_block300)) )   ){
		    $stringToPrint.=$p7_block300[ $arrayfield[1] - 301] ."\t".$p5_block[ $arrayfield[2] - 1]."\t".$arrayfield[0]."\n";
		}else{
		    $stringToPrint.=$p7_block[ $arrayfield[1] - 1] ."\t".$p5_block[ $arrayfield[2] - 1]."\t".$arrayfield[0]."\n";
		}
	    }

	    if($arrayfield[0] == "control"){
		$foundControl=1;
	    }
	}
	$indexOfLines++;
    } #for each field explode

    if($foundControl == 0){
	#if($indextype       == "single"){
	if($runinformation["cyclesread2"] == 0 ){
	    $stringToPrint.="TTGCCGC\tcontrol\n";
	}else{
	    $stringToPrint.="TTGCCGC\tAGATCTC\tcontrol\n";
	}
    }


    //////////////////////////////////
    //END checking step 5 variables //
    //////////////////////////////////
   # $runinformation["indicesseq"]=$stringToPrint;
   # $runinformation["indicesraw"]= $_POST["indextext"];

    #var_dump($runinformation);
    echo "<h3>Step 6: Verify indices</h3>";
    echo "<form action=\"runprocess.php\" method=\"post\">\n";

    echo "<input type=\"hidden\" name=\"step\" value=\"7\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    echo "<input type=\"hidden\" name=\"testindexorig\"   value=\"$indextext\" />\n";
    echo "The following indices will be used:<BR>\n";
    echo "<textarea readonly rows=\"20\" cols=\"100\" name=\"textindex\" wrap=\"physical\">".$stringToPrint."</textarea><br />\n";
    echo "<input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";
    #echo "<button onclick=\"history.go(-1);\">cancel</button>\n";
    

    echo "</form>\n";

}


function displayStep7() {
    ////////////////////////////////////
    //BEGIN checking step 6 variables //
    ////////////////////////////////////

    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));


    $runinformation["indicesseq"]= $_POST["textindex"];
    $runinformation["indicesraw"]= $_POST["testindexorig"];

    #var_dump($runinformation);


    //////////////////////////////////
    //END checking step 6 variables //
    //////////////////////////////////
    echo "<h3>Step 7: Mapping</h3>";
    echo "<form action=\"runprocess.php\" method=\"post\">\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"8\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    #echo "Mapping using BWA<BR>\n";
    echo "Mapping using BWA:  <input type=\"checkbox\" value=\"True\" name=\"usebwa\" checked ><br/>\n";
    $arrayofGenomes=getGenomes();
    natcasesort($arrayofGenomes);
    echo "Selecting the target genome:<BR><BR>\n";
    echo "<select name=\"genomebwa\" size=\"1\">\n";
    foreach($arrayofGenomes as $agenome){
	if($agenome == "hg19_1000g"){
	    echo "<option value=\"".$agenome."\" selected=\"selected\">".$agenome."</option>\n";
	}else{
	    echo "<option value=\"".$agenome."\" >".$agenome."</option>\n";
	}
    }
    echo "</select><BR><BR>\n";
    echo "Selecting the parameters for BWA:<BR><BR>\n";
    echo "<select name=\"parambwa\" size=\"1\">\n";
    echo "<option value=\"default\" >Default parameters (modern DNA)</option>\n";
    echo "<option value=\"ancient\" >Ancient parameters (ancient DNA)</option>\n";
    echo "</select><BR><BR>\n";
    echo "<input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";

    echo "</form>\n";
}


function displayStep8() {
    global $outputdirectory;
    ////////////////////////////////////
    //BEGIN checking step 7 variables //
    ////////////////////////////////////



    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));
    $runinformation["usebwa"]     = $_POST["usebwa"];
    $runinformation["genomebwa"]  = $_POST["genomebwa"];
    $runinformation["parambwa"]   = $_POST["parambwa"];


    //////////////////////////////////
    //END checking step 7 variables //
    //////////////////////////////////

    echo "<h3>Step 8: Summary</h3>";

    #echo "Please review the following information prior to pressing submit:<BR>\n";
    echo "<form action=\"runprocess.php\" method=\"post\">\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"9\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
#    var_dump($runinformation);
    echo "<BR>Please review the following information prior to submitting (submit buttom at the bottom of the page):<BR>\n";
    echo "<table   border=0>\n";
    echo "<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap> </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>General:      </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>Run ID      :</TD><TD> ".$runinformation["runid"]."</TD></TR>\n";
    echo "<TR><TD nowrap>Sequencer      :</TD><TD> ".$runinformation["sequencer"]."</TD></TR>\n";
    echo "<TR><TD nowrap>Your email      :</TD><TD> ".$runinformation["email"]."</TD></TR>\n";
    echo "<TR><TD nowrap>Cycle for read1 :</TD><TD> ".$runinformation["cyclesread1"]."</TD></TR>\n";
    echo "<TR><TD nowrap>Cycle for read2 :</TD><TD> ".$runinformation["cyclesread2"]."</TD></TR>\n";
    echo "<TR><TD nowrap>Cycle for index 1 :</TD><TD> ".$runinformation["cyclesindx1"]."</TD></TR>\n";
    echo "<TR><TD nowrap>Cycle for index 2 :</TD><TD> ".$runinformation["cyclesindx2"]."</TD></TR>\n";
    echo "<TR><TD nowrap>Lanes to analyze:</TD><TD> ".implode(",",$runinformation["lanes"])."</TD></TR>\n";
# echo "<TR><TD nowrap>Lanes to analyze:</TD><TD> ".implode(",",$runinformation["lanes"])."</TD></TR>\n";
    echo "<TR><TD nowrap> </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>Basecalling:      </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>Basecalling using Bustard      :</TD><TD> ".($runinformation["bustard"]?"yes":"no")."</TD></TR>\n";
    echo "<TR><TD nowrap>Basecalling using freeIbis     :</TD><TD> ".($runinformation["freeibis"]?"yes":"no")."</TD></TR>\n";
    echo "<TR><TD nowrap> </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>Merging/trimming:      </TD></TR>\n";
    echo "<TR><TD nowrap>Merge partially overlapping sequencing     :</TD><TD> ".($runinformation["mergeoverlap"]=="True"?"yes":"no")."</TD></TR>\n";
    echo "<TR><TD nowrap>Adapter 1     :</TD><TD> ".($runinformation["adapter1"])."</TD></TR>\n";
    echo "<TR><TD nowrap>Adapter 2     :</TD><TD> ".($runinformation["adapter2"])."</TD></TR>\n";
    echo "<TR><TD nowrap>Potential chimeras     :</TD><TD> ".($runinformation["chimeras"])."</TD></TR>\n";
    echo "<TR><TD nowrap>Protocol       :</TD><TD> ".($runinformation["protocol"])."</TD></TR>\n";
    echo "<TR><TD nowrap>Key read#1     :</TD><TD> ".(strlen($runinformation["key1"])==0?"none":$runinformation["key1"])."</TD></TR>\n";
    echo "<TR><TD nowrap>Key read#2     :</TD><TD> ".(strlen($runinformation["key2"])==0?"none":$runinformation["key2"])."</TD></TR>\n";
    echo "<TR><TD nowrap> </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>Quality filtering:      </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>Flag sequences with low likelihood  :</TD><TD> ".(($runinformation["filterseqlike"]=="1")?"yes":"no")."</TD></TR>\n";
    echo "<TR><TD nowrap>Likelihood cutoff                  :</TD><TD> ".($runinformation["seqlikecutoff"])."</TD></TR>\n";

    echo "<TR><TD nowrap>Flag sequences based on  entropy  :</TD><TD> ".(($runinformation["filterentropy"])?"yes":"no")."</TD></TR>\n";
    echo "<TR><TD nowrap>Entropy cutoff  :</TD><TD> ".($runinformation["entropycutoff"])."</TD></TR>\n";

    echo "<TR><TD nowrap>Flag sequences using frequency  :</TD><TD> ".(($runinformation["filterfrequency"])?"yes":"no")."</TD></TR>\n";
    echo "<TR><TD nowrap>Frequency cutoff  :</TD><TD> ".($runinformation["frequencycutoff"])."</TD></TR>\n";
    echo "<TR><TD nowrap> </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>BWA mapping:     </TD><TD></TD></TR>\n";
    echo "<TR><TD nowrap>Map using BWA  :</TD><TD> ".(($runinformation["usebwa"]=="True")?"yes":"no")."</TD></TR>\n";
    echo "<TR><TD nowrap>Genome to use  :</TD><TD> ".($runinformation["genomebwa"])."</TD></TR>\n";
    echo "<TR><TD nowrap>BWA parameters  :</TD><TD> ".($runinformation["parambwa"])."</TD></TR>\n";
    // echo "<TR><TD nowrap>Indices:     </TD><TD></TD></TR>\n";
    echo "</table>\n";

    echo "<BR>Indices to use  :<BR><PRE>\n".($runinformation["indicesseq"])."</PRE><BR>\n";
    echo "<BR>Indices raw     :<BR><PRE>\n".($runinformation["indicesraw"])."</PRE><BR>\n";

    #echo json_encode($runinformation);
    $runid=$runinformation["runid"];
    $targetfile=$outputdirectory."/".$runid."_".implode(",",$runinformation["lanes"]).".json";
    if(file_exists ( $targetfile )){
       echo "<font color=red size=+2>Warning:</font> file ".$targetfile." already exists <br><br>";
    }
    echo "<input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Submit\" />\n";

    echo "</form>\n";

}


function displayStep9() {
    global $outputdirectory;
    global $emailAddrToSend;


    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));
    $runid=$runinformation["runid"];
    $emailuser=$runinformation["email"];

    $targetfile=$outputdirectory."/".$runid."_".implode(",",$runinformation["lanes"]).".json";
    echo "printing to ".$targetfile."\n";
    echo "printing to ".json_encode($runinformation)."\n";
    $stringtoprint=json_encode($runinformation);

    $fh = fopen($targetfile, 'w') or die("can't open $targetfile");
    fwrite($fh, $stringtoprint);
    fclose($fh);

    
    
    $mail = new EMail;
    $mail->Username = 'sbsuser';
    $mail->Password = 'sbs123';
    
    $mail->SetFrom("sbsuser@eva.mpg.de","");  // Name is optional

    foreach($emailAddrToSend as $emailAddrTo){ //for each email in the array above
	$mail->AddTo($emailAddrTo,""); // Name is optional
    }

    $mail->Subject = "Analysis submitted for run ".$runid;
    $mail->Message = "User email ".$emailuser."\n".
	 "Created the following file : ".$targetfile."\n".
	 "with the following contents:\n------------------\n".json_encode($runinformation);
    $mail->ConnectTimeout = 30;  // Socket connect timeout (sec)
    $mail->ResponseTimeout = 8;  // CMD response timeout (sec)
    $success = $mail->Send();
    if($success != 1){
	$stringforscreen.="the admin notification email WAS NOT sent successfully, contact directly ".implode(",",$emailAddrToSend)."<br>\n";
    }else{
	$stringforscreen.="the admin notification email was sent successfully <br>\n";
    }
    

    if(0){

    }
    
    

	?>
    <h1>Thank You</h1>
	<p>Thank you, your application has been received.</p>

	<?php
	}
?>
</body>
</html>