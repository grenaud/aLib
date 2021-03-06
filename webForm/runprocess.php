<html><head>
<title>Illumina Form</title>
<body>
<?php

require "email.php";

ini_set('display_errors', 'On');
error_reporting(E_ALL);


if( ! file_exists(getcwd()."/config.json")  ){
    echo "Configuration file not found ".getcwd()."/config.json";
    exit(1);
 }

if( ! file_exists(getcwd()."/json2make.py")  ){
    echo "Python script does not exist: ".getcwd()."/json2make.py";
    exit(1);
}

$json2makePath=getcwd()."/json2make.py";

/* jsonconf = simplexml_load_file( getcwd()."/config.xml" ); */
$jsonconf = json_decode(file_get_contents( getcwd()."/config.json" ),true);

$emailAddrToSend  = $jsonconf["emailAddrToSend"];
$genomedirectory  = $jsonconf["genomedirectory"];
$illuminawritedir = $jsonconf["illuminawritedir"];
if(isset($jsonconf["couchdb"]) ){
    $couchdb=$jsonconf["couchdb"];
}

$protocol2chimera = array();
foreach($jsonconf["chimeras"]["chimera"] as $chimlem){
    $protocol2chimera[ (string)$chimlem["protocol"] ] = array((string)$chimlem["name"],
							      (string)$chimlem["adapter1"],
							      (string)$chimlem["adapter2"],
							      (string)$chimlem["chimera"] );
}
$ctrlindex=$jsonconf["controlindex"];
$ctrlindex2=$jsonconf["controlindex2"];

//first key is indexing scheme, second is number
$p7Indices=array();
$p5Indices=array();
/* $indexSchemes=array(); */


foreach($jsonconf["indices"]["p7indices"]["p7index"] as $p7ind){
    $p7Indices[ (string)$p7ind["id"] ] = (string)$p7ind["seq"];
}

foreach($jsonconf["indices"]["p5indices"]["p5index"] as $p5ind){
    $p5Indices[ (string)$p5ind["id"] ] = (string)$p5ind["seq"];
}

/* exit; */

$sequencers=array();

foreach($jsonconf["sequencers"]["sequencer"] as $seqelem){
    /* print "#".(string)$seqelem["id"]."#<BR>"; */
    /* print "#".(string)$seqelem["type"]."#<BR>"; */
    $sequencers[ (string)$seqelem["id"] ] = (string)$seqelem["type"];
}



function getGenomes(){
    global $genomedirectory;
    $myDirectory = opendir($genomedirectory);
    $filelist=array();
    
    while($entryName = readdir($myDirectory)) {
	//$item = $directoryToCheck."/".$entryName;
	if($entryName != "." and $entryName != ".."){
	    $filelist[] = $entryName; 
	}
    }   

    closedir($myDirectory);
    #$stringFound=$filelist[0]['name'];
    #$stringFound=substr($stringFound,1,strlen($stringFound)-3);
    return $filelist;
}

/* variables */
/* echo "-------------"; */
/* var_dump($_POST); */
/* echo "-------------"; */

if ( isset( $_POST["step"] ) and $_POST["step"] >= 1 and $_POST["step"]<= 9 ) {
    call_user_func( "processStep" . (int)$_POST["step"] );
} else {
    $runid=array_keys($_POST);
    $runid=$runid[0];
    if(date('jF') == "1April"){

?>
<!-- Progress bar holder -->
<div id="progress" style="width:500px;border:1px solid #ccc;"></div>
<!-- Progress information -->
<div id="information" style="width"></div>
<?
	    
	for($percent=1;$percent<=100;$percent++){
	    $numberFound=100;
	    $percentval = $percent."%";

	    echo '<script language="javascript">
    document.getElementById("progress").innerHTML="<div style=\"width:'.$percentval.';background-color:#ddd;\">&nbsp;</div>";
    document.getElementById("information").innerHTML="Deleting raw data from '.$runid.' to save space,  deleted '.($percent*10).' files so far";
    </script>';
	    echo str_repeat(' ',1024*64);
	    flush();
	    usleep(100000);
	    //sleep(0.1);
	}
	usleep(100000);
	echo "April fool's :-), here is the form:<BR>";
	usleep(200000);
	displayStep1($runid);
    }else{
	//exit;
	displayStep1($runid);
    }
    
}


function processStep1() {
    displayStep1();
}

function processStep2() {

    //if ( isset( $_POST["submitButton"] ) and $_POST["submitButton"] =="< Back" ) {
    displayStep2();
    /* } else { */
    /* 	displayThanks(); */
    /* } */
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
    //gathering run information

    //echo $illuminareaddir;
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
    $cannotGetRunInfo=False;
    //////////////////////////////////////////
    //GETTING SOME INFO FROM THE RunInfo.xml
    //////////////////////////////////////////
    if(file_exists($runxml)){
	
	$xmlp = simplexml_load_file($runxml);
	//getting number of cycles
	foreach($xmlp->Run->Reads->Read as $readelem){
	    
	    if($readelem["IsIndexedRead"]  == "N"){ //non index
		
		if($runinformation["cyclesread1"] == 0){//first read
		    $runinformation["cyclesread1"]=(int)$readelem["NumCycles"];
		}else{//second
		    $runinformation["cyclesread2"]=(int)$readelem["NumCycles"];
		}
		
	    }else{
		if($readelem["IsIndexedRead"]  == "Y"){//index
		    
		    if($runinformation["cyclesindx1"] == 0){//first indx
			$runinformation["cyclesindx1"]=(int)$readelem["NumCycles"];
		    }else{//second
			$runinformation["cyclesindx2"]=(int)$readelem["NumCycles"];
		    }
		    
		}else{
		    /* echo $readelem["FirstCycle"]."<BR>\n"; */

		    /* foreach($readelem->Index as $testindex){ */
		    /* 	echo $testindex."<BR>\n"; */
		    /* } */
		    /* echo "ERROR: cannot determine type of read in RunInfo.xml"; */
		    $cannotGetRunInfo=True;
		    /* exit(1);	  */
		}
	    }
	}
	
	//getting flowcell layout
	$runinformation["LaneCount"]    = (int)$xmlp->Run->FlowcellLayout["LaneCount"];
	$runinformation["SurfaceCount"] = (int)$xmlp->Run->FlowcellLayout["SurfaceCount"];
	$runinformation["SwathCount"]   = (int)$xmlp->Run->FlowcellLayout["SwathCount"];
	$runinformation["TileCount"]    = (int)$xmlp->Run->FlowcellLayout["TileCount"];
	if($cannotGetRunInfo){
	    echo "<font color=red size=+2>Warning:</font> cannot get run information from the RunInfo.xml, enter the information in the fields below <br><br>"; 
	    $runinformation["LaneCount"]    = 8;

	}

	

    }else{
	echo "Warning: The run has no RunInfo.xml file: ".$runxml;
    }

    //var_dump($runinformation);
    $runinformation["expname"]    = implode("_",array_slice(explode("_",$runid),1)); 



    echo "<h3>Step 1: Run Information</h3>";


    echo "This series of steps allows you to fill information about the basic computational processing of sequencing data. Please note that the processing will be identical for the lanes that you select. If you have multiple lanes on a single run and you want different processing (lanes have different read groups for example), please fill the form once for each set of lanes for which you want identical processing.<BR>";
    echo "<form action=\"runprocess.php\" method=\"post\">";
    echo "<input type=\"hidden\" name=\"runid\" value=\"".$runid."\" />\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"2\" />";
    

    echo "<label for=\"email\">Your email (put commas if multiple emails)</label>:\n";
    echo "<input type=\"text\" size=\"80\" maxlength=\"500\" name=\"email\"><br />";

    echo "<BR>Verify number of cycles: (Please note these values are used for the analysis hence, put the values observed not the ones that were planned i.e if a run had 2x76bp+double indices of 7bp but the machine was stopped after the first index, enter 0 for both the second index and read)<BR>\n";
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
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"LaneCount\" value=\"".$runinformation["LaneCount"]."\"";
    if(!$cannotGetRunInfo){
	echo  " readonly";
    }
    echo "><br />";

    echo "<label for=\"SurfaceCount\">Surface count</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"SurfaceCount\" value=\"".$runinformation["SurfaceCount"]."\"";
    if(!$cannotGetRunInfo){
	echo  " readonly";
    }
    echo "><br />";

    echo "<label for=\"SwathCount\">Swath count</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"SwathCount\" value=\"".$runinformation["SwathCount"]."\" ";
    if(!$cannotGetRunInfo){
	echo  " readonly";
    }
    echo "><br />";


    echo "<label for=\"TileCount\">Tile count</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"TileCount\" value=\"".$runinformation["TileCount"]."\" ";
    if(!$cannotGetRunInfo){
	echo  " readonly";
    }
    echo "><br />";

    echo "<label for=\"Experiment name\">Experiment name</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"100\" name=\"ExperimentName\" value=\"".$runinformation["expname"]."\" ";
    if(!$cannotGetRunInfo){
	echo  " readonly";
    }
    echo "><br />";


    if(!$cannotGetRunInfo){
	echo "<BR><b> Comment:</b>If you have an old Genome Analyser run, enter -1 for both the swath and tile and surface count.";

    }
    
    //echo "<label for=\"TileCount\">Lane count</label>:\n";

    if($runinformation["LaneCount"] != 1){
    echo "<br>Which lanes belong to you ? <br>\n";

    for ($i=1; $i<=$runinformation["LaneCount"]; $i++){
	    echo $i.":<input type=\"checkbox\" value=\"".$i."\" name=\"lanes".$i."\"><br/>\n";
	}
    }else{
	echo "<br>Which lane(s) belong to you ? <br>\n";
	echo "1:<input type=\"checkbox\" value=\"1\" name=\"lanes1\" checked ><br/>\n";
    }

    echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";
    echo "</form>\n";

}





function displayStep2() {
    global $ctrlindex;
    global $ctrlindex2;

    global $sequencers;

    echo "<form action=\"runprocess.php\" method=\"post\">\n";

    //////////////////////////////////
    //BEGIN checking step 1 variables
    ////////////////////////////////
    //lanes

    $arrayfield=explode("_",$_POST["runid"]);    
    $seqtype="unknown";


    /* echo  array_keys($sequencers); */
    if(isset($sequencers[ $arrayfield[1] ])) {
	$seqtype=$sequencers[ $arrayfield[1] ];
    }else{
	
	//exit;
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
	exit(1);
    }

    if($_POST["cyclesread1"] == 0 &&
       $_POST["cyclesread1"] == 0 &&
       $_POST["cyclesindx1"] == 0 &&
       $_POST["cyclesindx2"] == 0 ){
	echo "Cycles are all null";
	exit(1);
    }
       
    if(!isset($_POST["email"]) ||
       !strstr($_POST["email"],"@") ){
	echo "Please enter a valid email";
	exit(1);
    }
      

    $runinformation=array("runid"        => $_POST["runid"],
			  "expname"      => $_POST["ExperimentName"],
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
			  "lanes"        => $lanestoanalyze,
			  );

    //padding $ctrlindex and $ctrlindex2 for indices with 8 bp
    if($_POST["cyclesindx1"]==8){
	$ctrlindex .= "A";
    }

    if($_POST["cyclesindx2"]==8){
	$ctrlindex2 .= "A";
    }



    $runinformation["email"] = str_replace(' ','',$runinformation["email"]);

    /* echo "\nvar"; */
    /* echo var_dump($runinformation); */
    /* echo "\nend var\n"; */
    //////////////////////////////////
    //END checking step 1 variables
    ////////////////////////////////




    /* $arrayfield=explode("_",$entryName); */
    /* $seqtype="unknown"; */
    /* if(isset($sequencers[ $arrayfield[1] ])) { */
    /* $seqtype=$sequencers[ $arrayfield[1] ]; */
    /* } */

    echo "<form action=\"runprocess.php\" method=\"post\">\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"3\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    //echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";
    echo "<h3>Step 2: Basecalling (run id: ".$runinformation["runid"].") </h3><br>\n";   
    echo "<BR>Basecalling is the process that converts the raw intensities into sequences. The default Illumina basecaller is Bustard and our custom basecaller is freeIbis. freeIbis is more accurate in terms of sequence and quality score but required control sequences. We recommend using Bustard for small MiSeqs and freeIbis for GAs or HiSeqs<BR>";

    echo "<label for=\"LaneCount\">Detected platform</label>:\n";
    echo "<input type=\"text\" size=\"12\" maxlength=\"4\" name=\"platfType\" value=\"".$seqtype."\" readonly><br><br>";

    if($seqtype == "hiseq" || $seqtype == "ga" ){
	echo "Bustard:<input type=\"checkbox\" value=\"bustard\" name=\"bustard\"\"><br/>\n";
	echo "freeIbis:<input type=\"checkbox\" value=\"freeibis\" name=\"freeibis\"\" checked><br/>\n";
    }else{
	echo "Bustard:<input type=\"checkbox\" value=\"bustard\" name=\"bustard\"\" checked><br/>\n";
	echo "freeIbis:<input type=\"checkbox\" value=\"freeibis\" name=\"freeibis\"\"><br/>\n";
    }

    echo "Most runs spike-in PhiX DNA (for which we know the sequence in advance) used a control sequences to determine how successful a run was. freeIbis works by training on those control sequences<BR>"; 
    echo "<BR>If you picked freeIbis, how were control sequence (phiX) specified ?<BR>";



    /* echo "<input type=\"checkbox\" value=\"bustard\" name=\"bustard\"\"><br/>\n"; */
    /* echo "<input type=\"checkbox\" value=\"freeibis\" name=\"freeibis\"\" checked><br/>\n"; */
    
    echo "<input type=\"radio\" name=\"spikedin\" value=\"True\" checked>Spiked-in controls using P7 index:   <input type=\"text\" name=\"ctrlindex\" value=\"$ctrlindex\" size=\"".$_POST["cyclesindx1"]."\"> and P5 index:   <input type=\"text\" name=\"ctrlindex2\" value=\"$ctrlindex2\" size=\"".$_POST["cyclesindx2"]."\">  <BR>\n";
    echo "<input type=\"radio\" name=\"spikedin\" value=\"Mult\" >Spiked-in controls using the multiple spiked-in option<BR>\n";
    echo "<input type=\"radio\" name=\"spikedin\" value=\"False\">Dedicated lane (specify which below)\n";

    echo "<p style=\"padding-left:5em;\">";

    if($runinformation["LaneCount"] != 1){
	echo "If so, which lanes were used exclusively for phiX ? <br>\n";

	for ($i=1; $i<=$runinformation["LaneCount"]; $i++){
	    echo $i.":<input type=\"checkbox\" value=\"".$i."\" name=\"lanesdedicated".$i."\"><br/>\n";
	}
    }else{
	echo "<br>Which lanes were used exclusively for phiX    ? <br>\n";
	echo "1:<input type=\"checkbox\" value=\"1\" name=\"lanesdedicated1\" ><br/>\n";
    }

    echo "</p>";

    /* echo "<b>Number of separate processings</b><br/>\n"; */
    /* echo "Select number of different processings you wish to have:<br/>\n"; */
    /* echo "Use this option if different types of libraries were pooled on one lane<BR>e.g. ancient plus modern DNA on a single lane:<br/>\n"; */

    /* echo "<select name=\"numprocessing\">\n"; */
    /* for ($procnu = 1; $procnu <= 4; $procnu++) { */
    /* 	//echo "PROT ".var_dump($protocol)."<BR>\n"; */
    /* 	echo "<option value=".$procnu.">". $procnu ."</option>\n"; */
    /* } */
    /* echo "</select>\n"; */



    echo "<br><br><br><br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";

    echo "</form>\n";

	
}

function displayStep3() {
    global $protocol2chimera;
    global $illuminawritedir;

    //var_dump($_POST);
    //////////////////////////////////
    //BEGIN checking step 2 variables
    ////////////////////////////////



    if(isset($_POST["lanesdedicated1"]) && $_POST["spikedin"] == "True" ){ 	echo "Do not specify lane if the control was spiked in"; exit; }
    if(isset($_POST["lanesdedicated2"]) && $_POST["spikedin"] == "True" ){ 	echo "Do not specify lane if the control was spiked in"; exit; }
    if(isset($_POST["lanesdedicated3"]) && $_POST["spikedin"] == "True" ){ 	echo "Do not specify lane if the control was spiked in"; exit; }
    if(isset($_POST["lanesdedicated4"]) && $_POST["spikedin"] == "True" ){ 	echo "Do not specify lane if the control was spiked in"; exit; }
    if(isset($_POST["lanesdedicated5"]) && $_POST["spikedin"] == "True" ){ 	echo "Do not specify lane if the control was spiked in"; exit; }
    if(isset($_POST["lanesdedicated6"]) && $_POST["spikedin"] == "True" ){ 	echo "Do not specify lane if the control was spiked in"; exit; }
    if(isset($_POST["lanesdedicated7"]) && $_POST["spikedin"] == "True" ){ 	echo "Do not specify lane if the control was spiked in"; exit; }
    if(isset($_POST["lanesdedicated8"]) && $_POST["spikedin"] == "True" ){ 	echo "Do not specify lane if the control was spiked in"; exit; }

    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));
    $runinformation["freeibis"]       = isset($_POST["freeibis"]);
    $runinformation["bustard"]        = isset($_POST["bustard"]);
    $runinformation["spikedin"]       = ($_POST["spikedin"] == "True");
    $runinformation["spikedinmult"]   = False;

    if($runinformation["spikedin"]){ //if spiked-in with single index PhiX
	$runinformation["ctrlindex"]      = $_POST["ctrlindex"];
	$runinformation["ctrlindex2"]     = $_POST["ctrlindex2"];
	$runinformation["spikedinmult"]  = False;
    }else{
	if($_POST["spikedin"] == "Mult"){//multiple
	    $runinformation["spikedinmult"]  = True;
	}

    }
    /* $runinformation["numprocessing"]  = $_POST["numprocessing"]; */


    //Detecting # of processing
    for ($procnu = 1; $procnu <= 100 ; $procnu++) {
	//$targetfile=$illuminawritedir."/".$runinformation["runid"]."_".implode(",",$runinformation["lanes"])."-".$procnu.".json";    
	$targetfile=$illuminawritedir."/".$runinformation["runid"]."/build/".$runinformation["runid"]."_".implode(",",$runinformation["lanes"])."-".$procnu.".json";
	/* print $targetfile; */
	/* if(file_exists($targetfile)){ */
	/*     print "yes"; */
	/* }else{ */
	/*     print "no"; */
	/* } */
	//	exit;
	if(!file_exists ( $targetfile )){
	    $runinformation["numprocessingcurrent"] = $procnu;
	    break;
	}else{
	    /* print $procnu; */
	    /* exit; */
	}
    }


    $lanesdedicated=array();

    if(isset($_POST["lanesdedicated1"]) ){ array_push($lanesdedicated,1); }
    if(isset($_POST["lanesdedicated2"]) ){ array_push($lanesdedicated,2); }
    if(isset($_POST["lanesdedicated3"]) ){ array_push($lanesdedicated,3); }
    if(isset($_POST["lanesdedicated4"]) ){ array_push($lanesdedicated,4); }
    if(isset($_POST["lanesdedicated5"]) ){ array_push($lanesdedicated,5); }
    if(isset($_POST["lanesdedicated6"]) ){ array_push($lanesdedicated,6); }
    if(isset($_POST["lanesdedicated7"]) ){ array_push($lanesdedicated,7); }
    if(isset($_POST["lanesdedicated8"]) ){ array_push($lanesdedicated,8); }

    $runinformation["lanesdedicated"]  = $lanesdedicated;

    if(!$runinformation["spikedin"]     &&
       !$runinformation["spikedinmult"] &&
       count($lanesdedicated) == 0 ){
	echo "Error, you must specify which lanes were used for controls";
        exit(1);
    } 

    //////////////////////////////////
    //END checking step 2 variables
    ////////////////////////////////

    echo "<h3>Step 3: Merging/trimming adapters (run id: ".$runinformation["runid"].") </h3><br>\n";   
    
    

    // for ($procnu = 1; $procnu <= $runinformation["numprocessing"]; $procnu++) {
    echo "<form action=\"runprocess.php\" method=\"post\">\n";
    echo "<BR><b>Please note, you are filling the form for processing #".$runinformation["numprocessingcurrent"].":</b><BR><BR>";
    echo "<input type=\"hidden\" name=\"step\" value=\"4\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    /* echo "#"; */
    /* echo var_dump($runinformation); */
    /* echo var_dump($protocol2chimera); */
    echo "<b>Select the protocol used</b>:<br/>\n";
    echo "<select name=\"protocol".$runinformation["numprocessingcurrent"]."\">\n";
    foreach(array_keys($protocol2chimera) as $protocol){
	//echo "PROT ".var_dump($protocol)."<BR>\n";
	echo "<option value=".$protocol.">". $protocol2chimera[$protocol][0] ."</option>\n";
    }
    echo "</select>\n";

    echo "<BR><BR><b>Merging/trimming</b>:<br/>\n";

    echo "<BR><img src=\"images/diagramOverlapseq.gif\" alt=\"merge diagram\"  height=\"500\" width=\"600\">\n";
    echo "<BR>Merge paired reads if<BR>\n";
    echo "<input type=\"radio\" name=\"mergeoverlap".$runinformation["numprocessingcurrent"]."\" value=\"False\" checked>after adapter trimming, if they overlap completely<br>\n";
    echo "<input type=\"radio\" name=\"mergeoverlap".$runinformation["numprocessingcurrent"]."\" value=\"True\" >The above plus if they show partial overlap (recommended for ancient DNA)<br>\n";


    /* missing key */
    /* echo "<BR><BR>If a key is used, do you allow 1 mismatch in the key ?<BR>\n"; */
    /* echo "<input type=\"radio\" name=\"allowMissing\" value=\"True\" checked>Yes allow a mismatch  <input type=\"text\" name=\"ctrlindex\" value=\"$ctrlindex\" size=\"7\"><BR><BR>\n"; */
    /* echo "<input type=\"radio\" name=\"allowMissing\" value=\"False\"><br>\n"; */
    //    }

    echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";

    echo "</form>\n";
    //echo "#";   
    
    
}

function displayStep4() {
    global $protocol2chimera;

    ///////////////////////////////////
    //BEGIN checking step 3 variables//
    ///////////////////////////////////
    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));


    //    var_dump($_POST);
    //protocol
    $printkey=False;

    //    for ($procnu = 1; $procnu <= $runinformation["numprocessing"]; $procnu++) {
	
    if($_POST["protocol".$runinformation["numprocessingcurrent"] ]      == "IllMultiplexCR1"){
	$printkey=True;
    }elseif($_POST["protocol".$runinformation["numprocessingcurrent"]]  == "IllMultiplexCR2"){
	$printkey=True;
    }

	//}

    if($printkey){
	echo "<h3>Step 4: Key and Quality flagging  (run id: ".$runinformation["runid"].")</h3><br>\n";   
    }else{
	echo "<h3>Step 4: Quality flagging (run id: ".$runinformation["runid"].")</h3><br>\n";   
    }


    /* echo     $runinformation["numprocessing"]; */
    /* exit; */


    //    for ($procnu = 1; $procnu <= $runinformation["numprocessing"]; $procnu++) {
	///echo "<BR><b>Processing #".$procnu.":</b><BR><BR>";
    //$runinformation["numprocessingcurrent"]
    $runinformation["adapter1".$runinformation["numprocessingcurrent"]]     = $protocol2chimera[ $_POST["protocol".$runinformation["numprocessingcurrent"]] ][1];
    $runinformation["adapter2".$runinformation["numprocessingcurrent"]]     = $protocol2chimera[ $_POST["protocol".$runinformation["numprocessingcurrent"]] ][2];
    $runinformation["chimeras".$runinformation["numprocessingcurrent"]]     = $protocol2chimera[ $_POST["protocol".$runinformation["numprocessingcurrent"]] ][3];
    $runinformation["protocol".$runinformation["numprocessingcurrent"]]     = $_POST["protocol".$runinformation["numprocessingcurrent"]];
    $runinformation["mergeoverlap".$runinformation["numprocessingcurrent"]] = ($_POST["mergeoverlap".$runinformation["numprocessingcurrent"]]=="True");
    /* } */

    
    //////////////////////////////////
    //END checking step 3 variables //
    //////////////////////////////////
    

    echo "This section defines quality filters that can be applied on resulting sequences. Please note that the sequences will <b>not</b> be removed, they will simply be marked as QC fail in the BAM file (more information about BAM flags: http://samtools.sourceforge.net/samtools.shtml#5)<BR> We define two types of filtering  procedures:<BR>
<UL>
<LI> Standard filtering based on expected number of mismatches using quality scores 
<LI> Additional filtering based on the complexity of the sequence per se
</UL>
If you plan to genotype, the first one might improve calls for low-coverage data at the cost of a lesser amount of sequences<BR> If you plan to do <i>de novo</i> genome assembly, the second might help by removing low complexity sequences<BR>";
	
    echo "This step flags reads with an unusually high number of expected mismatches as failing the QC controls. (more verbose description <a href=\"#expmism\">below</a>)<BR><BR>";
    echo "<form action=\"runprocess.php\" method=\"post\">\n";       
    echo "<input type=\"hidden\" name=\"step\" value=\"5\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";

    //for ($procnu = 1; $procnu <= $runinformation["numprocessing"]; $procnu++) {
    echo "<BR><b>Processing #".$runinformation["numprocessingcurrent"].":</b><BR><BR>";


    $keypossible=False;
    $key="";
    if($runinformation["protocol".$runinformation["numprocessingcurrent"]] == "IllMultiplexCR1"){
	$key="ACTC";
	$keypossible=True;
    }elseif($runinformation["protocol".$runinformation["numprocessingcurrent"]] == "IllMultiplexCR2"){
	$key="GTCT";
	$keypossible=True;
    }



    if($keypossible){

	echo "The library protocol allows for reading a key in the beginning of a read, specify the key:<BR>\n";
	if($runinformation["cyclesread2"] == 0 ){ //single-end
	    echo "<input type=\"text\" name=\"key1".$runinformation["numprocessingcurrent"]."\" value=\"".$key."\" size=\"5\"> <BR>\n";
	    echo "<input type=\"hidden\" name=\"key2".$runinformation["numprocessingcurrent"]."\" value=\"\" />\n";
	}else{
	    echo "<input type=\"text\" name=\"key1".$runinformation["numprocessingcurrent"]."\" value=\"\"         size=\"5\"> <BR>\n";
	    echo "<input type=\"text\" name=\"key2".$runinformation["numprocessingcurrent"]."\" value=\"".$key."\" size=\"5\"> <BR>\n";
	}
	echo "<BR><BR>";
    }else{
	echo "<form action=\"runprocess.php\" method=\"post\">\n";
	echo "<input type=\"hidden\" name=\"key1".$runinformation["numprocessingcurrent"]."\" value=\"\" />\n";
	echo "<input type=\"hidden\" name=\"key2".$runinformation["numprocessingcurrent"]."\" value=\"\" />\n";
    }





    echo "Basic filtering:<BR>\n";



    echo "<input type=\"radio\" name=\"filterseqexp".$runinformation["numprocessingcurrent"]."\" value=\"False\" checked>Do not flag reads<BR>\n";
    echo "<input type=\"radio\" name=\"filterseqexp".$runinformation["numprocessingcurrent"]."\" value=\"True\">Flag reads with a high number of expected mismatches<br>\n";


    echo "Threshold for normalized number of expected mismatches  : <input type=\"text\" name=\"seqNormExpcutoff".$runinformation["numprocessingcurrent"]."\" value=\"0.01\" size=\"5\"><BR><BR>\n";
    echo "Additional flagging:<BR>\n";
    echo "<input type=\"radio\" name=\"addfilters".$runinformation["numprocessingcurrent"]."\" value=\"False\" checked>Do not use additional flagging<br>\n";
    echo "<input type=\"radio\" name=\"addfilters".$runinformation["numprocessingcurrent"]."\" value=\"entropy\">Apply sequence entropy [0.0-2.0]  flag at:  <input type=\"text\" name=\"entropycutoff".$runinformation["numprocessingcurrent"]."\" value=\"0.85\" size=\"5\"> <BR>\n";
    echo "<input type=\"radio\" name=\"addfilters".$runinformation["numprocessingcurrent"]."\" value=\"frequency\">Apply base frequency [0.0-1.0] flag at: <input type=\"text\" name=\"frequencycutoff".$runinformation["numprocessingcurrent"]."\" value=\"0.1\" size=\"5\"> <BR>\n";    
    echo "<BR>Do not fail reads due to poor sample assignment:\n";    
    echo "<input type=\"checkbox\" value=\"rgnofail".$runinformation["numprocessingcurrent"]."\" name=\"rgnofail".$runinformation["numprocessingcurrent"]."\"><br/>\n";

    //}
    echo "<br><input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n<BR><BR><hline>";
    echo "<a name=\"expmism\">";

    echo "<b>Explanation of the expected mismatches filter:</b><BR>
Every base that is sequenced has a quality score. This quality score indicates the probability that the sequenced base is different from the one on the flowcell. For example, if the quality score is 0.001, the probability of error is 1/1000 and the probability of correctness is 999/1000. For this given base, the average expectation of mismatches is also 0.001. When summing up the average expected mismatches over the sequence, we can compute the expected number of mismatches for this given sequence. If you want to filter sequences that have an expected number of mismatches greater than 1 mismatch over 100 bases, use the first filter and use the 0.01 threshold.</a>\n";

    echo "</form>\n";

}


function displayStep5() {
     global $illuminawritedir;

						     
    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));
    /* echo  var_dump($runinformation); */
    /* exit; */
    ///////////////////////////////////
    //BEGIN checking step 4 variables//
    ///////////////////////////////////
    //    for ($procnu = 1; $procnu <= $runinformation["numprocessing"]; $procnu++) {

    if((int)$_POST["seqNormExpcutoff".$runinformation["numprocessingcurrent"]]<0 || (int)$_POST["seqNormExpcutoff".$runinformation["numprocessingcurrent"]]>1){
	echo "Sequence likelihood cannot be less than 0 or more than 1";
	exit;
    }

    if($_POST["addfilters".$runinformation["numprocessingcurrent"]] == "entropy"){
	if((int)$_POST["entropycutoff"]<0 || (int)$_POST["entropycutoff"]>2){
	    echo "Sequence entropy cutoff cannot be less than 0 or more than 2";
	    exit;
	}
    }

    if($_POST["addfilters".$runinformation["numprocessingcurrent"]] == "frequency"){
	if((int)$_POST["frequencycutoff"]<0 || (int)$_POST["frequencycutoff"]>1){
	    echo "Sequence base frequency cutoff cannot be less than 0 or more than 1";
	    exit;
	}
    }

    $runinformation["rgnofail".$runinformation["numprocessingcurrent"]]           = isset($_POST["rgnofail".$runinformation["numprocessingcurrent"]]);
    $runinformation["key1".$runinformation["numprocessingcurrent"]]               = $_POST["key1".$runinformation["numprocessingcurrent"]];
    $runinformation["key2".$runinformation["numprocessingcurrent"]]               = $_POST["key2".$runinformation["numprocessingcurrent"]];

    $runinformation["filterseqexp".$runinformation["numprocessingcurrent"]]       = ($_POST["filterseqexp".$runinformation["numprocessingcurrent"]] == "True");
    $runinformation["seqNormExpcutoff".$runinformation["numprocessingcurrent"]]   = $_POST["seqNormExpcutoff".$runinformation["numprocessingcurrent"]];
    $runinformation["filterentropy".$runinformation["numprocessingcurrent"]]      = ($_POST["addfilters".$runinformation["numprocessingcurrent"]] == "entropy");
    $runinformation["filterfrequency".$runinformation["numprocessingcurrent"]]    = ($_POST["addfilters".$runinformation["numprocessingcurrent"]] == "frequency");		    
    $runinformation["entropycutoff".$runinformation["numprocessingcurrent"]]      = $_POST["entropycutoff".$runinformation["numprocessingcurrent"]];
    $runinformation["frequencycutoff".$runinformation["numprocessingcurrent"]]    = $_POST["frequencycutoff".$runinformation["numprocessingcurrent"]];

    if( !$runinformation["filterseqexp".$runinformation["numprocessingcurrent"]] &&
	($runinformation["filterentropy".$runinformation["numprocessingcurrent"]] || $runinformation["filterfrequency".$runinformation["numprocessingcurrent"]]) ){
	echo "Cannot filter on additional without using basic filtering";
	exit;
    }
    //    }

    /* var_dump($runinformation); */
    /* exit; */
    //////////////////////////////////
    //END checking step 4 variables //
    //////////////////////////////////

    echo "<h3>Step 5: Read group assignment (run id: ".$runinformation["runid"].") </h3><br>\n";   

    //if($runinformation["numprocessingcurrent"]==1){
    //$runinformation["foundIndices"]="False";
    $indicesAlreadyThere="";
    $indicesAlreadyThereAgree=True;
    $indicesFileAlreadyThere=False;
   
    foreach($runinformation["lanes"] as $lanetouse){	    
	if(!file_exists($illuminawritedir."/".$runinformation["runid"]."/build/lane".$lanetouse."/indices.raw.txt")){
	    /* echo "The file ".$illuminawritedir."/".$runinformation["runid"]."/build/lane".$lanetouse."/indices.txt does not exist<BR>The first processing for this lane should have created this file"; */
	    /* exit; */
	}else{
	    //$runinformation["foundIndices"]="True";
	    $indicesFileAlreadyThere=True;
	    $stringToPrint="";
	    $fh = fopen($illuminawritedir."/".$runinformation["runid"]."/build/lane".$lanetouse."/indices.raw.txt", 'r') or die("can't open ".$illuminawritedir."/".$runinformation["runid"]."/build/lane".$lanetouse."/indices.raw.txt");
	    while( !feof($fh) ){
		$stringToPrint.=fgets($fh);
	    }
	    fclose($fh);

	    if(strlen($indicesAlreadyThere)==0){
		$indicesAlreadyThere = $stringToPrint;
	    }else{
		if( $indicesAlreadyThere != $stringToPrint){
		    $indicesAlreadyThereAgree=False;
		}
	    }
	}

	break;
    }

    if($indicesAlreadyThereAgree == False){
	echo "<b>Warning</b>: There are already index files for the lanes you selected and they seem to not be equal<BR>";
	
    }


    echo "<form action=\"runprocess.php\" method=\"post\">\n";

    
    //    echo "<input type=\"radio\" name=\"spikedin\" value=\"True\" checked>Spiked-in controls using P7 index:   <input type=\"text\" name=\"ctrlindex\" value=\"$ctrlindex\" size=\"7\"><BR>\n";
    //echo "<input type=\"radio\" name=\"spikedin\" value=\"False\">Dedicated lane (specify which below)\n";
    /* foreach($indexSchemes as $indexscheme){ */
    /* 	echo "<input type=\"".$indexscheme."\" name=\"".$indexscheme."\" value=\"True\" checked><BR>\n"; */
    /* } */

    echo "<input type=\"hidden\" name=\"step\" value=\"6\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";

    /* reactivate if  */
    /* Maximum number of mismatches for lookup among the indices:   */
    /* 					    <select name="mmrgassign"> */
    /* 					    <option value="2">2</option> */
    /* 					    <option value="1">1</option> */
    /* 					    <option value="0">0</option> */
    /* 					    </select> */   
    /* 					    <BR> */

    /* if($runinformation["numprocessing"] > 1){ */
    /* 	echo "<BR>For runs with multiple processings: please make sure that you have entered all the possible read groups on the flowcell<BR><BR>"; */
    /* } */

    if($indicesFileAlreadyThere ){
	echo "<b>Warning</b>: There are already index files, they were copy-pasted in the text field below<BR>please review them, any modifications you make will overwrite the current content<BR>";       
    }
    ?>
    Put your indices here: (see format below)<br />
please enter the indices for everyone on the lane<br />
	<B>NOTE:</B> For TruSeq, enter a t before the number (e.g. : t4)<br />

	<br />
	<br />

	<textarea rows="20" cols="100" name="indextext" wrap="physical" placeholder="Paste your indices here, if a particular RG has only an index for the first adapter and is mixed with a multiplexed paired-end run, use is4 as the second index"><?php
					  if($indicesAlreadyThereAgree ){
					      echo $indicesAlreadyThere;
					  }
					  ?></textarea>

<br>
        <label <?php
				echo  "for=pool".$runinformation["numprocessingcurrent"];
	    ?>
    >Pool ID (optional):</label>

       <input type="text" size="12" maxlength="200" <?php
	    echo "name=pool".$runinformation["numprocessingcurrent"]
	    ?>>
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

	
    /* }else{ */
    /* 	/\* displayStep6(); *\/ */
	

    /* } */



}



function displayStep6() {
    //var_dump($_POST);
    global $ctrlindex;
    global $ctrlindex2;
    global $emailAddrToSend;
    global $p7Indices;
    global $p5Indices;
    global $couchdb;
    global $illuminawritedir;

    // echo "test6";
    

    ///////////////////////////////////
    //BEGIN checking step 5 variables//
    ///////////////////////////////////
    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));
    $stringToPrint="";

    /* if($runinformation["numprocessingcurrent"]!=1){ */
    /* 	//build/lane1/indices.txt */
    /* 	foreach($runinformation["lanes"] as $lanetouse){	     */
    /* 	    if(!file_exists($illuminawritedir."/".$runinformation["runid"]."/build/lane".$lanetouse."/indices.txt")){ */
    /* 		echo "The file ".$illuminawritedir."/".$runinformation["runid"]."/build/lane".$lanetouse."/indices.txt does not exist<BR>The first processing for this lane should have created this file"; */
    /* 		exit; */
    /* 	    } */

    /* 	    $fh = fopen($illuminawritedir."/".$runinformation["runid"]."/build/lane".$lanetouse."/indices.txt", 'r') or die("can't open ".$illuminawritedir."/".$runinformation["runid"]."/build/lane".$lanetouse."/indices.txt"); */
    /* 	    while( !feof($fh) ){ */
    /* 		$stringToPrint.=fgets($fh); */
    /* 	    } */
    /* 	    fclose($fh); */


    /* 	    break; */
    /* 	} */

    /* }else{ */
    
	$indextext      = $_POST["indextext"];
	$poolID         = $_POST["pool".$runinformation["numprocessingcurrent"]];

	//    $runinformation["mmrgassign"] = (int)$_POST["mmrgassign"];

	if($runinformation["cyclesindx2"] != 0 ){
	    $stringToPrint.="#Index1\tIndex2\tName\n";
	}else{
	    $stringToPrint.="#Index1\tName\n";
	}

	$indexOfLines=0;
	$foundControl=0;

	if(!$indextext ){
	    echo "ERROR: please enter your indices in the text field or enter an empty header if you do not want demultiplexing";
	    exit;	 
	}

    
	//CHECKING EACH LINE
	foreach(explode("\n",$indextext) as $line){
	
	    $line=trim($line);
	    //skip empty lines
	    if(strlen($line) == 0){
		continue 1;
	    }

	    //$arrayfield=explode("\t",$line);
	    $arrayfield=preg_split('/\s+/', $line);



	    if($runinformation["cyclesindx2"] == 0 ){
		if(count($arrayfield) != 2){
		    echo "ERROR: For single index (you entered a value of 0 for the second index at step 1), lines must have 2 fields check line \"".$line."\"";
		    exit;
		}
	    }else{
		if(count($arrayfield) != 3){
		    echo "ERROR: For double index (you entered a positive value for the second index at step 1), lines must have 3 fields check line \"".$line."\"";
		    exit;
		}
	    }


	    if($indexOfLines == 0){
		//CHECKING HEADER
		if($arrayfield[0] != "#index"){
		    echo "ERROR: the first field of the header must be #index (case sensitive)";
		    exit;
		}

		if($arrayfield[1] != "p7"){
		    echo "ERROR: the second field of the header must be p7 (case sensitive)";
		    exit;
		}

		//if($indextype  == "double"){
		if($runinformation["cyclesindx2"] != 0 ){

		    if($arrayfield[2] != "p5"){
			echo "ERROR: the first field of the header must be p5 (case sensitive)";
			exit;
		    }
		}

	    }else{
		//CHECKING REMAINING FIELDS

		if(strstr($arrayfield[0],"\"") ||
		   strstr($arrayfield[0],"'")  || 
		   strstr($arrayfield[0],"\\")  || 
		   strstr($arrayfield[0],"/")   ){	       
		    echo "ERROR: The first field cannot have quotes or (back)slashes characters for line ".$line;
		    exit;
		}

		if($runinformation["cyclesindx2"] == 0 ){
		    
		
		    if(!array_key_exists($arrayfield[1],$p7Indices)){
			echo "ERROR: index for p7 ".(string)$arrayfield[1]." is not within the expected range ".$line . "";
			exit;	    
		    }
		    $stringToPrint.=$p7Indices[ (string)$arrayfield[1] ] ."\t".$arrayfield[0]."\n";

		}else{



		    if(!array_key_exists((string)$arrayfield[1],$p7Indices)){
			echo "ERROR: index for p7 ".(string)$arrayfield[1]." is not within the expected range ".$line . "";
			var_dump($p7Indices);
			exit;	    
		    }

		    if(!array_key_exists((string)$arrayfield[2],$p5Indices)){
			echo "ERROR: index for p5 ".(string)$arrayfield[2]." is not within the expected range ".$line . "";
			exit;	    
		    }



		    $stringToPrint.=$p7Indices[ (string)$arrayfield[1] ] ."\t".$p5Indices[ (string)$arrayfield[2] ]."\t".$arrayfield[0]."\n";
																     
		}

		if($arrayfield[0] == "control"){
		    $foundControl=1;
		}
	    }
	    $indexOfLines++;
	} //for each field explode

	if($foundControl == 0){
	    //if($indextype       == "single"){
	    
	    if($runinformation["cyclesindx2"] == 0 ){ //single
		/* if($runinformation["cyclesindx1"] == 8 ){//pad with A ? */
		/*     $stringToPrint.=$ctrlindex."A"."\tcontrol\n"; */
		/* }else{ */
		if($runinformation["spikedinmult"] ){
		    $stringToPrint.= 
			( str_repeat("A", $runinformation["cyclesindx1"]) )."\tcontrolA\n".
			( str_repeat("C", $runinformation["cyclesindx1"]) )."\tcontrolC\n".
			( str_repeat("G", $runinformation["cyclesindx1"]) )."\tcontrolG\n".
			( str_repeat("T", $runinformation["cyclesindx1"]) )."\tcontrolT\n";

		}else{
		    $stringToPrint.=$ctrlindex."\tcontrol\n";
		}
		/* } */
	    }else{
		//first index
		/* if($runinformation["cyclesindx1"] == 8 ){//pad with A ? */
		if($runinformation["spikedinmult"] ){
		    $stringToPrint.= 
			( str_repeat("A", $runinformation["cyclesindx1"]) )."\t".( str_repeat("A", $runinformation["cyclesindx2"]) )."\tcontrolA\n".
			( str_repeat("C", $runinformation["cyclesindx1"]) )."\t".( str_repeat("C", $runinformation["cyclesindx2"]) )."\tcontrolC\n".
			( str_repeat("G", $runinformation["cyclesindx1"]) )."\t".( str_repeat("G", $runinformation["cyclesindx2"]) )."\tcontrolG\n".
			( str_repeat("T", $runinformation["cyclesindx1"]) )."\t".( str_repeat("T", $runinformation["cyclesindx2"]) )."\tcontrolT\n";
		}else{
		    $stringToPrint.=$ctrlindex."\t".$ctrlindex2."\tcontrol\n";
		}
	    }
	}

	//}
    //////////////////////////////////
    //END checking step 5 variables //
    //////////////////////////////////
     /* $runinformation["indicesseq"]=$stringToPrint; */
     /* $runinformation["indicesraw"]= $_POST["indextext"]; */

    //var_dump($runinformation);
    echo "<h3>Step 6: Verify indices (run id: ".$runinformation["runid"].")</h3><br>\n";   
    echo "<form action=\"runprocess.php\" method=\"post\">\n";

    echo "<input type=\"hidden\" name=\"step\" value=\"7\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    /* if($runinformation["numprocessingcurrent"]==1){     */
    echo "<input type=\"hidden\" name=\"testindexorig\"   value=\"$indextext\" />\n";
    echo "<input type=\"hidden\" name=\"poolid\"   value=\"$poolID\" />\n";

    /* }else{ */
    /* 	echo "A previous set of indices has been found, please verify that they are correct. If not, please email ".$emailAddrToSend."<BR><BR>\n"; */
    /* } */
    echo "The following indices will be used:<BR>\n";
    echo "<textarea readonly rows=\"20\" cols=\"100\" name=\"textindex\" wrap=\"physical\">".$stringToPrint."</textarea><br />\n";
    echo "<input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";
    //echo "<button onclick=\"history.go(-1);\">cancel</button>\n";
    

    echo "</form>\n";

}


function displayStep7() {
    ////////////////////////////////////
    //BEGIN checking step 6 variables //
    ////////////////////////////////////

    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));


    $runinformation["indicesseq"] = $_POST["textindex"];
    //if($runinformation["numprocessingcurrent"]==1){    
    $runinformation["indicesraw"] = $_POST["testindexorig"];
    $runinformation["poolid"]     = $_POST["poolid"];
   
	/* }else{ */
    /* 	$runinformation["indicesraw"]= ""; */
    /* } */
    /* echo "pool id :"; */
    /* var_dump($runinformation["poolid"]); */


    //////////////////////////////////
    //END checking step 6 variables //
    //////////////////////////////////
    echo "<h3>Step 7: Mapping (run id: ".$runinformation["runid"].")</h3><br>\n"; 

    echo "<form action=\"runprocess.php\" method=\"post\">\n";

    echo "Notes:<BR>\n";
    echo "<b>hg19_evan</b>: Is the 1000 genomes hg19 plus some decoy sequences (phix,herpes and unmapped potentially human sequences) .<BR><BR>\n";
    echo "<b>Please note</b>: If you have previous data aligned under a different version of a genome and pretend on merging the bam files for say genotyping, we strongly recommend to consistently use the same version of the genome.<BR><BR>\n";
    echo "Selecting the parameters for BWA:<BR><BR>\n";
    echo "Normally, the following are used:<BR> for modern DNA: \"-n 0.04 -o 1\"<BR>for ancient DNA: \"-n 0.01 -o 2 -l 16500\"<BR><BR>\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"8\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";

    //    for ($procnu = 1; $procnu <= $runinformation["numprocessing"]; $procnu++) {

    echo "<BR><b>Processing #".$runinformation["numprocessingcurrent"].":</b><BR><BR>";

    //echo "Mapping using BWA<BR>\n";
    echo "Mapping using BWA:  <input type=\"checkbox\" value=\"True\" name=\"usebwa".$runinformation["numprocessingcurrent"]."\" checked ><br/>\n";
    $arrayofGenomes=getGenomes();
    natcasesort($arrayofGenomes);
    echo "Selecting the target genome:<BR><BR>\n";
    echo "<select name=\"genomebwa".$runinformation["numprocessingcurrent"]."\" size=\"1\">\n";
    foreach($arrayofGenomes as $agenome){
	if($agenome == "hg19_evan"){
	    echo "<option value=\"".$agenome."\" selected=\"selected\">".$agenome."</option>\n";
	}else{
	    echo "<option value=\"".$agenome."\" >".$agenome."</option>\n";
	}
    }
    echo "</select><BR><BR>\n";

    
    echo "<select name=\"parambwa".$runinformation["numprocessingcurrent"]."\" size=\"1\">\n";
    echo "<option value=\"default\" >Default parameters (modern DNA)</option>\n";
    echo "<option value=\"ancient\" >Ancient parameters (ancient DNA)</option>\n";
    echo "</select><BR><BR>\n";
    //    }
    echo "<input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Next &gt;\" />\n";

    echo "</form>\n";
}


function displayStep8() {
    global $illuminawritedir;
    ////////////////////////////////////
    //BEGIN checking step 7 variables //
    ////////////////////////////////////




    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));
    //    for ($procnu = 1; $procnu <= $runinformation["numprocessing"]; $procnu++) {
    $runinformation["usebwa".$runinformation["numprocessingcurrent"]]     = isset($_POST["usebwa".$runinformation["numprocessingcurrent"]]);
    $runinformation["genomebwa".$runinformation["numprocessingcurrent"]]  = $_POST["genomebwa".$runinformation["numprocessingcurrent"]];
    $runinformation["parambwa".$runinformation["numprocessingcurrent"]]   = $_POST["parambwa".$runinformation["numprocessingcurrent"]];
    //}

    //////////////////////////////////
    //END checking step 7 variables //
    //////////////////////////////////

    echo "<h3>Step 8: Summary</h3>";
    /* echo "".var_dump($runinformation); */
    /* exit; */
    $indicesrawtext=$runinformation["indicesraw"];
    $indicesseqtext=$runinformation["indicesseq"];



    //raw
    $indicesrawjson=array();
    $firstline=True;
    foreach(explode("\n",$runinformation["indicesraw"]) as $line){

	/* if($firstline){ */
	/*     $firstline=False; */
	/*     continue 1; */
	/* } */

	$line=trim($line);
	if(strlen($line) == 0){
	    continue 1;
	}

	$temparray = preg_split('/\s+/', $line); //explode("\t",$line);

	if(      count($temparray) == 3){//double
	    /* $tempstr = "{name: ".$temparray[0]." , "." p7: "  .$temparray[1]." , "." p5: "  .$temparray[2]." } "; */
	    $tempstr = array("name"  => $temparray[0],
			     "p7"    => $temparray[1],
			     "p5"    => $temparray[2]);
	    array_push($indicesrawjson, $tempstr);
	    //echo $tempstr;
	}elseif( count($temparray) == 2){//single
	    /* $tempstr ="{name: ".$temparray[0]." , "." p7: "  .$temparray[1]." } "; */
	    $tempstr = array("name"  => $temparray[0],
			     "p7"    => $temparray[1]);
	    array_push($indicesrawjson, $tempstr);
	    //echo $tempstr;
	}else{
	    echo "ERROR: wrong number of fields (".count($temparray).") in line \"".$line."\" in raw indices ".$runinformation["indicesraw"];
	    exit;
	}
    }
    /* $runinformation["indicesraw"] = " [ ".implode(",",$indicesrawjson)." ] "; */
    $runinformation["indicesraw"] = $indicesrawjson;

    //seq
    $indicesseqjson=array();
    $firstline=True;
    foreach(explode("\n",$runinformation["indicesseq"]) as $line){
	/* if($firstline){ */
	/*     $firstline=False; */
	/*     continue 1; */
	/* } */

	$line=trim($line);
	if(strlen($line) == 0){
	    continue 1;
	}

	$temparray = preg_split('/\s+/', $line); //explode("\t",$line);

	if(      count($temparray) == 3){//double
	    //$tempstr = "{name: ".$temparray[0]." , "." p7: "  .$temparray[1]." , "." p5: "  .$temparray[2]." } ";
	    $tempstr = array("name" => $temparray[2],
			     "p7"    => $temparray[0],
			     "p5"    => $temparray[1]);
	    array_push($indicesseqjson, $tempstr);
	    //echo $tempstr;
	}elseif( count($temparray) == 2){//single
	    //$tempstr ="{name: ".$temparray[0]." , "." p7: "  .$temparray[1]." } ";
	    $tempstr = array("name" => $temparray[1],
			     "p7"    => $temparray[0]);
	    array_push($indicesseqjson, $tempstr);
	    //echo $tempstr;
	}else{
	    echo "ERROR: wrong number of fields (".count($temparray).") in line \"".$line."\" in raw indices ".$runinformation["indicesseq"];
	    exit;
	}
    }
    $runinformation["indicesseq"] = $indicesseqjson;
    //echo $runinformation["indicesseq"] ;


    //echo "Please review the following information prior to pressing submit:<BR>\n";
    echo "<form action=\"runprocess.php\" method=\"post\">\n";
    echo "<input type=\"hidden\" name=\"step\" value=\"9\" />\n";
    echo "<input type=\"hidden\" name=\"runinformation\" value=\"".htmlspecialchars(serialize($runinformation))."\" />\n";
    //    var_dump($runinformation);
    echo "<BR>Please review the following information prior to submitting (\"submit\" buttom at the bottom of the page):<BR>\n";

    $htmltable="";

    $htmltable.="<table   border=0>\n";
    $htmltable.="<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>General:      </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Run ID      :</TD><TD> ".$runinformation["runid"]."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Experiment name:</TD><TD> ".$runinformation["expname"]."</TD></TR>\n";

    $htmltable.="<TR><TD nowrap>Sequencer      :</TD><TD> ".$runinformation["sequencer"]."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Your email      :</TD><TD> ".$runinformation["email"]."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Cycle for read1 :</TD><TD> ".$runinformation["cyclesread1"]."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Cycle for read2 :</TD><TD> ".$runinformation["cyclesread2"]."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Cycle for index 1 :</TD><TD> ".$runinformation["cyclesindx1"]."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Cycle for index 2 :</TD><TD> ".$runinformation["cyclesindx2"]."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Lanes to analyze:</TD><TD> ".implode(",",$runinformation["lanes"])."</TD></TR>\n";

    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Basecalling:      </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Basecalling using Bustard      :</TD><TD> ".($runinformation["bustard"]?"yes":"no")."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Basecalling using freeIbis     :</TD><TD> ".($runinformation["freeibis"]?"yes":"no")."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    /* $htmltable.="<TR><TD nowrap>Number of independent processing     :</TD><TD> ".($runinformation["numprocessing"])."</TD></TR>\n"; */


    //    for ($procnu = 1; $procnu <= $runinformation["numprocessing"]; $procnu++) {
    //echo "<BR><b>Processing #".$procnu.":</b><BR><BR>";
    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";

    $htmltable.="<TR><TD nowrap><B>Processing #  ".$runinformation["numprocessingcurrent"]." </B> :      </TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Merging/trimming:      </TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Merge partially overlapping sequencing     :</TD><TD> ".($runinformation["mergeoverlap".$runinformation["numprocessingcurrent"]]?"yes":"no")."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Adapter 1     :</TD><TD> ".($runinformation["adapter1".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Adapter 2     :</TD><TD> ".($runinformation["adapter2".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Potential chimeras     :</TD><TD> ".($runinformation["chimeras".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Protocol       :</TD><TD> ".($runinformation["protocol".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Key read#1     :</TD><TD> ".(strlen($runinformation["key1".$runinformation["numprocessingcurrent"]])==0?"none":$runinformation["key1".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Key read#2     :</TD><TD> ".(strlen($runinformation["key2".$runinformation["numprocessingcurrent"]])==0?"none":$runinformation["key2".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    //$htmltable.="<TR><TD nowrap>Maximum number of mismatches for lookup in RG :</TD><TD> ".($runinformation["mmrgassign"])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Quality filtering:      </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Flag sequences with high exp. of mismatch  :</TD><TD> ".(($runinformation["filterseqexp".$runinformation["numprocessingcurrent"]]=="1")?"yes":"no")."</TD></TR>\n";
    if( ($runinformation["filterseqexp".$runinformation["numprocessingcurrent"]]=="1") ){
	$htmltable.="<TR><TD nowrap>Normalized expectency of mismatch cutoff    :</TD><TD> ".($runinformation["seqNormExpcutoff".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    }
    
    $htmltable.="<TR><TD nowrap>Flag sequences based on  entropy  :</TD><TD> ".(($runinformation["filterentropy".$runinformation["numprocessingcurrent"]])?"yes":"no")."</TD></TR>\n";
    if( ($runinformation["filterentropy".$runinformation["numprocessingcurrent"]]=="1") ){
	$htmltable.="<TR><TD nowrap>Entropy cutoff  :</TD><TD> ".($runinformation["entropycutoff".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    }
    $htmltable.="<TR><TD nowrap>Flag sequences using frequency  :</TD><TD> ".(($runinformation["filterfrequency".$runinformation["numprocessingcurrent"]])?"yes":"no")."</TD></TR>\n";

    if( ($runinformation["filterfrequency".$runinformation["numprocessingcurrent"]]=="1") ){
	$htmltable.="<TR><TD nowrap>Frequency cutoff  :</TD><TD> ".($runinformation["frequencycutoff".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    }

    $htmltable.="<TR><TD nowrap>Allow read groups to fail due to poor sample assignment  :</TD><TD> ".(($runinformation["rgnofail".$runinformation["numprocessingcurrent"]])?"no":"yes")."</TD></TR>\n";

    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>      </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>BWA mapping:     </TD><TD></TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Map using BWA   :</TD><TD> ".(($runinformation["usebwa".$runinformation["numprocessingcurrent"]])?"yes":"no")."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Genome to use   :</TD><TD> ".($runinformation["genomebwa".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>BWA parameters  :</TD><TD> ".($runinformation["parambwa".$runinformation["numprocessingcurrent"]])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap>Pool ID         :</TD><TD> ".( ($runinformation["poolid"]=="")?"none provided":$runinformation["poolid"])."</TD></TR>\n";
    $htmltable.="<TR><TD nowrap> </TD><TD></TD></TR>\n";
    // $htmltable.="<TR><TD nowrap>Indices:     </TD><TD></TD></TR>\n";
    //    }

    $htmltable.="</table>\n";


    /* if($runinformation["numprocessingcurrent"] == 1){ */
    $htmltable.="<BR>Indices to use  :<BR><PRE>\n".($indicesrawtext)."</PRE><BR>\n";
    /* } */
    $htmltable.="<BR>Indices raw     :<BR><PRE>\n".($indicesseqtext)."</PRE><BR>\n";

    echo $htmltable;
    //$runinformation["htmltable"] = $htmltable;
    echo "<input type=\"hidden\" name=\"htmltable\" value=\"".$htmltable."\" />\n";

    //echo json_encode($runinformation);

    /* $runid=$runinformation["runid"]; */
    /* $targetfile=$illuminawritedir."/".$runid."_".implode(",",$runinformation["lanes"]).".json"; */
    /* if(file_exists ( $targetfile )){ */
    /*    echo "<font color=red size=+2>Warning:</font> file:<BR>".$targetfile." <BR>already exists, pressing submit will overwrite the data <br><br>"; */
    /* } */
    echo "<input type=\"submit\" name=\"submitButton\" id=\"nextButton\" value=\"Submit\" />\n";

    echo "</form>\n";

}


function displayStep9() {
    global $illuminawritedir;
    global $emailAddrToSend;
    global $json2makePath;
    global $couchdb;


    $runinformation = unserialize(stripslashes(htmlspecialchars_decode($_POST["runinformation"])));

    /* echo var_dump($runinformation); */

    //build directory 
    $runid=$runinformation["runid"];

    if(!file_exists($illuminawritedir."/".$runid."/")){
	if(!mkdir($illuminawritedir."/".$runid."/")){
	    echo "ERROR: cannot create directory ".$illuminawritedir."/".$runid."\n";
	    exit;
	}
    }

    if(!file_exists($illuminawritedir."/".$runid."/build/")){
	if(!mkdir($illuminawritedir."/".$runid."/build/")){
	    echo "ERROR: cannot create directory ".$illuminawritedir."/".$runid."/build/"."\n";
	    exit;
	}
    }

    foreach($runinformation["lanes"] as $lanetouse){

	if(!file_exists($illuminawritedir."/".$runid."/build/lane".$lanetouse."/")){
	    if(!mkdir($illuminawritedir."/".$runid."/build/lane".$lanetouse."/")){
		echo "ERROR: cannot create directory ".$illuminawritedir."/".$runid."/build/lane".$lanetouse."/"."\n";
		exit;
	    }
	}
    
    }
    

    /* var_dump($runinformation); */
    $runid     = $runinformation["runid"];
    $emailuser = $runinformation["email"];
    $htmltable = $_POST["htmltable"];
    $htmltable = str_replace(array("<TR>","<TD>","</TD>","<TD nowrap>","<BR>","<PRE>","</PRE>","</table>","<table","border=0>") ,"" , $htmltable);
    $htmltable = str_replace(array("</TR>",) ,"" , $htmltable);
    /* echo $htmltable; */

    $targetfile=$illuminawritedir."/".$runid."/build/".$runid."_".implode(",",$runinformation["lanes"])."-".$runinformation["numprocessingcurrent"].".json";
    echo "printing to ".$targetfile."<BR>\n";
    
    //encode as json text fields

    $stringtoprint=json_encode($runinformation);
    
    

    
    //call json2make    
    //$json2makePath
    if(1){ //set to 0 to test
	$fh = fopen($targetfile, 'w') or die("can't open $targetfile");
	fwrite($fh, $stringtoprint);
	fclose($fh);

	foreach($runinformation["lanes"] as $lanetouse){	
	    $cmdToRun="python ".$json2makePath." -o ".$illuminawritedir."/".$runid."/build/ $targetfile";
	    //echo "<BR> Running command".$cmdToRun."<BR><BR>";
	    $outputStore="";
	    $returnCode=1;

	    exec($cmdToRun,$outputStore,$returnCode);
	    if($returnCode != 0){
		echo "Following command failed:<BR> ".$cmdToRun."<BR><BR>please contact directly ".$emailAddrToSend."<BR><BR>got the following output: ".var_dump($outputStore)." <br>\n";
		exit;
	    }else{
		echo "<BR>success!<BR>";
	    }

	}    
    }    
    //var_dump($couchdb);
    if(isset($couchdb)){
	
    	$url = $couchdb;

    	$file_headers = @get_headers($url);
    	if( strpos($file_headers[0], '404') ) {
    	    $exists = false;
    	}
    	else {
    	    $exists = true;
    	}
	
    	if(!$exists){
    	    echo "Unable to connect to couchDB:<BR> ".$url."<BR><BR>please contact directly ".$emailAddrToSend."<BR><BR>got the following output: ".var_dump($file_headers)." <br>\n";
    	}else{//
    	    echo "CouchDB detected, trying to post indices<br>";

    	    $url = $couchdb;
    	    if($runinformation["poolid"]!=""){//has pool ID
    		$url = $url."/".$runinformation["poolid"];
    	    }
	    
    	    $indicestopost="Library\tP7_index";
    	    if($runinformation["cyclesindx2"] != 0 ){
    		$indicestopost.="\tP5_index";
    	    }
    	    //$indicestopost.="";

	    foreach(array_slice($runinformation["indicesraw"],1) as $recordIndex){
		//var_dump($recordIndex);
		if(      count($recordIndex) == 3){//double
		    $indicestopost.="\n".$recordIndex["name"]."\t".$recordIndex["p7"]."\t".$recordIndex["p5"];
		}elseif( count($recordIndex) == 2){//single
		    $indicestopost.="\n".$recordIndex["name"]."\t".$recordIndex["p7"];
		}else{
		    echo "ERROR: wrong number of fields (".count($recordIndex).") in line \"".$recordIndex."\"";
		    exit;
		}

	    }
	    //var_dump($indicestopost);
    	    $options = array(
    			     'http' => array(
    					     'header'  => "Content-type: text/csv\r\n"
    					     ,'method'  => 'POST'
    					     ,'content' => $indicestopost	)
    			     );
	    
    	    $context  = stream_context_create($options);
    	    $result   = file_get_contents($url, false, $context);

    	    $couchdbresponse = json_decode($result,true);
    	    if(isset($couchdbresponse["status"])){
		echo "Status : ".$couchdbresponse["status"]."<BR>";
	    }else{
		echo "Status : error while posting to couchDB<BR>";
	    }
	    
   



    	}

    }else{
	echo "CouchDB was not set";

    }



    //send user email
    $mailu = new EMail;
    $mailu->Username = 'sbsuser';
    $mailu->Password = 'sbs123';

      
    $mailu->SetFrom("sbsuser@eva.mpg.de","");  // Name is optional

    foreach(explode(",",$emailuser) as $emailAddrTo){ //for each email in the array above
	$mailu->AddTo($emailAddrTo,""); // Name is optional
    }

    $mailu->Subject = "Analysis submitted for run ".$runid;
    $mailu->Message = "This is an automated message so please do not reply.\n------------------\nThe following user(s): ".$emailuser."\n".
	"Submitted a request for analysis for run ".$runid." processing #".$runinformation["numprocessingcurrent"]."\nPlease keep this email for your own personal record\n\n\n------------------\nThe following parameters will be used: ".$targetfile."\n------------------\n\n".$htmltable."\n\n\n------------------\njson file\n------------------\n".json_encode($runinformation);
    $mailu->ConnectTimeout = 30;  // Socket connect timeout (sec)
    $mailu->ResponseTimeout = 8;  // CMD response timeout (sec)
    $success = $mailu->Send();
    $stringforscreen="";

    if($success != 1){
	$stringforscreen.="the user notification email WAS NOT sent successfully, contact directly ".$emailAddrToSend."<br>\n";
    }else{
	$stringforscreen.="the user notification email was sent successfully <br>\n";
    }


    //send admin email
    $mail = new EMail;
    $mail->Username = 'sbsuser';
    $mail->Password = 'sbs123';
    
    $mail->SetFrom("sbsuser@eva.mpg.de","");  // Name is optional

    foreach(explode(",",$emailAddrToSend) as $emailAddrTo){ //for each email in the array above
	$mail->AddTo($emailAddrTo,""); // Name is optional
    }

    $mail->Subject = "User submitted request for run ".$runid;
    $mail->Message = "User email ".$emailuser."\n".
	"Processing #".$runinformation["numprocessingcurrent"]."\n".
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
    

    
    

	?>
    <h1>Thank You</h1>
	<p>Thank you, your application has been received.</p>

	<?php
	}
?>
</body>
</html>
