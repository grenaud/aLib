<?php
ini_set('display_errors', 'On');
error_reporting(E_ALL);

if( ! file_exists(getcwd()."/config.xml")  ){
    echo "Configuration file not found ".getcwd()."/config.xml";
    exit(1);
 }

$xmlconf = simplexml_load_file( getcwd()."/config.xml" );

#$illuminareaddir."/".$runid."/RunInfo.xml") ){ 

//CONFIG DATA
$illuminareaddir  = $xmlconf->illuminareaddir;
$illuminawritedir = $xmlconf->illuminawritedir;
//$illuminajson     = $xmlconf->illuminajson;
$sequencers=array();
foreach($xmlconf->sequencers->sequencer as $seqelem){
    $sequencers[ (string)$seqelem["id"] ]= (string)$seqelem["type"];
}
$runstodisplay= $xmlconf->runstodisplay;

//#exit(1);
//global vars
//$analysisrequests = array();


// $illuminareaddir="/mnt/solexa/";
// $illuminawritedir="/mnt/ngs_data/";
// $illuminajson="/mnt/solexa/tmp/";
// $analysisrequests=array();

// $sequencers=array("M00518"     =>   "miseq",
// 		  "SN7001204"  =>   "hiseq");
// $runstodisplay=30;




function cmp($a, $b) {
    if ($a['mtime'] == $b['mtime']) {
        return 0;
    }
    return ($a['mtime'] > $b['mtime']) ? -1 : 1;
} 

function lastestFile($directoryToCheck){
    $myDirectory = opendir($directoryToCheck);
    $filelist=array();
    
    while($entryName = readdir($myDirectory)) {
	$item = $directoryToCheck."/".$entryName;

	if($entryName != "." and $entryName != ".."){
	    $filelist[] = array('name'  => $entryName, 
				'mtime' => filemtime($item)); 
	}
    }
    
    usort($filelist, "cmp"); 

    closedir($myDirectory);
    $stringFound=$filelist[0]['name'];
    $stringFound=substr($stringFound,1,strlen($stringFound)-3);
    return $stringFound;
}

#check if run finished
function checkrun($runid) {

    if( file_exists($runid."/Run.completed") ){
	return "finished";
    }
    
    if(file_exists($runid."/Data/")                     &&
       file_exists($runid."/Data/Intensities/")         &&
       file_exists($runid."/Data/Intensities/L001/")    &&
       file_exists($runid."/Data/Intensities/BaseCalls/L001/") ){
	return "I=".lastestFile($runid."/Data/Intensities/L001/")."/B=".lastestFile($runid."/Data/Intensities/BaseCalls/L001/");
    }

    return "unknown";   
} 


 function getNumberLanes($runid) {
     global $illuminareaddir;
     if( file_exists($illuminareaddir."/".$runid."/RunInfo.xml") ){ 
	 $xmlp = simplexml_load_file($illuminareaddir."/".$runid."/RunInfo.xml");
	 return (int)$xmlp->Run->FlowcellLayout["LaneCount"];
     }// else{
     // 	 echo "file not found ".$illuminareaddir."/".$runid."/RunInfo.xml";
     // 	 exit(1);
     // }

    return 1;   //something failed
} 


/* function checkAnalysisRequest($runid,$numberLanes) { */
/*     global $illuminajson; */
/*     //    global $analysisrequests; */

/*     #$arrayfile=glob($illuminajson."/".); */
/*     #return var_dump($arrayfile); */
/*     $lanesSent=array(); */
/*     foreach($analysisrequests as $jsonfile){ */
/* 	if( !strncmp($jsonfile,$runid,strlen($runid)) ){ */
/* 	    $lanesSent=array_merge($lanesSent,explode(",",substr($jsonfile,strlen($runid)+1,-5))); */
/* 	} */
/*     } */

/*     if(count($lanesSent) != 0){ */
/* 	return "Submitted: ".implode(",",$lanesSent); */
/*     }else{ */
/* 	return "none"; */
/*     } */
/* } */

function checkAnalysisStatus($runid,$numberLanes) {
    global $illuminawritedir;
    if( file_exists($illuminawritedir."/".$runid ) ){
	$stringtoReturn = "";
	$bustardFinished=array();
	$freeibiFinished=array();

	#detect if bustard of freeibis
	if( file_exists($illuminawritedir."/".$runid."/Bustard/QC/") ){

	    for($lane=1;$lane<=$numberLanes;$lane++){
		$stringtoReturn.=" ".$lane.":";
		if( file_exists($illuminawritedir."/".$runid."/Bustard/QC/clusterTally_".$lane.".OK" ) ){
		    #array_push($bustardFinished,$lane);
		    $stringtoReturn.="done";
		    continue;
		}		

		if( file_exists($illuminawritedir."/".$runid."/Bustard/Final_Sequences/s_".$lane."_sequence.bam" ) ){
		    $stringtoReturn.="running";
		    continue;
		}

		if( file_exists($illuminawritedir."/".$runid."/Bustard/Raw_Sequences/s_".$lane."_sequence.bam" ) ){
		    $stringtoReturn.="basecall";
		    continue;
		}		
	    }   

	}elseif( file_exists($illuminawritedir."/".$runid."/Ibis/QC/") ){

	    for($lane=1;$lane<=$numberLanes;$lane++){
		$stringtoReturn.=" ".$lane.":";
		if( file_exists($illuminawritedir."/".$runid."/Ibis/QC/clusterTally_".$lane.".OK" ) ){
		    #array_push($bustardFinished,$lane);
		    $stringtoReturn.="done";
		    continue;
		}		

		if( file_exists($illuminawritedir."/".$runid."/Ibis/Final_Sequences/s_".$lane."_sequence.bam" ) ){
		    $stringtoReturn.="running";
		    continue;
		}

		if( file_exists($illuminawritedir."/".$runid."/Ibis/Raw_Sequences/s_".$lane."_sequence.bam" ) ){
		    $stringtoReturn.="basecall";
		    continue;
		}		
	    }   
	}else{
	    return "folder created";
	}
		

	
	return $stringtoReturn;

    }else{
	return "not started";
    }
}


function checkAnalysisRequest($runid,$numberLanes) {
    global $illuminawritedir;
    if( file_exists($illuminawritedir."/".$runid ) ){
	
	if( file_exists($illuminawritedir."/".$runid."/build/" ) ){
	    $lanesSent=array();
	    /* print $illuminawritedir."/".$runid."/build/"; */
	    //exit;
	    $myDirectory = opendir($illuminawritedir."/".$runid."/build/");
	    
	    while($entryName = readdir($myDirectory)) {
		if($entryName != "." and $entryName != ".."){		    

		    if(substr($entryName,-5) == ".json"){ 
			/* print $entryName; */
			/* exit; */
			$lanesSent=array_merge($lanesSent,explode(",",substr($entryName,strlen($runid)+1,-5)));
		    }
		}
	    }

	    if(count($lanesSent) != 0){
		return "Submitted: ".implode(",",$lanesSent);
	    }else{
		return "none";
	    }


	}

	return "none";
    }
    return "none";
}

#function fecho($string) {
#    echo $string;
#    ob_flush();
#}



// $myDirectory = opendir("/mnt/");

// while($entryName = readdir($myDirectory)) {
//     echo $entryName;
//  }
// exit;

#ob_flush();

######################
#DETECT CURRENT RUNS #
######################
// read dir and save name,size,date in array
$runlist = array(); 
$myDirectory = opendir($illuminareaddir);

while($entryName = readdir($myDirectory)) {

    if($entryName != "." and $entryName != ".."){
	if(substr($entryName,0,2) > "10" and substr($entryName,0,2) < "19"){ #should be good until 2019
	    $arrayfield=explode("_",$entryName);
	    $seqtype="unknown";
	    if(isset($sequencers[ $arrayfield[1] ])) {
		$seqtype=$sequencers[ $arrayfield[1] ];
	    }
	    $item = $illuminareaddir."/".$entryName;
	    
	    $runlist[] = array('dir' => $item, 
			       'name' => $entryName, 
			       'technology' => $seqtype,
			       'mtime' => filemtime($item)); 
	}
    }
 }

closedir($myDirectory);

// echo "ok2";

/* $myDirectory = opendir($illuminajson); */

/* while($entryName = readdir($myDirectory)) { */

/*     if($entryName != "." and $entryName != ".."){ */
/* 	#echo "#".$entryName."#\n#".substr($entryName,-4)."#\n"; */
/* 	if(substr($entryName,-5) == ".json"){  */
/* 	    array_push($analysisrequests,$entryName); */
/* 	} */
/*     } */
/* } */

/* closedir($myDirectory); */
// echo "ok3";
#

// <html>
// <head>
// <title>Illumina runs</title>
// </head>
// <body>
?>

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN">
<html lang="en">
<head>
    <title>Illumina runs</title>
</head>
<body>
<!-- Progress bar holder -->
<div id="progress" style="width:500px;border:1px solid #ccc;"></div>
<!-- Progress information -->
<div id="information" style="width"></div>

<?

usort($runlist, "cmp"); 



$stringWithTable="";

$stringWithTable.="<form action=\"runprocess.php\" method=\"post\">\n";
$stringWithTable.="<TABLE border=1>\n";

$stringWithTable.="<TR>\n";
$stringWithTable.="<TD>Run name</TD>\n";  
$stringWithTable.="<TD>Sequencer</TD>\n";  
$stringWithTable.="<TD>Sequencing<BR>progress<BR>I=Intensities<BR>B=basecalls</TD>\n";  
$stringWithTable.="<TD>Analysis</TD>\n";  

$stringWithTable.="<TD>Request</TD>\n";  

$stringWithTable.="<TD>Launch</TD>\n";  

#echo "<TD>status</TD>\n";  

$stringWithTable.="</TR>\n";


$runlist =  array_slice($runlist,0,$runstodisplay);


$numberFound=1;

foreach($runlist as $one_file) {
    $percent = intval($numberFound/count($runlist) *100)."%";

    // echo "p ".$percent."<BR>";
    // echo $numberFound."<BR>";
    // echo sizeof($runlist)."<BR>";
    // echo $numberFound/sizeof($runlist)."<BR>";

    echo '<script language="javascript">
    document.getElementById("progress").innerHTML="<div style=\"width:'.$percent.';background-color:#ddd;\">&nbsp;</div>";
    document.getElementById("information").innerHTML="'.$numberFound.' runs(s) found.";
    </script>';
    echo str_repeat(' ',1024*64);
    flush();
    //sleep(1);

    // store the information.
    $stringWithTable.="<TR>\n";
    $stringWithTable.="<TD>".$one_file['name']."</TD>\n";  
    $stringWithTable.="<TD>".$one_file['technology']."</TD>\n";  
    $stringWithTable.="<TD>".checkrun($one_file['dir'])."</TD>\n";  
    #$stringWithTable.=$one_file['technology']."<BR>\n";  
    #if(0){

    $numberLanes = getNumberLanes($one_file['name']);

    $stringWithTable.="<TD>".checkAnalysisStatus( $one_file['name'],$numberLanes)."</TD>\n";  
    $stringWithTable.="<TD>".checkAnalysisRequest($one_file['name'],$numberLanes)."</TD>\n";  

    // if($one_file['technology'] == "hiseq"){
    // 	$stringWithTable.="<TD>".checkAnalysisStatus( $one_file['name'],8)."</TD>\n";  
    // 	$stringWithTable.="<TD>".checkAnalysisRequest($one_file['name'],8)."</TD>\n";  
    // }elseif($one_file['technology'] == "miseq"){
    // 	$stringWithTable.="<TD>".checkAnalysisStatus( $one_file['name'],1)."</TD>\n";  
    // 	$stringWithTable.="<TD>".checkAnalysisRequest($one_file['name'],1)."</TD>\n";  
    // }else{
    // 	$stringWithTable.="<TD>".checkAnalysisStatus( $one_file['name'],1)."</TD>\n";  
    // 	$stringWithTable.="<TD>".checkAnalysisRequest($one_file['name'],1)."</TD>\n";  
    // }
    #}

    $stringWithTable.="<TD><input type=\"submit\" value=\"launch\" name=\"".$one_file['name']."\"><TD>\n";
    $stringWithTable.="</TR>\n";
    $numberFound++;
} 

echo '<script language="javascript">
    document.getElementById("progress").innerHTML="<div style=\"display:none\" >&nbsp;</div>";
    document.getElementById("information").innerHTML="'.$numberFound.' runs(s) found.";
    </script>';
echo str_repeat(' ',1024*64);
flush();
$stringWithTable.="</form>\n";

$stringWithTable.="</TABLE>\n";
$stringWithTable.="</HTML>";

echo $stringWithTable;
##########################
#end DETECT CURRENT RUNS #
##########################





?> 