<?php
ini_set('display_errors', 'On');
error_reporting(E_ALL);

if( ! file_exists(getcwd()."/config.json")  ){
    echo "Configuration file not found ".getcwd()."/config.json";
    exit(1);
}

function beginsWith($bigStr, $substrToFind){
    return !strncmp($bigStr, $substrToFind, strlen($substrToFind));
}


function endsWith($bigStr, $substrToFind){
    return (substr($bigStr, -strlen($substrToFind)) === $substrToFind);
}

function dealWithImage($imagePath,$basedirScript,$percent){

    if( file_exists( $imagePath )){  //$illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/duplication_levels.png"

	$tmpfname = tempnam($basedirScript, "TEMP").".png";

	//echo $tmpfname;
	//copying to temp dir

	//$illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/duplication_levels.png",
	copy($imagePath,
	     $tmpfname);
	
	$endname=implode("/",array_slice(explode("/",$tmpfname),-2));
			
	echo "<center><img src=\"".$endname."\" alt=\"image\" height=\"".$percent."%\" width=\"".$percent."%\"> </center> \n"; 
	flush();	
    }

}

function openAndPrint($fileToOpenAndPrint){
    echo "\n<pre>\n";
    $fh = fopen($fileToOpenAndPrint, 'r') or die("can't open $fileToOpenAndPrint");
    while( !feof($fh) ){
	echo fgets($fh) ;						
    }
    fclose($fh);
    echo "\n</pre>\n";
}

function convertToPng($imagePath,$basedirScript,$percent){
    if( file_exists( $imagePath )){  //$illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/duplication_levels.png"

	$tmpfname = tempnam($basedirScript, "TEMP").".png";

	//echo $tmpfname;
	//copying to temp dir

	//$illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/duplication_levels.png",
	$cmdToRun="/usr/bin/convert ".$imagePath." ".$tmpfname;
	$outputStore="";
	$returnCode=1;
	if(!file_exists("/usr/bin/convert")){
	    echo "Install convert on your webserver http://www.imagemagick.org/<br>\n";
	    exit;	    
	}
	exec($cmdToRun,$outputStore,$returnCode);
	if($returnCode != 0){
	    echo "Following command failed:<BR> ".$cmdToRun."<BR><BR>output : ".var_dump($outputStore)." please contact the webmaster <br>\n";
	    // exit;
	}
	
	$endname=implode("/",array_slice(explode("/",$tmpfname),-2));
			
	echo "<center><img src=\"".$endname."\" alt=\"image\" height=\"".$percent."%\" width=\"".$percent."%\"> </center> \n"; 
	flush();	
    }
}

//$xmlconf = simplexml_load_file( getcwd()."/config.xml" );
$jsonconf = json_decode(file_get_contents( getcwd()."/config.json" ),true);


//CONFIG DATA
$illuminareaddir  = $jsonconf["illuminareaddir"];
$illuminawritedir = $jsonconf["illuminawritedir"];
//$illuminajson     = $jsonconf->illuminajson;
$sequencers=array();
foreach($jsonconf["sequencers"]["sequencer"] as $seqelem){
    $sequencers[ (string)$seqelem["id"] ]= (string)$seqelem["type"];
}
$runstodisplay= $jsonconf["runstodisplay"];



$runid=array_keys($_POST);
$runid=$runid[0];

$percentScale = 60;

$basedirScript=dirname(__FILE__)."/images/"; 

				/* print $runid; */

				//exit(1);
?>

				<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN">
				<html lang="en">
    <head>
    <title>Postprocessing report</title>
    </head>
    <body>
    <!-- Progress bar holder -->
    <div id="progress" style="width:500px;border:1px solid #ccc;"></div>
    <!-- Progress information -->
    <div id="information" style="width"></div>

    <?

    //removing previous temp images
$myDirectory = opendir($basedirScript);
while($entryName = readdir($myDirectory)) {
    if(beginsWith($entryName,"TEMP")){
	unlink($basedirScript."/".$entryName);
    }
}

closedir($myDirectory);


/* //    echo     $basedirScript; */
if( !file_exists($basedirScript ) ){
    print "The directory to store the images ".$basedirScript." does not exist";
    exit(1);
}

if( !is_writable($basedirScript ) ){
    print "The directory to store the images ".$basedirScript." ist not writable";
    exit(1);
}


echo "<H1> Report for ".$runid."  </H1>\n<BR><BR>";


if( file_exists($illuminawritedir."/".$runid ) ){
    //print "exist";
    echo "<H2> RTA report </H2>\n<BR><BR>";
    
    

    if( file_exists($illuminawritedir."/".$runid."/Report/images/error_lanes.png")){
	//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";
	
	echo "<H3>PhiX error on a per lane basis</H3><BR>";
	dealWithImage($illuminawritedir."/".$runid."/Report/images/error_lanes.png",$basedirScript,$percentScale);
    }

    for($lane=1;$lane<=8;$lane++){//
	if( file_exists($illuminawritedir."/".$runid."/Report/images/focus_lane".$lane.".png")){
	    echo "<H3>Focus for lane = ".$lane."</H3><BR>";
	    dealWithImage($illuminawritedir."/".$runid."/Report/images/focus_lane".$lane.".png",$basedirScript,$percentScale);
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	if( file_exists($illuminawritedir."/".$runid."/Report/images/focus_lane".$lane.".png")){
	    echo "<H3>Intensity for lane = ".$lane."</H3><BR>";
	    dealWithImage($illuminawritedir."/".$runid."/Report/images/intensity_lane".$lane.".png",$basedirScript,$percentScale);
	}
    }

    echo "<H2> FastQC </H2>\n<BR><BR>";

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Duplication levels for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/duplication_levels.png",$basedirScript,$percentScale);
	    }
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Kmer profile for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/kmer_profiles.png",$basedirScript,$percentScale);
	    }
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Per base GC content for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/per_base_gc_content.png",$basedirScript,$percentScale);
	    }
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Per base n content for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/per_base_n_content.png",$basedirScript,$percentScale);
	    }
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Per base quality for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/per_base_quality.png",$basedirScript,$percentScale);
	    }
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Per base sequence content for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/per_base_sequence_content.png",$basedirScript,$percentScale);
	    }
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Per sequence GC content for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/per_sequence_gc_content.png",$basedirScript,$percentScale);
	    }
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Per sequence quality for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/per_sequence_quality.png",$basedirScript,$percentScale);
	    }
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/")){
		//print $illuminawritedir."/".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/";

		echo "<H3>Sequence length distribution for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		dealWithImage($illuminawritedir."".$runid."/".$basecaller."/FastQC/s_".$lane."_sequence_fastqc/Images/sequence_length_distribution.png",$basedirScript,$percentScale);
	    }
	}
    }



    echo "<H2> Sequence tally </H2>\n<BR><BR>\n";
     for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    $targetfileok=$illuminawritedir."/".$runid."/".$basecaller."/QC/clusterTally_".$lane.".OK";
	    $targetfileer=$illuminawritedir."/".$runid."/".$basecaller."/QC/clusterTally_".$lane.".ERROR";


	    if( file_exists($targetfileok) || file_exists($targetfileer)){		
		echo "<H3>Tally for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";
		if( file_exists($targetfileok) && file_exists($targetfileer)){
		    echo "\n<font color=red size=+2>Warning</font>: two files seem to exists, check which is more recent.<BR>";
		    openAndPrint($targetfileok);	
		}else{

		    if( file_exists($targetfileok) && !file_exists($targetfileer)){
			echo "Tally OK<BR>\n";
			openAndPrint($targetfileok);	
		    }else{
			if( !file_exists($targetfileok) && file_exists($targetfileer)){
			    echo "\n<font color=red size=+2>ERROR:there was an error in the processing for this lane.</font><BR>";
			    openAndPrint($targetfileer);	
			}else{
			    echo "Tally unavailable<BR>\n";
			}

		    }
		}
	    }


	}
     }

    echo "<H2> QC scores  </H2>\n<BR><BR>\n";

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){
	    
	    foreach(array("A","C","G","T") as $base){

		$targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/s_".$lane."_sequence".$base.".pdf";
		if( file_exists($targetfile)){
		    echo "<H3>QC scores for ".$basecaller." basecalls, lane = ".$lane." base = ".$base."</H3><BR>";

		    convertToPng($targetfile,$basedirScript,$percentScale);
		}
	    }
	}
    }


    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/qscores/s_".$lane."_control.bam.baseobspred.pdf";
	    if( file_exists($targetfile)){
		echo "<H3>RG quality score correlation ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";		
		convertToPng($targetfile,$basedirScript,$percentScale);
	    }	
	}
    }


    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/qscores/s_".$lane."_control.bam.baseobspred.dens.pdf";
	    if( file_exists($targetfile)){
		echo "<H3>RG quality score correlation ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";		
		convertToPng($targetfile,$basedirScript,$percentScale);
	    }	
	}
    }

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/qscores/s_".$lane."_control.type.pdf";
	    if( file_exists($targetfile)){
		echo "<H3>Error rate post-mapping on the PhiX ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";		
		convertToPng($targetfile,$basedirScript,$percentScale);
	    }

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/qscores/s_".$lane."_control.error.type.pdf";
	    if( file_exists($targetfile)){
		echo "<H3>Error rate post-mapping on the PhiX (per type) ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";		
		convertToPng($targetfile,$basedirScript,$percentScale);
	    }
	}
    }


    echo "<H2> Merger/trimming log  </H2>\n<BR><BR>\n";
    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/s_".$lane."_merge.log";
	    if( file_exists($targetfile)){
		echo "<H3>Merger Log ".$basecaller." basecalls, lane = ".$lane." </H3><BR>";
		openAndPrint($targetfile);	
	    }
	}
    }


    echo "<H2> Read group assignment  </H2>\n<BR><BR>\n";

    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/rg/s_".$lane."_rg_summary.txt";
	    if( file_exists($targetfile)){	    
		echo "<H3>RG assignment summary ".$basecaller." basecalls, lane = ".$lane." </H3><BR>";
		openAndPrint($targetfile);	
	    }
	}
    }

    //insert size
    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/insertsize/lane".$lane."/";
	    if( file_exists($targetfile)){
		$pdflist = array(); 
		$myDirectory = opendir($targetfile);
		
		while($entryName = readdir($myDirectory)) {

		    if($entryName != "." and $entryName != ".." ){
			if(substr($entryName,-3) == "pdf"){ 
			    $rgname=substr($entryName,3,-4);
			    echo "<H3>RG insert size distribution for ".$rgname.", lane = ".$lane."</H3><BR>";		
			    convertToPng($targetfile."/".$entryName,$basedirScript,$percentScale);
			}
		    }
		}
	    }

	    for($proc=1;$proc<=100;$proc++){
		
		$targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/insertsize/proc".$proc."lane".$lane."/";
		
		if( file_exists($targetfile)){
		    $pdflist = array(); 
		    $myDirectory = opendir($targetfile);
		
		    while($entryName = readdir($myDirectory)) {

			if($entryName != "." and $entryName != ".." ){
			    if(substr($entryName,-3) == "pdf"){ 
				$rgname=substr($entryName,3,-4);
				echo "<H3>RG insert size distribution for ".$rgname.", processing ".$proc.", lane = ".$lane."</H3><BR>";		
				convertToPng($targetfile."/".$entryName,$basedirScript,$percentScale);
			    }
			}
		    }

		}else{
		    break;
		}
		
	    }

	}
    }





    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/rg/s_".$lane."_unassigned.txt";
	    if( file_exists($targetfile)){	    
		echo "<H3>Unassigned indices ".$basecaller." basecalls, lane = ".$lane." </H3><BR>";
		openAndPrint($targetfile);	
	    }
	}
    }



    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/rg/s_".$lane."_rgqual.pdf";
	    if( file_exists($targetfile)){
		echo "<H3>RG quality score for ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";		
		convertToPng($targetfile,$basedirScript,$percentScale);
	    }	
	}
    }




    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/rg/s_".$lane."_ratio.pdf";
	    if( file_exists($targetfile)){
		echo "<H3>RG ratio ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";		
		convertToPng($targetfile,$basedirScript,$percentScale);
	    }	
	}
    }




    echo "<H2> Filtering  </H2>\n<BR><BR>\n";
    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/filter/s_".$lane."_filter.log";
	    if( file_exists($targetfile)){	    
		echo "<H3>Filtering log for ".$basecaller." basecalls, lane = ".$lane." </H3><BR>";
		openAndPrint($targetfile);	
	    }

	    for($proc=1;$proc<=100;$proc++){
		$targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/filter/proc".$proc."s_".$lane."_filter.log";
		if( file_exists($targetfile)){	    
		    echo "<H3>Filtering log for ".$basecaller." basecalls, processing ".$proc.", lane = ".$lane." </H3><BR>";
		    openAndPrint($targetfile);	
		}else{
		    break;
		}

	    }
	}
    }


    for($lane=1;$lane<=8;$lane++){//
	foreach(array("Bustard","Ibis") as $basecaller){

	    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/QC/filter/s_".$lane."_likelihood.pdf";
	    if( file_exists($targetfile)){
		echo "<H3>Sequence correctness likelihood ".$basecaller." basecalls, lane = ".$lane."</H3><BR>";		
		convertToPng($targetfile,$basedirScript,$percentScale);
	    }	
	}
    }

    echo "<H2> Mapping results  </H2>\n<BR><BR>\n";


    /* for($lane=1;$lane<=8;$lane++){// */
    /* 	foreach(array("Bustard","Ibis") as $basecaller){ */
    foreach(array("Bustard","Ibis") as $basecaller){

	if( file_exists($illuminawritedir."/".$runid."/".$basecaller."/BWA/")){
	    $myBWADirectory = opendir($illuminawritedir."/".$runid."/".$basecaller."/BWA/");
	    while($entryName = readdir($myBWADirectory)) {
		if(endsWith($entryName,"flgstx")){
		    $targetfile=$illuminawritedir."/".$runid."/".$basecaller."/BWA/".$entryName."";
		    //if( file_exists($targetfile)){
	    
		    echo "<H3>Flagstatx for file".$targetfile." </H3><BR>";
		    openAndPrint($targetfile);
		    //}
		}
	    }
	    closedir($myBWADirectory);
	}
    }
    /* 	} */
    /* } */

}else{
    print "<BR>The analysis of the run was not started<BR>";
}



/* $numberFound=1; */

/* /\* foreach($runlist as $one_file) { *\/ */
/* $percent = intval($numberFound/count($runlist) *100)."%"; */

/* $numberFound++; */
/* /\* }  *\/ */

/* echo '<script language="javascript"> */
/*     document.getElementById("progress").innerHTML="<div style=\"display:none\" >&nbsp;</div>"; */
/*     document.getElementById("information").innerHTML="'.$numberFound.' runs(s) found."; */
/*     </script>'; */
/* echo str_repeat(' ',1024*64); */


//$stringWithTable.="</HTML>";

//echo $stringWithTable;
echo "</HTML>";



?> 