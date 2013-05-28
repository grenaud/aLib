<?php

require "email.php";

ini_set('display_errors', 'On');
error_reporting(E_ALL);

//arguments: [email to send, comma separated] [runid] [runfolder] [lane]


$emailAddrToSend = $argv[1];
$runid           = $argv[2];
$runfolder       = $argv[3];
$lane            = $argv[4];

$mail = new EMail;
$mail->Username = 'sbsuser';
$mail->Password = 'sbs123';

foreach(explode(",",$emailAddrToSend) as $emailAddrTo){ //for each email in the array above
    $mail->AddTo($emailAddrTo,""); // Name is optional
}


    
$mail->SetFrom("sbsuser@eva.mpg.de","");  // Name is optional

$mail->Subject = "Processing finished ".$runid;
$mail->Message = "This is an automated message.\nThe following run :\n". 
    $runfolder."\n".
    "For the following lane: ".$lane."\n".
    "is done processing\nPlease refer to the webform for quality control information\n";
    
$mail->ConnectTimeout = 30;  // Socket connect timeout (sec)
$mail->ResponseTimeout = 8;  // CMD response timeout (sec)
$success = $mail->Send();
if($success != 1){
    echo "the admin notification email WAS NOT sent successfully\n";
}else{
    echo "the admin notification email was sent successfully\n";
}
    

    
    

?>
