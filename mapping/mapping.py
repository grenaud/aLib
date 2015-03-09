# File: simplehttpserver-example-1.py


try:
    # Python 2.x
    from SocketServer import ThreadingMixIn
except ImportError:
    # Python 3.x
    from socketserver import ThreadingMixIn

#from BaseHTTPServer import BaseHTTPRequestHandler,HTTPServer
import cgi;
import sys;
import thread
import threading;
from optparse import OptionParser
import sys,os
import operator
from threading import Thread, Lock
from  multiprocessing import Process

#from threading import Thread
import time;
#import zmq

#from socket import socket;
from socket import *
import json
import subprocess
#import socket
import tempfile;
import smtplib
from email.mime.text import MIMEText


def chomp(s):
    return s.rstrip('\n');

def timeString():
    return str(time.strftime("%d_%b_%Y_%H:%M:%S", time.localtime()));

def tprint(msg):
    sys.stdout.write(msg + '\n')
    sys.stdout.flush()



#BEGIN READ CONFIG FILE
pathToConfig=sys.argv[0];
pathToConfig=os.path.dirname(os.path.abspath(pathToConfig))+"/../webForm/config.json";
#alibdir=os.path.dirname(os.path.abspath( sys.argv[0]+"/../"))

#print pathToConfig;
#print alibdir;



try:
  fileHandleConf = open ( pathToConfig );
except IOError:
  print "Cannot open configuration json file "+pathToConfig
  sys.exit(1);

jsonstringConf="";
while 1:
  line = fileHandleConf.readline();
  if(not(line)):
    break
  line = chomp(line);
  jsonstringConf+=line;
fileHandleConf.close();

jsondataConf=json.loads(jsonstringConf);
BWAGENOMES       = jsondataConf["genomedirectory"];
alibdir          = jsondataConf["alibdir"];
sender           = jsondataConf["sender"];
smtpserver       = jsondataConf["smtpserver"];

#END READ CONFIG FILE



#############   BEGIN GLOBALS ########
PORT_NUMBERWEB = 8080
PORT_NUMBERCLI = 8088
#PORT_NUMBERZMQ = 8088
#PORT_NUMBERZMQ = range(5550,5650,2) #range(8088,8188,2) 
PORT_NUMBERBWA = 52690

PORT_NUMBERBWA = 52690
SLEEPTIME      = 10

qsubcmd  = "qsub";
qstatcmd = "qstat";
qdelcmd  = "qdel";


#tempDir = "/tmp/";

tempDirnetw = "/tmp/";
tempDirsort = "/tmp/";
tempDirbwa  = "/tmp/"

mutextocount   = Lock();
mutextomap     = Lock();
mutextosort    = Lock();
mutexdone      = Lock();

mutexismapping = Lock();

ismapping=False;

fileNtocount = "";
fileNtomap   = "";
fileNtosort  = "";
fileNdone    = "";

stringNtocount = "";
stringNtomap   = "";
stringNtosort  = "";

context=0;
sockzmq=0;
counterBAM            = "/tileCount/countClustersBAM";



counterBAM = alibdir+"/"+counterBAM;
if not os.path.exists(counterBAM):
  print "Required executable file not found "+counterBAM;
  sys.exit(1);


#############   END GLOBALS ########


                                                                             
def handle_job(cjob):
    #global options

  jobcreated=subprocess.Popen(cjob,shell=True,
                           stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE); 
  jobcreated.wait()

  out, err = jobcreated.communicate()
  errcode  = jobcreated.returncode;

  if(errcode != 0): #has finished but wrong code
    print "Job failed "+cjob+" failed";
    sys.exit(1);

  return out;


def handle_job_print(cjob):
    #global options

  jobcreated=subprocess.Popen(cjob,shell=True,
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE); 

  while jobcreated.poll() is None:
    l = jobcreated.stderr.readline(); # This blocks until it receives a newline.
    tprint(l);
  #jobcreated.wait()

  #out, err = jobcreated.communicate()
  errcode  = jobcreated.returncode;

  if(errcode != 0): #has finished but wrong code
    print "Job failed "+cjob+" failed";
    sys.exit(1);

  return out;
  


#This class will handles any incoming request from
#the browser 
class myHandler(BaseHTTPRequestHandler):
    def do_POST(self):
        form = cgi.FieldStorage(
            fp=self.rfile,
            headers=self.headers,
            environ={'REQUEST_METHOD':'POST',
                     'CONTENT_TYPE':self.headers['Content-Type'],
                     })
        #filename = form['file'].filename
        #data = form['file'].file.read()
        #open("/tmp/%s"%filename, "wb").write(data)        
        mygenomefile=form['genomebwa'].value;
        if(mygenomefile == "none"):
            mygenomefile=form['genomebwatext'].value;
            if( not(os.path.isfile(mygenomefile+".ann") ) ):
                self.respond("ERROR: the BWA file "+str(mygenomefile+".ann")+" should exist");
        else:
            mygenomefile=BWAGENOMES+"/"+form['genomebwa'].value+"/bwa-0.4.9";

        #for field in form.keys():
        #field_item = form[field];
        #self.wfile.write('\t%s=%s\n' % (field, form[field].value));

        if( len(str(form['input'].value)) == 0  ):
            self.respond("ERROR: the input file should be specified");

        jobnametouse="";
        if( len(str(form['jobname'].value)) == 0  ):
            #self.respond("ERROR: the input file should be specified");
            jobnametouse = str(int(time.time()));
        else:
            jobnametouse = str(form['jobname'].value);

        if( len(str(form['output'].value)) < 5  ):
            self.respond("ERROR: the output file should be specified and be at least 5 characters and end with .bam");

        if( str(form['output'].value)[-4:] != ".bam"  ):
            self.respond("ERROR: the output file should be specified and be at least 5 characters and end with .bam");

        if( not(os.path.isfile( str(form['input'].value) ) ) ):
            self.respond("ERROR: the input BAM "+str(form['input'].value)+" should exist");


        mutextocount.acquire()
        

        fileHandle = open (fileNtocount );
        

        fileHandle = open( fileNtocount , 'a' );
        fileHandle.write( jobnametouse+"\t"+str(form['input'].value) +"\t"+str(form['output'].value)+"\t"+str(mygenomefile) +"\t"+str(form['defaultparam'].value) +"\t"+str(form['email'].value)+"\t"+timeString()+"\n" );
        fileHandle.close();


        mutextocount.release()
        
        if( str(form['email'].value) == "none"):
            self.respond("thank you, please check again when you job is done.");
        else:
            self.respond("thank you, an email will be sent to "+str(form['email'].value)+" one processing is done.");

        #for field in form.keys():
        #    field_item = form[field];
        #    self.wfile.write('\t%s=%s\n' % (field, form[field].value));
        #self.respond("thank you, an email will be sent to "+str(form['email'].value))+" one processing is done";
#<form>
    def do_GET(self):
        global mutexdone;
        response = """
        <html><body>
<form enctype="multipart/form-data" method="post">
Jobname:<br>
<input type="text" name="jobname">
<br>
Email (leave to none if no email is required, put commas to separate multiple emails):<br>
<input type="text" name="email"  value="none">
<br>

Genome:<br>\n""";
        response = response+"<select name=\"genomebwa""\">\n";
        response = response+"<option value=\"none\"></option>\n";

        for filename in os.listdir(BWAGENOMES):
            filename=BWAGENOMES+"/"+filename;
            #response = response+str(filename);
            if(os.path.isfile(filename+"/bwa-0.4.9.ann")):
                #response = response+str(filename);
                response = response+"<option value=\""+str(filename[(len(BWAGENOMES)+1):])+"\">"+str(filename[(len(BWAGENOMES)+1):])+"</option>\n";
        response = response+"""
        </select>\n
        <BR>
or specify the genome file manually:<BR>
<input type="text" name="genomebwatext">
<br><br><br>
Input file:<BR>
<input type="text" name="input"><BR>
Output file<BR>
<input type="text" name="output"><BR>
<BR><BR>
<input type=\"radio\" name=\"defaultparam\" value=\"False\" checked>default \n
<input type=\"radio\" name=\"defaultparam\" value=\"True\" >ancient<br>\n
<br><br>
<input type="submit" value="Submit">

</form>""";

        mutextocount.acquire();
        
        fileHandle = open (fileNtocount );
        linesCount = fileHandle.readlines();
        fileHandle.close();

        mutextocount.release()


        mutextomap.acquire();
        
        fileHandle = open (fileNtomap );
        linesMAP = fileHandle.readlines();
        fileHandle.close();

        mutextomap.release();


        mutextosort.acquire();
        
        fileHandle = open (fileNtosort );
        linesSort = fileHandle.readlines();
        fileHandle.close();

        mutextosort.release();


        mutexdone.acquire();
        
        fileHandle = open (fileNdone );
        linesDone = fileHandle.readlines();
        fileHandle.close();

        mutexdone.release();


        response = response+"<h2>Jobs currently counting:</h2><BR><pre>"+str( stringNtocount )+"</pre><BR>"+"<h2>Jobs queued to count:</h2><BR><pre>"+str( '\n'.join( linesCount ) )+"</pre><BR>"+"<h2>Jobs currently mapping:</h2><BR><pre>"+str( stringNtomap  )+"</pre><BR>"+"<h2>Jobs queued to map:</h2><BR><pre>"+str( '\n'.join( linesMAP ) )+"</pre><BR>"+"<h2>Jobs currently sorting:</h2><BR><pre>"+str( stringNtosort  )+"</pre><BR>"+"<h2>Jobs queued to sort:</h2><BR><pre>"+str( '\n'.join( linesSort ) )+"</pre><BR>"+"<h2>Jobs done:</h2><BR><pre>"+str( '\n'.join( linesDone ) )+"</pre><BR>"+"</body></html>";

        self.respond(response)

    def respond(self, response, status=200):
        self.send_response(status)
        self.send_header("Content-type", "text/html")
        self.send_header("Content-length", len(response))
        self.end_headers()
        self.wfile.write(response)  
	
    #Handler for the GET requests
    #def do_GET(self):
        #self.send_response(200);
	#self.send_header('Content-type','text/html');
	#self.end_headers();
	# Send the html message
	#self.wfile.write("Hello World !");
        #self.wfile.write("""<html><body><form enctype="multipart/form-data" method="post"><p>File: <input type="file" name="file"></p><p><input type="submit" value="Upload"></p></form></body></html>""");
	#return

class ThreadingHTTPServer(ThreadingMixIn, HTTPServer):
    print """Handle requests in a separate thread."""
    pass

#
#def runserver():
#    try:
#        #Create a web server and define the handler to manage the
#        #incoming request
#        #server = HTTPServer(('', PORT_NUMBER), myHandler);
#        server = ThreadingHTTPServer(('', PORT_NUMBER), myHandler)
#
#        print 'Started httpserver on port ' , PORT_NUMBER;
#
#        #Wait forever for incoming htto requests
#        server.serve_forever();
#
#    except KeyboardInterrupt:
#        print '^C received, shutting down the web server'
#        server.socket.close();

server=0;
def runserver():
    global server;
        #Create a web server and define the handler to manage the
        #incoming request
        #server = HTTPServer(('', PORT_NUMBER), myHandler);
    server = ThreadingHTTPServer(('', PORT_NUMBERWEB), myHandler)

    tprint('Started httpserver on port ' +str( PORT_NUMBERWEB));

    #Wait forever for incoming http requests
    server.serve_forever();
    tprint( "done");
    #except KeyboardInterrupt:
    #    print '^C received, shutting down the web server'
    #    server.socket.close();


#runserver();

#try:
#    thread.start_new_thread( runserver , ()) 
#except :
#   print "Error: unable to start thread"
#   print "Unexpected error:", sys.exc_info()[0]
#   raise

#def zeromq( threadName ):


#mutex = Lock();



#def print_time( threadName, delay):
#   count = 0
#   mutex.acquire()
#   while count < 5:
#      time.sleep(delay)
#      count += 1
#      #print count;
#      print "%s %s %s %s" % ( count, delay, threadName, time.ctime(time.time()) )
#   mutex.release()

def returnFirstLines(fileToOpen ):

    fileHandle = open (fileToOpen );
    lines = fileHandle.readlines();
    fileHandle.close();

    lines=map(chomp,lines);

    if(len(lines)>0):
        return lines[0];
    else:
        return "";


def writeFileWithoutLine(fileToOpen,lineToRemove ):

    fileHandle = open (fileToOpen );
    lines = fileHandle.readlines();
    fileHandle.close();

    lineToRemove=chomp(lineToRemove);
    lines=map(chomp,lines);

    #if(len(lines)>0):
    foundLine=False;
    fileHandle = open( fileToOpen , 'w' );
    for l in lines:

        if(l==lineToRemove):
            foundLine=True;
        else:
            fileHandle.write( str(l)+"\n" );
    fileHandle.close();
    if(not(foundLine)):
        sys.stderr.write("ERROR line "+str(lineToRemove)+" was not found in file "+str(fileToOpen)+"\n");
        sys.exit(1);

    return True;
    #else:

def appendLineToFile(fileToOpen,lineToAdd ):
    lineToAdd = chomp(lineToAdd);
    fileHandle = open( fileToOpen , 'a' );
    fileHandle.write( lineToAdd+'\n' ); 
    fileHandle.close();

def isRecordAlreadyInFile(fileToOpen,lineToCheck,offsetInd,checkid ):
    fileHandle = open (fileToOpen );
    lines = fileHandle.readlines();
    fileHandle.close();

    listToCheck=lineToCheck.split("\t");
    for l in lines:
        lTemp=l.split("\t");
        if( (lTemp[offsetInd+1] == listToCheck[1] ) and #input  file
            (lTemp[offsetInd+2] == listToCheck[2] ) and #output file
            (lTemp[offsetInd+3] == listToCheck[3] ) and #genome file
            (lTemp[offsetInd+4] == listToCheck[4] ) ):  #ancient param
            if( checkid ):
                if( (lTemp[offsetInd+0] == listToCheck[0] ) ): #need to the be same ID
                    return True;
            else:
                return True;
    return False;
    
    

def handler(clientsocket, clientaddr):
    global mutextocount;
    global mutextomap;
    global mutextosort;


    while True :
        data = clientsocket.recv(1024)
        if not data:
            break
        else:

            tprint("Accepted connection from: "+str(clientaddr)+" with data: "+str(data));
            listTempFields=data.split("\t");
            cmdtoadd="\t".join( listTempFields[1:] );
            waitForJobs = (listTempFields[0] == "True");

            foundPreviously=False
            #check if already there in count
            mutextocount.acquire()
            foundPreviously = foundPreviously or isRecordAlreadyInFile(fileNtocount,cmdtoadd,0,False);
            mutextocount.release()

            #check if already there in map
            mutextomap.acquire()
            foundPreviously = foundPreviously or isRecordAlreadyInFile(fileNtomap,cmdtoadd,1,False);
            mutextomap.release()

            #check if already there in sort
            mutextosort.acquire()
            foundPreviously = foundPreviously or isRecordAlreadyInFile(fileNtosort,cmdtoadd,1,False);
            mutextosort.release()

            newlyAdded=True;
            if( not(foundPreviously) ):
                mutextocount.acquire();
                appendLineToFile(fileNtocount,cmdtoadd);
                mutextocount.release();
                newlyAdded=True;
            else:
                tprint("job "+str(data)+" already exists");
                newlyAdded=False;

            if( waitForJobs ):
                while True:
                    #tprint("newlyAdded "+str(newlyAdded));
                    mutexdone.acquire()
                    if( isRecordAlreadyInFile(fileNdone,cmdtoadd,1,newlyAdded) ): #must match the ID
                        mutexdone.release();
                        break;
                    else:
                        mutexdone.release();
                    time.sleep(300); #check back in 5 minutes
                msg = "Command "+data+" is done";
                try:
                    clientsocket.send(msg);
                except socket.timeout:
                    tprint("Timeout from : "+str(clientaddr)+" ");

            #time.sleep(10);
    try:
        clientsocket.close();
    except socket.timeout:
        tprint("Timeout from : "+str(clientaddr)+" ");

serversocket=0;
#if __name__ == "__main__":
def  threadedServer(threadName, delay):
    global serversocket;
    host = 'localhost'
    buf = 1024
    tprint("Socket server is listening for connections on port "+str(PORT_NUMBERCLI)+"");
    addr = (host, PORT_NUMBERCLI);

    try:
        serversocket = socket(AF_INET, SOCK_STREAM);
        serversocket.bind(addr);
        serversocket.listen(2);
    except socket.error:
        sys.stderr.write("ERROR socket already in use, check programs or wait a few minutes\n");
        sys.exit(1);

    while True:
        

        clientsocket, clientaddr = serversocket.accept()
        thread.start_new_thread(handler, (clientsocket, clientaddr))
    serversocket.close()


#def zeromqserver(port="5556"):
#    
#    context = zmq.Context();
#    sockzmq = context.socket(zmq.REP);
#    sockzmq.bind("tcp://*:%s" % port);
#
#    print "Running server on port: ", port
#    while True:
#        message = sockzmq.recv();
#        tprint("Received request: %s" % message)
##        #response = handle_message(message);
##        #check the queues
##
##
##        #if the mapping is done, 
#        #sockzmq.send("bye");
##
#    # serves only 5 request and dies
#    #for reqnum in range(5):
#        # Wait for next request from client
#     #   message = socket.recv()
#      #  print "Received request #%s: %s" % (reqnum, message)
#        #socket.send("World from %s" % port)
#         
#
#
#def zeromqserver( ):
##try:
#    global sockzmq;
#    global context;
#    global mutextocount;
#
#    tprint("starting zeromq on port "+str(PORT_NUMBERZMQ));
#    context = zmq.Context.instance()
#
#    sockzmq = context.socket(zmq.REP)
#    sockzmq.connect("tcp://localhost:"+str(PORT_NUMBERZMQ))
#
#    while True:
#        message = sockzmq.recv();
#        tprint("Received request: %s" % message)
#        #response = handle_message(message);
#        #check the queues
#
#
#        #if the mapping is done, 
#        sockzmq.send("bye");
#


def counter( ):
    global stringNtocount;
    global mutextocount;
    global mutextomap;

    tprint ("counter thread running");
    while True:
        #check for new mapping
        foundjob=False;
        linebamtocount=""


        #CHECKING FOR NEW JOBS IN COUNTER QUEUE
        mutextocount.acquire()
        #time.sleep(SLEEPTIME);

        linebamtocount=returnFirstLines(fileNtocount);
        #print "linebamtocount #"+linebamtocount+"#"+fileNtocount+"#";
        if(len(linebamtocount)>0):
            foundjob=True;

        fileHandle = open (fileNtocount );

        mutextocount.release()


        
        if(not(foundjob)):
            time.sleep(SLEEPTIME);
        else:
            #print "counter: #"+linebamtocount+"#";
            stringNtocount= linebamtocount;
            bamfiletocount=linebamtocount.split("\t")[1];
            #RUNNING
            cmd = counterBAM+" "+bamfiletocount;
            tprint ("counter: "+str(cmd)+"");

            countBAM=handle_job(cmd);
            countBAM = int ( countBAM );

            mutextocount.acquire()
            #UPDATING COUNTER QUEUE
            stringNtocount= "";
            writeFileWithoutLine(fileNtocount,linebamtocount);
            mutextocount.release()


            #UPDATING MAPPING QUEUE            
            mutextomap.acquire()

            #read previous
            #print "linebamtocount #"+str(linebamtocount)+"# #"+str(countBAM)+"#";
            mydictionary = {};
            mydictionary[linebamtocount]=countBAM; #adding new record to dictionary

            
            fileHandle = open ( fileNtomap );#reading  previous to dictionary
            while 1:
                line = fileHandle.readline();
                line = chomp(line);
                if(not(line)):
                    break                
                listemp=line.split("\t");
                mydictionary[ "\t".join( listemp[1:] ) ]=int(listemp[0]);
            fileHandle.close();
            #print "my dict #"+str(mydictionary)+"#";
            sorted_dict = sorted(mydictionary.items(), key=operator.itemgetter(1));
            
            fileHandle = open( fileNtomap , 'w' );
            for filesBAM in sorted_dict:
                fileHandle.write(  str(filesBAM[1])+"\t"+str(filesBAM[0])+"\n" );
            fileHandle.close();


            mutextomap.release()





def mapper( ):
    global stringNtomap;
    global ismapping;
    global mutexismapping;
    global mutextomap;
    global mutextosort;

    tprint ("mapper thread running");
    while True:
        #check for new mapping
        foundjob=False;
        jobbamtomap=""

        #CHECKING FOR NEW JOBS IN MAPPING QUEUE
        mutextomap.acquire()
        jobbamtomap=returnFirstLines(fileNtomap);
    
        if(len(jobbamtomap)>0):
            foundjob=True;

        mutextomap.release()

        
        if(not(foundjob)):
            time.sleep(SLEEPTIME);
        else:
            #ismapping for launcher
            mutexismapping.acquire();
            ismapping=True;
            mutexismapping.release();
            


            stringNtomap   = jobbamtomap;

            #BEGIN MAPPING
            listInfomapper=jobbamtomap.split("\t");

            cmd ="bwa bam2bam -t 3 ";
            cmd += " -p "+str(PORT_NUMBERBWA)+" ";
            if( str(listInfomapper[5]) == "True"):
                cmd += " -n 0.01 -o 2 -l 16500  ";     #bwa param
            cmd += " -g  "+str(listInfomapper[4])+" "; #genome
            cmd += " -f  "+str(listInfomapper[3])+"_temp_ "; #output
            cmd += " --temp-dir "+str(tempDirbwa)+" "; #temp dir
            cmd += "  "+str(listInfomapper[2])+" ";    #input
            tprint ("mapper running :"+str(cmd));
            handle_job(cmd);
            #end mapping

            stringNtomap   = "";

            mutexismapping.acquire();
            ismapping=False;
            mutexismapping.release();
            #deleting jobs
            deletemyjobs();


            #UPDATING MAPPING QUEUE            
            mutextomap.acquire()
            writeFileWithoutLine(fileNtomap,jobbamtomap);
            mutextomap.release()

            #UPDATE SORTING QUEUE
            mutextosort.acquire()
            appendLineToFile(fileNtosort,jobbamtomap);
            mutextosort.release()



def sorter( ):
    global stringNtosort;
    global mutextosort;
    global mutexdone;

    tprint ("sorter thread running");
    while True:
        #check for new mapping
        foundjob=False;
        bamfiletosort="";
        

        #READ SORTING QUEUE
        mutextosort.acquire()
        bamfiletosort=returnFirstLines(fileNtosort);

        if(len(bamfiletosort)>0):
            foundjob=True;

        mutextosort.release()
        
        if(not(foundjob)):
            time.sleep(SLEEPTIME);
        else:
            #tprint "cmd sort";

            #BEGIN SORTING
            stringNtosort = bamfiletosort;
            listInfosorter=bamfiletosort.split("\t");
            cmd = "sam sort -m 2G ";
            cmd+=" -O "+str(listInfosorter[3])[:-4]+" ";
            cmd+=" "+str(listInfosorter[3])+"_temp "+tempDirsort;
            tprint ("sorter "+str(cmd));
            handle_job(cmd);
            #if(str(listInfosorter[3]) 
            #print str(listInfosorter[3])[:-4]+".__bam";
            #move if wrote to another file
            if( not(os.path.isfile(  str(listInfosorter[3]) ) ) and
                (os.path.isfile(  str(listInfosorter[3])[:-4]+".__bam" ))  ):
                cmd = "mv "+str(listInfosorter[3])[:-4] + ".__bam " +  str(listInfosorter[3]);
                #print cmd;
                os.rename(str(listInfosorter[3])[:-4]+".__bam", str(listInfosorter[3]) );
                #handle_job(cmd);


            cmd="rm -fv "+str(listInfosorter[3])+"_temp ";
            #print cmd;
            os.remove( str(listInfosorter[3])+"_temp" );                
            #handle_job(cmd);
            stringNtosort = "";
    

            #UPDATE SORTING QUEUE

            mutextosort.acquire();
            writeFileWithoutLine(fileNtosort,bamfiletosort);
            mutextosort.release();
            #end sorting
            #write to done, send email
            #SEND email if not "none"
            if str(listInfosorter[6]) != "none":
                text =  "The mapping you requested is finished\n\n";
                text += "Input:\t"+str(listInfosorter[2])+"\n";
                text += "Output:\t"+str(listInfosorter[3])+"\n";
                text += "Genome:\t"+str(listInfosorter[4])+"\n";
                text += "Ancient parameters:\t"+str(listInfosorter[5])+"\n";
                text += "Requested:\t"+str(listInfosorter[7])+"\n";
                text += "Done:\t"+str(timeString())+"\n";

                msg            = MIMEText(text);
                msg['Subject'] = 'Mapping finished';
                msg['From']    = sender;
                msg['To']      = str(listInfosorter[6]); 

                s = smtplib.SMTP(smtpserver);
                s.set_debuglevel(1);
                s.sendmail(sender, str(listInfosorter[6]).split(','), msg.as_string());
                s.quit();

            #UPDATE DONE QUEUE
            mutexdone.acquire()
            tprint ("Finished: "+str(listInfosorter[3]));
            appendLineToFile(fileNdone, bamfiletosort+"\t"+timeString() );            
            mutexdone.release()



def jobsAreAllRunning():
    cmd = str(qstatcmd)+" ";
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (stdout, stderr) = p.communicate()
    pstat = p.wait();
    #print "done";
    if(pstat != 0):
        print "Error cmd returned a non zero code :"+cmd;
        sys.exit(1);
        #print stdout;
    found=False;
    runcode = "done";
    for line in stdout.split('\n'):
        fields= line.split();
        if(len(fields) < 3):
            continue;

        if(fields[0] == "job-ID"):
            continue;

        if(fields[2] == "alib"):
            if(fields[4] == "qw"): #code for queued
                return False;
            
    return True;

def deletemyjobs():
    cmd = str(qstatcmd)+" ";
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (stdout, stderr) = p.communicate()
    pstat = p.wait();
    #print "done";
    if(pstat != 0):
        print "Error cmd returned a non zero code :"+cmd;
        sys.exit(1);
        #print stdout;
    found=False;
    runcode = "done";
    for line in stdout.split('\n'):
        fields= line.split();
        if(len(fields) < 3):
            continue;

        if(fields[0] == "job-ID"):
            continue;

        if(fields[2] == "alib"):
            cmddel = str(qdelcmd)+" "+str(fields[0]);
            #print "killing "+str(fields[0]);
            handle_job(cmddel);

            

def runQsub():    
    global tempDirnetw;
    # -o /dev/null  -e /dev/null
    #temp = tempfile.NamedTemporaryFile(prefix=options.tmp+"script",suffix=".sge",delete=False)
    cmd= " echo \"/home/public/usr/bin/bwa worker -T 20000 -t \$NSLOTS -p "+str(PORT_NUMBERBWA)+" -h "+str( (gethostname()) )+"\" | "+str(qsubcmd)+"   -N alib -S /bin/bash -l \"class=*,h_vmem=6.8G,s_vmem=6.8G,virtual_free=6.8G \" -V  -pe smp 1- -e "+str(tempDirnetw)+" -o "+str(tempDirnetw)+" ";
    #tprint( cmd);
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (stdout, stderr) = p.communicate()
    pstat = p.wait();

    if(pstat != 0):
        print "Error cmd returned a non zero code :"+cmd;
        sys.exit(1);

    fields= stdout.split();
    if(len(fields) < 3):
        print "Error cmd did not have the expected format : "+cmd;
        sys.exit(1);

    return fields[2];

def launcher( ):
    global ismapping;
    global mutexismapping;

    tprint ("launcher thread running");

    while True:
        #check for new mapping
        #imappingLocal=False;
        mutexismapping.acquire();


        #print "ismapping "+str(ismapping);
        if(ismapping):
            mutexismapping.release();
            #print "trying to launch";
            #TODO add qstat and qsub commands
            if( jobsAreAllRunning() ):
                #launch more
                #print "All jobs running, launching more";
                runQsub();
            else:
                #print "All jobs queued, sleeping";
                time.sleep(SLEEPTIME);
                continue;
            #if qstat says all running
            #   launch jobs SGE (socket.gethostname())
            #print (socket.gethostname())

        else:            
            mutexismapping.release();
            #not mapping, go to sleep
            time.sleep(SLEEPTIME);



#3 files: tocount, tomap, tosort

#1 thread webserver writes to tocount
#2 thread mapping
#   reads from tomap
#   runs 
#   stores in file which runs to sort
#3 thread launching jobs
#4 sorts reads from tocount
#5 zeromq writes to tocount
#6 counts, reads from tocount write to tomap


parser = OptionParser(usage="usage: %prog [options]");
parser.add_option("-t", "--temp",       dest="tempnetwork", help="Temporary directory for the network",default="/tmp/");
parser.add_option("-l", "--templocal",  dest="tempsamsort", help="Temporary directory on local server for sam tools sort",default="/tmp/");
parser.add_option("-d", "--dirlocal",   dest="tempbam2bam",  help="Temporary directory on local server for bwa bam2bam",default="/tmp/");

parser.add_option("--qstat",  dest="qstat", help="qstat location [default: %default]",default="/opt/sge/bin/lx-amd64/qstat");
parser.add_option("--qsub",  dest="qsub", help="qsub location [default: %default]",default="/opt/sge/bin/lx-amd64/qsub");
parser.add_option("--qdel",  dest="qdel", help="qdel location [default: %default]",default="/opt/sge/bin/lx-amd64/qdel");

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)


(options, args) = parser.parse_args();

tempDirnetw = options.tempnetwork;
tempDirsort = options.tempsamsort;
tempDirbwa  = options.tempbam2bam;
#print tempDirnetw;
#print tempDirsort;
#print tempDirbwa;

#sys.exit(1);
qsubcmd  = options.qsub;
qstatcmd = options.qstat;
qdelcmd  = options.qdel;

fileNtocount = tempDirbwa+"/alib.tocount"
fileNtomap   = tempDirbwa+"/alib.tomap"
fileNtosort  = tempDirbwa+"/alib.tosort"
fileNdone    = tempDirbwa+"/alib.done"

if not os.path.exists(fileNtocount):
    fileHandleWrite = open ( fileNtocount, 'w' ) ;
    fileHandleWrite.write("");
    fileHandleWrite.close();

if not os.path.exists(fileNtomap):
    fileHandleWrite = open ( fileNtomap, 'w' ) ;
    fileHandleWrite.write("");
    fileHandleWrite.close();

if not os.path.exists(fileNtosort):
    fileHandleWrite = open ( fileNtosort, 'w' ) ;
    fileHandleWrite.write("");
    fileHandleWrite.close();

if not os.path.exists(fileNdone):
    fileHandleWrite = open ( fileNdone, 'w' ) ;
    fileHandleWrite.write("");
    fileHandleWrite.close();



#MAIN
try:

    

    tweb = threading.Thread(target = runserver);
    tweb.daemon = True
    tweb.start()

    #if False:
    tcounter = threading.Thread(target = counter );
    tcounter.daemon = True
    tcounter.start()        

    tmap = threading.Thread(target = mapper );
    tmap.daemon = True
    tmap.start()    

    tlauncher = threading.Thread(target = launcher );
    tlauncher.daemon = True
    tlauncher.start()    

    tsorter = threading.Thread(target = sorter );
    tsorter.daemon = True
    tsorter.start()    


    thread.start_new_thread(threadedServer,("Thread-2", 1, ))

    #tzeromqserver = threading.Thread(target = zeromqserver );
    #tzeromqserver.daemon = True
    #tzeromqserver.start()    

#    for server_port in PORT_NUMBERZMQ:
 #       Process(target=zeromqserver, args=(server_port,)).start()
    #tc2 = threading.Thread(target = print_time,args=("Thread-2", 1, ) );
    #tc2.daemon = True
    #tc2.start()    

    while True: 
        tprint ("program alive at "+timeString());
        time.sleep(3600)
except (KeyboardInterrupt, SystemExit):
    tprint ("\n! Received keyboard interrupt, quitting threads.\n")
    tprint ("deleting jobs");
    deletemyjobs();
    tprint ("closing webserver server");
    server.socket.close();
    #tprint ("closing zmq socket");
    #sockzmq.close();
    #tprint ("closing zmq context");
    #context.term();
    tprint ("closing socket server");
    serversocket
    tprint ("shutdown server done");
    
