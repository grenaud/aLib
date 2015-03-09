# File: simplehttpserver-example-1.py


try:
    # Python 2.x
    from SocketServer import ThreadingMixIn
except ImportError:
    # Python 3.x
    from socketserver import ThreadingMixIn

from BaseHTTPServer import BaseHTTPRequestHandler,HTTPServer
import cgi;
import sys;
import threading;
from optparse import OptionParser
import sys,os
import operator
from threading import Thread, Lock
#from threading import Thread
import time;
#import zmq
import json
import subprocess
import socket
import tempfile;
import smtplib
from email.mime.text import MIMEText


def chomp(s):
    return s.rstrip('\n');

def timeString():
    return str(time.strftime("%d_%b_%Y_%H:%M:%S", time.localtime()));





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
PORT_NUMBERZMQ = 8080
PORT_NUMBERBWA = 52690
SLEEPTIME      = 10

qsubcmd  = "qsub";
qstatcmd = "qstat";
qdelcmd  = "qdel";


tempDir = "/tmp/";

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
    print l;
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

    print 'Started httpserver on port ' , PORT_NUMBERWEB;

    #Wait forever for incoming http requests
    server.serve_forever();
    print "done";
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

def zeromqserver( ):
#try:
    context = zmq.Context.instance()

    sockzmq = context.socket(zmq.REP)
    sockzmq.connect("tcp://localhost:"+str(PORT_NUMBERZMQ))

    while True:
        message = sockzmq.recv();
        print("Received request: %s" % message)
        #response = handle_message(message);
        sockzmq.send("bye");



def counter( ):
    global stringNtocount;
    print "counter thread running";
    while True:
        #check for new mapping
        foundjob=False;
        linebamtocount=""
        mutextocount.acquire()
        #time.sleep(SLEEPTIME);


        fileHandle = open (fileNtocount );

        lines = fileHandle.readlines();
        fileHandle.close();
        #print "lines #"+str(lines)+"#";
        lines=map(chomp,lines);
        #print "lines #"+str(lines)+"#";        
                  
        if(len(lines)>0):
            foundjob=True;
            linebamtocount=chomp(lines[0]);

            #print "line 0 "+str(linebamtocount);
            #print "line 1 "+str( '\n'.join( lines[1:] ) );

            fileHandle = open( fileNtocount , 'w' );
            fileHandle.write( ('\n'.join( lines[1:] )) );
            fileHandle.close();

        #print "length #"+str(len(lines))+"#";
        mutextocount.release()
        
        if(not(foundjob)):
            time.sleep(SLEEPTIME);
        else:
            #print "counter: #"+linebamtocount+"#";
            stringNtocount= linebamtocount;
            bamfiletocount=linebamtocount.split("\t")[1];

            cmd = counterBAM+" "+bamfiletocount;
            print "counter: "+str(cmd)+"";

            countBAM=handle_job(cmd);
            countBAM = int ( countBAM );

            stringNtocount= "";
            mutextomap.acquire()

            #read previous

            mydictionary = {};
            mydictionary[linebamtocount]=countBAM; #adding new record
            #print "dict1 "+str(mydictionary);

            fileHandle = open ( fileNtomap );#read previous to map
            while 1:
                line = fileHandle.readline();
                if(not(line)):
                    break
                line = chomp(line);
                #print "line map #"+line+"#";
                listemp=line.split("\t");
                mydictionary[ "\t".join( listemp[1:] ) ]=int(listemp[0]);
            fileHandle.close();
            #print mydictionary;
            #print "dict2 "+str(mydictionary);

            sorted_dict = sorted(mydictionary.items(), key=operator.itemgetter(1));
            #print "counter "+str(sorted_dict);
            fileHandle = open( fileNtomap , 'w' );
            for filesBAM in sorted_dict:
                fileHandle.write(  str(filesBAM[1])+"\t"+str(filesBAM[0])+"\n" );
            fileHandle.close();

            #print str(tuple1)+"\t"+str(mydictionary[tuple1]);

            mutextomap.release()





def mapper( ):
    global stringNtomap;
    global ismapping;

    print "mapper thread running";
    while True:
        #check for new mapping
        foundjob=False;
        jobbamtomap=""

        mutextomap.acquire()
        
        #print "mapping reading file"

        fileHandle = open (fileNtomap );
        lines = fileHandle.readlines();
        fileHandle.close();
        lines=map(chomp,lines);

        if(len(lines)>0):
            foundjob=True;
            jobbamtomap=chomp(lines[0]);
            
            fileHandle = open( fileNtomap , 'w' );
            fileHandle.write( ('\n'.join( lines[1:] )) );
            fileHandle.close();


        mutextomap.release()
        
        if(not(foundjob)):
            time.sleep(SLEEPTIME);
        else:
            mutexismapping.acquire();
            ismapping=True;
            mutexismapping.release();
            


            stringNtomap   = jobbamtomap;

            #Begin mapping
            listInfomapper=jobbamtomap.split("\t");
            #print "mapper :#"+str(jobbamtomap)+"#";

            cmd ="bwa bam2bam -t 3 ";
            cmd += " -p "+str(PORT_NUMBERBWA)+" ";
            if( str(listInfomapper[5]) == "True"):
                cmd += " -n 0.01 -o 2 -l 16500  ";     #bwa param
            cmd += " -g  "+str(listInfomapper[4])+" "; #genome
            cmd += " -f  "+str(listInfomapper[3])+"_temp_ "; #output
            cmd += " --temp-dir "+str(tempDir)+" "; #temp dir
            cmd += "  "+str(listInfomapper[2])+" ";    #input
            print "mapper running :"+str(cmd);
            handle_job(cmd);
            #end mapping

            stringNtomap   = "";

            mutexismapping.acquire();
            ismapping=False;
            mutexismapping.release();

            deletemyjobs();

            

            mutextosort.acquire()
            fileHandle = open( fileNtosort , 'a' );
            fileHandle.write( jobbamtomap+'\n' ); 
            fileHandle.close();
            mutextosort.release()



def sorter( ):
    global stringNtosort;
    print "sorter thread running";
    while True:
        #check for new mapping
        foundjob=False;
        bamfiletosort=""

        mutextosort.acquire()
        
        #print "sorter reading file"

        fileHandle = open (fileNtosort );
        lines = fileHandle.readlines();
        fileHandle.close();

        lines=map(chomp,lines);

        if(len(lines)>0):
            foundjob=True;
            bamfiletosort=chomp(lines[0]);
            
            fileHandle = open( fileNtosort , 'w' );
            fileHandle.write( ('\n'.join( lines[1:] )) );
            fileHandle.close();


        mutextosort.release()
        
        if(not(foundjob)):
            time.sleep(SLEEPTIME);
        else:
            #print "cmd sort";

            #Begin sorting
            stringNtosort = bamfiletosort;
            listInfosorter=bamfiletosort.split("\t");
            cmd = "sam sort -m 2G ";
            cmd+=" -O "+str(listInfosorter[3])[:-4]+" ";

            cmd+=" "+str(listInfosorter[3])+"_temp "+tempDir;
            print "sorter "+str(cmd);
            handle_job(cmd);
            stringNtosort = "";
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

            #end sorting
            #write to done, send email

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

            mutexdone.acquire()
            print "Finished: "+str(listInfosorter[3]);
            fileHandle = open( fileNdone , 'a' );
            fileHandle.write( bamfiletosort+"\t"+timeString()+"\n" );
            fileHandle.close();
            
            
            mutexdone.release()

            #TODO email if not "none"


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
    # -o /dev/null  -e /dev/null
    #temp = tempfile.NamedTemporaryFile(prefix=options.tmp+"script",suffix=".sge",delete=False)
    cmd= " echo \"/home/public/usr/bin/bwa worker -T 20000 -t \$NSLOTS -p "+str(PORT_NUMBERBWA)+" -h "+str( (socket.gethostname()) )+"\" | "+str(qsubcmd)+"   -N alib -S /bin/bash -l \"class=*,h_vmem=6.8G,s_vmem=6.8G,virtual_free=6.8G \" -V  -pe smp 1- ";
    #print cmd;
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
    print "launcher thread running";

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
parser.add_option("-t", "--temp",  dest="temp", help="Temporary directory",default="/tmp/");
parser.add_option("--qstat",  dest="qstat", help="qstat location [default: %default]",default="/opt/sge/bin/lx-amd64/qstat");
parser.add_option("--qsub",  dest="qsub", help="qsub location [default: %default]",default="/opt/sge/bin/lx-amd64/qsub");
parser.add_option("--qdel",  dest="qdel", help="qdel location [default: %default]",default="/opt/sge/bin/lx-amd64/qdel");

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)


(options, args) = parser.parse_args();

tempDir = options.temp;

qsubcmd  = options.qsub;
qstatcmd = options.qstat;
qdelcmd  = options.qdel;

fileNtocount = tempDir+"/alib.tocount"
fileNtomap   = tempDir+"/alib.tomap"
fileNtosort  = tempDir+"/alib.tosort"
fileNdone    = tempDir+"/alib.done"

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


    tcounter = threading.Thread(target = counter );
    tcounter.daemon = True
    tcounter.start()    

    #if False:
    

    tmap = threading.Thread(target = mapper );
    tmap.daemon = True
    tmap.start()    

    tlauncher = threading.Thread(target = launcher );
    tlauncher.daemon = True
    tlauncher.start()    

    tsorter = threading.Thread(target = sorter );
    tsorter.daemon = True
    tsorter.start()    
    
    #tlauncher = threading.Thread(target = launcher );
    #tlauncher.daemon = True
    #tlauncher.start()    

    #tc2 = threading.Thread(target = print_time,args=("Thread-2", 1, ) );
    #tc2.daemon = True
    #tc2.start()    

    while True: 
        print "program alive at "+timeString();
        time.sleep(3600)
except (KeyboardInterrupt, SystemExit):
    print '\n! Received keyboard interrupt, quitting threads.\n'
    deletemyjobs();
    server.socket.close();
    sockzmq.close();
    context.term();

    print "shutdown server";

