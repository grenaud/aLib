#!/usr/bin/python




import sys,os
import json
import time
import datetime
import subprocess

#import xml.parsers.expat
#import xml.etree.ElementTree as ET
import pprint;
from optparse import OptionParser
from optparse import OptionGroup
from random import randint
from os import listdir
from os.path import isfile, isdir, join

illuminareaddir="";
illuminawritedir=""

parser = OptionParser(usage="usage: %prog [options]");
parser.add_option("--mock", dest="mock", help="Do a mock run for testing",default=False,action="store_true")
parser.add_option("--run", dest="torunnow", help="Run these makefiles now (comma delimited)",default="")


if(len(sys.argv) == 0):
  parser.print_help()
  sys.exit(1)


(options, args) = parser.parse_args();

arrayOfJobs=[]; #array of 3-upple  (cmd string,Popen object,runid) of jobs currently running


def chomp(s):
  return s.rstrip('\n');

def getMakefilesstack():
  toreturn = [];
  for fileIll in listdir(illuminawritedir):
    if(not isfile(fileIll)):
      #print illuminawritedir+"/"+fileIll;
      for lane in range(1,8):
        if( isdir(illuminawritedir+"/"+fileIll+"/build/lane"+str(lane)) ):
          #print illuminawritedir+"/"+fileIll+"/build/lane"+str(lane);
          for proc in range(1,100):
            if( isdir(illuminawritedir+"/"+fileIll+"/build/lane"+str(lane)+"/proc"+str(proc)) ):
              #print illuminawritedir+"/"+fileIll+"/build/lane"+str(lane)+"/proc"+str(proc);
              if(isfile(illuminawritedir+"/"+fileIll+"/build/lane"+str(lane)+"/proc"+str(proc)+"/Makefile")):
                toreturn.append(illuminawritedir+"/"+fileIll+"/build/lane"+str(lane)+"/proc"+str(proc)+"/Makefile");
  return toreturn;



  
def handle_jobs(cjob):
  global options
  #  global jobs
  #print "launching "+str(cjob);
  if options.mock:
    return None

  jobcreated=subprocess.Popen(cjob,shell=True);

  return jobcreated;

def checkFinishedJobs():
  global arrayOfJobs;      
  #checking for finished jobs
  arrayOfJobsCopy=[];#arrayOfJobs[]#[:];
  for toverify in arrayOfJobs:
    if(toverify[1].poll() != None): #has finished      

      if(toverify[1].returncode != 0): #has finished but wrong code
        print "WARNING: process "+str(toverify[0])+" failed, relaunch it manually";

    else: #is finished
      arrayOfJobsCopy.append(toverify);

  arrayOfJobs=arrayOfJobsCopy;

def handleListOfjobs(alljobs):
  global arrayOfJobs;  

  while(len(alljobs) != 0):
    print "currently "+str(len(arrayOfJobs))+ " jobs are running\n";

    checkFinishedJobs();
    unresolved=alljobs[:]; #copy

    for jobToAdd in unresolved:
      if(len(arrayOfJobs) < 2):
        print "launching "+str(jobToAdd[0]);
        arrayOfJobs.append([jobToAdd[0],handle_jobs(jobToAdd[0]),jobToAdd[1]]); #the cmd, the object, runid
        alljobs.remove(jobToAdd);

    if(len(alljobs) != 0):
      runidtoprint=[];
      for toverify in arrayOfJobs:
        runidtoprint.append(toverify[2]);
      
      print "Could not launch all jobs, sleeping for 15 mins, currently the following runs are running: "+(','.join(runidtoprint))+"\n";
      time.sleep(900);


pathToConfig=sys.argv[0];
pathToConfig=os.path.dirname(os.path.abspath(pathToConfig))+"/config.json";


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



illuminawritedir = jsondataConf["illuminawritedir"];
illuminareaddir  = jsondataConf["illuminareaddir"];

illuminawritedir=illuminawritedir.rstrip('/');
#print illuminawritedir;

#startup 
initialList = getMakefilesstack();
#initialList=initialList[:-1];

#print "Initial makefiles found:";
#for makefilesfound in initialList:
#  print "Found: "+str(makefilesfound);

#for runtocheck in (options.torunnow).split(","):
#  print "remove "+str(runtocheck);

initialListCopy = [];

if(len(options.torunnow)>1):
  for runtocheck in ( (options.torunnow).split(",") ):
    if not(runtocheck in initialList):
      print "The run provided "+runtocheck+" was not found among the initial runs";
      sys.exit(1);


for runtocheck in initialList:
  if not(runtocheck in (options.torunnow).split(",") ):
    initialListCopy.append(runtocheck);

initialList = initialListCopy;
print "Initial makefiles found:";
for makefilesfound in initialList:
  print "Found: "+str(makefilesfound);

#sys.exit(1);

while True:

  newList = getMakefilesstack();
  
  runsToProcess={};
  listNextIteration = [];

  #find unique runs
  tolaunch=[]; 
  runid2makefiles={};

  for makef in newList:
    boolFound = False;

    for oldmakef in initialList:
      if(oldmakef == makef):
        boolFound = True;
    
    arrayfields=makef.split("/")[:-3];
    runid='/'.join(arrayfields[:-1]);
    runid=runid[len(illuminawritedir):]
    #print runid;
    runidpath='/'.join(arrayfields);
    
    if(not boolFound):
      print "New makefile found "+makef;
      print "Checking for .completed file: "+str(illuminareaddir)+"/"+str(runid)+"/Run.completed";
      if(isfile(illuminareaddir+"/"+runid+"/Run.completed")):
        runalreadyrunning=False;
        for toverify in arrayOfJobs:
          if toverify[2] == runid:
            runalreadyrunning=True;
            
        if not(runalreadyrunning):
          runsToProcess[runid]=True;
          print "Ready to process";
          listNextIteration.append(makef); #append to next iteration because it will be processed

          if(runid in runid2makefiles):
            print "This makefile belongs to a previously discovered run : "+str(runid);
            runid2makefiles[runid].append(makef);
          else:
            print "This makefile is for a newly discovered run : "+str(runid);
            runid2makefiles[runid]= [makef];
        else:
          print "Another instance of run "+str(runid)+" is already running";

      else:
        print "not ready";
    else:
      listNextIteration.append(makef);

  for runidtoadd in runsToProcess.keys():
    cmd = "cd "+illuminawritedir+"/"+runidtoadd+"/build/ && make -j 4 > "+illuminawritedir+"/"+runidtoadd+"/build/Makefile.stdout 2> "+illuminawritedir+"/"+runidtoadd+"/build/Makefile.stderr ";#&& make sendemail"
    for makefileassociated in runid2makefiles[runidtoadd]:
      cmd += " && make -f "+str(makefileassociated)+" sendemail ";
    tolaunch.append( [cmd,runidtoadd] );

  handleListOfjobs(tolaunch);
  initialList = listNextIteration; #list of current make files

  timeNowRaw = time.time()
  timeNow    = datetime.datetime.fromtimestamp(timeNowRaw).strftime('%Y-%m-%d_%H:%M:%S');
  runidtoprint=[];
  for toverify in arrayOfJobs:
    runidtoprint.append(toverify[2]);

  print "Sleeping for 1 hour "+str(timeNow)+" currently the following runs are running: "+(','.join(runidtoprint))+"\n";
  time.sleep(3600);
  checkFinishedJobs();
