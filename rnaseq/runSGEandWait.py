#!/usr/bin/python

import sys,os
from optparse import OptionParser
from optparse import OptionGroup
import subprocess
import random;

#qsub ="/opt/sge/bin/lx-amd64/qsub"
#qstat="/opt/sge/bin/lx-amd64/qstat"
import time;
#param = 
import tempfile;

def runQsub(cmd):    
#  print cmd;
  p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  (stdout, stderr) = p.communicate()
  pstat = p.wait();
 # print "done";
  if(pstat != 0):
    print "Error cmd returned a non zero code :"+cmd;
    sys.exit(1);

  fields= stdout.split();
  if(len(fields) < 3):
    print "Error cmd did not have the expected format : "+cmd;
    sys.exit(1);
  return fields[2];
      #while(True):
      #retcode = p.poll() #returns None while subprocess is running

      #stdout = stdout+p.stdout.readline();

def runQstat(cmd,code):    
  #print cmd;
  cmd = cmd +" -j "+str(code);
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
  for line in  stdout.split('\n'):
    fields= line.split();
    if(len(fields) < 3):
      continue;
    #if(fields[0][0] == "-"):
    #  continue;    
    if(fields[0] == "job-ID"):
      continue;

    if(fields[0] == code):
      runcode = fields[4];
      found=True;
      #print "line "+line;
      #if(not found ):
  return runcode;
  
  #if(len(fields) < 3):
  #  print "Error cmd did not have the expected format : "+cmd;
  #  sys.exit(1);
  #return fields[2];

def runQacct(cmd):    
  #print cmd;
  p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  (stdout, stderr) = p.communicate()
  pstat = p.wait();

  for line in  stdout.split('\n'):
    #print "out "+line;

    fields= line.split();
    if(len(fields) < 2):
      continue;
    if(fields[0] == "exit_status"):
      return fields[1];
  return "NA";

parser = OptionParser(usage="usage: %prog [options]");

parser.add_option("-i", "--infile",  dest="infile", help="File with commands or specify:");
parser.add_option("-c", "--cmd",  dest="command", help="Command to run");
parser.add_option("--tmp",  dest="tmp", help="Temp directory [default: %default]",default="/mnt/scratch/tmp/");
parser.add_option("--time",  dest="checktime", help="Delay for checking on the job in seconds [default: %default]" ,type="int",default=600);

#parser.add_option("-i", "--infile",  dest="qsub", help="qsub location",default="/opt/sge/bin/lx-amd64/qsub");
parser.add_option("--qsub",  dest="qsub", help="qsub location [default: %default]",default="/opt/sge/bin/lx-amd64/qsub");
parser.add_option("--qacct",  dest="qacct", help="qacct location [default: %default]",default="/opt/sge/bin/lx-amd64/qacct");
parser.add_option("--qstat",  dest="qstat", help="qstat location [default: %default]",default="/opt/sge/bin/lx-amd64/qstat");
parser.add_option("--param",  dest="param", help="Parameters for qsub [default: \"%default\"]",default=" -S /bin/bash -l \"class=*,h_vmem=2500M,virtual_free=2500M\" -V -R y -pe smp 1 ");

#group = OptionGroup(parser, "Group1","First Group of parameters");
#group.add_option("-o", "--outfile", dest="outfile", help="output");

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)

#parser.add_option_group(group);
(options, args) = parser.parse_args();

#check for either infile or cmd
#if cmd make temp file
if((options.infile)  and
   (options.command)  ):
  print "Specify either -i or -c but not both ";
  sys.exit(1);


temp = tempfile.NamedTemporaryFile(prefix=options.tmp+"script",suffix=".sge",delete=False)
#print 'temp.name:', temp.name
if( (options.command) ):
  
  #temp = tempfile.NamedTemporaryFile(prefix=options.tmp,suffix=".sge")

  #print 'temp:', temp

  temp.write("#!/bin/bash\n");
  temp.write("\n");
  temp.write("cd "+os.getcwd()+"\n");
  temp.write(""+options.command+"\n");
  #sys.exit(1);
  options.infile=temp.name;
  temp.close();
  #print 'temp.name:', temp.name

#print 'temp.name2:', temp.name

if not os.path.exists(options.qsub):
  print "Required executable file not found "+options.qsub;
  sys.exit(1);

if not os.path.exists(options.qstat):
  print "Required executable file not found "+options.qstat;
  sys.exit(1);


if not os.path.exists(options.infile):
  print "Required input file not found "+options.infile;
  sys.exit(1);
#else:
  #print "exist "+options.infile;

if os.path.exists(options.infile+".o"):
  sys.stderr.write("removing "+options.infile+".o");
  os.remove(options.infile+".o");

if os.path.exists(options.infile+".e"):
  sys.stderr.write("removing "+options.infile+".e");
  os.remove(options.infile+".e");

cmd = options.qsub+" "+options.param+" -e "+options.infile+".e"+" -o "+options.infile+".o "+options.infile;
#print cmd;
code = runQsub(cmd);
#removing

while True:
  time.sleep(  max(int(random.normalvariate(options.checktime, 60)),1) );

  codeqstat = runQstat(options.qstat,code)
  #print codeqstat;
  if(codeqstat == "done"):
    break;

#time.sleep( 15 );
#print code;
checks=0;
while True:
  acctcode = runQacct(options.qacct+" -j "+code);
  checks+=1;
  if(checks>1000):
    print "cannot validate code for job "+code;
    acctcode="1";
    break;
  if(acctcode == "NA"):
    time.sleep(  max(int(random.normalvariate( options.checktime , 2)),1) );
  else:
    break;

if(acctcode != "0" ):
  print "Error cmd returned a non zero code ("+str(acctcode)+") in the qacct, code :"+code;
  sys.exit(1);
  
if( (options.command)  ):

  if os.path.exists(options.infile):
    #sys.stderr.write("removing "+options.infile+".o");
    os.remove(options.infile);

  if os.path.exists(options.infile+".o"):
    #sys.stderr.write("removing "+options.infile+".o");
    os.remove(options.infile+".o");

  if os.path.exists(options.infile+".e"):
    #sys.stderr.write("removing "+options.infile+".e");
    os.remove(options.infile+".e");
