#!/usr/bin/python

import sys,os
import glob
import subprocess

from optparse import OptionParser
from optparse import OptionGroup
from struct import unpack
                                                                                                                                     

def clean_path(cpath):
  npath = cpath.replace('/./','/').replace('//','/')
  while (cpath != npath):
    cpath = npath
    npath = cpath.replace('/./','/').replace('//','/')
  return npath
if (len(sys.argv[0]) > 0)  and os.path.isdir(os.getcwd()+'/'+os.path.dirname(sys.argv[0])):
  progPath = os.path.dirname(os.getcwd()+'/'+sys.argv[0])+'/'                                                                                                    
elif (len(sys.argv[0]) > 0)  and os.path.isdir(os.path.dirname(sys.argv[0])): 
  progPath = os.path.dirname(sys.argv[0])+'/'                                                                                                 
else:
  progPath = os.getcwd()+'/'                                                                                                         
progPath=clean_path(progPath);
clocsProg=progPath+"/count_clocs";
countBAM =progPath+"/countClustersBAM";   

if not os.path.isfile(clocsProg):
  print "File "+clocsProg+" not found";   
  sys.exit()      

if not os.path.isfile(countBAM):
  print "File "+countBAM+" not found";   
  sys.exit()      

parser = OptionParser(usage="This program counts the number of expected clusters in the position files\nand checks if the number of sequences is equal in the resulting BAM files.\nThis script is used for quality control.\n usage: %prog [options]");
parser.add_option("-o", "--outdir",  dest="outdir", help="Where to store the report file");
parser.add_option("-b", "--bamfile", dest="bamfile", help="Comma-separated list of BAM files to analyze");
parser.add_option("-r", "--run", dest="runfolder", help="Runfolder containing the positions files in pos.txt, locs, clocs");
parser.add_option("-l", "--lane", dest="lane", help="Lane to use");

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)

#parser.add_option_group(group);
(options, args) = parser.parse_args();

if(options.outdir == None):
  print "Output directory must be defined";
  sys.exit(1);

if(options.bamfile == None):
  print "BAM files must be defined";
  sys.exit(1);

if(options.runfolder == None):
  print "Runfolder must be defined";
  sys.exit(1);

if(options.lane == None):
  print " must be defined";
  sys.exit(1);
  


totalClusters=0;
foundPos=0;
foundLOCS=0;

#pos
for name in glob.glob(str(options.runfolder)+"/L00"+str(options.lane)+"/s_"+str(options.lane)+"_*_pos.txt"):
  print name;
  foundPos=1;
  cmd="wc -l "+str(name);
  jobcreated=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE);
  out = str(jobcreated.communicate()[0]);  

  if(jobcreated.returncode != 0):
    print "Wrong return code for cmd = "+str(cmd);
    sys.exit(1);

  totalClusters+=int(out.split()[0]);


#locs, mostly 
for name in glob.glob(str(options.runfolder)+"/L00"+str(options.lane)+"/s_"+str(options.lane)+"_*.locs"):
  if(foundPos != 0):
    print "Cannot have two types of position files";
    sys.exit(1);
  foundLOCS=1
  f = open(name, "rb");
  nclus=(f.read(8));
  nclus=unpack('I',f.read(4))[0];
  totalClusters+=nclus;


#clocs mostly Hiseq
for name in glob.glob(str(options.runfolder)+"/L00"+str(options.lane)+"/s_"+str(options.lane)+"_*.clocs"):
  if(foundPos != 0 or foundLOCS != 0):
    print "Cannot have two types of position files";
    sys.exit(1);

  cmd=clocsProg+" "+str(name);
  jobcreated=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE);
  out = str(jobcreated.communicate()[0]);  

  if(jobcreated.returncode != 0):
    print "Wrong return code for cmd = "+str(cmd);
    sys.exit(1);

  totalClusters+=int(out.split()[0]);

if(totalClusters == 0):
  print "Could not compute the number of clusters";
  sys.exit(1);

#put bam


tilesAreFine=1;
stringToPrint=str(totalClusters)+" initial positions\n";

for name in options.bamfile.split(","):
#  print name;
#countBAM
  cmd=countBAM+" "+str(name);
  jobcreated=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE);
  out = str(jobcreated.communicate()[0]);  

  if(jobcreated.returncode != 0):
    print "Wrong return code for cmd = "+str(cmd);
    sys.exit(1);

  foundClusters=int(out.split()[0]);

  if(foundClusters != totalClusters):
    tilesAreFine=0;
    stringToPrint+=str(foundClusters)+" ERROR BAM:"+str(name)+"\n";
  else:
    stringToPrint+=str(foundClusters)+" OK BAM:"+str(name)+"\n";

if(tilesAreFine ):
  fileHandleWrite = open ( str(options.outdir)+"clusterTally_"+str(options.lane)+".OK", 'w' ) ;
else:
  fileHandleWrite = open ( str(options.outdir)+"clusterTally_"+str(options.lane)+".ERROR", 'w' ) ;


fileHandleWrite.write(stringToPrint);

fileHandleWrite.close();

