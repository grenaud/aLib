#!/usr/bin/python

import sys,os
import json
import time
import datetime

from optparse import OptionParser
from optparse import OptionGroup

splitByRG            = "pipeline/splitByRG";
bam2fastq            = "BCL2BAM2FASTQ/bam2fastq";

def chomp(s):
  return s.rstrip('\n');


parser = OptionParser(usage="usage: %prog [options]");
#parser.add_option("-i", "--infile",  dest="infile", help="input");
input = OptionGroup(parser, "Input1","Input files");
input.add_option("-i", "--infile", dest="infile", help="Input BAM file");
input.add_option("-g", "--bowtie", dest="bowtieIndex", help="Bowtie index");


parser.add_option_group(input);
(options, args) = parser.parse_args();

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)



pathToConfig=sys.argv[0];
pathToConfig=os.path.dirname(os.path.abspath(pathToConfig))+"/../webForm/config.json";

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

alibdir          = jsondataConf["alibdir"];

splitByRG = alibdir+"/"+splitByRG;
if not os.path.exists(splitByRG):
  print "Required executable file not found "+splitByRG;
  sys.exit(1);

bam2fastq = alibdir+"/"+bam2fastq;
if not os.path.exists(bam2fastq):
  print "Required executable file not found "+bam2fastq;
  sys.exit(1);


#check for bowtie index
