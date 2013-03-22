#!/usr/bin/python


# Date: Mar-20-2013 
# Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com


import sys,os
import json
#import xml.parsers.expat
import xml.etree.ElementTree as ET
from optparse import OptionParser
from optparse import OptionGroup

# directory variables

alibdir="";
illuminareaddir="";
illuminawritedir=""


BCL2BAM              = "BCL2BAM/bcl2bam"

MergeReads           = "pipeline/mergeTrimReadsBAM"
IndexReassign        = "pipeline/assignRG"
IndexReassignRGR     = "pipeline/rg.R"
IndexReassignRATIOR  = "pipeline/ratio.R"

BAMFilter            = "pipeline/filterReads"
BAMFilterLIKER       = "pipeline/likeli.R"
ERRORPERCYCLE        = "pipeline/errorRatePerCycle"
ERRORPERCYCLER       = "pipeline/errorRatePerCycle.R"

Tilecounter          = "tileCount/tileCount.py";
Qualplotter          = "plotQualScores/plotQualScores.py";

Ctrlextract          = "extractControlReadsBam/getCtrlReadsBAM";
Predvsobs            = "qualScoreC++/qualScoresObsVsPred";
PredvsobsR           = "qualScoreC++/generateplot.R";





def chomp(s):
    return s.rstrip('\n');


def readConfig(pathToConfig):
  global alibdir;
  global illuminawritedir;
  global illuminareaddir;

  pathToConfig=os.path.dirname(os.path.abspath(pathToConfig))+"/config.xml";
  try:
    tree =ET.parse(pathToConfig);
  except IOError:
    print "Cannot open XML config file "+pathToConfig;
    sys.exit(1);

  root = tree.getroot()
  alibdir          = root.find("alibdir").text;
  illuminawritedir = root.find("illuminawritedir").text;
  illuminareaddir  = root.find("illuminareaddir").text;


def checkPrograms():
  global alibdir;
  

parser = OptionParser(usage="usage: %prog [options] <json doc>");
#parser.add_option("-i", "--infile",  dest="infile", help="input");
#group = OptionGroup(parser, "Group1","First Group of parameters");
#group.add_option("-o", "--outfile", dest="outfile", help="output");

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)

#parser.add_option_group(group);
(options, args) = parser.parse_args();


readConfig(sys.argv[0]);



#check if execs are there

BCL2BAM = alibdir+"/"+BCL2BAM;
if not os.path.exists(BCL2BAM):
  print "Required executable file not found "+BCL2BAM;
  sys.exit(1);


MergeReads = alibdir+"/"+MergeReads;
if not os.path.exists(MergeReads):
  print "Required executable file not found "+MergeReads;
  sys.exit(1);

IndexReassign = alibdir+"/"+IndexReassign;
if not os.path.exists(IndexReassign):
  print "Required executable file not found "+IndexReassign;
  sys.exit(1);

IndexReassignRGR = alibdir+"/"+IndexReassignRGR;
if not os.path.exists(IndexReassignRGR):
  print "Required executable file not found "+IndexReassignRGR;
  sys.exit(1);

IndexReassignRATIOR = alibdir+"/"+IndexReassignRATIOR;
if not os.path.exists(IndexReassignRATIOR):
  print "Required executable file not found "+IndexReassignRATIOR;
  sys.exit(1);

BAMFilter = alibdir+"/"+BAMFilter;
if not os.path.exists(BAMFilter):
  print "Required executable file not found "+BAMFilter;
  sys.exit(1);

BAMFilterLIKER = alibdir+"/"+BAMFilterLIKER;
if not os.path.exists(BAMFilterLIKER):
  print "Required executable file not found "+BAMFilterLIKER;
  sys.exit(1);

ERRORPERCYCLE = ALIBDIR+"/"+ERRORPERCYCLE;
if not os.path.exists(ERRORPERCYCLE):
  print "Required executable file not found "+ERRORPERCYCLE;
  sys.exit(1);

ERRORPERCYCLER = ALIBDIR+"/"+ERRORPERCYCLER;
if not os.path.exists(ERRORPERCYCLER):
  print "Required executable file not found "+ERRORPERCYCLER;
  sys.exit(1);

Tilecounter = alibdir+"/"+Tilecounter;
if not os.path.exists(Tilecounter):
  print "Required executable file not found "+Tilecounter;
  sys.exit(1);

Qualplotter = alibdir+"/"+Qualplotter;
if not os.path.exists(Qualplotter):
  print "Required executable file not found "+Qualplotter;
  sys.exit(1);

Ctrlextract = alibdir+"/"+Ctrlextract;
if not os.path.exists(Ctrlextract):
  print "Required executable file not found "+Ctrlextract;
  sys.exit(1);

Predvsobs = alibdir+"/"+Predvsobs;
if not os.path.exists(Predvsobs):
  print "Required executable file not found "+Predvsobs;
  sys.exit(1);

PredvsobsR = alibdir+"/"+PredvsobsR;
if not os.path.exists(PredvsobsR):
  print "Required executable file not found "+PredvsobsR;
  sys.exit(1);












print illuminareaddir;

try:
  fileHandle = open ( sys.argv[-1] );
except IOError:
    print "Cannot open json file "+sys.argv[-1];
    sys.exit(1);

jsonstring="";
while 1:
    line = fileHandle.readline();
    if(not(line)):
        break
    line = chomp(line);
    jsonstring+=line;
fileHandle.close();

#print jsonstring;
jsondata=json.loads(jsonstring);

print jsondata["genomebwa"];






#
#  try:
#    fileHandleXML = open ( pathToConfig );
#  except IOError:
#    print "Cannot open XML config file "+pathToConfig;
#    sys.exit(1);
#
#  xmlstring="";
#  while 1:
#    line = fileHandleXML.readline();
#    if(not(line)):
#      break
#    #    print line;
#    xmlstring+=line;
#    
#  fileHandleXML.close();
#  print xmlstring;
#  #p = xml.parsers.expat.ParserCreate()
#
