#!/usr/bin/python


# Date: Mar-20-2013 
# Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
#TODO
# dedicated lane for control
# error plot freeIbis
# allow one mismatch in merger program
# number of bases to trim in single ended reads


import sys,os
import json
import time
import datetime
#import xml.parsers.expat
#import xml.etree.ElementTree as ET
import pprint;
from optparse import OptionParser
from optparse import OptionGroup
from random import randint

# directory variables
def makedirs(path):
  global options
  if os.path.isdir(path): return False
  if not options.mock:
    try: os.makedirs(path)
    except: pass
    if not os.path.isdir(path):
      print "Error: Could not create",path
      return False
    else:
      return True
  else:
    print "os.makedirs(%s)"%path
    return True

#yes, copy pasted
def parse_rangestr(rangestr):
  res = []
  fields = rangestr.split(',')
  for elem in fields:
    if "-" in elem:
      se = elem.split('-')
      if len(se) == 2:
        try:
          start = int(se[0])
          end = int(se[1])
          res.extend(range(start,end+1))
        except: return None
      else: return None
    else:
      try:
        pos = int(elem)
        res.append(pos)
      except: return None
  if len(res) == 0: return None
  else:
    res = list(set(res))
    res.sort()
    return map(str,res)

alibdir="";
illuminareaddir="";
illuminawritedir=""

timeNowRaw = time.time()
timeNow    = datetime.datetime.fromtimestamp(timeNowRaw).strftime('%Y-%m-%d_%H:%M:%S');

BCL2BAM              = "BCL2BAM2FASTQ/BCL2BAM/bcl2bam"
FREEIBIS             = ""

rtaReportHask        = "pipeline/generate_report"
#rtaReportPython      = "pipeline/generate_report.py"

MergeReads           = "leehom/src/leeHom"
IndexReassign        = "deML/src/deML"
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
FastQCreport         = "fastqc"

insertsizebi         = "insertSize/insertSizeRG"
insertsizeR          = "insertSize/insertSizeRG.R"

sendemail            = "webForm/sendemail.php"



def chomp(s):
  return s.rstrip('\n');





#BEGIN READ CONFIG FILE
pathToConfig=sys.argv[0];
pathToConfig=os.path.dirname(os.path.abspath(pathToConfig))+"/config.json";
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

#try:
#  tree =ET.parse(pathToConfig);
#except IOError:
#  print "Cannot open JSON config file "+pathToConfig;
#  sys.exit(1);



#XMLconfig        = tree.getroot()


illuminawritedir = jsondataConf["illuminawritedir"];
illuminareaddir  = jsondataConf["illuminareaddir"];
tempdir          = jsondataConf["tempdirectory"];
alibdir          = jsondataConf["alibdir"];
#alibdir          = "/home/gabriel_renaud/projects/aLib/aLib/";

FREEIBIS         = jsondataConf["freeibispath"];
BWAGENOMES       = jsondataConf["genomedirectory"];

FastQCreport     = jsondataConf["fastqcdir"]+"/"+FastQCreport;


def checkPrograms():
  global alibdir;
  

parser = OptionParser(usage="usage: %prog [options] <json doc>");
parser.add_option("-o", "--outdir", dest="outdir", help="Out directory to put the makefile, will create directories like [outdir]/../Report/ etc");
parser.add_option("--trainlanes",dest="trainlanes", help="Consider only a subset of lanes for training freeIbis (default use data from json)",default="1-8")
parser.add_option("--tempdir",dest="tempdir", help="Use this directory as temporary directory",default=tempdir)

parser.add_option("--lanes",     dest="lanes",help="Consider only a subset of lanes for the Makefile (default use all from json)")

parser.add_option("-c", "--cores", dest="cores", help="Maximum number of CPU cores to be used for freeIbis(default 8)",default=8,type="int")
parser.add_option("--mock", dest="mock", help="Just print command",action="store_true",default=False)

#group = OptionGroup(parser, "Output","Output options");
#group.add_option("-o", "--outdir", dest="outdir", help="Out directory to put the makefile");

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)

#parser.add_option_group(group);
(options, args) = parser.parse_args();



if options.outdir == None or not os.path.isdir(options.outdir):
  print "Error: Need valid path for output directory";
  sys.exit()
  

options.outdir=os.path.abspath(options.outdir);
outBaseDirectory=os.path.abspath(options.outdir+"/../");


#readConfig(sys.argv[0]);



#check if execs are there

#/mnt/solexa/Genomes/

CTRLGENOME = BWAGENOMES +"/"+"phiX/control/whole_genome.fa";
if not os.path.exists(CTRLGENOME):
  print "Control whole genome file not found "+CTRLGENOME;
  sys.exit(1);

CTRLGENOMEBWA = BWAGENOMES +"/"+"phiX/control/bwa-0.4.9.amb";
if not os.path.exists(CTRLGENOMEBWA):
  print "Control BWA index file not found "+CTRLGENOMEBWA;
  sys.exit(1);

CTRLGENOMEBWA = BWAGENOMES +"/"+"phiX/control/bwa-0.4.9";





BCL2BAM = alibdir+"/"+BCL2BAM;
if not os.path.exists(BCL2BAM):
  print "Required executable file not found "+BCL2BAM;
  sys.exit(1);


MergeReads = alibdir+"/"+MergeReads;
if not os.path.exists(MergeReads):
  print "Required executable file not found "+MergeReads;
  sys.exit(1);

rtaReportHask = alibdir+"/"+rtaReportHask;
if not os.path.exists(rtaReportHask):
  print "Required executable file not found ->"+rtaReportHask+"<-";
  sys.exit(1);

#rtaReportPython = alibdir+"/"+rtaReportPython;
#if not os.path.exists(rtaReportPython):
#  print "Required executable file not found "+rtaReportPython;
#  sys.exit(1);

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

ERRORPERCYCLE = alibdir+"/"+ERRORPERCYCLE;
if not os.path.exists(ERRORPERCYCLE):
  print "Required executable file not found "+ERRORPERCYCLE;
  sys.exit(1);

ERRORPERCYCLER = alibdir+"/"+ERRORPERCYCLER;
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

insertsizebi = alibdir+insertsizebi;
if not os.path.exists(insertsizebi):
  print "Required executable file not found "+insertsizebi;
  sys.exit(1);

insertsizeR = alibdir+insertsizeR;
if not os.path.exists(insertsizeR):
  print "Required executable file not found "+insertsizeR;
  sys.exit(1);


FastQCreport = FastQCreport;
if not os.path.exists(FastQCreport):
  print "Required executable file not found "+FastQCreport;
  sys.exit(1);


if not os.path.exists(FREEIBIS+"/runBaseCalling.py"):
  print "Required executable file not found "+FREEIBIS+"/runBaseCalling.py";
  sys.exit(1);


#print illuminareaddir;

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

#pprint.pprint(jsondata);
lanesToUse=[];

if(options.lanes):
  lanesToUse = parse_rangestr(options.lanes);
else:
  lanesToUse = jsondata["lanes"]

#print "LANE "+str(lanesToUse);
#sys.exit(1);

numprocessingcurrent=jsondata["numprocessingcurrent"];


for lanetopredict in lanesToUse:
  #allsubdir.append("QC/insertsize/lane"+str(lanetopredict) );
  targetdirbuild=options.outdir+"/build/lane"+str(lanetopredict)+"/proc"+str(numprocessingcurrent)+"/";    
  makedirs(targetdirbuild);


if( jsondata["usebwa"+str(numprocessingcurrent)] ):
  TARGETGENOME = BWAGENOMES +"/"+jsondata["genomebwa"+str(numprocessingcurrent)]+"/whole_genome.fa";
  if not os.path.exists(TARGETGENOME):
    print "Target whole genome file not found "+TARGETGENOME;
    sys.exit(1);

  TARGETGENOMEBWA = BWAGENOMES +"/"+jsondata["genomebwa"+str(numprocessingcurrent)]+"/bwa-0.4.9.amb";
  if not os.path.exists(TARGETGENOMEBWA):
    print "Target BWA index file not found "+TARGETGENOMEBWA;
    sys.exit(1);





#BEGIN WRITING TARGET DIRECTORIES    
targetdir=outBaseDirectory+"/Report/";    
#makeWrite.write(targetdir+"\n\tmkdir "+targetdir);
makedirs(outBaseDirectory+"/Report/");


allsubdir=["FastQC/",
           "QC/",
           "QC/qscores/",
           "QC/rg/",
           "QC/filter/",
           "QC/insertsize/",
           "Raw_Sequences/",
           "Final_Sequences/",
           "BWA",
           "FastQ"];

for lanetopredict in lanesToUse:
  allsubdir.append("QC/insertsize/lane"+str(lanetopredict) );

allsubdir.append("Final_Sequences/proc"+str(numprocessingcurrent)+"/");
allsubdir.append("BWA/proc"+str(numprocessingcurrent)+"/");
allsubdir.append("QC/insertsize/proc"+str(numprocessingcurrent)+"/");
allsubdir.append("QC/filter/proc"+str(numprocessingcurrent)+"/");


if(jsondata["bustard"]):  
  for subdir in allsubdir:
    targetdir=outBaseDirectory+"/Bustard/"+subdir;    
    makedirs(targetdir);

if(jsondata["freeibis"]):  
  for subdir in allsubdir:
    targetdir=outBaseDirectory+"/Ibis/"+subdir;    
    makedirs(targetdir);

#bam files for tile count
listOfBAMfilesToCheck={};
#flagstatx files for tile count
listOfFLGfilesToCheck={};


#all files from basecall
listOfFilesBasecall={};
#all files from final_sequence
listOfFilesFinal={};
#all files from BWA
listOfFilesBWA={};
#all files from QC
listOfFilesQC={};
#all files from QC
listOfTargetFiles={};
#all files from QC
listOfFilesWithIndices={};



if(jsondata["bustard"]):
  for lanetopredict in lanesToUse:
    listOfBAMfilesToCheck[lanetopredict] = [];
    listOfFLGfilesToCheck[lanetopredict] = [];
    listOfFilesBasecall[lanetopredict]   = [];
    listOfFilesFinal[lanetopredict]      = [];
    listOfFilesBWA[lanetopredict]        = [];
    listOfFilesQC[lanetopredict]         = [];
    listOfTargetFiles[lanetopredict]     = [];

if(jsondata["freeibis"]):
  for lanetopredict in lanesToUse:
    listOfBAMfilesToCheck[lanetopredict] = [];
    listOfFLGfilesToCheck[lanetopredict] = [];
    listOfFilesBasecall[lanetopredict]   = [];
    listOfFilesFinal[lanetopredict]      = [];
    listOfFilesBWA[lanetopredict]        = [];
    listOfFilesQC[lanetopredict]         = [];
    listOfTargetFiles[lanetopredict]     = [];


#END WRITING TARGET DIRECTORIES    


#################################################
#                                               #
#              BEGIN MAKEFILE                   #
#                                               #
#################################################


makeWrite=[]
for l in range(int(jsondata["LaneCount"])+1):
  #print "l"+str(l)+"\n";
  makeWrite.append(None);



#for l in range(1,int(jsondata["LaneCount"])+1):
#  print "l2 "+str(l)+"\n";
#  makeWrite.append(None);

for lanetopredict in lanesToUse:
  makefilePath=options.outdir+"/lane"+str(lanetopredict)+"/proc"+str(numprocessingcurrent)+"/Makefile";
  makeWrite[int(lanetopredict)] = open (makefilePath , 'w' ) ;
#BEGIN Writing indices
  #testing = json.loads(jsondata["indicesseq"]);
  #print testing;

  fileWithIndices=options.outdir+"/lane"+str(lanetopredict)+"/indices.txt";
  if(os.path.exists(fileWithIndices)):
    #move
    if not options.mock:
      os.rename(fileWithIndices,  fileWithIndices+"_"+timeNow);
    else:
      print "mv "+fileWithIndices +  " " +fileWithIndices+"_"+timeNow;

  indicesWrite = open (fileWithIndices , 'w' ) ;

  for jsonrow in jsondata["indicesseq"]:

    indicesWrite.write( jsonrow["p7"]+"\t" );
    if("p5" in jsonrow):
      indicesWrite.write( jsonrow["p5"]+"\t" );

    indicesWrite.write( jsonrow["name"] ); 
    indicesWrite.write( "\n" );

    #indicesWrite.write( jsondata["indicesseq"] ); 
  indicesWrite.close();

  fileWithIndicesRaw=options.outdir+"/lane"+str(lanetopredict)+"/indices.raw.txt";
  if(os.path.exists(fileWithIndices)):
    #move
    #os.rename(fileWithIndicesRaw,  fileWithIndicesRaw+"_"+timeNow);
    if not options.mock:
      os.rename(fileWithIndicesRaw,  fileWithIndicesRaw+"_"+timeNow);
    else:
      print "mv "+fileWithIndicesRaw +  " " +fileWithIndicesRaw+"_"+timeNow;

  indicesWriteRaw = open (fileWithIndicesRaw , 'w' ) ;

  for jsonrow in jsondata["indicesraw"]:
    indicesWriteRaw.write( jsonrow["name"] ); 

    indicesWriteRaw.write( jsonrow["p7"]+"\t" );
    if("p5" in jsonrow):
      indicesWriteRaw.write( jsonrow["p5"] );

    indicesWriteRaw.write( "\n" );

    #indicesWrite.write( jsondata["indicesseq"] ); 
  indicesWriteRaw.close();




  listOfFilesWithIndices[int(lanetopredict)] =  fileWithIndices;

#END Writing indices



#################################################
#GENERATE REPORT
#################################################



for lanetopredict in lanesToUse:
  makeWrite[int(lanetopredict)].write("SHELL := /bin/bash\n\nDefault:\tall\n\n");  
  if(not listOfFilesBasecall[lanetopredict]):#empty
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/Report/Summary.html");

  #if(jsondata["sequencer"] == "miseq"):
  makeWrite[int(lanetopredict)].write(outBaseDirectory+"/Report/Summary.html:\n\t"+rtaReportHask+" -o "+outBaseDirectory+"/Report/ "+illuminareaddir+"/"+jsondata["runid"]+"\n\n");
  #else:
  #  makeWrite[int(lanetopredict)].write(outBaseDirectory+"/Report/Summary.html:\n\t"+rtaReportPython+" -o "+outBaseDirectory+"/Report/ "+illuminareaddir+"/"+jsondata["runid"]+"\n\n");




#################################################
#BASECALLING
#################################################
#BUSTARD
if(jsondata["bustard"]):
  for lanetopredict in lanesToUse:
    #print lanetopredict;
    listOfBAMfilesToCheck[lanetopredict].append(outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished");
    listOfFilesBasecall[lanetopredict].append(outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam");
    listOfFilesBasecall[lanetopredict].append(outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished");
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished");

    makeWrite[int(lanetopredict)].write(outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished: "+illuminareaddir+"/"+jsondata["runid"]+"/Run.completed \n\t"+
                    BCL2BAM+
                    " -f "+jsondata["cyclesread1"]+
                    " -r "+jsondata["cyclesread2"]+
                    " -i "+jsondata["cyclesindx1"]+
                    " -j "+jsondata["cyclesindx2"]+
                    " -p "+illuminareaddir+"/"+jsondata["runid"]+"/Data/Intensities/"+
                    " -b "+illuminareaddir+"/"+jsondata["runid"]+"/Data/Intensities/BaseCalls/"+
                    " -o "+outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam"+
                    " -e "+jsondata["expname"]+
                    " -l "+str(lanetopredict)+"\n\n");
#FREEIBIS    
if(jsondata["freeibis"]):
  finishedFiles=[];
  for lanetopredict in lanesToUse:
    #print lanetopredict;
    listOfBAMfilesToCheck[lanetopredict].append(outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished");
    listOfFilesBasecall[lanetopredict].append(outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam");
    listOfFilesBasecall[lanetopredict].append(outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished");

    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished");

    finishedFiles.append(outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished");
    makeWrite[int(lanetopredict)].write(outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished:  "+illuminareaddir+"/"+jsondata["runid"]+"/Run.completed \n\t"+"sleep "+str(randint(10,100))+"\n\t"+"ulimit  -s 8192\n\t"+ #used to appease SOAP
                    FREEIBIS+"/runBaseCalling.py "            
                    " --NoFinishCheck -c "+str(options.cores)+" "+
                    " -o "+outBaseDirectory+"/Ibis/Raw_Sequences/ "+
                    " -b "+illuminareaddir+"/"+jsondata["runid"]+"/Data/Intensities/BaseCalls/ "+
                    " --temp="+options.tempdir+" "
                    " -e "+jsondata["expname"]+" "
                    " --indexlength="+str(jsondata["cyclesindx1"] )+" --2nd_indexlength="+str( jsondata["cyclesindx2"] )+" ");
    start="";
    end  ="";
    clen = len(jsondata["key1"+str(numprocessingcurrent)]);
    clen2= len(jsondata["key2"+str(numprocessingcurrent)]);

    if(int(jsondata["cyclesread2"]) > 0): #double index
      start = str(clen+1)+","+str(int(jsondata["cyclesread1"])+int(jsondata["cyclesindx1"])+1+clen2);
      end   = str(int(jsondata["cyclesread1"]))+","+str(int(jsondata["cyclesread1"])+int(jsondata["cyclesindx1"])+int(jsondata["cyclesread2"]));
#      if (clen > 0) and (clen2 == 0):  
#        start=str(clen+2)+","+str(int(jsondata["cyclesread1"])+int(jsondata["cyclesindx1"])+1+clen2);
#      elif (clen == 0) and (clen2 > 0): 
#        start=str(clen+1)+","+str(int(jsondata["cyclesread1"])+int(jsondata["cyclesindx1"])+2+clen2);
#      else:  #both keys are zero
#        start=str(clen+2)+","+str(int(jsondata["cyclesread1"])+int(jsondata["cyclesindx1"])+2+clen2);
#        end=str(int(jsondata["cyclesread1"]))+","+str(int(jsondata["cyclesread1"])+int(jsondata["cyclesindx1"])+int(jsondata["cyclesread2"]));
    else:
      start = str(clen+1);
      end   = str(int(jsondata["cyclesread1"]));

    #checking lanes
    lanesToUseTrain=[];
        
    if(jsondata["spikedin"]):
      makeWrite[int(lanetopredict)].write(" --control_index="+jsondata["ctrlindex"] +" ");
      if(options.trainlanes):
        lanesToUseTrain = parse_rangestr(options.trainlanes);
      else:
        lanesToUseTrain=",".join( range(1,jsondata["LaneCount"]) );
    else:
      lanesToUseTrain = jsondata["lanesdedicated"]

    if(len(lanesToUseTrain) == 0):
      print "Error: no lanes selected for training";
      sys.exit(1);

    #checking tiles
    tilesToUse=[];

    for tile in range(1,int(jsondata["TileCount"])+1):
        for swath in range(1,int(jsondata["SwathCount"])+1):
          for surface in range(1,int(jsondata["SurfaceCount"])+1):
            tilesToUse.append("%d%d%02d" % (surface, swath, tile));

    
    makeWrite[int(lanetopredict)].write(" --start="+start+" --end="+end+"  -l '"+(",".join(lanesToUseTrain))+"'    -t '"+(",".join(tilesToUse))+"' --recalibration --plotqual --lock ");
    makeWrite[int(lanetopredict)].write(" --reference="+jsondataConf["phixref"]+" ");
    makeWrite[int(lanetopredict)].write("\n\n");

    #ERROR PROFILE
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/Ibis/error_profile.pdf");

    makeWrite[int(lanetopredict)].write(outBaseDirectory+"/Ibis/error_profile.pdf:\t"+outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\n");
    makeWrite[int(lanetopredict)].write("\t"+FREEIBIS+"/plot_error.cmd.R "+outBaseDirectory+"/Ibis/Raw_Sequences/Models/SVMlight_models.index "+outBaseDirectory+"/Ibis/error_profile.pdf\n\n");
    
    
#makeWrite[int(lanetopredict)].write("\n");






















BasecallersUsed=[];

if(jsondata["bustard"]):
  BasecallersUsed.append("Bustard");

if(jsondata["freeibis"]):
  BasecallersUsed.append("Ibis");


#FastQC:



for baseCaller in BasecallersUsed:

#################################################
#FASTQC
#################################################

  for lanetopredict in lanesToUse:
    listOfFilesBasecall[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/FastQC/s_"+str(lanetopredict)+"_sequence_fastqc.zip");
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/FastQC/s_"+str(lanetopredict)+"_sequence_fastqc.zip");

    makeWrite[int(lanetopredict)].write(outBaseDirectory+"/"+baseCaller+"/FastQC/s_"+str(lanetopredict)+"_sequence_fastqc.zip:\t"+outBaseDirectory+"/"+baseCaller+"/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\n");
    makeWrite[int(lanetopredict)].write("\t"+FastQCreport+" -f bam -q -o "+outBaseDirectory+"/"+baseCaller+"/FastQC/ "+outBaseDirectory+"/"+baseCaller+"/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam\n\n");
    







































#################################################
#MERGE/TRIM + QC FILTER + RGASSIGN
#################################################

  for lanetopredict in lanesToUse:
    listOfBAMfilesToCheck[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished");
    listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam");
    listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished");

    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished");

    makeWrite[int(lanetopredict)].write(outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished:\t"+outBaseDirectory+"/"+baseCaller+"/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\n");



    conversion_str = MergeReads + " ";

    if(int(jsondata["cyclesread2"]) > 0):
      conversion_str += " -k '%s,%s' -f %s -s %s -c %s "%( jsondata["key1"+str(numprocessingcurrent)] ,jsondata["key2"+str(numprocessingcurrent)],jsondata["adapter1"+str(numprocessingcurrent)] ,jsondata["adapter2"+str(numprocessingcurrent)],jsondata["chimeras"+str(numprocessingcurrent)])
      if jsondata["mergeoverlap"+str(numprocessingcurrent)] :
        conversion_str += "--ancientdna "
    else:

      conversion_str += " -k '%s' -f '%s'  -c '%s' "%( jsondata["key1"+str(numprocessingcurrent)] ,jsondata["adapter1"+str(numprocessingcurrent)] ,jsondata["chimeras"+str(numprocessingcurrent)])
    #    if lane in oneErrorKey: conversion_str += "--allowMissing "

    conversion_str += " --log  "+outBaseDirectory+"/"+baseCaller+"/QC/s_"+str(lanetopredict)+"_merge.log";
    listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/s_"+str(lanetopredict)+"_merge.log");

    #conversion_str += " -t %d "%(adapterTrim[lane]);
    conversion_str += " -u  -o /dev/stdout "; #to send to pipe
    conversion_str += outBaseDirectory+"/"+baseCaller+"/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam "; #actual input

    if(jsondata["filterseqexp"+str(numprocessingcurrent)]):
       conversion_str += "| "+BAMFilter+" -u -o /dev/stdout ";
    
       if( jsondata["filterfrequency"+str(numprocessingcurrent)]):
          conversion_str += "--frequency --comp_cutoff=%.4f "%(jsondata["frequencycutoff"+str(numprocessingcurrent)])
          conversion_str += "--freq  "+outBaseDirectory+"/"+baseCaller+"/QC/filter/s_"+str(lanetopredict)+"_freq.dat"            
          listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/filter/s_"+str(lanetopredict)+"_freq.dat");

       if( jsondata["filterentropy"+str(numprocessingcurrent)]):
          conversion_str += "--entropy --comp_cutoff=%.4f "%(jsondata["entropycutoff"+str(numprocessingcurrent)])
          conversion_str += "--ent  "+outBaseDirectory+"/"+baseCaller+"/QC/filter/s_"+str(lanetopredict)+"_ent.dat";     
          listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/filter/s_"+str(lanetopredict)+"_ent.dat");

       conversion_str += " --like "+outBaseDirectory+"/"+baseCaller+"/QC/filter/s_"+str(lanetopredict)+"_likelihood.dat ";
       conversion_str += " --log "+outBaseDirectory+"/"+baseCaller+"/QC/filter/s_"+str(lanetopredict)+"_filter.log "; 
       listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/filter/s_"+str(lanetopredict)+"_likelihood.dat ");
       listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/filter/s_"+str(lanetopredict)+"_filter.log ");

       conversion_str += " /dev/stdin ";
    conversion_str += "| "+IndexReassign+" --summary "+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rg_summary.txt" 
    listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rg_summary.txt"); 


    conversion_str +=  " -i "+listOfFilesWithIndices[lanetopredict]+ " ";
    #conversion_str +=  " --mm "+str(jsondata["mmrgassign"])+ " ";

    conversion_str += " --error "+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_unassigned.txt";
    listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_unassigned.txt");
    conversion_str += " --rgval "+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rgqual.dat";
    listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rgqual.dat");
    conversion_str += " --ratio "+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_ratio.dat";
    listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_ratio.dat");
    conversion_str += " -o "+outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam"


    conversion_str += " /dev/stdin "



    makeWrite[int(lanetopredict)].write("\t"+conversion_str+"\n\ttouch "+outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished\n\n");





























#################################################
#BWA ALIGN
#################################################

  max_threads=3;
  if(jsondata["usebwa"+str(numprocessingcurrent) ]  and str(jsondata["sequencer"]) == "miseq"):

    for lanetopredict in lanesToUse:
      listOfBAMfilesToCheck[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam");
      listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam");
      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam");


      makeWrite[int(lanetopredict)].write(outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".flgstx:\t"+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam\n");

      listOfFilesFinal[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".flgstx");
      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".flgstx");
      listOfFLGfilesToCheck[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".flgstx");

      conversion_str = "mv "+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam "+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".temp.bam\n\t"+"sam sort -m 16G -o "+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".temp.bam -i "+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam.bai -x "+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".flgstx -c "+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".covstat "+str(tempdir)+" > "+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam\n\t"+"rm -f "+outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".temp.bam";


      makeWrite[int(lanetopredict)].write("\t"+conversion_str+"\n\n");


      makeWrite[int(lanetopredict)].write(outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam:\t"+outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished\n");
      bwaparameter="";
      if(str(jsondata["parambwa"+str(numprocessingcurrent)]) == "default"):
        bwaparameter="";
      elif(str(jsondata["parambwa"+str(numprocessingcurrent)]) == "ancient") :
        bwaparameter=" -n 0.01 -o 2 -l 16500 ";
      else:
        print "unexpected bwa param";
        sys.exit(1);



      conversion_str = "bwa bam2bam -p 52690 -t %d -f %s -g  %s %s %s "%(max_threads,
                                                                outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam_",
                                                                BWAGENOMES+"/"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+"/bwa-0.4.9",
                                                                bwaparameter,
                                                                outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam");
      
      #conversion_str +=  ("\n\n"+outBaseDirectory+"/"+baseCaller+"/BWA/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".flgstx: "+outBaseDirectory+"/"+baseCaller+"/BWA/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam"+"\n\tsam flagstatx "+outBaseDirectory+"/"+baseCaller+"/BWA/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam  > "+outBaseDirectory+"/"+baseCaller+"/BWA/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".flgstx\n");

      #listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/BWA/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".flgstx");


      makeWrite[int(lanetopredict)].write("\t"+conversion_str+"\n\n");


















#################################################
#QC CONTROL
#################################################
  #listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/*");
  #listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/*");
  #listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/*");
  #listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/filter/*");
  #listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/insertsize/lane*/");
  #listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/*");

  for lanetopredict in lanesToUse:
#INSERT SIZE
    targetbaminsert=outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam";
    if(jsondata["usebwa"+str(numprocessingcurrent)]  and str(jsondata["sequencer"]) == "miseq"):
      targetbaminsert=outBaseDirectory+"/"+baseCaller+"/BWA/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence"+"_"+str(jsondata["parambwa"+str(numprocessingcurrent)])+"_"+str(jsondata["genomebwa"+str(numprocessingcurrent)])+".bam";

    makeWrite[int(lanetopredict)].write(outBaseDirectory+"/"+baseCaller+"/QC/insertsize/proc"+str(numprocessingcurrent)+"/lane"+str(lanetopredict)+"/insert_"+str(lanetopredict)+".dat:\t"+outBaseDirectory+"/"+baseCaller+"/QC/clusterTally_"+str(lanetopredict)+".OK"+"\n");
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/insertsize/proc"+str(numprocessingcurrent)+"/lane"+str(lanetopredict)+"/insert_"+str(lanetopredict)+".dat");


    cmdInsertSize   =  insertsizebi;
    cmdInsertSize  +=  " "+str(targetbaminsert);
    cmdInsertSize  +=  " > "+outBaseDirectory+"/"+baseCaller+"/QC/insertsize/proc"+str(numprocessingcurrent)+"/lane"+str(lanetopredict)+"/insert_"+str(lanetopredict)+".dat\n";
    cmdInsertSize  +=  "\t"+insertsizeR;
    cmdInsertSize  +=  "  "+outBaseDirectory+"/"+baseCaller+"/QC/insertsize/proc"+str(numprocessingcurrent)+"/lane"+str(lanetopredict)+"/insert_"+str(lanetopredict)+".dat ";
    cmdInsertSize  +=  "  "+outBaseDirectory+"/"+baseCaller+"/QC/insertsize/proc"+str(numprocessingcurrent)+"/lane"+str(lanetopredict)+"/\n";
    makeWrite[int(lanetopredict)].write("\t"+cmdInsertSize+"\n");
    listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/insertsize/proc"+str(numprocessingcurrent)+"/lane"+str(lanetopredict)+"/*");



    
#CLUSTER COUNT
    makeWrite[int(lanetopredict)].write(outBaseDirectory+"/"+baseCaller+"/QC/clusterTally_"+str(lanetopredict)+".OK:\t"+" ".join(listOfBAMfilesToCheck[lanetopredict])+" "+" ".join(listOfFLGfilesToCheck[lanetopredict])+"\n");
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/clusterTally_"+str(lanetopredict)+".OK");

    cmdTileCounter =  Tilecounter;
    cmdTileCounter += " -o "+str(outBaseDirectory+"/"+baseCaller+"/QC/");
    cmdTileCounter += " -b "+str(",".join(listOfBAMfilesToCheck[lanetopredict])).replace(".finished","");
    cmdTileCounter += " -r "+str(illuminareaddir+"/"+str(jsondata["runid"])+"/Data/Intensities/")
    cmdTileCounter += " -l "+str(lanetopredict);
    makeWrite[int(lanetopredict)].write("\t"+cmdTileCounter+"\n");
    listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/clusterTally_"+str(lanetopredict)+".OK");

#QC SCORE PLOT
    makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/s_"+str(lanetopredict)+"_sequenceA.pdf:\t"+outBaseDirectory+"/"+baseCaller+"/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\n");  
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/s_"+str(lanetopredict)+"_sequenceA.pdf");

    cmdQualPlot = "\t"+Qualplotter;
    cmdQualPlot += " -o "+outBaseDirectory+"/"+baseCaller+"/QC/";
    cmdQualPlot += " -i "+outBaseDirectory+"/"+baseCaller+"/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam";
    cmdQualPlot += " --outprefix="+"s_"+str(lanetopredict)+"_sequence";
    cmdQualPlot += " -t "+str(tempdir);
    makeWrite[int(lanetopredict)].write(cmdQualPlot+"\n");
    listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/s_"+str(lanetopredict)+"_sequenceA.pdf");
    listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/s_"+str(lanetopredict)+"_sequenceC.pdf");
    listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/s_"+str(lanetopredict)+"_sequenceG.pdf");
    listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/s_"+str(lanetopredict)+"_sequenceT.pdf");

#RG QUAL PLOT
    makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rgqual.pdf:\t"+outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished\n");  
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rgqual.pdf");

    cmdRGQual  = "\t"+IndexReassignRGR;
    cmdRGQual += " "+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rgqual.dat";
    cmdRGQual += " "+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rgqual.pdf";
    makeWrite[int(lanetopredict)].write(cmdRGQual+"\n");
    listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_rgqual.pdf");


#RG QUAL PLOT (2)
    makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_ratio.pdf: "+outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished\n\t");  
    listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_ratio.pdf");

    cmdRatioPlot  = IndexReassignRATIOR;
    cmdRatioPlot += " "+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_ratio.dat";
    cmdRatioPlot += " "+outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_ratio.pdf";
    makeWrite[int(lanetopredict)].write(cmdRatioPlot+"\n");
    listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/rg/s_"+str(lanetopredict)+"_ratio.pdf");

#SEQ LIKELIHOOD
    if(jsondata["filterseqexp"+str(numprocessingcurrent)]):
      makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/filter/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_likelihood.pdf: "+outBaseDirectory+"/"+baseCaller+"/Final_Sequences/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_sequence.bam.finished\n\t");  

      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/filter/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_likelihood.pdf");
      cmdQualPlot  = BAMFilterLIKER;
      cmdQualPlot += " "+outBaseDirectory+"/"+baseCaller+"/QC/filter/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_likelihood.dat";
      cmdQualPlot += " "+outBaseDirectory+"/"+baseCaller+"/QC/filter/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_likelihood.pdf";
      makeWrite[int(lanetopredict)].write(cmdQualPlot+"\n");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/filter/proc"+str(numprocessingcurrent)+"/s_"+str(lanetopredict)+"_likelihood.pdf");      




#################################
#EXTRACTING AND MAPPING CONTROLS#
#################################
    if(jsondata["spikedin"]):
#EXTRACTING
      makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_ctrl.bam:\t"+outBaseDirectory+"/"+baseCaller+"/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\n"); 

      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_ctrl.bam");
      cmdCtrlExt = Ctrlextract;
      cmdCtrlExt += " "+jsondata["ctrlindex"]+" ";
      cmdCtrlExt += " "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_ctrl.bam";
      cmdCtrlExt += " "+outBaseDirectory+"/"+baseCaller+"/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam";
      makeWrite[int(lanetopredict)].write("\t"+cmdCtrlExt+"\n");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_ctrl.bam");

#MAPPING
      makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam:\t"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_ctrl.bam\n");
      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam");

      cmdBWACTRL = "bwa bam2bam -t 1 -g "+CTRLGENOMEBWA+" -f "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam_"+" "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_ctrl.bam";
      makeWrite[int(lanetopredict)].write("\t"+cmdBWACTRL+"\n");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam");


#COMPUTING ERROR RATES
      makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.dat:\t"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam\n");    

      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.dat");
      cmdBWACTRLC = ERRORPERCYCLE+" "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam"+" > "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.dat";
      # cmdBWACTRLC = ERRORPERCYCLER+" "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.dat")+" "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.all.pdf")+" "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.type.pdf");
      makeWrite[int(lanetopredict)].write("\t"+cmdBWACTRLC+"\n");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.dat");



      makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.type.pdf:\t"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.dat\n");
      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.type.pdf");

      cmdBWACTRLCR = ERRORPERCYCLER+" "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.dat "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.pdf "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.type.pdf";
      makeWrite[int(lanetopredict)].write("\t"+cmdBWACTRLCR+"\n");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.pdf");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.error.type.pdf");


#QC SCORES, OBS VS PREDICTED
      readlengths=str(jsondata["cyclesread1"]);
      if(int(jsondata["cyclesread2"]) != 0):
        readlengths=readlengths+","+str(jsondata["cyclesread2"]);


      makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred6:\t"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam\n");    
      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred6 ");

      cmdDataOBSPRED = Predvsobs;
      cmdDataOBSPRED += " "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam"+" /dev/null "+readlengths;#replace with actual mask files for freeIbis
      makeWrite[int(lanetopredict)].write("\t"+cmdDataOBSPRED+"\n");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred1");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred2");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred3");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred4");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred5");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred6");

      makeWrite[int(lanetopredict)].write("\n"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred.dens.pdf:\t"+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred6\n");
      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred.pdf ");
      listOfTargetFiles[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspred.dens.pdf ");

      cmdDataOBSPREDR = "R CMD BATCH --vanilla --quiet ";
      cmdDataOBSPREDR += "  '--args "+outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam"+" "+jsondata["expname"]+" ' ";
      cmdDataOBSPREDR += PredvsobsR;
      cmdDataOBSPREDR += " /dev/null ";                           
      makeWrite[int(lanetopredict)].write("\t"+cmdDataOBSPREDR+"\n");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspredcycle1cycle.pdf");
      listOfFilesQC[lanetopredict].append(outBaseDirectory+"/"+baseCaller+"/QC/qscores/s_"+str(lanetopredict)+"_control.bam.baseobspredcycle2cycle.pdf");


#handle_jobs("chmod a+rX -R -f "+options.outpath+options.folderID)


typesofclean=[];

firstPassClean=True;

for lanetopredict in lanesToUse: 

  makeWrite[int(lanetopredict)].write("\n"+"all:\t");

  makeWrite[int(lanetopredict)].write((" ".join(listOfTargetFiles[lanetopredict]))+"\n\n");
  

  if(firstPassClean):
    typesofclean.append("cleanall");
  makeWrite[int(lanetopredict)].write("\n"+"cleanall:\n"+"\trm -fv ");

  makeWrite[int(lanetopredict)].write((" ".join(listOfFilesBasecall[lanetopredict]))+" "+
                                      (" ".join(listOfFilesFinal[lanetopredict]))+" "+
                                      (" ".join(listOfFilesBWA[lanetopredict]))+" "+
                                      (" ".join(listOfFilesQC[lanetopredict]))+"\n\n");

  
  if(firstPassClean):
    typesofclean.append("cleanfromfinal");
  makeWrite[int(lanetopredict)].write("\n"+"cleanfromfinal:\n"+"\trm -fv ");

  makeWrite[int(lanetopredict)].write(#(" ".join(listOfFilesBasecall[lanetopredict]))+" "+
                                      (" ".join(listOfFilesFinal[lanetopredict]))+" "+
                                      (" ".join(listOfFilesBWA[lanetopredict]))+" "+
                                      (" ".join(listOfFilesQC[lanetopredict]))+"\n\n");

  if(firstPassClean):
    typesofclean.append("cleanfrombwa");
  makeWrite[int(lanetopredict)].write("\n"+"cleanfrombwa:\n"+"\trm -fv ");

  makeWrite[int(lanetopredict)].write(#(" ".join(listOfFilesBasecall[lanetopredict]))+" "+
                                      #(" ".join(listOfFilesFinal[lanetopredict]))+" "+
                                      (" ".join(listOfFilesBWA[lanetopredict]))+" "+
                                      (" ".join(listOfFilesQC[lanetopredict]))+"\n\n");

  if(firstPassClean):
    typesofclean.append("cleanfromqc");
  makeWrite[int(lanetopredict)].write("\n"+"cleanfromqc:\n"+"\trm -fv ");

  makeWrite[int(lanetopredict)].write(#(" ".join(listOfFilesBasecall[lanetopredict]))+" "+
                                      #(" ".join(listOfFilesFinal[lanetopredict]))+" "+
                                      #(" ".join(listOfFilesBWA[lanetopredict]))+" "+
                                      (" ".join(listOfFilesQC[lanetopredict]))+"\n\n");

  
  makeWrite[int(lanetopredict)].write("\n"+"sendemail:\t"+outBaseDirectory+"/"+baseCaller+"/QC/clusterTally_"+str(lanetopredict)+".OK\n");
  makeWrite[int(lanetopredict)].write("\tphp "+alibdir+"/"+sendemail+" "+jsondata["email"]+" "+jsondata["runid"]+" "+illuminawritedir+"/"+jsondata["runid"]+" "+str(lanetopredict)+"\n");

  makeWrite[int(lanetopredict)].close();

  firstPassClean=False;

#################################################
#                                               #
#               END  MAKEFILE                   #
#                                               #
#################################################

#general makefile

makefilePath=options.outdir+"/Makefile";
makeGeneralWrite= open (makefilePath , 'w' ) ;


lanelist=[];
for lanepred in range(1,int(jsondata["LaneCount"])+1):
  lanelist.append("lane"+str(lanepred));
  for procnum in range(1,100):
    if( os.path.exists("lane"+str(lanepred)+"/proc"+str(procnum)) ):
      lanelist.append( "lane"+str(lanepred)+"/proc"+str(procnum) );
    else:
      break;

makeGeneralWrite.write("SUBDIRS = "+(" ".join(lanelist))+"\n\n");

makeGeneralWrite.write(".PHONY: subdirs $(SUBDIRS)\n\n");

makeGeneralWrite.write("subdirs: $(SUBDIRS)\n\n");

makeGeneralWrite.write("$(SUBDIRS):\n");
makeGeneralWrite.write("\t$(MAKE) -j 1 -C $@ all\n\n");

makeGeneralWrite.write("sendemail:\n");
makeGeneralWrite.write("\tfor dir in $(SUBDIRS); do ");
makeGeneralWrite.write(" $(MAKE) -C $$dir sendemail; ");
makeGeneralWrite.write(" done\n\n");




#I need to fix this:
for cleantarg in typesofclean:
  makeGeneralWrite.write(cleantarg+":\n");
  makeGeneralWrite.write("\tfor dir in $(SUBDIRS); do ");
  makeGeneralWrite.write(" $(MAKE) -C $$dir "+cleantarg+"; ");
  makeGeneralWrite.write(" done \n\n");
#
#

makeGeneralWrite.close();
