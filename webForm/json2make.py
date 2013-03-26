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
#import xml.parsers.expat
import xml.etree.ElementTree as ET
from optparse import OptionParser
from optparse import OptionGroup

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



alibdir="";
illuminareaddir="";
illuminawritedir=""


BCL2BAM              = "BCL2BAM/bcl2bam"
FREEIBIS             = ""

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
FastQCreport         = "FastQC/fastqc"




def chomp(s):
  return s.rstrip('\n');





#BEGIN READ CONFIG FILE
pathToConfig=sys.argv[0];
pathToConfig=os.path.dirname(os.path.abspath(pathToConfig))+"/config.xml";
alibdir=os.path.dirname(os.path.abspath( sys.argv[0]+"/../"))

print pathToConfig;
print alibdir;
try:
  tree =ET.parse(pathToConfig);
except IOError:
  print "Cannot open XML config file "+pathToConfig;
  sys.exit(1);

XMLconfig        = tree.getroot()


illuminawritedir = XMLconfig.find("illuminawritedir").text;
illuminareaddir  = XMLconfig.find("illuminareaddir").text;
FREEIBIS         = XMLconfig.find("freeibispath").text;


def checkPrograms():
  global alibdir;
  

parser = OptionParser(usage="usage: %prog [options] <json doc>");
parser.add_option("-o", "--outdir", dest="outdir", help="Out directory to put the makefile");
parser.add_option("--lanes", help="Consider only a subset of lanes for training freeIbis (default '1-8')",default="1-8")

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
  


outBaseDirectory=os.path.abspath(options.outdir+"/../");


#readConfig(sys.argv[0]);



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

FastQCreport = alibdir+"/"+FastQCreport;
if not os.path.exists(FastQCreport):
  print "Required executable file not found "+FastQCreport;
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


#BEGIN Writing indices

fileWithIndices=options.outdir+"indices.txt";
indicesWrite = open (fileWithIndices , 'w' ) ;
indicesWrite.write( jsondata["indicesseq"] ); 
indicesWrite.close();

#END Writing indices






#BEGIN WRITING TARGET DIRECTORIES    
targetdir=outBaseDirectory+"/Report/";    
#makeWrite.write(targetdir+"\n\tmkdir "+targetdir);
makedirs(outBaseDirectory+"/Report/");


allsubdir=["FastQC/",
           "QC/",
           "QC/qscores/",
           "QC/rg/",
           "QC/filter/",
           "Raw_Sequences/",
           "Final_Sequences/",
           "BWA",
           "FastQ"];


if(jsondata["bustard"]):
  
  for subdir in allsubdir:
    targetdir=outBaseDirectory+"/Bustard/"+subdir;    
    makedirs(targetdir);

if(jsondata["freeibis"]):
  
  for subdir in allsubdir:
    targetdir=outBaseDirectory+"/Ibis/"+subdir;    
    makedirs(targetdir);


listOfBAMfilesToCheck={};

if(jsondata["bustard"]):
  for lanetopredict in jsondata["lanes"]:
    listOfBAMfilesToCheck[lanetopredict]=[];

if(jsondata["freeibis"]):
  for lanetopredict in jsondata["lanes"]:
    listOfBAMfilesToCheck[lanetopredict]=[];

#END WRITING TARGET DIRECTORIES    


#################################################
#                                               #
#              BEGIN MAKEFILE                   #
#                                               #
#################################################

makefilePath=options.outdir+"Makefile";
makeWrite = open (makefilePath , 'w' ) ;

#GENERATE REPORT
makeWrite.write(outBaseDirectory+"/Report/Summary.htm:\n\tgenerate -o "+outBaseDirectory+"/Report/ "+illuminareaddir+"/"+jsondata["runid"]+"\n\n");


#BASECALLING
if(jsondata["bustard"]):
  for lanetopredict in jsondata["lanes"]:
    print lanetopredict;
    makeWrite.write(outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam:\n\t"+
                    BCL2BAM+
                    " -f "+jsondata["cyclesread1"]+" -r "+jsondata["cyclesread2"]+
                    " -i "+jsondata["cyclesindx1"]+
                    " -j "+jsondata["cyclesindx2"]+
                    " -p "+illuminareaddir+"/"+jsondata["runid"]+"/"+
                    " -b "+illuminareaddir+"/"+jsondata["runid"]+"/Data/Intensities/BaseCalls/"+
                    " -o "+outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam"+
                    #" -e "+str(options.expID)+
                    " -l "+str(lanetopredict)+"\n\n");
    
if(jsondata["freeibis"]):
  for lanetopredict in jsondata["lanes"]:
    print lanetopredict;
    makeWrite.write(outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam:\n\t"+
                    FREEIBIS+             
                    " --NoFinishCheck -c "+str(options.cores)+" "+

                    #" -e "+options.expID+"
                    " --indexlength="+str(jsondata["cyclesindx1"] )+" --2nd_indexlength="+str( jsondata["cyclesindx2"] ));
    start="";
    end  ="";
    clen = len(jsondata["key1"]);
    clen2= len(jsondata["key2"]);

    if(jsondata["cyclesread2"] > 0):
      if (clen > 0) and (clen2 == 0): 
        start=str(clen+2)+","+str(jsondata["cyclesread1"]+jsondata["cyclesindx1"]+1+clen2);
      elif (clen == 0) and (clen2 > 0): 
        start=str(clen+1)+","+str(jsondata["cyclesread1"]+jsondata["cyclesindx1"]+2+clen2);
      else:  #both keys are zero
        start=str(clen+2)+","+str(jsondata["cyclesread1"]+jsondata["cyclesindx1"]+2+clen2);
        end=str(jsondata["cyclesread1"])+","+str(jsondata["cyclesread1"]+jsondata["cyclesindx1"]+jsondata["cyclesread2"]);
    else:
      if clen > 0: 
        start=str(clen+2)
      else: 
        start=str(clen+1)
        end=str(readlength)

        makeWrite.write("--start="+start+" --end="+end+"  -l "+options.train_lanes+" --recalibration --plotqual --lock ");
        makeWrite.write(" --control_index="+XMLconfig.find("controlindex").text+" ");
        makeWrite.write(" -r="+XMLconfig.find("phixref").text+" ");

    #ERROR PROFILE
    #makeWrite.write(outBaseDirectory+"/Ibis/error_profile.pdf/:\t");
    #for lanetopredict in jsondata["lanes"]:
    #    makeWrite.write(outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\t");
    #makeWrite.write("\n");
    #makeWrite.write("cd "+ibisfolder+"Models; R --vanilla --quiet < "+IbisInstallFolder+"plot_error.R; cp error_profile.pdf "+options.outpath+options.folderID+"/Ibis/    
    #--start="+start+" --end="+end+" --indexlength="+str(maxlengthindex)+" --2nd_indexlength="+str(maxlengthindex2)+" -l "+options.train_lanes+"
    
    
makeWrite.write("\n");

#FastQC:

if(jsondata["bustard"]):
  for lanetopredict in jsondata["lanes"]:
    makeWrite.write(outBaseDirectory+"/Bustard/FastQC/s_"+str(lanetopredict)+"_sequence_fastqc.zip:\t"+outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\n");
    makeWrite.write("\t"+FastQCreport+" -f bam -q -o "+outBaseDirectory+"/Bustard/FastQC/ "+outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam\n");

if(jsondata["freeibis"]):
  for lanetopredict in jsondata["lanes"]:
    makeWrite.write(outBaseDirectory+"/Ibis/FastQC/s_"+str(lanetopredict)+"_sequence_fastqc.zip:\t"+outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\n");
    makeWrite.write("\t"+FastQCreport+" -f bam -q -o "+outBaseDirectory+"/Ibis/FastQC/ "+outBaseDirectory+"/Ibis/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam\n");

#merging/trimming
makeWrite.write("\n");


if(jsondata["bustard"]):
  for lanetopredict in jsondata["lanes"]:
    listOfBAMfilesToCheck[lanetopredict].append(outBaseDirectory+"/Bustard/Final_Sequences/s_"+str(lanetopredict)+"_sequence.bam");
    makeWrite.write(outBaseDirectory+"/Bustard/Final_Sequences/s_"+str(lanetopredict)+"_sequencebam:\t"+outBaseDirectory+"/Bustard/Raw_Sequences/s_"+str(lanetopredict)+"_sequence.bam.finished\n");



    conversion_str = MergeReads + " ";
    if ispaired:
      conversion_str += " -k '%s,%s' -f '%s' -s '%s' -c '%s' "%( jsondata["key1"] ,jsondata["key2"],jsondata["adapter1"] ,jsondata["adapter2"],jsondata["chimeras"])
      if jsondata["mergereads"] == "true"
        conversion_str += "--mergeoverlap "
    else:
      #conversion_str += " -k '%s' -f '%s' -c '%s' "%(lane_keys[lane],lane_adapter[lane],",".join(lane_chimera[lane]));
      conversion_str += " -k '%s' -f '%s'  -c '%s' "%( jsondata["key1"] ,jsondata["adapter1"] ,jsondata["chimeras"])
#    if lane in oneErrorKey: conversion_str += "--allowMissing "

    conversion_str += " --log  "+outBaseDirectory+"/Bustard/QC/s_"+lane+"_merge.log";
    #conversion_str += " -t %d "%(adapterTrim[lane]);
    conversion_str += " -u  -o /dev/stdout "; #to send to pipe
    conversion_str += outBaseDirectory+"/Bustard/Raw_Sequences/s_"+lane+"_sequence.bam "; #actual input

      
#    if not (qualTrim[lane] == 0 and trimLastN[lane] == 0 and qualCutoff[lane] == -1 and compCutoff[lane] == -1 and lengthCutoff[lane] == 0):
#      conversion_str += "| "+BAMFilter+" -u -o /dev/stdout "
#
#      if compCutoff[lane] != -1:
#        if lane in cMethod: 
#          conversion_str += "--frequency --comp_cutoff=%.4f "%(compCutoff[lane])            
#          conversion_str += "--freq  "+options.outpath+options.folderID+"/Bustard/QC/filter/s_"+lane+"_freq.dat"            
#        else: 
#          conversion_str += "--entropy --comp_cutoff=%.4f "%(compCutoff[lane])
#          conversion_str += "--ent  "+options.outpath+options.folderID+"/Bustard/QC/filter/s_"+lane+"_ent.dat"            
#
#      if lengthCutoff[lane] != 0:
#        if lane in qualTrim and qualTrim[lane] != 0:
#          if lane in lMethod: 
#            conversion_str += " --min_length %d --max_length %d "%(qualTrim[lane],lengthCutoff[lane])
#          else: 
#            conversion_str += " --min_length %d --max_length -1 "%(max(lengthCutoff[lane],qualTrim[lane]))
#        else:
#          if lane in lMethod: 
#            conversion_str += " --min_length 0 --max_length %d "%(lengthCutoff[lane])
#          else: 
#            conversion_str += " --min_length %d --max_length -1 "%(lengthCutoff[lane])
#
#      if qualCutoff[lane] != -1:
#        if qualNumberCutoff[lane] == -1: 
#          conversion_str += " --average --qual_cutoff=%d "%(qualCutoff[lane]) #why would we want to use that ?  Feb 05--GR
#          print "--average for filtering not implemented "
#          sys.exit()
#
#        elif lane in qualTrim and qualTrim[lane] != 0: 
#          conversion_str += " --trim --qual_cutoff=%d --qual_number=%d "%(qualCutoff[lane],qualNumberCutoff[lane]) #not implemented yet Feb 05 --GR
#          print "--trim for filtering not implemented "
#          sys.exit()
#        else: 
#          #conversion_str += " --quality -c=%f  "%(qualCutoff[lane],qualNumberCutoff[lane])
#          conversion_str += " --like "+options.outpath+options.folderID+"/Bustard/QC/filter/s_"+lane+"_likelihood.dat" 
#
#      if trimLastN[lane] != 0:
#        print "--clip for filtering not implemented "
#        sys.exit()
#        conversion_str += " --clip=-%d "%abs(trimLastN[lane]) #why would we want to use that ?  Feb 05--GR
#      
#      conversion_str += " /dev/stdin | "+IndexReassign+" --summary "+options.outpath+options.folderID+"/Bustard/QC/rg/s_"+lane+"_rg_summary.txt" 
#
#      if (lengthindex[lane] > 0):
#        conversion_str += " -i %s "%hindexfiles[lane]
#        if lane in options.noindexdist:
#          print " option no longer supported "                                                                                                                                
#          sys.exit()   
#          conversion_str += "--no_skip_first_base --no_mutants --no_Ns "
##      if lane in indexQualCutoff:
##        conversion_str += " --quality=%i "%indexQualCutoff[lane]
#      conversion_str += " --error "+options.outpath+options.folderID+"/Bustard/QC/rg/s_"+lane+"_unassigned.txt";
#      conversion_str += " --rgval "+options.outpath+options.folderID+"/Bustard/QC/rg/s_"+lane+"_rgqual.dat";
#      conversion_str += " --ratio "+options.outpath+options.folderID+"/Bustard/QC/rg/s_"+lane+"_ratio.dat";
#      conversion_str += " -o "+options.outpath+options.folderID+"/Bustard/Final_Sequences/s_"+lane+"_sequence.bam"
#      conversion_str += " /dev/stdin "
#      handle_jobs(conversion_str)




makeWrite.close();

#################################################
#                                               #
#               END  MAKEFILE                   #
#                                               #
#################################################
