#!/usr/bin/python

import sys,os
import json
import time
import datetime
import subprocess
import tempfile

from optparse import OptionParser
from optparse import OptionGroup

splitByRG            = "pipeline/splitByRG";
bam2fastq            = "BCL2BAM2FASTQ/bam2fastq/bam2fastq";
runSGEandWait        = "rnaseq/runSGEandWait.py";

                           
def chomp(s):
  return s.rstrip('\n');

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

parser = OptionParser(usage="This program splits a BAM file into RG and generates a makefile to run on the SGE.\n\ncommand: %prog [options]");
#parser.add_option("-i", "--infile",  dest="infile", help="input");
input = OptionGroup(parser, "Input","Input files");
input.add_option("-i", "--infile", dest="infile", help="Input BAM file");
input.add_option("-g", "--bowtie",    dest="bowtieIndex", help="Bowtie index");
input.add_option("-c", "--cufflinks", dest="cufflinksgtf", help="Cufflinks gtf file");

output = OptionGroup(parser, "Output","Output options");
output.add_option("-o", "--outdir", dest="outdir", help="Output directory");

misc = OptionGroup(parser, "Miscellaneous","Miscellaneous");
misc.add_option("--tmp", dest="tempdir", default=jsondataConf["tempdirectory"], help="Temporary directory");
misc.add_option("--mock", dest="mock", help="Do a mock run for testing ",default=False,action="store_true");
misc.add_option("--nosge", dest="nosge", help="Do not use the SGE, run on the node ",default=False,action="store_true");
misc.add_option("--nosplit", dest="nosplit", help="Do not split the BAM file, use the current BAM file in the directory specified");


parser.add_option_group(input);
parser.add_option_group(misc);
parser.add_option_group(output);

(options, args) = parser.parse_args();

                                                                                                                                         
def handle_jobs(cjob):
  global options
  if options.mock:
    return None

  jobcreated=subprocess.Popen(cjob,shell=True); 
  jobcreated.wait()

  if(jobcreated.returncode != 0): #has finished but wrong code
    print "Job failed "+cjob+" failed";
    sys.exit(1);
     
  #return jobcreated; 


if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)


if (options.infile == None) :
  print "Need to specify input file ";
  sys.exit(1);


if (options.bowtieIndex == None) :
  print "Need to specify bowtie index file ";
  sys.exit(1);
  
if (options.cufflinksgtf == None) :
  print "Need to specify cufflinks gtf file ";
  sys.exit(1);

alibdir          = jsondataConf["alibdir"];
tophat           = jsondataConf["tophat"];
cufflinks        = jsondataConf["cufflinks"];

splitByRG = alibdir+"/"+splitByRG;
if not os.path.exists(splitByRG):
  print "Required executable file not found "+splitByRG;
  sys.exit(1);

bam2fastq = alibdir+"/"+bam2fastq;
if not os.path.exists(bam2fastq):
  print "Required executable file not found "+bam2fastq;
  sys.exit(1);

runSGEandWait = alibdir+"/"+runSGEandWait;
if not os.path.exists(runSGEandWait):
  print "Required executable file not found "+runSGEandWait;
  sys.exit(1);

#check for bowtie index
if not os.path.exists(options.bowtieIndex+".1.bt2"):
  print "Bowtie index file not found "+options.bowtieIndex+".1.bt2";
  sys.exit(1);

if not os.path.exists(options.bowtieIndex+".rev.1.bt2"):
  print "Bowtie index file not found "+options.bowtieIndex+".rev.1.bt2";
  sys.exit(1);

#check for gtf files
if not os.path.exists(options.cufflinksgtf):
  print "Bowtie index file not found "+options.cufflinksgtf+"";
  sys.exit(1);

########################################

#make temp dir

#split RG
if (len(options.nosplit)==0) :
  tempdirname=tempfile.mkdtemp(prefix=options.tempdir);
  #print tempdirname;
  if(tempdirname[-1:] != "/"):
    tempdirname= tempdirname+"/";

  cmd = splitByRG+" "+options.infile+" "+tempdirname+"/rna";
  sys.stderr.write("spliting BAM file using command : "+str(cmd)+"\nWait this may take a while");
  handle_jobs(cmd);
  sys.stderr.write("done spliting\n");
  #write makefile

else:
  tempdirname=options.nosplit;
  if(tempdirname[-1:] != "/"):
    tempdirname= tempdirname+"/";

print "SHELL := /bin/bash\n\nDefault:\tall\n\n";  

listOfTargetFiles=[];
listOfTargetFilesTH=[];

fileFoundIndex=0;
for filefound in os.listdir(tempdirname):
  #print "file "+filefound[-3:];
  if(filefound[-3:] == "bam"):
   # print "file2 "+filefound;
    filefound = tempdirname+filefound;
    sampleid=filefound[ len(str(tempdirname+"/rna")) :-4];
    #print "sample #"+sampleid+"#";
    #os.
    if(sampleid.startswith("control")):
      continue;
    if(sampleid.startswith("controlA")):
      continue;
    if(sampleid.startswith("controlC")):
      continue;
    if(sampleid.startswith("controlG")):
      continue;
    if(sampleid.startswith("controlT")):
      continue;

    if(sampleid.startswith("unknown")):
      continue;

    #define target
    listOfTargetFiles.append(filefound[:-4]+"_r1.fq.gz");  
    listOfTargetFilesTH.append(filefound[:-4]+"_r1.fq.gz");  

    print "\n\n"+filefound[:-4]+"_r1.fq.gz:\n\t"+bam2fastq+" "+filefound+" "+filefound[:-4];
    
    #os.mkdir(filefound[:-4]+"/tophat/");
    
    listOfTargetFiles.append(filefound[:-4]+"/tophat/accepted_hits.bam");  
    listOfTargetFilesTH.append(filefound[:-4]+"/tophat/accepted_hits.bam");  

    #TODO add runSGE
    print "\n\n"+filefound[:-4]+"/tophat/accepted_hits.bam: "+filefound[:-4]+"_r1.fq.gz\n\tmkdir -p "+filefound[:-4]+"/tophat/"+"";
    if options.nosge:
      print "\t"+tophat+"   -o "+filefound[:-4]+"/tophat/ "+options.bowtieIndex+" "+filefound[:-4]+"_r1.fq.gz "+" "+filefound[:-4]+"_r2.fq.gz\"";
    else:
      print "\t"+runSGEandWait+" --param=\" -S /bin/bash -l \\\"h_vmem=5500M,virtual_free=5500M\\\" -V -R y -pe smp 1- \" --tmp="+tempdirname+"tophat_"+str(++fileFoundIndex)+" -c \""+tophat+" -p \\$$NSLOTS  -o "+filefound[:-4]+"/tophat/ "+options.bowtieIndex+" "+filefound[:-4]+"_r1.fq.gz "+" "+filefound[:-4]+"_r2.fq.gz\"";
#bam2fastq+" "+filefound+" "+filefound[-4:];
    
    listOfTargetFiles.append(filefound[:-4]+"/cuffdiff_output/genes.fpkm_tracking");  
    #/cuffdiff_output/genes.fpkm_tracking
    print "\n\n"+filefound[:-4]+"/cuffdiff_output/genes.fpkm_tracking: "+filefound[:-4]+"/tophat/accepted_hits.bam\n\tmkdir -p "+filefound[:-4]+"/cuffdiff_output/"+"";
    if options.nosge:
      print "\t"+cufflinks+"   -o "+filefound[:-4]+"/cuffdiff_output/ --GTF "+options.cufflinksgtf+" "+filefound[:-4]+"/tophat/accepted_hits.bam 2> "+filefound[:-4]+"/cuffdiff_output/error.log";
    else:
      print "\t"+runSGEandWait+" --param=\" -S /bin/bash -l \\\"h_vmem=5500M,virtual_free=5500M\\\" -V -R y -pe smp 1- \" --tmp="+tempdirname+"cufflinks_"+str(fileFoundIndex)+" -c \""+cufflinks+" -p \\$$NSLOTS  -o "+filefound[:-4]+"/cuffdiff_output/ --GTF "+options.cufflinksgtf+" "+filefound[:-4]+"/tophat/accepted_hits.bam 2> "+filefound[:-4]+"/cuffdiff_output/error.log\"";

print "\n"+"all:\t"+(" ".join(listOfTargetFiles))+"\n\n";
    
#move output to output dir
#os.removedirs(tempdirname);
