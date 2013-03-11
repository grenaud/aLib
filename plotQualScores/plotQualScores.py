#!/usr/bin/python

import sys,os
import time
from optparse import OptionParser
from optparse import OptionGroup
import subprocess;


timestamp = int(round(time.time() * 1000))



parser = OptionParser(usage="Script to create pdfs to represent the quality scores for each nucleotide and each cycle\nusage: %prog [options]");
group = OptionGroup(parser, "General","General parameters");
group.add_option("-i", "--infile", dest="infile", help="Input file in BAM format");
group.add_option("-o", "--outdir", dest="outdir", help="Output directory");
group.add_option("--outprefix", dest="outprefix", help="Prefix for output files");
group.add_option("-t", "--tmpdir", dest="tmpdir", help="Temp directory",default="/tmp/");
parser.add_option_group(group);
(options, args) = parser.parse_args();


if (options.infile == None or not os.path.isfile(options.infile) ):
  print "Need to specify input file";
  sys.exit(1);
if (options.outdir == None or not os.path.isdir(options.outdir) ):
  print "Need to specify out directory";
  sys.exit(1);

fileprefix=options.outprefix;
if (options.outprefix == None ):
    fileprefix=os.path.splitext( options.infile )[0].split("/")[-1];


filenametemp = "/"+options.tmpdir+"/%s.qdat" % timestamp;
#print filename;

cmd =str(sys.path[0])+"/plotQualScores "+str(options.infile)+" > "+str(filenametemp);
proc = subprocess.Popen(cmd, shell=True);
proc.wait();
proc.poll();
if(proc.returncode != 0):
  print "Wrong return code "+str(proc.returncode)+" for command "+str(cmd);
  sys.exit(1);


cmd ="R CMD BATCH --no-save --no-restore '--args "+str(filenametemp)+" "+str(options.outdir)+" "+str(fileprefix)+"' "+str(sys.path[0])+"/qualScore2graph.R /dev/null"

proc = subprocess.Popen(cmd, shell=True);
proc.wait();
proc.poll();
if(proc.returncode != 0):
  print "Wrong return code "+str(proc.returncode)+" for command "+str(cmd);
  sys.exit(1);



os.remove(filenametemp);
