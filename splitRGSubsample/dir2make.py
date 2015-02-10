#!/usr/bin/python

import sys,os
import time
from optparse import OptionParser
from optparse import OptionGroup
import subprocess;


timestamp = int(round(time.time() * 1000))



parser = OptionParser(usage="Script to create a makefile given a directory listing\nusage: %prog [options]");
group = OptionGroup(parser, "General","General parameters");

group.add_option("-d",                dest="dir",          help="Directory containing the BAMs");
group.add_option("-c",                dest="cont",         help="Estimate contamination",  action="store_true", default=False);
group.add_option("--schmutzi",        dest="schmutzi",     help="command line for schmutzi");
group.add_option("--ref",             dest="reference",    help="reference used for mapping");
group.add_option("--bam2prof",        dest="bam2prof",     help="command line for bam2prof");
group.add_option("--deamProf2pdf",    dest="deamProf2pdf", help="command line for deamProf2pdf.R");
group.add_option("--length",          dest="length",       help="Length for bam2prof",      default=20);

parser.add_option_group(group);
(options, args) = parser.parse_args();


if (options.dir == None or not os.path.isdir(options.dir) ):
  print "Need to specify directory";
  sys.exit(1);

if ( options.cont  ):

  if ( options.schmutzi == None ):
    print "Need to specify schmutzi command";
    sys.exit(1);

if ( options.cont  ):
 if ( options.reference == None ):
    print "Need to specify reference used for mapping";
    sys.exit(1);

if ( options.bam2prof == None ):
  print "Need to specify bam2prof command";
  sys.exit(1);

if ( options.deamProf2pdf == None ):
  print "Need to specify deamProf2pdf command";
  sys.exit(1);


print "SHELL := /bin/bash\n\nDefault:\tall\n\n";  

listOfTargetFiles=[];

for filefound in os.listdir(options.dir):
  if(filefound[-3:] == "bam"):
    sampleid=filefound[(filefound.find(".")+1):-4];

    if(sampleid.startswith("control")):
      continue;

    listOfTargetFiles.append(options.dir+filefound+".5p.prof");  
    print "\n"+options.dir+filefound+".5p.prof:\n\t"+options.bam2prof+" -length "+str(options.length)+" -5p "+options.dir+filefound+".5p.prof "+" -3p "+options.dir+filefound+".3p.prof "+options.dir+filefound;

    listOfTargetFiles.append(options.dir+filefound+".5p.prof.pdf");
    print "\n"+options.dir+filefound+".5p.prof.pdf: "+options.dir+filefound+".5p.prof"+"\n\t"+options.deamProf2pdf+"  "+options.dir+filefound+".5p.prof "+options.dir+filefound+".5p.prof.pdf \"5p deamination patterns for "+str(sampleid)+"\" ";

    listOfTargetFiles.append(options.dir+filefound+".3p.prof.pdf");
    print "\n"+options.dir+filefound+".3p.prof.pdf: "+options.dir+filefound+".5p.prof"+"\n\t"+options.deamProf2pdf+"  "+options.dir+filefound+".3p.prof "+options.dir+filefound+".3p.prof.pdf \"3p deamination patterns for "+str(sampleid)+"\" ";

    if ( options.cont  ):
      listOfTargetFiles.append(options.dir+filefound[:-4]+".cont.est");
      print "\n"+options.dir+filefound[:-4]+".cont.est:\n\t"+options.schmutzi+"  --library single "+" --out "+options.dir+filefound[:-4]+"  --ref "+options.reference+"  --title \"Posterior on contamination using deamination for "+sampleid+"\" "+options.dir+filefound;


print "\n\n"+"all:\t"+(" ".join(listOfTargetFiles))+"\n\t"+"touch "+options.dir+"finished\n\n";
print "\n"+"clean:\n"+"\trm -fv "+(" ".join(listOfTargetFiles))+"\n\n";

#  makeWrite[int(lanetopredict)].write((" ".join(listOfFilesBasecall[lanetopredict]))+" "+
#  makeWrite[int(lanetopredict)].write((" ".+"\n\n");
  

#make touch finish when all are done



