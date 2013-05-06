#!/usr/bin/python

import sys,os
import xml.etree.ElementTree as ET
import shutil

from datetime import date
from optparse import OptionParser
from optparse import OptionGroup

def isInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

#define #days for deleting
deleteALLafterdays    = 3650;
deleteBCLCIFafterdays = 180;



alibdir="";
illuminareaddir="";
todaysdate=date.today();


#BEGIN READ CONFIG FILE
pathToConfig=sys.argv[0];
pathToConfig=os.path.dirname(os.path.abspath(pathToConfig))+"/config.xml";
#alibdir=os.path.dirname(os.path.abspath( sys.argv[0]+"/../"))

try:
  tree =ET.parse(pathToConfig);
except IOError:
  print "Cannot open XML config file "+pathToConfig;
  sys.exit(1);

XMLconfig        = tree.getroot()



illuminareaddir  = XMLconfig.find("illuminareaddir").text;



parser = OptionParser(usage="usage: %prog [options]\nThis program will delete the runs that are more than ");
parser.add_option("--mock",  dest="mock", help="Do not delete anything, just print what you would delete",action="store_true",default=False);
parser.add_option("--list",  dest="list", help="Only delete the runs on this list");

#group = OptionGroup(parser, "Group1","First Group of parameters");
#group.add_option("-o", "--outfile", dest="outfile", help="output");

#if(len(sys.argv) == 1):
#  parser.print_help()
#  sys.exit(1)


#parser.add_option_group(group);
(options, args) = parser.parse_args();

allowedRuns = {};

if(options.list != None):
    fileHandle = open ( options.list );
    while 1:
        line = fileHandle.readline();
        if(not(line)):
            break
        line = line.rstrip('\n');
        allowedRuns[line]=1;
    fileHandle.close();


#print options.mock;

def delDirectory(dirtodel):   
    if(options.mock):
        print "deleting "+dirtodel;
    else:
        print "deleting "+dirtodel;
        #os.removedirs(dirtodel);
        shutil.rmtree(dirtodel);

for filename in os.listdir(illuminareaddir):
    originaldirname=filename;
    filename=illuminareaddir+"/"+filename;
    if(os.path.isdir(filename)):
      if("_" in filename):
        firstfield=filename.split("_")[0];
        firstfield=firstfield.split("/")[-1];
        if(len(firstfield) == 6    and
           isInt(firstfield[0:2]) and
           isInt(firstfield[2:4]) and
           isInt(firstfield[4:6]) ):        
          year  = int(firstfield[0:2])+2000; #will break in 2019 !
          month = int(firstfield[2:4]);
          day   = int(firstfield[4:6]);
          if(month>0 and month<13 and day>0 and day<32):#it's a run
            dateofrun = date(year,month,day);


            if(options.list != None):
                if(not(originaldirname in allowedRuns)):
                    print "keeping not in list "+filename;
                    continue;

            datediff=todaysdate-dateofrun;

            if(datediff.days>deleteALLafterdays):
              delDirectory(filename+"/");
            else:
              if(datediff.days>deleteBCLCIFafterdays):
                  if(os.path.isdir(filename+"/Data/Intensities/")):
                      for todel in os.listdir(filename+"/Data/Intensities/"):
                          #if(todel      == "BaseCalls" ):
                          #    delDirectory(filename+"/Data/Intensities/"+todel+"/");
                          if(todel[0:3] == "L00"):
                              delDirectory(filename+"/Data/Intensities/"+todel+"/");
                                  
              else:
                print "keeping "+filename;
