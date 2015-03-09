#!/usr/bin/python

from socket import *

import time;
import sys,os;
from optparse import OptionParser
from optparse import OptionGroup



def timeString():
    return str(time.strftime("%d_%b_%Y_%H:%M:%S", time.localtime()));

parser = OptionParser(usage="usage: %prog [options]");
parser.add_option("-n", "--name",     dest="name",    help="Name of Job");
parser.add_option("-i", "--in",       dest="infile",  help="Input BAM file");
parser.add_option("-o", "--out",      dest="outfile", help="Output BAM file");
parser.add_option("-g", "--genome",   dest="genome",  help="BWA indexed genome file");
parser.add_option("-a", "--ancient",  dest="ancient", help="Use ancient parameters",default=False,action="store_true");
parser.add_option("-e", "--email",    dest="email",   help="Your email, use commas for multiple ones");
parser.add_option("-w", "--wait",     dest="wait",    help="Wait for done signal from server",default=False,action="store_true");
parser.add_option("-f", "--force",    dest="force",   help="Ignore the existing output file if it exists",default=False,action="store_true");

if(len(sys.argv) == 1):
  parser.print_help()
  sys.exit(1)

(options, args) = parser.parse_args();

PORT_NUMBERCLI = 8088

if(not (options.name ) ):
    options.name = str(int(time.time()));



if(not (options.infile ) ):
  sys.stderr.write("Must specify input file");
  sys.exit(1)
    
if( not(os.path.isfile(options.infile) ) ):
    sys.stderr.write("ERROR: the input file "+str(options.infile)+" should exist");
    sys.exit(1);
    
if(not (options.outfile ) ):
  sys.stderr.write("Must specify output file");
  sys.exit(1)

if(not (options.force ) ):
    if( os.path.isfile(options.outfile) ) :
        sys.stderr.write("WARNING: the output file "+str(options.outfile)+" already exist, remove it for the program to send to the server");
        sys.exit(0); #let's exit nicely

if(not (options.genome ) ):
  sys.stderr.write("Must specify genome file");
  sys.exit(1);

if( not(os.path.isfile(options.genome+".ann") ) ):
    sys.stderr.write("ERROR: the BWA file "+str(options.genome+".ann")+" should exist");
    sys.exit(1);

if(not (options.email ) ):    
    options.email = "none";

cmd = str(options.wait)+"\t"+str(options.name)+"\t"+str(options.infile)+"\t"+str(options.outfile)+"\t"+str(options.genome)+"\t"+str(options.ancient)+"\t"+str(options.email)+"\t"+timeString();
print "launching "+str(cmd)+"\n";
#sys.exit(1)

host = 'localhost'
buf = 1024

addr = (host, PORT_NUMBERCLI);
    
clientsocket = socket(AF_INET, SOCK_STREAM);

clientsocket.connect(addr);

clientsocket.send(cmd);


if(options.wait):
    data = clientsocket.recv(buf);
    if  data:        
        #else:
        print data;
        clientsocket.close();
