#!/usr/bin/python



import sys,os
import json


def chomp(s):
  return s.rstrip('\n');



if(len(sys.argv) != 2):
  sys.stderr.write("This python script puts the sequences in a json file without the sequences\nUsage: [json file in] > [json file out]");
  sys.exit(1);
  


#print sys.argv[1];

pathToConfig=sys.argv[0];
pathToConfig=os.path.dirname(os.path.abspath(pathToConfig))+"/config.json";


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




try:
  fileHandle = open ( sys.argv[1] );
except IOError:
  print "Cannot open json file "+sys.argv[1];
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


#print "TESTI "+str(jsondata);
#print "TEST2 "+str(jsondataConf);
#print "TEST2 "+str(jsondataConf["indices"]["p7indices"]["p7index"]);
#print "TEST2 "+str(jsondataConf["indices"]["p5indices"]["p5index"]);
p7indices2seq={};
p5indices2seq={};


for jsonrow in jsondataConf["indices"]["p7indices"]["p7index"]:
  if(jsonrow["id"] in p7indices2seq):
    sys.stderr.write("Found id twice "+str(jsonrow["id"])+" in "+str(pathToConfig));
    sys.exit(1); 
  else:
    p7indices2seq[ jsonrow["id"] ] = jsonrow["seq"];

for jsonrow in jsondataConf["indices"]["p5indices"]["p5index"]:
  if(jsonrow["id"] in p5indices2seq):
    sys.stderr.write("Found id twice "+str(jsonrow["id"])+" in "+str(pathToConfig));
    sys.exit(1); 
  else:
    p5indices2seq[ jsonrow["id"] ] = jsonrow["seq"];
    





#print "TEST "+str(type(jsondata["indicesraw"]));


indicesseq = [];
for jsonrow in jsondata["indicesraw"]:
  #print "TEST "+str(type(jsonrow));
  newrow = {};
  newrow["name"] = jsonrow["name"];


  #print jsonrow["name"];
  #print jsonrow["p7"];
  #if("p5" in jsonrow):
  if( jsonrow["name"] == "#index" ):
    newrow["p7"] = jsonrow["p7"];
    if("p5" in jsonrow):
      newrow["p5"] = jsonrow["p5"];
  else:
    if not(jsonrow["p7"] in p7indices2seq):
      sys.stderr.write("Could not find id : "+str(jsonrow["id"])+" in "+str(pathToConfig));
      sys.exit(1); 
    newrow["p7"] = p7indices2seq[jsonrow["p7"]];

    if("p5" in jsonrow):
      newrow["p5"] = p5indices2seq[jsonrow["p5"]];

  indicesseq.append(newrow);


jsondata["indicesseq"] = indicesseq;




#print "TESTF "+str(jsondata);
print json.dumps(jsondata);
