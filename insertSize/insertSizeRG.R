#!/usr/bin/env Rscript


args=(commandArgs(TRUE))

data <- read.table(args[1],header=FALSE,comment.char="");
outdir<-args[2];

inserts<-c();
insertc<-c();
rg<-"";

for(i in seq(1,length(data[,1]))){
  #print(i);
  if(grepl("^#",data[i,1])){
    #print(rg);
    if(length(inserts) > 10){
   #   print(rg);
      pdf(paste(outdir,"rg_",rg,".pdf",sep=""));
      barplot(insertc,names.arg=c(inserts),main=paste("RG = ",rg),col="cyan");
      dev.off();
    }
    rg<-substring(as.character(data[i,1]), 2)
    inserts<-c();
    insertc<-c();
  }else{
    inserts<-c(inserts, as.numeric(as.character(data[i,1])));
    insertc<-c(insertc, as.numeric(as.character(data[i,2])));
    #print( as.numeric(as.character(data[i,1]) ) );
    #print( as.numeric(as.character(data[i,2]) ) );          
                                        #print(insertc);    
  }
}
#print(inserts);
#print(insertc);
if(length(inserts) > 10){
   pdf(paste(outdir,"rg_",rg,".pdf",sep=""));
   barplot(insertc,names.arg=c(inserts),main=paste("RG = ",rg),col="cyan");
   dev.off();
 }
