library(Hmisc)
args=(commandArgs(TRUE))


data<-read.table(args[1],header=TRUE)
nucarray<-c("A","C","G","T");

maxq<-ncol(data)-3;
for(n in 1:4){  #for each nuc
x<-c();
y<-c();
yminus<-c();
yplus<-c();
ym1<-c();
yp1<-c();
ym5<-c();
yp5<-c();
ym10<-c();
yp10<-c();
ym25<-c();
yp25<-c();



for(c in 1:max(data$cycle)) { #for each cycle

quals<-seq(maxq+1)-1;

print(c);
arraytoconsider<-data[seq(n, nrow(data), by=4),][c,3:(maxq+3)];
quals<-quals[arraytoconsider!=0]
arraytoconsider<-arraytoconsider[arraytoconsider!=0]
q<-wtd.quantile(quals,weights=arraytoconsider,probs=c(0.01,0.05, .10, .25,0.75, .9, 0.95, 0.99));

x<-c(x,c);
y<-c(y,wtd.mean(quals,weights=arraytoconsider));
yminus<-c(yminus,min(quals[arraytoconsider>0]));
yplus <-c(yplus, max(quals[arraytoconsider>0]));
ym1<- c(ym1,as.numeric(q[1]));
ym5<- c(ym5,as.numeric(q[2]));
ym10<-c(ym10,as.numeric(q[3]));
ym25<-c(ym25,as.numeric(q[4]));
yp25<-c(yp25,as.numeric(q[5]));
yp10<-c(yp10,as.numeric(q[6]));
yp5<- c(yp5,as.numeric(q[7]));
yp1<- c(yp1,as.numeric(q[8]));
       
}


pdf(paste(paste(args[2],"/",sep=""),args[3],nucarray[n],".pdf",sep=""), width=16, height=9.6);

par(fg = "grey")
errbar(x, y,yplus,yminus,ylim=c(0,maxq),main=paste("Distribution of quality score for ",nucarray[n],sep=""),xlab="Cycle",ylab="Quality score (phred scale)",lwd=0.1)
par(fg = rgb(1,1,0))
#errbar(x, y, yp1,   ym1     ,ylim=c(0,maxq),add=TRUE,lwd=0.2)

errbar(x, y, yp5,   ym5     ,ylim=c(0,maxq),add=TRUE,lwd=0.3)
par(fg = rgb(0.9,0.9,0))
errbar(x, y, yp10,  ym10    ,ylim=c(0,maxq),add=TRUE,lwd=0.4)
par(fg = rgb(0.9,0.5,0))
#par(fg = rgb(1,0,0))
errbar(x, y, yp25,  ym25    ,ylim=c(0,maxq),add=TRUE,lwd=0.5)
par(fg = "black")
errbar(x, y,0,0,ylim=c(0,maxq),add=TRUE,lwd=0)


dev.off();


}
