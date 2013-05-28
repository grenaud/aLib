library(hydroGOF)
library(shape)

args=(commandArgs(TRUE))

#############################
#   PRED vs OBS nucleotide  #
#############################


run<-args[1];
title<-args[2];

filetoopenA<-paste(run,".baseobspred3",sep="");
dataA<-read.table(filetoopenA);
xA<-dataA$V1
yA<-dataA$V2+1
zA<-dataA$V3

filetoopenC<-paste(run,".baseobspred4",sep="");
dataC<-read.table(filetoopenC);
xC<-dataC$V1
yC<-dataC$V2+1
zC<-dataC$V3

filetoopenG<-paste(run,".baseobspred5",sep="");
dataG<-read.table(filetoopenG);
xG<-dataG$V1
yG<-dataG$V2+1
zG<-dataG$V3


filetoopenT<-paste(run,".baseobspred6",sep="");
dataT<-read.table(filetoopenT);
xT<-dataT$V1
yT<-dataT$V2+1
zT<-dataT$V3

filetoopen<-paste(run,".baseobspred",sep="");

pdf(paste(filetoopen,".pdf",sep=""));
limit<-41
                                        #limit<-max(c(xA,xC,xG,xT))+5;

plot(xA,-10*log10(zA/(zA+yA)),xlim=c(0,limit),ylim=c(0,limit+6) ,main=title, xlab="Predicted score",ylab="Obsverved score",col=2,type="n");
lines(xA,-10*log10(zA/(zA+yA)),col=3,type="o",pch=20);
lines(xC,-10*log10(zC/(zC+yC)),col=4,type="o",pch=22);
lines(xG,-10*log10(zG/(zG+yG)),col=1,type="o",pch=23);
lines(xT,-10*log10(zT/(zT+yT)),col=2,type="o",pch=24);
valA<-(-10*log10(zA/(zA+yA)))
valC<-(-10*log10(zC/(zC+yC)))
valG<-(-10*log10(zG/(zG+yG)))
valT<-(-10*log10(zT/(zT+yT)))
rmseA<-rmse(xA[valA!=Inf],valA[valA!=Inf],na.rm=TRUE);
rmseC<-rmse(xC[valC!=Inf],valC[valC!=Inf],na.rm=TRUE);
rmseG<-rmse(xG[valG!=Inf],valG[valG!=Inf],na.rm=TRUE);
rmseT<-rmse(xT[valT!=Inf],valT[valT!=Inf],na.rm=TRUE);

                                        #legend(5,35,paste(sum(wrongloglike[okay]*z[okay],correctloglike[okay]*y[okay]),rmse(x[val!=Inf],val[val!=Inf],na.rm=TRUE),sep="\n")  )

legend(0, 41,  cex=0.8, c(paste("A",rmseA,sep=" "),paste("C",rmseC,sep=" "),paste("G",rmseG,sep=" "),paste("T",rmseT,sep=" ")),col=c(3,4,1,2),lty=c(1,1,1,1),   pch=c(20,22,23,24),  title="Predicted Nucleotide")

abline(0,1);
dev.off()

bwtouse=1;

ylimmax<-max(density(xA,weights=(yA+zA), bw = bwtouse)$y,
             density(xC,weights=(yC+zC), bw = bwtouse)$y,
             density(xG,weights=(yG+zG), bw = bwtouse)$y,
             density(xT,weights=(yT+zT), bw = bwtouse)$y);

pdf(paste(filetoopen,".dens.pdf",sep=""));
                                        #main="Distribution of predicted quality score",
plot(density(xA,weights=(yA+zA), bw = bwtouse),main=title,xlim=c(0,42),ylim=c(0,ylimmax),xlab="Predicted score",ylab="Distribution",col=2,type="n");
lines(density(xA,weights=(yA+zA), bw = bwtouse),col=3,type="l",pch=20);
lines(density(xC,weights=(yC+zC), bw = bwtouse),col=4,type="l",pch=20);
lines(density(xG,weights=(yG+zG), bw = bwtouse),col=1,type="l",pch=20);
lines(density(xT,weights=(yT+zT), bw = bwtouse),col=2,type="l",pch=20);
legend(0, ylimmax/2,  cex=0.8, c("A","C","G","T"),col=c(3,4,1,2),lty=c(1,1,1,1),   pch=c(20,22,23,24),  title="Predicted Nucleotide");
dev.off();





########################
#   PRED vs OBS cycle  #
########################
# run<-".";
for(file in c(paste(run,".baseobspredcycle1",sep=""),
              paste(run,".baseobspredcycle2",sep=""))){

  
  #data<-cbind(read.table( paste(dir,"/mappingr1/",file,sep="") ),read.table( paste(dir,"/mappingr2/",file,sep="") )[,-1] );
  data<-read.table(file);
  oddvals <- seq(3, ncol(data), by=2);
  evenvals <- seq(2, ncol(data), by=2);
  numberofcs<-(length(data[1,])-1)/2
  maxqual<-max(data$V1)
  minqual<-min(data$V1)
  steps<-c()
 
  steps<-c(6,6,6,6,6,5,5);
  minqual<-1;


  pdf(paste(file,"cycle.pdf",sep=""));
  plot(seq(1,numberofcs,by=1),-10*log(data[23,oddvals]/(data[23,oddvals]+data[23,evenvals]))/log(10),ylim=c(0,maxqual+5),type="n", main="Predicted vs observed quality score per cycle",xlab="Cycle",ylab="Obsverved score")   
  xpos<-5;
  color<-1
  colortouse=rainbow(length(steps),v=0.8)[color]
  step<-1

  steptouse<-steps[step];

  i<-1;

  while(step<=length(steps)){
    xvalues<-seq(1,numberofcs,by=1);
    numerator<-data[i:(i+steptouse-1),oddvals];
    denominator<-data[i:(i+steptouse-1),evenvals];
    valuesok<-numerator+denominator>200;
    numerator[!valuesok]<-0
    denominator[!valuesok]<-0
    yvalues<-apply(rbind(-10*log(numerator/(numerator+denominator))/log(10),(numerator+denominator)),
                   2,
                   function(x) weighted.mean(x[1:steptouse],x[(steptouse+1):(2*steptouse)],na.rm=TRUE) );
    sumvalues<-apply( (data[i:(i+steptouse-1),oddvals]+data[i:(i+steptouse-1),evenvals]),2, function(x) sum(x))
    okay <- !is.infinite(xvalues) & !is.infinite(yvalues) ;
    if(i!=1){
      lines( na.omit(data.frame(xvalues[okay],yvalues[okay]))   ,col=rainbow(length(steps),v=0.8)[color],type="l")
      Arrows(-0.5,i+minqual-1-0.5,
             -0.5,(i+minqual-1+steptouse-1)+0.5,
             col=rainbow(length(steps),v=0.8)[color],code = 3,lend=1,arr.length=0.1,);
                                        #rmsetoprint<- mean( ( mean(yvalues[okay])-yvalues[okay])^2);
                                        #text(xpos,mean(na.omit( yvalues[okay]) )+1,paste(paste(i+minqual-1,(i+minqual-1+steptouse-1),sep="-"),rmsetoprint,sep=" ") ,col=rainbow(length(steps),v=0.8)[color] );
      text(xpos,mean(na.omit( yvalues[okay]) )+1,paste(i+minqual-1,(i+minqual-1+steptouse-1),sep="-") ,col=rainbow(length(steps),v=0.8)[color] );
    }

    xpos<-xpos+5;
    color<-color+1;
    step<-step+1;
    i<-i+steptouse;
    steptouse<-steps[step];
  }
  
  dev.off();

}







