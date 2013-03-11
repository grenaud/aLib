#!/usr/bin/env Rscript

library(gplots)
library(ggplot2)
library(reshape)


args=(commandArgs(TRUE))


toread <- read.table(args[1],header=TRUE);


pdf(args[2]);
ggplot(toread) + geom_line(aes(x=cycle, y=mismatches.),colour="red")+xlab("Sequencing Cycle #") + ggtitle("PhiX error rate")+ ylab("Error rate (%)") 
dev.off()
                

melttoread<-melt.data.frame(toread,id.vars="cycle",measure.vars=c(6,8,10,12,14,16,18,20,22,24,26,28));
melttoread$variable<-sub(".","",sub(".","->",melttoread$variable,fixed=TRUE),fixed=TRUE)              

pdf(args[3]);
ggplot(melttoread) + geom_line(aes(x=cycle, y=value,colour=variable))      +xlab("Sequencing Cycle #") + ggtitle("PhiX error rate")+ ylab("Error rate (%)") + scale_colour_discrete(name = "Type of error")           
dev.off()


