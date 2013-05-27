#!/usr/bin/env Rscript



args=(commandArgs(TRUE))

toread = file(args[1], "rb")
d<-readBin(toread, double(), n=100000000,endian = "little")
d[d>0.01]<- 0.01;

pdf(args[2]);
plot(density(d),main="Distribution of the sequence likelihood")
dev.off()


