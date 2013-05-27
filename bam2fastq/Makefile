include ../makefile.defs

CXX      = g++   
LIBGAB   = /home/gabriel_renaud/lib/
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}utils.o -lm -lz

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB} -I${BAMTOOLS}/include/ -I${LIBGAB}/VCFparser/gzstream/  -c
LDFLAGS  = -lz


all: bam2fastq 

bam2fastq.o:	bam2fastq.cpp
	${CXX} ${CXXFLAGS} bam2fastq.cpp


bam2fastq:	bam2fastq.o ${LIBGAB}utils.o  ${LIBGAB}VCFparser/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f bam2fastq.o bam2fastq

