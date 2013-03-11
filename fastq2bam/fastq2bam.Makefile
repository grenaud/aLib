
CXX      = g++   
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB}  -c
LDFLAGS  = -lz


all: fastq2bam 

fastq2bam.o:	fastq2bam.cpp
	${CXX} ${CXXFLAGS} fastq2bam.cpp


fastq2bam:	fastq2bam.o ${LIBGAB}utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f fastq2bam.o fastq2bam

