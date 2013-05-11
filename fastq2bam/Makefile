include ../makefile.defs


CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}VCFparser/ -I${LIBGAB}VCFparser/gzstream/ -I${BAMTOOLS}/include/ -c
LDFLAGS  = ${BAMTOOLS}/lib/libbamtools.a -lz


all: fastq2bam 

fastq2bam.o:	fastq2bam.cpp
	${CXX} ${CXXFLAGS} fastq2bam.cpp


fastq2bam:	fastq2bam.o ${LIBGAB}utils.o ${LIBGAB}PutProgramInHeader.o ${LIBGAB}/VCFparser/gzstream/gzstream.o ${LIBGAB}VCFparser/FastQObj.o  ${LIBGAB}VCFparser/FastQParser.o
	${CXX} -o $@ $^ $(LDFLAGS)

clean :
	rm -f fastq2bam.o fastq2bam

