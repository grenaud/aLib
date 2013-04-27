
CXX      = g++
BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2/
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB}   -I${BAMTOOLS}/include/ -c
LDLIBS   += ${BAMTOOLS}/lib/libbamtools.a -lm -lz

all: plotQualScores 


plotQualScores.o:	plotQualScores.cpp
	${CXX} ${CXXFLAGS} plotQualScores.cpp

plotQualScores:	plotQualScores.o ${LIBGAB}utils.o  
	${CXX} -o $@ $^ $(LDLIBS)

clean:
	rm -f plotQualScores
