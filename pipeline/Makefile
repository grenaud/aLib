#BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2


include ../makefile.defs

CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I../SimpleJSON/src/
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz


all:     filterReads errorRatePerCycle generate_report 

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

filterReads: filterReads.o FilterBAMal.o ${LIBGAB}/PutProgramInHeader.o 
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

errorRatePerCycle: errorRatePerCycle.o   ${LIBGAB}/ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS)



generate_report: generate_report.hs
	cabal update
	cabal install xml
	cabal install blaze-html
	ghc --make generate_report.hs

clean:
	rm -f filterReads errorRatePerCycle generate_report *.o

