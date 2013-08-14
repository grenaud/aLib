#BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2
#LIBGAB   = /home/gabriel_renaud/lib/

include ../makefile.defs


CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I../SimpleJSON/src/
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz


all:    mergeTrimReadsBAM assignRG filterReads errorRatePerCycle generate_report 

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

assignRG: assignRG.o PrefixTree.o RGAssign.o ${LIBGAB}/PutProgramInHeader.o ../SimpleJSON/obj/JSON.o ../SimpleJSON/obj/JSONValue.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

mergeTrimReadsBAM: mergeTrimReadsBAM.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o 
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

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
	rm -f mergeTrimReadsBAM  assignRG filterReads errorRatePerCycle generate_report *.o

