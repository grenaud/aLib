BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB}
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a -lm -lz


all:    mergeTrimReadsBAM assignRG filterReads errorRatePerCycle
#MergeTrimReads.o: 	MergeTrimReads.h
#mergeTrimReadsBAM.o: MergeTrimReads.h
%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

assignRG: assignRG.o PrefixTree.o RGAssign.o ${LIBGAB}/PutProgramInHeader.o 
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

mergeTrimReadsBAM: mergeTrimReadsBAM.o  MergeTrimReads.o ${LIBGAB}/PutProgramInHeader.o 
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

filterReads: filterReads.o FilterBAMal.o ${LIBGAB}/PutProgramInHeader.o 
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

errorRatePerCycle: errorRatePerCycle.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS)

clean:
	rm -f mergeTrimReadsBAM MergeTrimReads.o PrefixTree.o RGAssign.o assignRG filterReads errorRatePerCycle errorRatePerCycle.o

