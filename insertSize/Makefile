
include ../makefile.defs

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/VCFparser/ -I${LIBGAB}/VCFparser/gzstream/  -I/home/gabriel_renaud/Software/tabix-0.2.6/ -I${BAMTOOLS}/include/ -c

LDLIBS   += ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/PutProgramInHeader.o  -lm -lz


all: insertSizeRG


%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@



insertSizeRG: insertSizeRG.o  ${LIBGAB}/utils.o
	${CXX} -o $@ $^ $(LDLIBS) 




clean :
	rm -f  insertSizeRG *.o
