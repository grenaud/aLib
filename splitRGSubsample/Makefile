include ../makefile.defs


CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB} -I${LIBGAB}/gzstream/
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz


all:	splitByRG 

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

splitByRG:	splitByRG.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)


clean:
	rm -f splitByRG  splitByRG.o
