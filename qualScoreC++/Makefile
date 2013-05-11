
include ../makefile.defs

CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB}
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}utils.o -lm -lz

all:	qualScoresObsVsPred

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

qualScoresObsVsPred:	qualScoresObsVsPred.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)


clean:
	rm -f qualScoresObsVsPred qualScoresObsVsPred.o
