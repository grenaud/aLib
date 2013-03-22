BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB}
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}utils.o -lm -lz

all:	qualScoresObsVsPred

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

qualScoresObsVsPred:	qualScoresObsVsPred.o
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)


clean:
	rm -f qualScoresObsVsPred qualScoresObsVsPred.o
