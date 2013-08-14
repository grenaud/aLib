include ../makefile.defs


CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB}
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz


all:	countClustersBAM count_clocs	

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

countClustersBAM:	countClustersBAM.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

count_clocs:	count_clocs.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)


clean:
	rm -f countClustersBAM count_clocs
