BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB}
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}utils.o -lm -lz


all:	countClustersBAM count_clocs	

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

countClustersBAM:	countClustersBAM.o
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

count_clocs:	count_clocs.o
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)


clean:
	rm -f countClustersBAM count_clocs
