BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB}
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}utils.o -lm -lz

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

all:	getCtrlReadsBAM	getCtrlReadsFASTQ

getCtrlReadsBAM:	getCtrlReadsBAM.o
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

getCtrlReadsFASTQ:	getCtrlReadsFASTQ.o
	g++ $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f getCtrlReadsBAM getCtrlReadsFASTQ *.o
