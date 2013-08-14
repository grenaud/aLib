include ../makefile.defs



CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB}
LDLIBS+= ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/utils.o -lm -lz

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

all:	getCtrlReadsBAM	getCtrlReadsFASTQ

getCtrlReadsBAM:	getCtrlReadsBAM.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

getCtrlReadsFASTQ:	getCtrlReadsFASTQ.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f getCtrlReadsBAM getCtrlReadsFASTQ *.o
