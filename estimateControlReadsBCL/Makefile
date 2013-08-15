
include ../makefile.defs


CXXFLAGS+= -O3 -Wall -I${BAMTOOLS}/include/ -I${LIBGAB}
LDLIBS+= -lm

all:	estimateControlReadsBCL	#printIndicesBCL

%.o: %.cpp
	${CXX} -c ${CXXFLAGS} $^ 

estimateControlReadsBCL:	estimateControlReadsBCL.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

printIndicesBCL:	printIndicesBCL.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)


clean:
	rm -f estimateControlReadsBCL  *.o #printIndicesBCL
