CPPFLAGS := -lm -I/mnt/solexa/bin/bamtools-1.0.2/include/ -L/mnt/solexa/bin/bamtools-1.0.2/lib/ -lbamtools 
all:	getCtrlReadsBAM	getCtrlReadsFASTQ

getCtrlReadsBAM:	getCtrlReadsBAM.cpp
	g++ $(CPPFLAGS) -o getCtrlReadsBAM getCtrlReadsBAM.cpp

getCtrlReadsFASTQ:	getCtrlReadsFASTQ.cpp
	g++ $(CPPFLAGS) -o getCtrlReadsFASTQ getCtrlReadsFASTQ.cpp

clean:
	rm -f getCtrlReadsBAM getCtrlReadsFASTQ
