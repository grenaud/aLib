BAMTOOLS= /mnt/solexa/bin/bamtools-1.0.2



CPPFLAGS := -O3 -lm -I${BAMTOOLS}/include/  -L${BAMTOOLS}/lib/ -lz -lbamtools


getCtrlReadsBAM:	plotQualScores.cpp
	g++ $(CPPFLAGS) -o plotQualScores plotQualScores.cpp 

clean:
	rm -f plotQualScores
