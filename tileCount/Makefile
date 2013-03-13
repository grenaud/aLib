BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2



CPPFLAGS := -O3 -lm -I${BAMTOOLS}/include/  -L${BAMTOOLS}/lib/ -lz -lbamtools

all:	countClustersBAM count_clocs	

countClustersBAM:	countClustersBAM.cpp
	g++ $(CPPFLAGS) -o countClustersBAM countClustersBAM.cpp

count_clocs:	count_clocs.cpp
	g++ -O3 -lm -o count_clocs count_clocs.cpp


clean:
	rm -f countClustersBAM count_clocs
