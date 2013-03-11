CPPFLAGS :=  -I/mnt/solexa/bin/bamtools-1.0.2/include/ -I/home/gabriel_renaud/lib/
LDFLAGS := -lm -I/mnt/solexa/bin/bamtools-1.0.2/include/ -L/mnt/solexa/bin/bamtools-1.0.2/lib/ -lbamtools 

all:	qualScoresObsVsPred


qualScoresObsVsPred.o:	qualScoresObsVsPred.cpp
	g++ $(CPPFLAGS) -c  qualScoresObsVsPred.cpp 

qualScoresObsVsPred:	qualScoresObsVsPred.o
	g++ $(LDFLAGS) -o qualScoresObsVsPred qualScoresObsVsPred.o /home/gabriel_renaud/lib/utils.o

clean:
	rm -f qualScoresObsVsPred qualScoresObsVsPred.o
