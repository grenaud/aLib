all: 
	make -C SimpleJSON
	make -C libgab
	make -C pipeline
	make -C estimateControlReadsBCL
	make -C extractControlReadsBam
	make -C plotQualScores
	make -C BCL2BAM2FASTQ
	make -C qualScoreC++
	make -C tileCount
	make -C insertSize

	cd biohazard ; [ -d ~/.cabal ] || cabal update ; cabal install
	make -C network-aware-bwa all


clean:
	make -C SimpleJSON clean
	make -C libgab clean
	make -C pipeline clean
	make -C estimateControlReadsBCL clean
	make -C extractControlReadsBam clean
	make -C plotQualScores clean
	make -C qualScoreC++ clean
	make -C tileCount clean
	make -C BCL2BAM2FASTQ clean
	make -C insertSize clean

.PHONY: all
