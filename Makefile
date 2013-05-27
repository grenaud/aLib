all: 
	make -C pipeline
	make -C bam2fastq
	make -C BCL2BAM
	make -C estimateControlReadsBCL
	make -C extractControlReadsBam
	make -C fastq2bam
	make -C plotQualScores
	make -C qualScoreC++
	make -C tileCount
	cd biohazard ; [ -d ~/.cabal ] || cabal update ; cabal configure ; cabal build
	make -C network-aware-bwa all

.PHONY: all
