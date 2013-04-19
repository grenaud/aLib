all: 
	make -C network-aware-bwa all
	make -C pipeline
	cd biohazard ; [ -d ~/.cabal ] || cabal update ; cabal configure ; cabal build


.PHONY: all
