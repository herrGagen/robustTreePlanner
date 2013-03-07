all: rtp

.PHONY: rtp

rtp:
	make -C src
	cp src/rtp .
