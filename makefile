all: rtp

.PHONY: rtp

rtp:
	make -C src rtp
	cp src/rtp .

clean:
	rm -f rtp
	make -C src clean