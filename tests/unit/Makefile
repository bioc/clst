TOP=../..
RUNIT_DIR=$(shell pwd)
INST_PATH=$(shell dirname $(RUNIT_DIR))
PKG_PATH=$(shell dirname $(INST_PATH))
R=R $(R_ARGS)

all: install test

full: build install test

install:
	cd ${TOP}; \
	make -f devmakefile install; \
	cd ${RUNIT_DIR}

build:
	cd ${TOP}; \
	make build; \
	cd ${RUNIT_DIR}

test:
	cd ..; \
	${R} --slave < runalltests.R

clean:
	rm -rf ../unit_output; \
	cd ${RUNIT_DIR}; \
	rm -f *.pdf *.Rdata *.rda

%.R: FORCE
	cd ..; \
	./runfile.sh $@ ${TARGET}

FORCE:
