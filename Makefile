# this makefile requires GNU make
# pass options using FC and FFLAGS
# e.g. make FC=gfortran FFLAGS=-O2

FTARGETS=ram1.5 ramclr rams0.5 ramsurf1.5 ramsurfclr2.0 tests/ramcmp
CTARGETS=ramsurf
TARGETS=$(FTARGETS) $(CTARGETS)
SOURCES=$(patsubst %,%.f, $(FTARGETS)) $(patsubst %,%.c, $(CTARGETS))

override CFLAGS+= -std=c99
override LDFLAGS+= -fPIC
LDLIBS=-lm

EXTRA_DIST= ramsurf.h readme.orig README.rst tests Makefile ramsurflib.c
DIST_NAME=ramsurf

ramsurf:ramsurf.o ramsurflib.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

all:$(TARGETS)

clean:
	rm -rf *.o $(TARGETS) $(DIST_NAME) $(DIST_NAME).tar.gz
	find -name 'tl.*' -o -name ref.grid -delete

dist: clean Makefile $(SOURCES) $(EXTRA_DIST)
	mkdir $(DIST_NAME)
	cp -r $(SOURCES) $(EXTRA_DIST) $(DIST_NAME)/
	tar czf $(DIST_NAME).tar.gz $(DIST_NAME)
	rm -rf $(DIST_NAME)

distcheck: dist
	tar xzf $(DIST_NAME).tar.gz
	# test both normal and accelerated versions
	$(MAKE) -C $(DIST_NAME) check
	$(MAKE) -C $(DIST_NAME) clean
	$(MAKE) CFLAGS=-march=native -C $(DIST_NAME) check
	rm -rf $(DIST_NAME)

check: all
	for d in `find tests -mindepth 1 -type d` ; do \
		cd $$d ; \
		../../ramsurf1.5 ;\
		mv tl.grid ref.grid ;\
		../../ramsurf ;\
		printf "$$d " ; \
		if `../ramcmp | grep -q T` ; then \
			printf 'KO\n' ;\
		else\
			printf 'OK\n' ;\
		fi ;\
		cd - 2>&1 1>/dev/null ;\
	done
