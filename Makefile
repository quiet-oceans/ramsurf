# this makefile requires GNU make
# pass options using FC and FFLAGS
# e.g. make FC=gfortran FFLAGS=-O2

TARGETS=ram1.5 ramclr rams0.5 ramsurf1.5 ramsurfclr2.0
SOURCES=$(patsubst %,%.f, $(TARGETS))

EXTRA_DIST= readme.orig README.rst tests
DIST_NAME=ram

all:$(TARGETS)

clean:
	rm -rf *.o $(TARGETS) $(DIST_NAME) $(DIST_NAME).tar.gz
	find -name 'tl.*' -delete

dist: clean Makefile $(SOURCES) $(EXTRA_DIST)
	mkdir $(DIST_NAME)
	cp -r $(SOURCES) $(EXTRA_DIST) $(DIST_NAME)/
	tar czf $(DIST_NAME).tar.gz $(DIST_NAME)

check: all
	for d in tests/flat* ; do cd $$d ; { ../../ramsurf1.5 && echo "$$d: OK" ; } || echo "$$d: KO" ; cd - 2>&1 1>/dev/null ; done
