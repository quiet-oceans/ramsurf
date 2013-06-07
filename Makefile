TARGETS=ram1.5 ramclr rams0.5 ramsurf1.5 ramsurfclr2.0
SOURCES=$(patsubst %,%.f, $(TARGETS))
EXTRA_DIST= readme.orig README.rst
NAME=ram

all:$(TARGETS)

clean:
	rm -rf *.o $(TARGETS) $(NAME) $(NAME).tar.gz

dist: Makefile $(SOURCES) $(EXTRA_DIST)
	rm -rf $(NAME)
	mkdir $(NAME)
	cp $(SOURCES) $(EXTRA_DIST) $(NAME)/
	tar czf $(NAME).tar.gz $(NAME)
	rm -rf $(NAME)
