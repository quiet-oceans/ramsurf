check_PROGRAMS=ramcmp ramdump

ramcmp_SOURCES=ramcmp.f
ramdump_SOURCES=ramdump.f

TESTS=decreasing_slop.sh\
		  flat.sh\
		  flat_stability.sh\
		  deep_flat.sh\
		  flat_pade.sh\
		  flat_limit.sh\
		  increasing_slop.sh\
		  flat_2ssps.sh\
		  intensive.sh

EXTRA_DIST=decreasing_slop.test\
		  flat.test\
		  flat_stability.test\
		  deep_flat.test\
		  flat_pade.test\
		  flat_limit.test\
		  increasing_slop.test\
		  flat_2ssps.test\
		  intensive.test

SUFFIXES=.test .sh

.test.sh:$(srcdir)/Makefile.am
	echo '#!/bin/sh' > $@
	echo 'set -e' 	>> $@
	echo 'TARGET="`basename $<`"' >> $@
	echo 'test ! -d $$TARGET.d || chmod -R u+rwx $$TARGET.d' >> $@
	echo 'cp -r $(srcdir)/$$TARGET $$TARGET.d' >> $@
	echo 'chmod -R u+rwx $$TARGET.d' >> $@
	echo 'cd $$TARGET.d' >> $@
	echo '$(abs_top_builddir)/ramsurf1.5' >> $@
	echo 'mv tl.grid ref.grid' >> $@
	echo 'mv tl.line ref.line' >> $@
	echo '$(abs_top_builddir)/ramsurf' >> $@
	echo 'if ../ramcmp | grep -q T ; then false ; fi' >> $@
	chmod u+x $@

CLEANFILES=$(TESTS)

clean-local:
	rm -rf `echo $(EXTRA_DIST) | sed -e 's/test/test.d/g'`
