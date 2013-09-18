SUBDIRS=. tests

bin_PROGRAMS=ram1.5  ramclr  rams0.5  ramsurf1.5  ramsurfclr2.0 ramsurf
lib_LIBRARIES=libramsurflib.a
dist_include_HEADERS=ramsurf.h

EXTRA_DIST=readme.orig README.rst

ram1_5_SOURCES=ram1.5.f
ramclr_SOURCES=ramclr.f
rams0_5_SOURCES=rams0.5.f
ramsurf1_5_SOURCES=ramsurf1.5.f
ramsurfclr2_0_SOURCES=ramsurfclr2.0.f

libramsurflib_a_SOURCES=ramsurflib.c
ramsurf_SOURCES=ramsurf.c 
ramsurf_LDADD=libramsurflib.a
