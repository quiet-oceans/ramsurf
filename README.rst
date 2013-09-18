RAM PARABOLIC EQUATION
======================

This source code is a translation in optimized C of the original ram code from
Dr. David C. Calvo. See `readme.orig` for original instructions.

The original code shipped five sources:

- `ram1.5.f`
- `ramclr.f`
- `rams0.5.f`
- `ramsurf1.5.f`
- `ramsurfclr2.0.f`

The aim of this project is to make these program run faster by using several
optimization techniques, while keeping the original behavior. This implies a
move to the C language.

Building
--------

Simply run::

    ./configure
    make

You can install all the ramsurf variant as well as the headers using::

    sudo make install

To turn on the vectorized C code, run::


    make CFLAGS="-O2 -march=native"

All the traditional autotool flags work, including not building in the sources.

Testing
-------

Compare the output produced by the FORTRAN code and the C code using::

    make check
    
Using the Ramsurf Program
--------------------------

The C version mimics the FORTRAN interface and assumes all inputs are described
in `ram.in`. It outputs two files: `lt.line` as plain text and `tl.grid` using
the binary format used by FORTRAN.

The input format is *the same* as the one described in `readme.orig`::
    
    # sample config. The first line is a header and is skipped. 
    freq zs zr          #frequency in Hz, source depth, receiver depth (data at this depth goes into tl.line)
    rmax dr ndr         #maximum range (rmax), range step (delta r or dr), ndr (writes output every ndr range steps)
    zmax dz ndz zmplt   #maximum depth of the whole computational domain (zmax), depth step (dz), ndz (output data at every ndz vertical point at a given range), vertical extent of the domain to output in tl.grid (zmplt).
    c0 np ns rs         # reference sound speed, number of Pade terms (np), # of stability constraints (ns), radius of stability constraint (rs)
    surface_range0 surface_depth0 # surface modeling
    ...
    surface_rangeN surface_depthN 
    -1 -1               # END of SECTION marker
    bathymetry_range0 bathymetry_depth0 #bathymetry modeling
    ...
    bathymetry_rangeM bathymetry_depthM
    -1 -1
    ssp_depth0 ssp_value0       # sound speed profile in water
    ...
    -1 -1
    bssp_depth0 bssp_value0    # sound speed profile in bottom
    ...
    -1 -1
    bd_depth0 bd_value0    # density of bottom
    ...
    -1 -1
    ba_depth0 ba_value0    # attenuation of bottom
    ...
    -1 -1
    optional_range # if given, starts a new ssp description starting from this range, otherwise the previous values are valid for all remaining ranges
    [...]


Using the Ramsurf Library
-------------------------

The C version can be used as a C library. Link with `-lramsurflib` and include
`ramsurf.h`. The interface documentation is contained in the header file.
