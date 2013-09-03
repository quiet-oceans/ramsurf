RAM PARABOLIC EQUATION
======================

This source code is a translation in optimized C of the orginal ram code from
Dr. David C. Calvo. See `readme.orig` for original instructions.

The original code shipped five sources:

- `ram1.5.f`
- `ramclr.f`
- `rams0.5.f`
- `ramsurf1.5.f`
- `ramsurfclr2.0.f`

The aim of this project is to make these program run faster by using vector
instructions and multi-threading, while keeping the original beahavior. This
implies a move to the C language.

Building
--------

Simply run `make` from the top directory, eventually passing additionnal options such as::

    make FCFLAGS=-Ofast FC=gfortran

GNU make is assumed.

To turn on the vectorized C code, run::

    make CFLAGS="-O2 -march=native"

Testing
-------

Compare the output produced by the FORTRAN code and the C code using::

    make check
    
Usage
-----

The C version mimics the FORTRAN interface.
Refer to the original documentation in `readme.orig`.

