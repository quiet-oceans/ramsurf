#ifndef RAMSURF_H
#define RAMSURF_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct {
        size_t mr; // max number of elements between rb and rsrf
        float freq, zs, zr; //frequency in Hz, source depth, receiver depth (data at this depth goes into tl.line)
        float rmax, dr; int ndr; // maximum range (rmax), range step (delta r or dr), ndr (writes output every ndr range steps)
        float zmax, dz; int ndz; float zmplt; // maximum depth of the whole computational domain (zmax), depth step (dz), ndz (output data at every ndz vertical point at a given range), vertical extent of the domain to output in tl.grid (zmplt).

        float c0; int np, ns; float rs; // reference sound speed, number of Pade terms (np), # of stability constraints (ns), radius of stability constraint (rs)
        float *rsrf, *zsrf; // surface depth: rsrf contains ranges and zsrf the depths. Both terminate by -1
        float *rb, *zb; // bathymetry: rb contains ranges , zb contains bathymetry values. Both teminate by -1
        // array of float arrays. Each subarray contains a profile, as an array of floats containing depth and informations alternatively.
        float **cw, **cb, **rhob, **attn; // sound speed in water, sound speed in bottom, density of bottom, attenuation of bottom, 
        float *rp; // new profile boundaries (use new profile when a value is reached)
    } ramsurf_t;


    int ramsurf(ramsurf_t const* rsurf, //ram parameters
                FILE *tl_line, // transmission loss versus range at a specified receiver depth
                FILE *tl_grid  // the entire range-depth transmission loss field (binary file)
                );

#ifdef __cplusplus
}
#endif

#endif
