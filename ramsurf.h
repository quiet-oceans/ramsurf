#ifndef RAMSURF_H
#define RAMSURF_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct {
        size_t mr, mz, mp;
        float freq, zs, zr;
        float rmax, dr; int ndr;
        float zmax, dz; int ndz; float zmplt;
        float c0; int np, ns; float rs;
        float *rsrf, *zsrf;
        float *rb, *zb;
        float **cw, **attw, **cb, **rhob, **attn;
        float *rp;
    } ramsurf_t;


    int ramsurf(ramsurf_t const* rsurf, FILE *tl_line, FILE *tl_grid);

#ifdef __cplusplus
}
#endif

#endif
