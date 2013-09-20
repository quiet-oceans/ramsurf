/*
 * Copyright (c) 2006, Dr. David C. Calvo <calvo@ccs.nrl.navy.mil> Copyright
 * (c) 2012, Serge Guelton <serge.guelton@quiet-oceans.com>
 * 
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 *     Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *     Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *     Neither the name of the <ORGANIZATION> nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#include "ramsurf.h"

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <setjmp.h>

enum { BUFSIZE = 1024 };
static jmp_buf exception_env;
#ifndef max
   #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

enum { PARSE_ERROR=1 };

static
void read_dimensions( FILE* fs1, size_t *mr, size_t *mz, size_t *mp)
{

    float zr,zmax,zmplt;
    float dr, dz, rs;
    int ndr, ndz, np, ns;

    static char tmp[BUFSIZE];
    //
    fgets(tmp,BUFSIZE,fs1);
    fscanf(fs1,"%*f %*f %f",&zr);
    fgets(tmp,BUFSIZE,fs1);
    fscanf(fs1,"%*f %f %d",&dr,&ndr);
    fgets(tmp,BUFSIZE,fs1);
    fscanf(fs1,"%f %f %d %f",&zmax,&dz,&ndz,&zmplt);
    fgets(tmp,BUFSIZE,fs1);
    fscanf(fs1,"%*f %d %d %f",&np,&ns,&rs);
    fgets(tmp,BUFSIZE,fs1);
    //
    int i=0;
    float rsrf, zsrf;
    while(1) {
        if (fscanf(fs1,"%f %f",&rsrf,&zsrf)!=2)
    	    longjmp(exception_env,PARSE_ERROR);
        fgets(tmp,BUFSIZE,fs1);
        if( rsrf < 0.0) break;
        i=i+1;
    }

    //
    float rb, zb;
    int j=0;
    while(1) {
        if (fscanf(fs1,"%f %f",&rb,&zb)!=2)
            longjmp(exception_env,PARSE_ERROR);
        fgets(tmp,BUFSIZE,fs1);
        if( rb < 0.0) break;
        j=j+1;
    }

    //
    *mr = 1 + max(i,j);
    *mz = zmax/dz + 1.5;
    *mp = np;
    rewind(fs1);
}

static
void raw_read(float ***data, FILE* fs1, size_t n)
{
    static char tmp[BUFSIZE];
    *data = realloc(*data, sizeof(float*)*(n+1));
    size_t max_m = 2;
    float **curr = (*data) + n ;
    *curr = malloc(sizeof(float) * max_m);

    float profi,zi = 0.;
    for(size_t m=0; zi>=0. && !feof(fs1); m+=2) {
        if( m == max_m) {
            max_m *= 2;
            *curr = realloc(*curr, sizeof(float) * max_m);
        }
        fscanf(fs1,"%f %f",&zi,&profi);
        fgets(tmp,BUFSIZE,fs1);
        (*curr)[m] = zi ;
        (*curr)[m+1] = profi ;
    }
    if(feof(fs1)) {
        free(*curr);
        *curr = NULL;
    }
}

static
void rsurf_init(ramsurf_t *rsurf, FILE* fs1)
{
    /* grab dimensions */
    size_t mr, mz, mp;
    read_dimensions(fs1, &mr, &mz, &mp);

    rsurf->mr = mr;

    /* allocate memory */
    rsurf->rsrf = malloc(sizeof(float)*mr);
    rsurf->zsrf = malloc(sizeof(float)*mr);
    rsurf->rb = malloc(sizeof(float)*mr);
    rsurf->zb = malloc(sizeof(float)*mr);

    rsurf->cw = malloc(sizeof(float*));   *rsurf->cw = NULL;
    rsurf->cb = malloc(sizeof(float*));   *rsurf->cb = NULL;
    rsurf->rhob = malloc(sizeof(float*)); *rsurf->rhob = NULL;
    rsurf->attn = malloc(sizeof(float*)); *rsurf->attn = NULL;
    rsurf->rp = NULL;

    /* parse parameters */
    static char tmp[BUFSIZE];
    //
    fgets(tmp,BUFSIZE,fs1);
    fscanf(fs1,"%f %f %f",&rsurf->freq,&rsurf->zs,&rsurf->zr);
    fgets(tmp,BUFSIZE,fs1);
    fscanf(fs1,"%f %f %d",&rsurf->rmax,&rsurf->dr,&rsurf->ndr);
    fgets(tmp,BUFSIZE,fs1);
    fscanf(fs1,"%f %f %d %f",&rsurf->zmax,&rsurf->dz,&rsurf->ndz,&rsurf->zmplt);
    fgets(tmp,BUFSIZE,fs1);
    fscanf(fs1,"%f %d %d %f",&rsurf->c0,&rsurf->np,&rsurf->ns,&rsurf->rs);
    fgets(tmp,BUFSIZE,fs1);
    
    //
    int i=0;
    while(1) {
        if (fscanf(fs1,"%f %f", &rsurf->rsrf[i], &rsurf->zsrf[i])!=2)
    	    longjmp(exception_env,PARSE_ERROR);
        fgets(tmp,BUFSIZE,fs1);
        if( rsurf->rsrf[i] < 0.0) break;
        i=i+1;
    }

    //
    i=0;
    while(1) {
        if (fscanf(fs1,"%f %f",&rsurf->rb[i],&rsurf->zb[i])!=2)
            longjmp(exception_env,PARSE_ERROR);
        fgets(tmp,BUFSIZE,fs1);
        if( rsurf->rb[i] < 0.0) break;
        i=i+1;
    }

    // read profiles
    for(size_t step = 0; !feof(fs1); ++step) {
        raw_read(&rsurf->cw, fs1, step);
        raw_read(&rsurf->cb, fs1, step);
        raw_read(&rsurf->rhob, fs1, step);
        raw_read(&rsurf->attn, fs1, step);

        rsurf->rp = realloc(rsurf->rp, sizeof(float)*(1 + step));
        if(fscanf(fs1, "%f", rsurf->rp + step) != 1)
            rsurf->rp[step] = 2.0 * rsurf->rmax;
    }
}

static
void rsurf_del(ramsurf_t *rsurf)
{
    free(rsurf->rsrf);
    free(rsurf->zsrf);
    free(rsurf->rb);
    free(rsurf->zb);

#define freeall(f) \
    for(float **iter = f; *iter; ++iter)\
        free(*iter);\
    free(f);

    freeall(rsurf->cw);
    freeall(rsurf->cb);
    freeall(rsurf->rhob);
    freeall(rsurf->attn);

#undef freeall

    free(rsurf->rp);
}

//
//     This version of ram handles variable surface height, including
//     rough surfaces and propagation in beach environments [J. Acoust.
//     Soc. Am. 97, 2767-2770 (1995)]. The input file is similar to the 
//     input file for ram but contains a block of data for the location 
//     of the surface. The computational grid extends from z=0 to z=zmax. 
//     The pressure release surface is located at z=zsrf(r), which is 
//     linearly interpolated between input points just like the bathymetry. 
//     The inputs zsrf must be greater than or equal to zero. As in ramgeo, 
//     the layering in the sediment parallels the bathymetry. 
//
//     ******************************************************************
//     ***** Range-dependent Acoustic Model, Version 1.5b, 19-Oct-00 ****
//     ******************************************************************
//
int main(int argc __attribute__((unused)), char *argv[] __attribute__((unused))) {
    const char ram_in[] = "ram.in";
    const char tl_line[] = "tl.line";
    const char tl_grid[] = "tl.grid";
    FILE* fs1 = fopen(ram_in,"r");
    if(!fs1) {
        perror(ram_in);
        return 1;
    }
    FILE* fs2 = fopen(tl_line,"w+");
    if(!fs2) {
        fclose(fs1);
        perror(tl_line);
        return 1;
    }
    FILE* fs3 = fopen(tl_grid,"w+");
    if(!fs3) {
        fclose(fs1);
        fclose(fs2);
        perror(tl_grid);
        return 1;
    }

    int errorCode;
    if(!(errorCode=setjmp(exception_env)))
    {
        ramsurf_t rsurf;
        rsurf_init(&rsurf, fs1);
        ramsurf(&rsurf, fs2, fs3);
        rsurf_del(&rsurf);
    }

    fclose(fs1);
    fclose(fs2);
    fclose(fs3);
    return 0;
}
