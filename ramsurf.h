#ifndef RAMSURF_H
#define RAMSURF_H

#ifdef __cplusplus
extern "C" {
#endif

enum ERRCODE{
    PARSE_ERROR=1,
    //mz need to be increase
    INC_MZ,
    //mp need to be increase
    INC_MP,
    //mr need to be increase
    INC_MR,
    //laguerre is not converging
    LAG_NOT_CON
};

/* M_PI dropped in c99 */
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

int ramsurf(char const *input, char const *resultat, char const *output);

#ifdef __cplusplus
}
#endif

#endif
