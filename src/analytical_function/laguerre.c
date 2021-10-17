#include <stdio.h>
#include <gsl/gsl_sf_laguerre.h>


void laguerre_(int *l, int *p, double *x, double *res, int *size)
{
    for(int i = 0; i < *size; i++){
        res[i] = gsl_sf_laguerre_n(*l, (double) *p, x[i]);
    }
    return;
}

void laguerre_scalar_(int *l, int *p, double *x, double *res)
{
    *res = gsl_sf_laguerre_n(*l, (double) *p, *x);
    return;
}
