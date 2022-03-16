#include <stdio.h>

#ifdef CPP_GSL

#include <gsl/gsl_sf_laguerre.h>


void laguerre_(int *l, int *p, double *x, double *res, int *size)
{
    for(int i = 0; i < *size; i++){
        res[i] = gsl_sf_laguerre_n(*p, (double) *l, x[i]);
    }
    return;
}

void laguerre_scalar_(int *l, int *p, double *x, double *res)
{
    // printf("l: %d, p: %lf\n", *p, (double) *l);
    *res = gsl_sf_laguerre_n(*p, (double) *l, *x);
    return;
}

#else

void laguerre_(int *l, int *p, double *x, double *res, int *size)
{
	printf("Laguerre requires GSL library");
	exit();
}

void laguerre_scalar_(int *l, int *p, double *x, double *res)
{
    printf("Laguerre requires GSL library");
	exit();
}

#endif
