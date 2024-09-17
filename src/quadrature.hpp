#ifndef SWD_QUADRATURE
#define SWD_QUADRATURE
#include <cmath>

//GLL
void gauss_legendre_lobatto(double* knots, double* weights, int length);
void lagrange_poly(double xi,int nctrl,const double *xctrl,
                double *h,double*  hprime);

// GRL
void gauss_radau_laguerre(double *xgrl,double *wgrl,size_t length);
double laguerre_func(size_t n, double x);

#endif