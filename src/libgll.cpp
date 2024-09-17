/**
 * 
 * Created by Nanqiao Du and Bin He on 11/26/2021
 * Copyright @ 2021 Nanqiao Du and Bin He. All rights reserved.
 */

/*
 *  Adapted from Fortran Standard Library :
 *  https://github.com/fortran-lang/stdlib/blob/master/src/stdlib_quadrature_gauss.f90
 */
#include <limits>
#include <cmath>

/**
 * @brief compute the Lagrange interpolants based upon the interpolation points   
 * 
 * @param xi current location,
 * @param nctrl no. of control points 
 * @param xctrl control nodes, shape (nctrl) 
 * @param h polynomial value, shape (nctrl)
 * @param hprime derivative, shape (nctrl)
 */
void lagrange_poly(double xi,int nctrl,const double *xctrl,
                double * h,double*  hprime)
{
    //! note: this routine is hit pretty hard by the mesher, optimizing the loops here will be beneficial
    for(int dgr = 0;dgr<nctrl;dgr++){
        double prod1 = 1., prod2 = 1.;

        // lagrangian interpolants
        double x0 = xctrl[dgr];
        for(int i=0;i< nctrl;i++){
            if(i != dgr){
                double x = xctrl[i];
                prod1 = prod1*(xi-x);
                prod2 = prod2*(x0-x);
            }
        }

        //! takes inverse to avoid additional divisions
        //! (multiplications are cheaper than divisions)
        double prod2_inv = 1. / prod2;
        h[dgr] = prod1 * prod2_inv;

        // first derivatives
        double sum = 0.0;
        for(int i=0;i<nctrl;i++){
            if (i != dgr){
                double prod3 = 1.0;
                for(int j=0;j<nctrl;j++){
                    if (j != dgr && j != i) prod3 = prod3*(xi-xctrl[j]);
                }
                sum = sum + prod3;
            }
        }

        hprime[dgr] = sum * prod2_inv;

    }

}

/**
 * @brief compute the derivative of Legendre Polynomial of order n
 * 
 * @param n order 
 * @param x current location
 */
static double dlegendre(int n,double x){
    double dleg{};
    if(n == 0){
        dleg = 0.;
    }
    else if (n==1){
        dleg = 1.;
    }
    else{
        double leg_down1 = x, leg_down2 = 1., leg;
        double dleg_down1 = 1., dleg_down2 = 0.;
        for(int i=2;i<=n;i++){
            leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i;
            dleg = dleg_down2 + (2*i-1)*leg_down1;
            leg_down2 = leg_down1;
            leg_down1 = leg;
            dleg_down2 = dleg_down1;
            dleg_down1 = dleg;
        }
    }

    return dleg;
}

/**
 * @brief compute the value of Legendre Polynomial of order n
 * 
 * @param n order 
 * @param x current location
 */
static double legendre(int n,double x){
    double leg{};
    if(n == 0){
        leg = 1;
    }
    else if(n==1){
        leg = x;
    } 
    else{
        double leg_down1 = x, leg_down2 = 1.;
        for(int i=2;i<=n;i++){
            leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i;
            leg_down2 = leg_down1;
            leg_down1 = leg;
        }
    }

    return leg;
}


/**
 * @brief compute Gauss-Legendre-Lobatto Nodes and Weights, Newton's iteration is used
 * 
 * @param knots nodes, shape (length)
 * @param weights weights shape (length)
 * @param length length of the array or the order+1 of the polynomial
 */
void gauss_legendre_lobatto(double* knots, double* weights, int length)
{
    // special length
    if(length == 5) {
        const double v = std::sqrt(21.) / 7.;
        double xgll[5] = {-1.,-v,0.,v,1.};
        double wgll[5] = {0.1,49./90.,32./45.,49./90.,0.1};
        for(int i = 0; i < 5; i ++) {
            knots[i] = xgll[i];
            weights[i] = wgll[i]; 
        }
        return;
    }

    // define some constants
    const double tolerance = 4.0 * std::numeric_limits<double>::epsilon();
    const int nnewton_iter = 100;
    const double pi = std::atan(1.0) * 4.0;

    // local arrays
    double *x = new double[length], *w = new double[length];

    int n = length - 1; // order of polynomial
    if(n == 1){
        x[0] = -1.; x[1] = 1.;
       w[0] = 1.; w[1] = 1.;
    }
    else{
        double leg = 0.,dleg = 0.,delta = 0.;
        //set end points
        x[0]   = -1.0; x[n] = 1.;
        w[0]   =  2./(n*(n+1.)); w[n] =  2./(n*(n+1.));
        for(int i=1;i<=(n+1)/2-1;i++){
            // initial guess from an approximate form given by SV Parter (1999)
            x[i] = -std::cos( (i+0.25)*pi/n  - 3/(8*n*pi*(i+0.25)));

            // newton iteration
            for(int j=0;j<nnewton_iter;j++){
                leg = legendre(n+1,x[i]) - legendre(n-1,x[i]);
                dleg = dlegendre(n+1,x[i]) - dlegendre(n-1,x[i]);
                delta = -leg/dleg;
                x[i] += delta;
                if ( std::abs(delta) <= tolerance * std::abs(x[i]) ) break;
            }
            x[n-i] = - x[i];
            leg = legendre(n, x[i]);
            w[i] = 2./(n*(n+1.)*leg*leg);
            w[n-i] = w[i];
        }

        if (n %2 == 0) {
            x[n/2] = 0.;
            leg = legendre(n, 0.0);
            w[n/2]  = 2./(n*(n+1.)*leg*leg); 
        }
    }

    // copy to knots/weights
    for(int i = 0; i < length; i ++) {
        knots[i] = (double) x[i];
        weights[i] = (double) w[i];
    }

    // free space
    delete[] w; delete[] x;
}