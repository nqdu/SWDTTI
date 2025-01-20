
#include <cmath>
#include <vector>
#include <unsupported/Eigen/Polynomials>

/**
 * @brief Laguerre polynomials coefs
 * 
 * @param p order of polynomails
 * @param coefs shape(p+1)  P_n(x) = sum_{i=0}^p coefs[i] * x^i
 */
static void
laguerre_poly(double *coefs,size_t p) {
    for(size_t k = 0; k <=p; k ++) {
        double s = 1.;
        for(size_t i = 0; i< k; i ++ ){
            s *= (p - i + 0.0) / std::pow((k * 1. - i),2);
        }
        coefs[k] = s * std::pow(-1,k);;
    }
}

/**
 * @brief Laguerre function L_{n}(x) = e^{-x/2} * l_n(x)
 * 
 * @param n degree of laguerre polynomail
 * @param x 
 */
double
laguerre_func(size_t n, double x) {
    std::vector<double> coefs(n+1);
    laguerre_poly(coefs.data(),n);
    double s = 0.;
    double x0 = 1.;
    for(size_t i = 0; i < n + 1; i ++) {
        s += x0 * coefs[i];
        x0 *= x;
    }
    s = s * std::exp(-x * 0.5);

    return s;
}

/**
 * @brief compute gauss radau laguerre knots and weights
 * 
 * @param xgrl 
 * @param wgrl 
 * @param length size of xgrl
 */
void gauss_radau_laguerre(double *xgrl,double *wgrl,size_t length) 
{

    // special length
    if (length == 20) {
        double  xx[20] = {0.,                  0.1836651730961082,  0.6168163821266257,
                1.3007994029666046,  2.239949399008818,   3.4404719428812083,
                4.910645879674723,   6.661160612058037,   8.70559717712797,
                11.061108683576073,  13.74939412009166,   16.798123005945243,
                20.243084963302437,  24.131567267602627,  28.527947971314145,
                33.52361761728886,   39.25629429078257,   45.95295064863841,
                54.047093024701944,  64.64971243781604};
       double ww[20] = {
            0.05               , 0.3087051288288945  ,0.5579967647079931,
            0.8106921971976067 , 1.068638425565268   ,1.333780906686964,
            1.608340646098338  , 1.8949411986015714  ,2.1967749086342163,
            2.5178414711577393 , 2.863306272240217   ,3.240062283392964,
            3.6576540837679574 , 4.129885111104641   ,4.677813187132266,
            5.335849599667416  , 6.165741520073857   ,7.294555062758996,
            9.049660483553179  ,12.76854717545136};

        for(size_t i = 0; i < length; i ++) {
            xgrl[i] = xx[i];
            wgrl[i] = ww[i];
        }

        return;
    }
    
    size_t n = length - 1;
    std::vector<double> coefs(n+2),coefs_deriv(n+1);
    laguerre_poly(coefs.data(),n+1); // l_{n+1}(x)

    // derivative l_{n+1}'(x)
    for(size_t i = 0; i < n+1 ;i ++) {
        coefs_deriv[i] =  (double)(i + 1) * coefs[i+1];
    }

    // construct a polynomial
    Eigen::PolynomialSolver<double,-1> sol;
    Eigen::Map<Eigen::VectorXd> poly(coefs_deriv.data(),n+1);
    sol.compute(poly);
    Eigen::VectorXd z = sol.roots().real();
    std::sort(z.begin(),z.end());

    xgrl[0] = 0.;
    for(size_t i = 0; i < n; i ++ ) {
        xgrl[i+1] = z[i];
    }
    for(size_t i = 0; i < n + 1; i ++) {
        wgrl[i] = 1. / (n + 1) / std::pow(laguerre_func(n,xgrl[i]),2);
    }
}