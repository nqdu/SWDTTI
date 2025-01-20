#include "vti/swdlayervti.hpp"

#include <iostream>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

typedef std::complex<double> dcmplx;

void 
filter_swd(double vmin, double vmax,double om, const Eigen::MatrixXcd &displ_all,
            const Eigen::Array<std::complex<double>,-1,1> &k,std::vector<double> &c,
           std::vector<double> &displ);


/**
 * @brief compute Love wave dispersion and eigenfunctions
 * 
 * @param freq current frequency
 * @param vmin,vmax min/max velocity for your model 
 * @param c dispersion, shape(nc)
 * @param displ eigen functions(displ at y direction), shape(nc,nglob)
 */
void  LayerModelVTI::
compute_slegn(double freq,std::vector<double> &c,
              std::vector<double> &displ) const
{
    // map M/E/K matrices
    typedef Eigen::Matrix<double,-1,-1,1> dmat2;
    Eigen::Map<const Eigen::VectorXd> M(Mmat.data(),nglob);
    Eigen::Map<const Eigen::VectorXd> K(Kmat.data(),nglob);
    Eigen::Map<const dmat2> E(Emat.data(),nglob,nglob);
    
    // construct A  = K^{-1}(om^2 M - E)
    double om = 2. * M_PI * freq;
    double omega2 = std::pow(om,2);
    Eigen::MatrixXd A = dmat2((1. / K.array()).matrix().asDiagonal()) * (-E + dmat2(M.asDiagonal()) * omega2);

    // get eigen values/eigen vector
    Eigen::EigenSolver<Eigen::MatrixXd> sol(A);
    Eigen::MatrixXcd displ_all = sol.eigenvectors();
    Eigen::Array<dcmplx,-1,1> k = sol.eigenvalues().array().sqrt();

    // filter swd
    filter_swd(PHASE_VELOC_MIN,PHASE_VELOC_MAX,om,displ_all,k,c,displ);
}

/**
 * @brief compute Love wave dispersion and eigenfunctions
 * 
 * @param freq current frequency
 * @param vmin,vmax min/max velocity for your model 
 * @param c dispersion, shape(nc)
 * @param displ eigen functions(displ at x/z direction), shape(nc,2,nglob)
 */
void LayerModelVTI:: 
compute_sregn(double freq,std::vector<double> &c,
                std::vector<double> &displ) const 
{
    // map M/E/K matrices
    typedef Eigen::Matrix<double,-1,-1,1> dmat2;
    Eigen::Map<const Eigen::VectorXd> M(Mmat.data(),nglob * 2);
    Eigen::Map<const dmat2> E(Emat.data(),nglob * 2,nglob * 2);
    Eigen::Map<const dmat2> K(Kmat.data(),nglob * 2,nglob * 2);

    // prepare matrix A = om^2 M -E
    double om = 2. * M_PI * freq;
    double omega2 = om * om;
    Eigen::MatrixXf A = (omega2 * dmat2(M.asDiagonal()) - E).cast<float>();
    Eigen::MatrixXf K1 = K.cast<float>().inverse();
    A = K1 * A; 

    // solve generalized eigenvalue problem A x = k^2 K x
    //Eigen::GeneralizedEigenSolver<Eigen::MatrixXf> sol;
    Eigen::EigenSolver<Eigen::MatrixXf> sol;
    sol.compute(A);
    Eigen::MatrixXcd displ_all = sol.eigenvectors().cast<dcmplx>();
    Eigen::Array<dcmplx,-1,1> k = sol.eigenvalues().array().sqrt().cast<dcmplx>();

    // filter swd
    filter_swd(PHASE_VELOC_MIN,PHASE_VELOC_MAX,om,displ_all,k,c,displ);

    // note we obtained is [U,V^{\bar}]
    int nc = c.size();
    for(int ic = 0; ic < nc; ic ++) {
        double k0 = om / c[ic];
        for(int iglob = nglob; iglob < nglob * 2; iglob ++) {
            displ[ic * 2 * nglob + iglob] /= k0;
        }
    }
}