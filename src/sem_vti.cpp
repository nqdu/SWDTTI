#include "swdlayervti.hpp"

#include <iostream>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

typedef std::complex<double> dcmplx;

/**
 * @brief filter swd in the [ vmin,vmax] and Rek >= Imk
 * 
 * @param vmin/vmax min/max velocity in this region 
 * @param om angular frequency
 * @param displ_all eigen functions, shape(size,nc_all)
 * @param k eigenvalues, wavenumber
 * @param c filtered phase velocity, shape(nc)
 * @param displ filtered eigen function, shape(nc,size)
 */
void 
filter_swd(double vmin, double vmax,double om, const Eigen::MatrixXcd &displ_all,
            const Eigen::Array<std::complex<double>,-1,1> &k,std::vector<double> &c,
           std::vector<double> &displ)
{
    Eigen::Array<double,-1,1> c_all = (om / k).real();
    //std::cout << c_all.transpose() << "\n";

    // filter swd in [vmin * 0.85,vmax] region
    Eigen::Array<double,-1,1> c_filt;
    Eigen::MatrixXd displ_filt;
    using Eigen::all;
    //auto mask = ((c_all >= 0.85 * vmin) && (c_all <= vmax)) && (k.real().abs() >= k.imag().abs());
    auto mask = ((c_all >= vmin)&& (c_all <= vmax)) && (k.real().abs() >= k.imag().abs());
    // std::cout << mask << "\n";
    std::vector<int> idx0; idx0.reserve(mask.cast<int>().sum());
    for(int i = 0; i < c_all.size(); i ++) {
        if(mask[i]) {
            idx0.push_back(i);
        }
    }

    int nc = idx0.size();
    int size = displ_all.rows();
    c_filt.resize(nc); displ_filt.resize(size,nc);
    for(int i = 0; i < nc; i ++) {
        c_filt[i] = c_all[idx0[i]];
        displ_filt(all,i) = displ_all(all,idx0[i]).real();
    }
    //printf("%g %g\n",c_filt.minCoeff(),c_filt.maxCoeff());

    // sort according to ascending order 
    std::vector<int> idx;
    idx.resize(nc);
    for(int i = 0; i < nc; i ++ ) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(),
        [&c_filt](size_t i1, size_t i2) {return c_filt[i1] < c_filt[i2];}); 

    // copy to c/displ
    c.resize(nc); displ.resize(nc * size);
    for(int ic = 0; ic < nc; ic ++) {
        c[ic] = c_filt[idx[ic]];
        for(int i = 0; i < size; i ++) {
            displ[ic * size + i] = displ_filt(i,idx[ic]);
        }
    } 
}


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
    double omega2 = std::pow(om,2);
    dmat2 A = omega2 * dmat2(M.asDiagonal()) - E;

    // solve generalized eigenvalue problem A x = k^2 K x
    Eigen::GeneralizedEigenSolver<dmat2> sol;
    sol.compute(A,K);
    Eigen::MatrixXcd displ_all = sol.eigenvectors();
    Eigen::Array<dcmplx,-1,1> k = sol.eigenvalues().array().sqrt();

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