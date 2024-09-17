#include "swdlayertti.hpp"
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>


/**
 * @brief compute phase velocity and eigen displacements for a given direction
 * 
 * @param freq current frequency
 * @param phi current direction, in deg
 * @param c phase velocity
 * @param displ displacement
 */
void LayerModelTTI::
compute_egnfun(double freq, double phi, std::vector<double> &c, std::vector<dcmplx> &displ) const 
{
    typedef Eigen::Matrix<dcmplx,-1,-1,1> dcmat;

    // angular frequency
    double om = 2. * M_PI * freq;
    double omega2 = std::pow(om,2);

    // mapping memory
    Eigen::Map<const Eigen::VectorXcd> M(Mmat.data(),nglob*3);
    Eigen::Map<const dcmat> E(Emat.data(),nglob*3,nglob*3);
    Eigen::Map<const dcmat> K1(K1mat.data(),nglob*3,nglob*3);
    Eigen::Map<const dcmat> K2(K2mat.data(),nglob*3,nglob*3);

    // initialize auxiliary matrices A,B
    dcmat A(nglob * 6,nglob * 6), B(nglob * 6,nglob * 6);
    A.setZero(); B.setZero();

    // construct auxiliary matrices A,B
    using Eigen::seq;
    auto idx1 = seq(0,nglob*3-1), idx2 = seq(nglob*3,nglob*6-1);
    A(idx1,idx2).setIdentity(); // A[:nglob*3,nglob*3:] = I;
    A(idx2,idx1) = omega2 * dcmat(M.asDiagonal()) - E; // A[nglob*3:,:nglob*3] = om^2 * M - E
    A(idx2,idx2) =  -K1;
    B(idx1,idx1).setIdentity();
    B(idx2,idx2) = K2.inverse(); // B = [[I 0],[0,M]]
    A = B * A;

    // get eigen values/eigen vector
    Eigen::ComplexEigenSolver<dcmat> sol(A);
    Eigen::MatrixXcd displ_all = sol.eigenvectors();;
    Eigen::Array<dcmplx,-1,1> k = sol.eigenvalues().array();

   // filter swd in [vmin * 0.85,vmax] region
    double vmin = PHASE_VELOC_MIN, vmax = PHASE_VELOC_MAX;
    Eigen::Array<double,-1,1> c_filt;
    Eigen::MatrixXcd displ_filt;
    Eigen::Array<double,-1,1> c_all = (om / k).real();
    auto mask = ((c_all >= vmin)&& (c_all <= vmax)) && (k.real().abs() >= k.imag().abs());
    std::vector<int> fidx; fidx.reserve(mask.cast<int>().sum()); // filterd indexes
    for(int i = 0; i < c_all.size(); i ++) {
        if(mask[i]) {
            fidx.push_back(i);
        }
    }

    // copy filtered eigenvalues/funcs to a temporary array
    using Eigen::all;
    int nc = fidx.size();
    c_filt.resize(nc); displ_filt.resize(nglob*3,nc);
    for(int i = 0; i < nc; i ++) {
        c_filt[i] = c_all[fidx[i]];
        displ_filt(all,i) = displ_all(idx1,fidx[i]);
    }

    // sort according to ascending order 
    std::vector<int> sidx;
    sidx.resize(nc);
    for(int i = 0; i < nc; i ++ ) sidx[i] = i;
    std::stable_sort(sidx.begin(), sidx.end(),
        [&c_filt](size_t i1, size_t i2) {return c_filt[i1] < c_filt[i2];}); 

    // copy to c/displ
    int size = nglob * 3;
    c.resize(nc); displ.resize(nc * size);
    for(int ic = 0; ic < nc; ic ++) {
        c[ic] = c_filt[sidx[ic]];
        for(int i = 0; i < size; i ++) {
            displ[ic * size + i] = displ_filt(i,sidx[ic]);
        }
    } 
}