#include "multiphysics/vti_acoustic.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
typedef std::complex<double> dcmplx;

/**
 * @brief prepare M/K/E matrices
 * 
 */
void LayerModelMultiPhyVTI:: 
prepare_matrices(double freq)
{
    // allocate space and set zero
    int ng = nglob_ac + nglob_el * 2;
    Mmat.resize(ng);
    Emat.resize(ng * ng);
    Kmat.resize(ng * ng);
    std::fill(Mmat.begin(),Mmat.end(),0.);
    std::fill(Emat.begin(),Emat.end(),0.);
    std::fill(Kmat.begin(),Kmat.end(),0.);

    // compute M/K/E for gll/grl layer, elastic
    std::array<double,NGRL> sumC,sumL;
    for(int ispec = 0; ispec < nspec_el + nspec_el_grl; ispec ++) {
        int iel = el_elmnts[ispec];
        int id = ispec * NGLL;
        const double *weight = wgll.data();
        const double *hpT = hprimeT.data();
        const double *hp = hprime.data();
        int NGL = NGLL;

        // grl case
        if(ispec == nspec_el) {
            weight = wgrl.data();
            hpT = hprimeT_grl.data();
            hp = hprime_grl.data();
            NGL = NGRL;
        }   
        // cache temporary arrays
        for(int i = 0; i < NGL; i ++) {
            sumC[i] = xC[id + i] * weight[i] / jaco[iel];
            sumL[i] = xL[id + i] * weight[i] / jaco[iel];
        }

        // compute M/K/E
        for(int i = 0; i < NGL; i ++) {
            int iglob = ibool_el[id + i];
            double temp = weight[i] * jaco[iel];

            // element wise M/K1/K3
            double M0 = temp * xrho_el[id + i];
            double K1 = temp * xA[id + i];
            double K3 = temp * xL[id + i];

            // assemble
            Mmat[iglob] += M0;
            Mmat[iglob + nglob_el] += M0;
            Kmat[iglob * ng + iglob] += K1;
            Kmat[(nglob_el + iglob) * ng + (nglob_el + iglob)] += K3;

            // other matrices
            for(int j = 0; j < NGL; j ++) {
                int iglob1 = ibool_el[id + j];
                double E1{},E3{};
                for(int m = 0; m < NGL; m ++) {
                    E1 += sumL[m] * hpT[i * NGL + m] * hpT[j * NGL + m];
                    E3 += sumC[m] * hpT[i * NGL + m] * hpT[j * NGL + m];
                }
                Emat[iglob * ng + iglob1] += E1;
                Emat[(iglob + nglob_el) * ng + (iglob1 + nglob_el)] += E3;

                // K2/E2
                double K2 = weight[j] * xF[id + j] * hpT[i * NGL + j] - 
                            weight[i] * xL[id + i] * hp[i * NGL + j];
                double E2 = weight[i] * xF[id + i] * hp[i * NGL + j] - 
                            weight[j] * xL[id + j] * hpT[i * NGL + j];
                Kmat[(nglob_el + iglob) * ng + iglob1] += K2;
                Emat[iglob * ng + nglob_el +  iglob1] += E2;
            }
        }
    }

    // acoustic case
    for(int ispec = 0; ispec < nspec_ac + nspec_ac_grl; ispec ++) {
        int iel = ac_elmnts[ispec];
        int id = ispec * NGLL;
        const double *weight = wgll.data();
        const double *hpT = hprimeT.data();
        const double *hp = hprime.data();
        int NGL = NGLL;

        // grl case
        if(ispec == nspec_ac) {
            weight = wgrl.data();
            hpT = hprimeT_grl.data();
            hp = hprime_grl.data();
            NGL = NGRL;
        }   
        // cache temporary arrays
        for(int i = 0; i < NGL; i ++) {
            sumL[i] =  weight[i] / jaco[iel] / xrho_ac[id+i];
        }

        // compute M/K/E
        for(int i = 0; i < NGL; i ++) {
            int ig0 = ibool_ac[id + i];
            if(ig0 == -1) continue;
            int iglob = ig0 + nglob_el * 2;
            double temp = weight[i] * jaco[iel];

            // assemble M and K
            Mmat[iglob] += temp / xkappa_ac[id + i];
            Kmat[iglob * ng + iglob] += temp / xrho_ac[id + i];

            // assemble E
            for(int j = 0; j < NGL; j ++) {
                int ig1 = ibool_ac[id + j];
                if(ig1 == -1) continue;
                int iglob1 = ig1 + nglob_el * 2;
                double s{};
                for(int m = 0; m < NGL; m ++) {
                    s += sumL[m] * hpT[i * NGL + m] * hpT[j * NGL + m];
                }
                Emat[iglob * ng + iglob1] += s;
            }
        }
    }

    // acoustic free boundary condition
    if(ac_elmnts.size() < 0 &&  ac_elmnts[0] == 0) { 
        double coef = 1. / (jaco[0] * xrho_ac[0]); 
        std::array<double,NGLL> sumE;
        for(int i = 1; i < NGLL; i ++) {
            double s{};
            for(int m = 0; m < NGLL; m ++) {
                s += wgll[m] / jaco[0] / xrho_ac[m] * hprimeT[m] * hprimeT[i * NGLL + m];
            }
            sumE[i] = s;
        }

        for(int i = 1; i < NGLL; i ++) {
            int iglob = ibool_ac[i] + nglob_el * 2;
            Emat[nglob_el * 2 * ng + iglob] += coef * hprime[i] + sumE[i];
        }
    }

    // acoustic-elastic boundary
    double om = M_PI * 2 * freq;
    for(int iface = 0; iface < nfaces_bdry; iface ++) {
        int ispec_ac = ispec_bdry_loc[iface * 2 + 0];
        int ispec_el = ispec_bdry_loc[iface * 2 + 1];
        double norm = -1.;
        int igll_el = 0;
        int igll_ac = NGLL - 1;
        if(! is_top_ac_bdry[iface]) {
            norm = 1.;
            igll_ac = 0;
            igll_el = NGLL - 1;
        }

        // get ac/el global loc
        int iglob_el = ibool_el[ispec_el * NGLL + igll_el];
        int iglob_ac = ibool_ac[ispec_ac * NGLL + igll_ac];

        // add contribution to E mat, elastic case
        // E(nglob_el + iglob_el, nglob_el*2 + iglob_ac) += 
        Emat[(nglob_el + iglob_el) * ng + (nglob_el * 2 + iglob_ac)] += om * om * norm;
        
        // acoustic case
        // E(nglob_el*2 + iglob_ac, nglob_el + iglob_el) += norm
        Emat[(nglob_el*2 + iglob_ac) * ng + (nglob_el + iglob_el)] += norm;
    }
}

void
filter_swd(double vmin, double vmax,double om, const Eigen::MatrixXcd &displ_all,
            const Eigen::Array<std::complex<double>,-1,1> &k,std::vector<double> &c,
           std::vector<double> &displ);

/**
 * @brief compute dispersion curves and eigenfunction at frequency freq
 * 
 * @param freq target frequency
 * @param c output dispersion curves, all modes, shape(nc)
 * @param egn eigenfunctions (U,V,chi), shape(nc,nglob_el*2 + nglob_ac)
 */
void LayerModelMultiPhyVTI:: 
compute_egnfun(double freq, std::vector<double> &c, std::vector<double> &egn) const
{
    // map M/E/K matrices
    typedef Eigen::Matrix<double,-1,-1,1> dmat2;
    int ng = nglob_el * 2 + nglob_ac;
    Eigen::Map<const Eigen::VectorXd> M(Mmat.data(),ng);
    Eigen::Map<const dmat2> E(Emat.data(),ng,ng);
    Eigen::Map<const dmat2> K(Kmat.data(),ng,ng);

    // prepare matrix A = om^2 M -E
    double om = 2. * M_PI * freq;
    double omega2 = om * om;
    Eigen::MatrixXf A = (omega2 * dmat2(M.asDiagonal()) - E).cast<float>();
    Eigen::MatrixXf K1 = K.cast<float>().inverse();
    A = K1 * A;

    // solve eigenvalue problem A x = k^2 K x
    Eigen::EigenSolver<Eigen::MatrixXf> sol;
    sol.compute(A);
    Eigen::MatrixXcd egn_all = sol.eigenvectors().cast<dcmplx>();
    Eigen::Array<dcmplx,-1,1> k = sol.eigenvalues().array().sqrt().cast<dcmplx>();

    // filter swd
    filter_swd(PHASE_VELOC_MIN,PHASE_VELOC_MAX,om,egn_all,k,c,egn);

    // note we obtained is [U,V^{\bar}，char^{\bar}]
    int nc = c.size();
    for(int ic = 0; ic < nc; ic ++) {
        double k0 = om / c[ic];
        for(int iglob = nglob_el; iglob < ng; iglob ++) {
            egn[ic * ng + iglob] /= k0;
        }
    }
}
