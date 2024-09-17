#include "swdlayertti.hpp"
typedef std::complex<double> dcmplx;

extern "C" {
/**
 * @brief compute tti sensitivity kernels
 * 
 * @param NGL shape of input parameters 
 * @param k wavenumber
 * @param A,C,F,L,N,theta0,dphi TTI parameters, shape (NGL) 
 * @param U,V,W,dU,dV,dW displacements and their z derivative, shape (NGL) 
 * @param Kwvnm dL/dk shape (NGL)
 * @param KA,KC,KF,KL,KN,KT,KP dL/dparam, shape(NGL) 
 * @param dL_dkv  dL/ d(kx,ky), shape(NGL)
 */
void get_kernels_(int NGL,double k,const double *A,const double *C,const double *F,
                const double *L,const double *N,const double *theta0,
                const double *dphi, const dcmplx *U,const dcmplx *V,const dcmplx *W,
                const dcmplx *dU, const dcmplx *dV, const dcmplx *dW,
                double *Kwvnm, double *KA, double *KC, double *KF, double *KL,
                double *KN, double *KT, double *KP, double *dL_dkv);
}

/**
 * @brief compute group velocity and kernels for tti model
 * 
 * @param freq current frequency
 * @param c  phase velocity at this frequency
 * @param phi azimuthal angle of c
 * @param displ eigen function, shape(nglob * 3)
 * @param frekl Frechet kernels A/C/F/L/N/T/P/rho_kl kernels for elastic parameters, shape(8,nspec*NGLL + NGRL) 
 * @return double u group velocity and it's azimthual angle
 */
std::array<double,2> LayerModelTTI ::
compute_kernels(double freq, double c,double phi,
                const std::vector<dcmplx> &displ,
                std::vector<double> &frekl) const 
{
   // first allocate element wise displ
    std::array<dcmplx,NGRL> U,dU,V,dV,W,dW;
    std::array<double,NGRL> dphi,Kwvnm;
    std::array<double,NGRL * 2> dL_dkv;

    // resize
    size_t size = xrho.size();
    frekl.resize(8 * size);
    double *rho_kl = &frekl[7 * size];

    // consts
    double om = 2 * M_PI * freq;
    double k = om / c;

    // loop every element
    double I1{}, I2{};
    double I3x{},I3y{}; 
    for(int ispec = 0; ispec < nspec + 1; ispec += 1) {
        const double *hp = &hprime[0];
        const double *w = &wgll[0];
        int NGL = NGLL;
        int id = ispec * NGLL;

        // GRL layer
        if(ispec == nspec) {
            hp = &hprime_grl[0];
            w = &wgrl[0];
            NGL = NGRL;
        }

        // cache displ in a element
        for(int i = 0; i < NGL; i ++) {
            int iglob = ibool[id + i];
            U[i] = displ[iglob];
            V[i] = displ[nglob + iglob];
            dphi[i] = xP[id + i] - M_PI / 180. * phi;
        }

        // compute derivative 
        for(int i = 0; i < NGL; i ++) {
            dcmplx sx{},sy{},sz{};
            for(int j = 0; j < NGL; j ++) {
                sx += U[j] * hp[i * NGL + j];
                sy += W[j] * hp[i * NGL + j];
                sz += V[j] * hp[i * NGL + j];
            }
            dU[i] = sx / jaco[ispec];
            dW[i] = sy / jaco[ispec];
            dV[i] = sz / jaco[ispec];
        }

        // compute A/C/F/L/N/T/P kernels
        get_kernels_(NGL,k,&xA[id],&xC[id],&xF[id],&xL[id],&xN[id],&xT[id],dphi.data(),
                    U.data(),V.data(),W.data(),dU.data(),dV.data(),dW.data(),Kwvnm.data(),
                    &frekl[id],&frekl[id+size],&frekl[id+ 2*size],&frekl[id+3*size],
                    &frekl[id+4*size],&frekl[id+5*size],&frekl[id+6*size],dL_dkv.data());

        // compute kernel/energy integral 
        // I1 = int dL/dk dz 
        // I2 = int dL/dom dz
        // I3x = int dL/d_kx dz  I3y = int dL/d_ky dz
        for(int i = 0; i < NGL; i ++) {
            // rho kernel
            double rho = xrho[id + i];
            double disp_sq = (U[i] * std::conj(U[i]) + V[i] * std::conj(V[i]) + W[i] * std::conj(W[i])).real();
            rho_kl[id + i] = om * om * disp_sq;

            // accumulate I1/I2/I3x/I3y
            I1 += Kwvnm[i] * w[i] * jaco[ispec];
            I2 += 2 * om * rho * w[i] * disp_sq * jaco[ispec];
            I3x += dL_dkv[i] * w[i] * jaco[ispec];
            I3y += dL_dkv[i + NGL] * w[i] * jaco[ispec];
        }
    }

    // compute group velocity
    std::array<double,2> cg;
    cg[0] = -I3x / I2; cg[1] = -I3y / I2;

    // rescale kernels with I1
    double coef1 = c / k / I1;
    for(int iker = 0; iker < 8; iker ++) {
        for(int i = 0; i < nspec * NGLL + NGRL; i ++) {
            frekl[iker * 8 + i] *= coef1;
        }
    }

    return cg;
}