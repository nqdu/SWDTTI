#include "vti/swdlayervti.hpp"

/**
 * @brief compute group velocity and kernels for love wave
 * 
 * @param freq current frequency
 * @param c  current phase velocity
 * @param displ eigen function, shape(nglob)
 * @param frekl Frechet kernels (N/L/rho) for elastic parameters, shape(3,nspec*NGLL + NGRL) 
 * @return double u group velocity
 */
double LayerModelVTI:: 
compute_love_kl(double freq,double c,const double *displ, 
                std::vector<double> &frekl) const
{
   // first allocate element wise displ
    std::array<double,NGRL> W,dW;

    // resize
    frekl.resize(3 * xrho.size());
    double *N_kl = &frekl[0];
    double *L_kl = &frekl[xrho.size()];
    double *rho_kl = &frekl[xrho.size() * 2];

    // init I1/I2
    // I1 = \int rho W^2 dz 
    // I2 = int N W^2 dz 
    double I1{}, I2{};

    // kernels for N/L/rho (Aki and Richards, 2002, (7.71))
    double om = 2 * M_PI * freq;
    double k = om / c;
    for(int ispec = 0; ispec < nspec + 1; ispec ++) {
        const double *hp = &hprime[0];
        const double *w = &wgll[0];
        int NGL = NGLL;
        int id0 = ispec * NGLL;

        // GRL layer
        if(ispec == nspec) {
            hp = &hprime_grl[0];
            w = &wgrl[0];
            NGL = NGRL;
        }

        // cache displ in a element
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            int iglob = ibool[id];
            W[i] = displ[iglob];
        }

        // compute derivative 
        for(int i = 0; i < NGL; i ++) {
            double s{};
            for(int j = 0; j < NGL; j ++) {
                s += W[j] * hp[i * NGL + j];
            }
            dW[i] = s  / jaco[ispec];
        }

        // compute kernel/energy integral 
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            double rho = xrho[id];
            double N = xN[id];

            // N/L/rho kernel
            N_kl[id] = std::pow(k * W[i],2);
            L_kl[id] = std::pow(dW[i],2);
            rho_kl[id] = -std::pow(om * W[i],2);

            // accumulate I1/I2 (Aki and Richards, 2002, (7.66))
            auto tmp = rho * w[i] * W[i] * W[i] * jaco[ispec];
            auto tmp1 = N * w[i] * W[i] * W[i] * jaco[ispec];
            I1 += tmp;
            I2 += tmp1;
        }
    }

    // group velocity
    double u = I2 / (c * I1);

    // rescale kernels
    const double coef = c / (2. * k * k * I2);
    for(int i = 0; i < nspec * NGLL + NGRL; i ++) {
        N_kl[i] *= coef; L_kl[i] *= coef;
        rho_kl[i] *= coef;
    }

    return u;
}

/**
 * @brief compute group velocity and kernels for love wave
 * 
 * @param freq current frequency
 * @param c  current phase velocity
 * @param displ eigen function, shape(nglob * 2)
 * @param frekl Frechet kernels A/C/L/F/rho_kl kernels for elastic parameters, shape(5,nspec*NGLL + NGRL) 
 * @return double u group velocity
 */
double LayerModelVTI:: 
compute_rayl_kl(double freq,double c,const double *displ, 
                std::vector<double> &frekl) const
{
   // first allocate element wise displ
    std::array<double,NGRL> U,dU,V,dV;

    // resize
    size_t size = xrho.size();
    frekl.resize(5 * size);
    double *A_kl = &frekl[0], *C_kl = &frekl[size];
    double *L_kl = &frekl[2*size], *F_kl = &frekl[size * 3];
    double *rho_kl = &frekl[4 * size];

    // consts
    double om = 2 * M_PI * freq;
    double k = om / c;

    // init energy integral
    // I1 = \int rho (v^2 + u^2) dz
    // I2 = \int (Au^2 + LV^2) dz 
    // I3 = \int (-LV \dot{U} + FU\dot{V}) dz 
    double I1{},I2{},I3{};

    // loop every element
    for(int ispec = 0; ispec < nspec + 1; ispec += 1) {
        const double *hp = &hprime[0];
        const double *w = &wgll[0];
        int NGL = NGLL;
        int id0 = ispec * NGLL;

        // GRL layer
        if(ispec == nspec) {
            hp = &hprime_grl[0];
            w = &wgrl[0];
            NGL = NGRL;
        }

        // cache displ in a element
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            int iglob = ibool[id];
            U[i] = displ[iglob];
            V[i] = displ[nglob + iglob];
        }

        // compute derivative 
        for(int i = 0; i < NGL; i ++) {
            double sx{},sz{};
            for(int j = 0; j < NGL; j ++) {
                sx += U[j] * hp[i * NGL + j];
                sz += V[j] * hp[i * NGL + j];
            }
            dU[i] = sx  / jaco[ispec];
            dV[i] = sz  / jaco[ispec];
        }

        // compute kernel/energy integral 
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            double rho = xrho[id];
            double A = xA[id];
            double L = xL[id], F = xF[id];

            // kernels
            A_kl[id] = k * k * U[i] * U[i];
            C_kl[id] = std::pow(dV[i],2);
            F_kl[id] = 2 * k * U[i] * dV[i];
            L_kl[id] = std::pow(dU[i],2) + std::pow(k * V[i],2) - 2. * k * V[i] * dU[i];
            rho_kl[id] = -om * om * (U[i] * U[i] + V[i] * V[i]);

            // accumulate I1/I2 (Aki and Richards, 2002, (7.66))
            auto tmp1 = rho * w[i] * (U[i] * U[i] + V[i] * V[i]) * jaco[ispec];
            auto tmp2 = (A * U[i] * U[i] + L * V[i] * V[i]) * w[i] * jaco[ispec];
            auto tmp3 = (-L * V[i] * dU[i] + F * U[i] * dV[i]) * w[i] * jaco[ispec];
            I1 += tmp1;
            I2 += tmp2;
            I3 += tmp3;
        }
    }

    // group velocity
    double u = (I2 + I3 / k) /(c * I1);

    // rescale kernels
    double coef = 1. / (2. * k * k * u * I1);
    for(int i = 0; i < nspec * NGLL + NGRL; i ++) {
        A_kl[i] *= coef; L_kl[i] *= coef;
        C_kl[i] *= coef; F_kl[i] *= coef;
        rho_kl[i] *= coef;
    }

    return u;
}

/**
 * @brief transform kernels from base to rho/vsv/vsh(love) and rho/vpv/vph/vsv/eta (rayleigh)
 * 
 * @param frekl base Frechet kernels, shape(3/5,nspec*NGLL+NGRL)
 */
void LayerModelVTI:: 
transform_kernels(std::vector<double> &frekl) const
{
    // get # of kernels
    int npts = nspec * NGLL + NGRL;
    int nker = frekl.size() / npts;

    // loop every point to transform kernels
    for(int ipt = 0; ipt < npts; ipt ++) {
        // set zero of all base kernels
        double A_kl{},C_kl{},F_kl{}, L_kl{}, N_kl{}, rho_kl{};

        // fetch kernels from global array
        switch (nker)
        {
        case 3:
            N_kl = frekl[0 * npts + ipt];
            L_kl = frekl[1 * npts + ipt];
            break;
        
        default:
            A_kl = frekl[0 * npts + ipt];
            C_kl = frekl[1 * npts + ipt];
            L_kl = frekl[2 * npts + ipt];
            F_kl = frekl[3 * npts + ipt];
            break;
        }
        rho_kl = frekl[(nker - 1) * npts + ipt];

        // get variables
        double A = xA[ipt], C = xC[ipt], L = xL[ipt];
        double F = xF[ipt], N = xN[ipt], rho = xrho[ipt];

        // compute vph/vpv/vsh/vsv/eta
        double vph = std::sqrt(A / rho), vpv = std::sqrt(C / rho);
        double vsh = std::sqrt(N / rho), vsv = std::sqrt(L / rho);
        //double eta = F / (A - 2. * L);

        // transform kernels
        double vph_kl = 2. * rho * vph * A_kl, vpv_kl = 2. * rho * vpv * C_kl;
        double vsh_kl = 2. * rho * vsh * N_kl, vsv_kl = 2. * rho * vsv * L_kl;
        double eta_kl = (A - 2. * L) * F_kl;
        double r_kl = vph * vph * A_kl + vpv * vpv * C_kl + 
                      vsh * vsh * N_kl + vsv * vsv * L_kl +
                      rho_kl;

        // copy back to frekl array
        switch (nker)
        {
        case 3:
            frekl[1 * npts + ipt] = vsv_kl;
            frekl[2 * npts + ipt] = vsh_kl;
            break;
        
        default:
            frekl[1 * npts + ipt] = vpv_kl;
            frekl[2 * npts + ipt] = vph_kl;
            frekl[3 * npts + ipt] = vsv_kl;
            frekl[4 * npts + ipt] = eta_kl;
            break;
        }
        frekl[0 * npts + ipt] = r_kl;
    }       
}