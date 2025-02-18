#include "multiphysics/vti_acoustic.hpp"

/**
 * @brief compute group velocity and kernels for love wave
 * 
 * @param freq current frequency
 * @param c  current phase velocity
 * @param egn eigen function, shape(nglob_el*2 + nglob_ac), (U,V,chi)
 * @param frekl Frechet kernels A/C/L/F/rho_kl or kappa/rho_kl (5,npts_el),in fluid C/L/F_kl = 0
 * @return double u group velocity
 */
double LayerModelMultiPhyVTI:: 
compute_kernels(double freq, double c,const double *egn,
                std::vector<double> &frekl) const
{
   // first allocate element wise egn
    std::array<double,NGRL> U,dU,V,dV;
    std::array<double,NGRL> chi,dchi;

    // resize
    size_t size_all = ibool.size();
    frekl.resize(5 * size_all);
    std::fill(frekl.begin(),frekl.end(),0.);

    // consts
    double om = 2 * M_PI * freq;
    double k = om / c;
    int ng = nglob_ac + nglob_el * 2;

    // init energy integral
    // I1 = \int rho (v^2 + u^2) dz = \partial_{\omega} L / om
    //  k * I2 + I3 = \int (-\partial_k L) dz
    double I1{},I2{},I3{};

    // loop every elastic element
    for(int ispec = 0; ispec < nspec_el + nspec_el_grl; ispec += 1) {
        int iel = el_elmnts[ispec];
        const double *hp = &hprime[0];
        const double *w = &wgll[0];
        int NGL = NGLL;
        int id0 = ispec * NGLL;
        const double J = jaco[iel];

        // GRL layer
        if(ispec == nspec_el) {
            hp = &hprime_grl[0];
            w = &wgrl[0];
            NGL = NGRL;
        }

        // cache egn in a element
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            int iglob = ibool_el[id];
            U[i] = egn[iglob];
            V[i] = egn[nglob_el + iglob];
        }

        // compute derivative 
        for(int i = 0; i < NGL; i ++) {
            double sx{},sz{};
            for(int j = 0; j < NGL; j ++) {
                sx += U[j] * hp[i * NGL + j];
                sz += V[j] * hp[i * NGL + j];
            }
            dU[i] = sx  / J;
            dV[i] = sz  / J;
        }

        // compute kernel/energy integral 
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            int id1 = iel * NGLL + i;
            double rho = xrho_el[id];
            double A = xA[id];
            double L = xL[id], F = xF[id];
            double U2 = U[i] * U[i], V2 = V[i] * V[i];
            double dU2 = dU[i] * dU[i], dV2 = dV[i] * dV[i];

            // kernels dc = \int (-2 KL) /(2 k^2 U I1)
            // A/C/L/F/rho
            frekl[size_all * 0 + id1] = k * k * U2;
            frekl[size_all * 1 + id1] = dV2;
            frekl[size_all * 2 + id1] = dU2 + std::pow(k,2) * V2 - 2. * k * V[i] * dU[i];
            frekl[size_all * 3 + id1] = 2 * k * U[i] * dV[i];
            frekl[size_all * 3 + id1] = -om * om * (U2 + V2);

            // accumulate I1/I2 (Aki and Richards, 2002, (7.66))
            auto tmp1 = rho * (U2 + V2) * w[i] * J;
            auto tmp2 = (A * U2 + L * V2) * w[i] * J;
            auto tmp3 = (-L * V[i] * dU[i] + F * U[i] * dV[i]) * w[i] * J;
            I1 += tmp1;
            I2 += tmp2;
            I3 += tmp3;
        }
    }

    // acoustic elment
    for(int ispec = 0; ispec < nspec_ac + nspec_ac_grl; ispec += 1) {
    // for(int ispec = 0; ispec <0; ispec ++) {
        int iel = ac_elmnts[ispec];
        const double *hp = &hprime[0];
        const double *w = &wgll[0];
        int NGL = NGLL;
        int id0 = ispec * NGLL;
        const double J = jaco[iel];

        // GRL layer
        if(ispec == nspec_ac) {
            hp = &hprime_grl[0];
            w = &wgrl[0];
            NGL = NGRL;
        }

        // cache chi in an element
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            int iglob = ibool_ac[id];
            chi[i] = (iglob == -1) ? 0.: egn[nglob_el * 2 + iglob];
        }
    
        // compute derivative 
        for(int i = 0; i < NGL; i ++) {
            double s{};
            for(int j = 0; j < NGL; j ++) {
                s += chi[j] * hp[i * NGL + j];
            }
            dchi[i] = s / J;
        }

        // compute kernel/energy integral 
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            double rhoinv = 1. / xrho_ac[id];
            double kapainv = 1. / xkappa_ac[id];
            double dchi2 = dchi[i] * dchi[i];
            double chi2 = chi[i] * chi[i];
            int id1 = iel * NGLL + i;

            // note: dc = \int (-2 KL) /(2 k^2 U I1)
            frekl[size_all * 0 + id1] = std::pow(kapainv * om * om,2) * chi2; // kappa_kl
            frekl[size_all * 4 + id1] = -std::pow(rhoinv * om,2) * (dchi2 + k*k * chi2);; // rho_kl
            // rhoa_kl[id] = -std::pow(rhoinv * om,2) * (dchi2 + k*k * chi2);
            // kappa_kl[id] = std::pow(kapainv * om * om,2) * chi2;

            // energy integral
            // I1 = \partial_{\omega} L / om
            //  k * I2 + I3 = \int (-\partial_k L) dz
            auto tmp1 = rhoinv * (dchi2 + k*k * chi2); //- 2. * std::pow(om,2) * kapainv * chi2;
            auto tmp2 = om * om * rhoinv * chi2; // (-\partial_k L) / k
            I1 += w[i] * J * tmp1;
            I2 += w[i] * J * tmp2;
        }
    }

    // group velocity
    double u = (I2 + I3 / k) /(c * I1);

    // rescale kernels
    double coef = 1. / (2. * k * k * u * I1);
    for(size_t i = 0; i < frekl.size(); i ++) {
        frekl[i] *= coef;
    }

    return u;
}


/**
 * @brief transform kernels from base to rho/vpv/vph/vsv/eta, for acoustic only rho/vpv
 * 
 * @param frekl base frechet kernels, shape(5,nspec*NGLL+NGRL),A/C/L/F/rho
 */
void LayerModelMultiPhyVTI:: 
transform_kernels(std::vector<double> &frekl) const
{
    int npts = ibool.size();

    // loop every elastic element
    for(int ispec = 0; ispec < nspec_el + nspec_el_grl; ispec += 1) {
        int iel = el_elmnts[ispec];
        int NGL = NGLL;
        int id0 = ispec * NGLL;

        // GRL layer
        if(ispec == nspec_el) {
            NGL = NGRL;
        }

        // loop every point in this element
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            int ipt = iel * NGLL + i;

            // kernels
            double A_kl{},C_kl{},F_kl{}, L_kl{}, N_kl{}, rho_kl{};

            // kernels
            A_kl = frekl[0 * npts + ipt];
            C_kl = frekl[1 * npts + ipt];
            L_kl = frekl[2 * npts + ipt];
            F_kl = frekl[3 * npts + ipt];
            rho_kl = frekl[4 * npts + ipt];

            // get variables
            double A = xA[id], C = xC[id], L = xL[id];
            double F = xF[id], rho = xrho_el[id];

            // compute vph/vpv/vsh/vsv/eta
            double vph = std::sqrt(A / rho), vpv = std::sqrt(C / rho);
            double vsv = std::sqrt(L / rho);
            //double eta = F / (A - 2. * L);

             // transform kernels
            double vph_kl = 2. * rho * vph * A_kl, vpv_kl = 2. * rho * vpv * C_kl;
            double vsv_kl = 2. * rho * vsv * L_kl;
            double eta_kl = (A - 2. * L) * F_kl;
            double r_kl = vph * vph * A_kl + vpv * vpv * C_kl + 
                         + vsv * vsv * L_kl + rho_kl;

            // copy back
            frekl[0 * npts + ipt] = rho_kl;
            frekl[1 * npts + ipt] = vpv_kl;
            frekl[2 * npts + ipt] = vph_kl;
            frekl[3 * npts + ipt] = vsv_kl;
            frekl[4 * npts + ipt] = eta_kl;
        }
    }

    // acoustic domain
    for(int ispec = 0; ispec < nspec_ac + nspec_ac_grl; ispec += 1) {
        int iel = ac_elmnts[ispec];
        const double *hp = &hprime[0];
        const double *w = &wgll[0];
        int NGL = NGLL;
        int id0 = ispec * NGLL;
        const double J = jaco[iel];

        // GRL layer
        if(ispec == nspec_ac) {
            NGL = NGRL;
        }

        // loop every point in this element
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            int ipt = iel * NGLL + i;

            // kernels
            double kappa_kl{}, rho_kl{};

            // kernels
            kappa_kl = frekl[0 * npts + ipt];
            rho_kl = frekl[4 * npts + ipt];

            // get variables
            double kappa = xkappa_ac[id], rho = xrho_ac[id];

            // compute vp
            double vp = std::sqrt(kappa / rho);

            // transform kernels
            double vp_kl = 2 * rho * vp * kappa_kl; 
            double r_kl = vp * vp * kappa_kl + rho_kl;

            // copy back
            frekl[0 * npts + ipt] = rho_kl;
            frekl[1 * npts + ipt] = vp_kl;
            frekl[2 * npts + ipt] = 0.;
            frekl[3 * npts + ipt] = 0.;
            frekl[4 * npts + ipt] = 0.;
        }
    }
}

/**
 * @brief convert eigen functions to displacement 
 * 
 * @param egn shape(nglob_ac + nglob_el * 2)
 * @param displ shape(2,ibool.size())
 */
void LayerModelMultiPhyVTI:: 
transform_egn2disp(double freq,double c,const double *egn, double *__restrict displ) const
{
    // consts
    int npts = ibool.size();
    std::array<double,NGRL> chi;
    double k = 2 * M_PI * freq / c;

    // loop every elastic element
    for(int ispec = 0; ispec < nspec_el + nspec_el_grl; ispec += 1) {
        int iel = el_elmnts[ispec];
        int NGL = NGLL;
        int id0 = ispec * NGLL;
        int id1 = iel * NGLL;

        // GRL layer
        if(ispec == nspec_el) {
            NGL = NGRL;
        }

      // cache egn in a element
        for(int i = 0; i < NGL; i ++) {
            int iglob = ibool_el[id0 + i];
            displ[0 * npts + id1 + i] = egn[iglob];
            displ[1 * npts + id1 + i] = egn[iglob + nglob_el];
        }
    }

    // loop each acoustic element
    for(int ispec = 0; ispec < nspec_ac + nspec_ac_grl; ispec += 1) {
        int iel = ac_elmnts[ispec];
        int NGL = NGLL;
        int id0 = ispec * NGLL;
        int id1 = iel * NGLL;
        const double *hp = &hprime[0];
        const double *w = &wgll[0];
        const double J = jaco[iel];

        // GRL layer
        if(ispec == nspec_ac) {
            NGL = NGRL;
            hp = &hprime_grl[0];
            w = &wgrl[0];
        }

        // cache chi in an element
        for(int i = 0; i < NGL; i ++) {
            int id = id0 + i;
            int iglob = ibool_ac[id];
            chi[i] = (iglob == -1) ? 0.: egn[nglob_el * 2 + iglob];
        }

        // compute derivative  dchi / dz
        for(int i = 0; i < NGL; i ++) {
            double dchi{};
            for(int j = 0; j < NGL; j ++) {
                dchi += chi[j] * hp[i * NGL + j];
            }
            dchi /= J;

            // set value to displ
            displ[0 * npts + id1 + i] = k / xrho_ac[id0 + i] * chi[i];
            displ[1 * npts + id1 + i] = dchi / xrho_ac[id0 + i];
        }

    }
}