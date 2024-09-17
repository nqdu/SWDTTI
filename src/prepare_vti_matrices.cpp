#include "swdlayer.hpp"

/**
 * @brief prepare M/K/E matrices
 * 
 * @param wavetype = 1 for Love = 2 for Rayleigh
 */
void LayerModel:: 
prepare_matrices(int wavetype)
{
    switch (wavetype)
    {
    case 1:
        this -> prepare_matrices_love();
        break;
    case 2:
        this -> prepare_matrices_rayl();
        break;
    default:
        printf("wavetype should be one of [1,2]\n");
        exit(1);
        break;
    }
}

/**
 * @brief prepare M/K/E matrices for Love wave
 * 
 */
void LayerModel:: 
prepare_matrices_love()
{
    // allocate space and set zero
    Mmat.resize(nglob);
    Kmat.resize(nglob);
    Emat.resize(nglob * nglob);
    for(int i = 0; i < nglob; i ++) {
        Mmat[i] = 0.;
        Kmat[i] = 0.;
        for(int j = 0; j < nglob; j ++) {
            Emat[i * nglob + j] = 0.;
        }
    }

    // compute M/K/E for gll/grl layer
    std::array<double,NGRL> sum_terms;
    for(int ispec = 0; ispec < nspec + 1; ispec ++) {
        int id = ispec * NGLL;
        const double *weight = wgll.data();
        const double *hpT = hprimeT.data();
        int NGL = NGLL;

        // grl case
        if(ispec == nspec) {
            weight = wgrl.data();
            hpT = hprimeT_grl.data();
            NGL = NGRL;
        }   

        // cache temporary arrays
        for(int i = 0; i < NGL; i ++) {
            sum_terms[i] = xL[id + i] * weight[i] / jaco[ispec];
        }

        // compute M/K/E
        for(int i = 0; i < NGL; i ++) {
            int iglob = ibool[id + i];
            double temp = weight[i] * jaco[ispec];
            Mmat[iglob] += temp * xrho[id + i];
            Kmat[iglob] += temp * xN[id + i];

            for(int j = 0; j < NGL; j ++) {
                int iglob1 = ibool[id + j];
                double s{};
                for(int m = 0; m < NGL; m ++) {
                    s += sum_terms[m] * hpT[i * NGL + m] * hpT[j * NGL + m];
                }
                Emat[iglob * nglob + iglob1] += s;
            }
        }
    }
}

/**
 * @brief prepare M/K/E matrices for Rayleigh wave
 * 
 */
void LayerModel:: 
prepare_matrices_rayl()
{
    // allocate space and set zero
    Mmat.resize(nglob * 2);
    Emat.resize(nglob * nglob * 4);
    Kmat.resize(nglob * nglob * 4);
    for(int i = 0; i < nglob * 2; i ++) {
        Mmat[i] = 0.;
        for(int j = 0; j < nglob * 2; j ++) {
            Emat[i * nglob * 2 + j] = 0.;
            Kmat[i * nglob * 2 + j] = 0.;
        }
    }

    // size of Kmat.rows()
    int row = nglob * 2;

    // compute M/K/E for gll/grl layer
    std::array<double,NGRL> sumC,sumL;
    for(int ispec = 0; ispec < nspec + 1; ispec ++) {
        int id = ispec * NGLL;
        const double *weight = wgll.data();
        const double *hpT = hprimeT.data();
        const double *hp = hprime.data();
        int NGL = NGLL;

        // grl case
        if(ispec == nspec) {
            weight = wgrl.data();
            hpT = hprimeT_grl.data();
            hp = hprime_grl.data();
            NGL = NGRL;
        }   

        // cache temporary arrays
        for(int i = 0; i < NGL; i ++) {
            sumC[i] = xC[id + i] * weight[i] / jaco[ispec];
            sumL[i] = xL[id + i] * weight[i] / jaco[ispec];
        }

        // compute M/K/E
        for(int i = 0; i < NGL; i ++) {
            int iglob = ibool[id + i];
            double temp = weight[i] * jaco[ispec];

            // element wise M/K1/K3
            double M0 = temp * xrho[id + i];
            double K1 = temp * xA[id + i];
            double K3 = temp * xL[id + i];

            // assemble
            Mmat[iglob] += M0;
            Kmat[iglob * row + iglob] += K1;
            Kmat[(nglob + iglob) * row + (nglob + iglob)] += K3;

            // other matrices
            for(int j = 0; j < NGL; j ++) {
                int iglob1 = ibool[id + j];
                double E1{},E3{};
                for(int m = 0; m < NGL; m ++) {
                    E1 += sumL[m] * hpT[i * NGL + m] * hpT[j * NGL + m];
                    E3 += sumC[m] * hpT[i * NGL + m] * hpT[j * NGL + m];
                }
                Emat[iglob * row + iglob1] += E1;
                Emat[(iglob + nglob) * row + (iglob1 + nglob)] += E3;

                // K2/E2
                double K2 = weight[j] * xF[id + j] * hpT[i * NGL + j] - 
                            weight[i] * xL[id + i] * hp[i * NGL + j];
                double E2 = weight[i] * xF[id + i] * hp[i * NGL + j] - 
                            weight[j] * xL[id + j] * hpT[i * NGL + j];
                Kmat[(nglob + iglob) * row + iglob1] += K2;
                Emat[iglob * row + nglob +  iglob1] += E2;
            }
        }
    }

    // M[nglob:-1] = M[0:nglob]
    for(int i = 0; i < nglob; i ++) {
        Mmat[i + nglob] = Mmat[i];
    }
}