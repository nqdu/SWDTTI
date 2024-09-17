#include "swdlayertti.hpp"

typedef std::complex<double> dcmplx;

// fortran functions
extern "C" {

void get_comp1_(int NGL,const double *A,const double *C,const double *F,
                const double *L,const double *N,const double *theta0,
                const double *dphi,double jaco, const double *weight,
                const double *hp,const double *hpT, dcmplx *K0U,dcmplx *K0V,
                dcmplx *K0W,dcmplx *K1U,dcmplx *K1V,dcmplx *K1W,dcmplx *K2U,
                dcmplx *K2V,dcmplx *K2W);
void get_comp2_(int NGL,const double *A,const double *C,const double *F,
                const double *L,const double *N,const double *theta0,
                const double *dphi,double jaco, const double *weight,
                const double *hp,const double *hpT, dcmplx *K0U,dcmplx *K0V,
                dcmplx *K0W,dcmplx *K1U,dcmplx *K1V,dcmplx *K1W,dcmplx *K2U,
                dcmplx *K2V,dcmplx *K2W);
void get_comp3_(int NGL,const double *A,const double *C,const double *F,
                const double *L,const double *N,const double *theta0,
                const double *dphi,double jaco, const double *weight,
                const double *hp,const double *hpT, dcmplx *K0U,dcmplx *K0V,
                dcmplx *K0W,dcmplx *K1U,dcmplx *K1V,dcmplx *K1W,dcmplx *K2U,
                dcmplx *K2V,dcmplx *K2W);
}

/**
 * @brief wrapper function to compute weak form matrices in one element
 * 
 * @param comp = 1,2,3
 * @param NGL 
 * @param A,C,F,L,N,theta0,dphi parameters, shape(NGL)
 * @param jaco jacobian for this element
 * @param weight GLL/GRL weights, shape (NGL)
 * @param hp,hpT GLL/GRL derivatives, shape(NGL,NGL), note it's column major!
 * @param K(0/1/2)/(U/V/W)  matrices related to k^p (p=0,1,2)
 */
static void 
get_comp_wrapper(int comp,
    int NGL,const double *A,const double *C,const double *F,
    const double *L,const double *N,const double *theta0,
    const double *dphi,double jaco, const double *weight,
    const double *hp,const double *hpT, dcmplx *K0U,dcmplx *K0V,
    dcmplx *K0W,dcmplx *K1U,dcmplx *K1V,dcmplx *K1W,dcmplx *K2U,
    dcmplx *K2V,dcmplx *K2W
)
{
    switch (comp + 1)
    {
    case 1:
        get_comp1_(NGL,A,C,F,L,N,theta0,dphi,jaco,weight,hp,hpT,
                   K0U,K0V,K0W,K1U,K1V,K1W,K2U,K2V,K2W);
        break;
    case 2:
        get_comp2_(NGL,A,C,F,L,N,theta0,dphi,jaco,weight,hp,hpT,
                   K0U,K0V,K0W,K1U,K1V,K1W,K2U,K2V,K2W);
        break;
    case 3:
        get_comp3_(NGL,A,C,F,L,N,theta0,dphi,jaco,weight,hp,hpT,
                   K0U,K0V,K0W,K1U,K1V,K1W,K2U,K2V,K2W);
        break;
    default:
        printf("comp should be one of [1,2,3]!\n");
        exit(1);
        break;
    }
}

/**
 * @brief prepare M/K1/K2/E matrices for TTI model
 * 
 * @param phi polar angle of k vector, in deg
 */
void LayerModelTTI :: 
prepare_matrices(double phi)
{
    // allocate matrices
    size_t msize = nglob * 3;
    Mmat.resize(msize);
    K1mat.resize(msize * msize);
    K2mat.resize(msize * msize);
    Emat.resize(msize * msize);
    std::fill(Mmat.begin(),Mmat.end(),(0.,0.));
    std::fill(K1mat.begin(),K1mat.end(),(0.,0.));
    std::fill(K2mat.begin(),K2mat.end(),(0.,0.));
    std::fill(Emat.begin(),Emat.end(),(0.,0.));

    // compute M/K/E for gll/grl layer
    std::array<double,NGRL> dphi;
    std::array<dcmplx,NGRL*NGRL> EU,EV,EW,K1U,K1V,K1W,K2U,K2V,K2W;
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

        // compute temporary polar angle
        for(int i = 0; i < NGL; i ++) {
            dphi[i] = xP[id + i] - M_PI / 180. * phi;
        }

        // get K1/K2/E matrix
        for(int ic = 0; ic < 3; ic ++ ) {
            get_comp_wrapper(ic,NGL,&xA[id],&xC[id],&xF[id],&xL[id],&xN[id],&xT[id],
                            dphi.data(),jaco[ispec],weight,hpT,hp,EU.data(),EV.data(),
                            EW.data(),K1U.data(),K1V.data(),K1W.data(),K2U.data(),
                            K2V.data(),K2W.data());

            // if(ic ==0 && ispec == 0) {
            //     printf("Eu = \n");
            //     for(int i = 0; i < NGL; i ++) {
            //         for(int j = 0; j < NGL; j ++) {
            //             dcmplx a = K1U[i * NGL + j];
            //             printf("%g + %gI, ",a.real(),a.imag());
            //         }
            //         printf("\n");
            //     }

            //     printf("EW = \n");
            //     for(int i = 0; i < NGL; i ++) {
            //         for(int j = 0; j < NGL; j ++) {
            //             dcmplx a = K1W[i * NGL + j];
            //             printf("%g + %gI, ",a.real(),a.imag());
            //         }
            //         printf("\n");
            //     }
            //     printf("EV = \n");
            //     for(int i = 0; i < NGL; i ++) {
            //         for(int j = 0; j < NGL; j ++) {
            //             dcmplx a = K1V[i * NGL + j];
            //             printf("%g + %gI, ",a.real(),a.imag());
            //         }
            //         printf("\n");
            //     }

            //     exit(1);
            // }
            
            for(int i = 0; i < NGL; i ++) {
                int iglob = ibool[id + i];
                
                // mass matrix
                double M0 = weight[i] * jaco[ispec] * xrho[id + i];

                // assemble
                Mmat[ic * nglob + iglob] += M0;

                for(int j = 0; j < NGL; j ++) {
                    int iglob1 = ibool[id + j];
                    
                    size_t row = (ic * nglob + iglob) * msize;
                    Emat[row + nglob * 0 + iglob1] += EU[i * NGL + j]; 
                    Emat[row + nglob * 1 + iglob1] += EW[i * NGL + j]; 
                    Emat[row + nglob * 2 + iglob1] += EV[i * NGL + j]; 
                    K1mat[row + nglob * 0 + iglob1] += K1U[i * NGL + j]; 
                    K1mat[row + nglob * 1 + iglob1] += K1W[i * NGL + j]; 
                    K1mat[row + nglob * 2 + iglob1] += K1V[i * NGL + j]; 
                    K2mat[row + nglob * 0 + iglob1] += K2U[i * NGL + j]; 
                    K2mat[row + nglob * 1 + iglob1] += K2W[i * NGL + j]; 
                    K2mat[row + nglob * 2 + iglob1] += K2V[i * NGL + j]; 
                }
                
            }
        }
    }
}