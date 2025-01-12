#include "swdlayertti.hpp"
typedef std::complex<double> dcmplx;

/**
 * @brief initalize SEM mesh and create a TTI database from a layered model
 * 
 * @param freq current frequency
 * @param nlayer # of nlayers
 * @param vpv/vph/vsv/vsh/eta/rho layer model vti parameters, shape(nlayer) 
 * @param theta0/phi0 axis direction
 */
void LayerModelTTI:: 
create_database(double freq,int nlayer, const float *rho,
                const float *vpv, const float* vph,
                const float *vsv, const float *vsh, const float *eta,
                const float *theta0, const float *phi0,
                const float *thk,bool is_layer)
{
    // find min/max thickness, min veloc and create depth list
    std::vector<float> zlist(nlayer);
    zlist[0] = 0.;
    for(int i = 0; i < nlayer - 1; i ++) {
        zlist[i + 1] = zlist[i] + thk[i];
    }

    // find min/max vs
    PHASE_VELOC_MAX = -1.;
    PHASE_VELOC_MIN = 1.0e20;
    for(int i = 0; i < nlayer; i ++) {
        double vmin = std::min(vsv[i],vsh[i]);
        double vmax = std::max(vsv[i],vsh[i]);
        PHASE_VELOC_MAX = std::max(vmax,PHASE_VELOC_MAX);
        PHASE_VELOC_MIN = std::min(vmin,PHASE_VELOC_MIN);
    }
    PHASE_VELOC_MIN *= 0.85;

    // copy is_layer
    IS_DICON_MODEL = is_layer;

    //double scale =  Tmin /(2.0 * M_PI * std::sqrt(1. / (4.99 * 4.99) - 1 / (vmax * vmax)));
    double scale = PHASE_VELOC_MAX / freq / xgrl[NGRL-1] * 10;

    // determine no. elements in each layer
    std::vector<int> nel;
    if(this -> IS_DICON_MODEL) {
        nel.resize(nlayer - 1);
        for(int i = 0; i < nlayer - 1; i ++) {
            float v0 = std::min(vsv[i],vsh[i]) * 0.85;
            nel[i] = (int)(thk[i] * freq / v0 * 1.5 ) + 1;
        }
    }
    else { // continuous model, constructed with min velocity
        nel.resize(1);
        float maxdepth = zlist[nlayer - 1] - zlist[0];
        nel[0] = 1.5 * (maxdepth * freq / PHASE_VELOC_MIN) + 1;
    }

    // create mesh
    this -> create_mesh(nel.data(),thk,zlist.data(),nlayer,scale);

    // allocate space
    size_t size = nspec * NGLL + NGRL;

    // set value to modulus
    xrho.resize(size); xA.resize(size); xC.resize(size);
    xL.resize(size); xN.resize(size); xF.resize(size);
    xT.resize(size); xP.resize(size);

    // interpolate
    std::vector<double> temp_eta(size);
    this -> interp_model(zlist.data(),rho,xrho);
    this -> interp_model(zlist.data(),vph,xA);
    this -> interp_model(zlist.data(),vpv,xC);
    this -> interp_model(zlist.data(),vsh,xN);
    this -> interp_model(zlist.data(),vsv,xL);
    this -> interp_model(zlist.data(),eta,temp_eta);
    this -> interp_model(zlist.data(),theta0,xT);
    this -> interp_model(zlist.data(),phi0,xP);

    // convert to elastic parameters
    for(int i = 0; i < nspec * NGLL + NGRL; i ++) {
        xA[i] = xA[i] * xA[i] * xrho[i];
        xC[i] = xC[i] * xC[i] * xrho[i];
        xN[i] = xN[i] * xN[i] * xrho[i];
        xL[i] = xL[i] * xL[i] * xrho[i];
        xF[i] = temp_eta[i] * (xA[i] - 2. * xL[i]);
    }
}



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
    const dcmplx ZERO{0.,0.};
    std::fill(Mmat.begin(),Mmat.end(),ZERO);
    std::fill(K1mat.begin(),K1mat.end(),ZERO);
    std::fill(K2mat.begin(),K2mat.end(),ZERO);
    std::fill(Emat.begin(),Emat.end(),ZERO);

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