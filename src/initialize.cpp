#include "swdlayertti.hpp"
#include "swdlayervti.hpp"
#include "quadrature.hpp"

/**
 * @brief initialize GLL/GRL nodes/weights
 * 
 */
void LayerModel:: initialize()
{   
    // GLL nodes/weights
    gauss_legendre_lobatto(xgll.data(),wgll.data(),NGLL);

    // compute hprime and hprimeT
    double poly[NGLL];
    for(int i = 0; i < NGLL; i ++) {
        double xi = xgll[i];
        lagrange_poly(xi,NGLL,xgll.data(),poly,&hprime[i*NGLL]);
    }
    for(int i = 0; i < NGLL; i ++) {
    for(int j = 0; j < NGLL; j ++) {
        hprimeT[i * NGLL + j] = hprime[j * NGLL + i];
    }}


    // GRL nodes/weights
    gauss_radau_laguerre(xgrl.data(),wgrl.data(),NGRL);

    // compute hprime_grl and hprimeT_grl
    for(int i = 0; i < NGRL; i ++) {
    for(int j = 0; j < NGRL; j ++) {
        int id = i * NGRL + j;
        if(i != j) {
            hprimeT_grl[id] = laguerre_func(NGRL,xgrl[j]) / laguerre_func(NGRL,xgrl[i]) / (xgrl[j] - xgrl[i]);
        }
        else if(i == j && i == 0) {
            hprimeT_grl[id] = -NGRL / 2.;
        }
        else {
            hprimeT_grl[id] = 0.;
        }
    }} 
    for(int i = 0; i < NGRL; i ++) {
    for(int j = 0; j < NGRL; j ++) {
        int id = i * NGRL + j;
        hprime_grl[j * NGRL + i] = hprimeT_grl[id];
    }}
}


/**
 * @brief interpolate a model by using coordinates
 * 
 * @param z input model z coordinates, shape(nlayer)
 * @param param input model parameter, shape(nlayer)
 * @param md model required to interpolate, shape(nspec*NGLL + NGRL)
 * @return float 
 */
void LayerModel::
interp_model(const float *z,const float *param,std::vector<double> &md) const
{
    // get ilayer
    int nlay = ilayer_flag[nspec] + 1;

    for(int ispec = 0; ispec < nspec + 1; ispec ++) {
        int NGL = NGLL;
        if(ispec == nspec) {
            NGL = NGRL;
        }

        // get layer index
        int ilay = ilayer_flag[ispec];

        for(int i = 0; i < NGL; i ++) {
            int id = ispec * NGLL + i;
            double z0 = znodes[id];

            // set model
            if(IS_DICON_MODEL) {
                md[id] = param[ilay];
            }
            else {
                if(ilay == nlay - 1) {
                    md[id] = param[nlay - 1];
                }
                else {
                    float dz = z[ilay + 1] - z[ilay];
                    md[id] = param[ilay] + (param[ilay + 1] - param[ilay]) / dz * (z[ilay + 1] - z0);
                }
            }
        }
    }
}

/**
 * @brief Create SEM mesh by using input model info
 * 
 * @param nel no. of elements for each layer, shape(nlayer - 1)
 * @param thk thickness of each layer, shape(nlayer)
 * @param zlist cumsum(thk)
 * @param nlayer no. of layers
 * @param scale scale factor for GRL layer, zbot = sum(thk) + xgrl[-1] * scale
 */
void LayerModel:: 
create_mesh(const int *nel, const float *thk,const float *zlist,int nlayer,double scale)
{
    // count # of elements
    nspec = 0;
    for(int i = 0; i < nlayer - 1; i ++) {
        nspec += nel[i];
    }
    nspec_grl = 1;

    // allocate space
    size_t size = nspec * NGLL + NGRL;
    ibool.resize(size); znodes.resize(size);
    jaco.resize(nspec + 1); skel.resize(nspec * 2 + 2);
    ilayer_flag.resize(nspec + 1);

    // connectivity matrix
    int idx = 0;
    for(int ispec = 0; ispec < nspec; ispec ++) {
        for(int igll = 0; igll < NGLL; igll ++) {
            ibool[ispec * NGLL + igll] = idx;
            idx += 1;
        }
        idx -= 1;
    }
    for(int i = 0; i < NGRL; i ++) {
        ibool[nspec * NGLL + i] = idx;
        idx += 1;
    }
    nglob = ibool[nspec * NGLL + NGRL - 1] + 1;

    // compute skeleton coordinates in GLL layer
    int id = 0;
    for(int i = 0; i < nlayer - 1; i ++) {
        float h = thk[i] / nel[i];
        for(int j = 0; j < nel[i]; j ++ ) {
            skel[id * 2 + 0] = zlist[i] + h  * j;
            skel[id * 2 + 1] = zlist[i] + h  * (j+1);
            ilayer_flag[id] = i;
            id += 1;
        }
    }
    ilayer_flag[nspec] = nlayer - 1;

    // compute coordinates and jaco in GLL layer
    for(int ispec = 0; ispec < nspec; ispec ++) {
        double h = skel[ispec * 2 + 1] - skel[ispec * 2 + 0];
        jaco[ispec] = h / 2.;
        for(int i = 0; i < NGLL; i ++) {
            double xi = xgll[i];
            znodes[ispec * NGLL + i] = skel[ispec * 2] + h * 0.5 * (xi + 1);
        }
    }

    // min/max z for GRL layer
    float zmax = zlist[nlayer - 1];
    skel[nspec * 2 + 0] = zmax;
    skel[nspec * 2 + 1] = zmax + xgrl[NGRL-1] * scale;

    // GRL layer
    for(int ispec = nspec; ispec < nspec + 1; ispec ++) {
        jaco[ispec] = scale;
        for(int i = 0; i < NGRL; i ++) {
            double xi = xgrl[i];
            znodes[ispec * NGLL + i] = skel[ispec * 2] + xi * scale;
        }
    }

    // UNIQUE coordinates
    zstore.resize(nglob);
    for(int i = 0; i < nspec * NGLL + NGRL; i ++) {
        int iglob = ibool[i];
        zstore[iglob] = znodes[i];
    }
}

/**
 * @brief project SEM-type kernels to original model
 * 
 * @param z z coordiantes of previous model, shape(nlayer)
 * @param param_kl SEM typed kernel, shape(nspec * NGLL + NGRL)
 * @param kl_out kernels in original model, shape(nlayer)
 */
void LayerModel :: 
project_kl(const float *z, const double *param_kl, double *kl_out) const
{
    // zero out kl_out 
    int nlayer = ilayer_flag[nspec] + 1;
    for(int i = 0; i < nlayer; i ++) kl_out[i] = 0.;

    for(int ispec = 0; ispec < nspec + 1; ispec ++) {
        int NGL = NGLL;
        if(ispec == nspec) {
            NGL = NGRL;
        }

        // get layer index
        int ilay = ilayer_flag[ispec];

        for(int i = 0; i < NGL; i ++) {
            int id = ispec * NGLL + i;
            double z0 = znodes[id];

            // interp
            if(IS_DICON_MODEL) {
                kl_out[ilay] += param_kl[id];
            }
            else {
                if(ilay == nlayer - 1) {
                    kl_out[ilay] += param_kl[id];
                }
                else {
                    double dz = z[ilay + 1] - z[ilay];
                    double coef = (z[ilay + 1] - z0) / dz;
                    kl_out[ilay] += 1 - coef;
                    kl_out[ilay + 1] += coef;
                }
            }
        }
    }
}