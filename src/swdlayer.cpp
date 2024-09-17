#include "swdlayer.hpp"
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
 * @brief initalize SEM mesh and create a VTI database from a layered model
 * 
 * @param freq current frequency
 * @param nlayer # of nlayers
 * @param vpv/vph/vsv/vsh/eta/rho layer model vti parameters, shape(nlayer) 
 */
void LayerModel:: 
create_database(double freq,int nlayer, const float *vph, const float* vpv,
                const float *vsh, const float *vsv, const float *eta,
                const float *rho,const float *thk)
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

    // determine element size
    std::vector<int> nel(nlayer - 1);
    nspec = 0;
    for(int i = 0; i < nlayer - 1; i ++) {
        float v0 = std::min(vsv[i],vsh[i]);
        nel[i] = int(thk[i] * freq / v0 * 1) + 1;
        nspec += nel[i];
    }
    nspec_grl = 1;
    //printf("# of elements = %d\n",nspec);

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

    // skeleton coordinates in GRL layer
    //double scale =  Tmin /(2.0 * M_PI * std::sqrt(1. / (4.99 * 4.99) - 1 / (vmax * vmax)));
    double v0 = std::max(vsh[nlayer -1],vsv[nlayer - 1]);
    double scale = v0 / freq / xgrl[NGRL-1] * 5;
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
    z.resize(nglob);
    for(int i = 0; i < nspec * NGLL + NGRL; i ++) {
        int iglob = ibool[i];
        z[iglob] = znodes[i];
    }

    // set value to modulus
    xrho.resize(size); xA.resize(size); xC.resize(size);
    xL.resize(size); xN.resize(size); xF.resize(size);

    // loop every elements to set values
    for(int ispec = 0; ispec < nspec + 1; ispec ++) {
        int n = NGLL;
        if(ispec == nspec) n = NGRL;

        // get layer flag
        int ilay = ilayer_flag[ispec];

        // set value
        for (int i = 0; i < n; i++) {
            int id = ispec * NGLL + i;
            xrho[id] = rho[ilay];
            double A = vph[ilay] * vph[ilay] * rho[ilay];
            double C = vpv[ilay] * vpv[ilay] * rho[ilay];
            double N = vsh[ilay] * vsh[ilay] * rho[ilay];
            double L = vsv[ilay] * vsv[ilay] * rho[ilay];
            double F = eta[ilay] * (A - 2. * L);

            xA[id] = A; xC[id] = C; xN[id] = N;
            xL[id] = L; xF[id] = F;
        }
        
    }
}