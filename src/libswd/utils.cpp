#include "swdlayertti.hpp"
#include "swdlayervti.hpp"

#include <iostream>
#include <vector>

/**
 * @brief compute vti surface wave dispersion
 * 
 * @param thkm thickness of each layer, shape(nlayer)
 * @param vphm/vpvm/vshm/vsvm/etam vti 5 parameters, shape(nlayer)
 * @param rhom density, shape(nlayer)
 * @param nlayer no. of layers
 * @param nt no. of period points used
 * @param mode which mode you want, start from 0 = fundamental
 * @param wavetype =1 for love and =2 for rayleigh
 * @param only_phase = false only for phase, or for both phase/group
 * @param t period vector
 * @param cp/cg output dispersion,phase/group
 */
void 
__surfvti(const float *thkm,const float *vphm,const float *vpvm,
            const float *vshm,const float *vsvm,const float *etam,
            const float *rhom,int nlayer,int nt,int mode, int wavetype,
            bool only_phase,const double *t,double *cp,double *cg)
{
    LayerModelVTI model;
    model.initialize();
    int nglob = model.nglob;

    // allocate space for displ/velocity
    std::vector<double> c,displ,u,frekl;
    int eig_size = nglob * wavetype;

    // loop every period 
    for(int it = 0; it < nt; it ++) {
        double freq = 1. / t[it];
        model.create_database(freq,nlayer,vphm,vpvm,
                                vshm,vsvm,etam,rhom,thkm,
                                true);
        model.prepare_matrices(wavetype);
        switch (wavetype)
        {
        case 1:
            model.compute_slegn(freq,c,displ);
            break;
        
        default:
            model.compute_sregn(freq,c,displ);
            break;
        }

        // cache the one you want
        cp[it] = c[mode];
        if(!only_phase) return;

        // compute group velocity
        if(wavetype == 1) {
            cg[it] = model.compute_love_kl(freq,c[mode],&displ[mode * eig_size],frekl);
        }
        else {
            cg[it] = model.compute_rayl_kl(freq,c[mode],&displ[mode * eig_size],frekl);
        }
    }
}

/**
 * @brief compute vti surface wave dispersion
 * 
 * @param thkm thickness of each layer, shape(nlayer)
 * @param vphm/vpvm/vshm/vsvm/etam vti 5 parameters, shape(nlayer)
 * @param rhom density, shape(nlayer)
 * @param nlayer no. of layers
 * @param nt no. of period points used
 * @param mode which mode you want, start from 0 = fundamental
 * @param wavetype =1 for love and =2 for rayleigh
 * @param t period vector
 * @param cp/cg output dispersion,phase/group
 * @returns kernels: shape(nt,nker,nlayer) sensitivity kernels
 */
std::vector<double> 
__surfvti_kl(const float *thkm,const float *vphm,const float *vpvm,
            const float *vshm,const float *vsvm,const float *etam,
            const float *rhom,int nlayer,int nt,int mode, int wavetype,
            const double *t,double *cp,double *cg)
{
    LayerModelVTI model;
    model.initialize();
    int nglob = model.nglob;

    // allocate space for displ/velocity
    std::vector<double> c,displ,u,frekl;
    int eig_size = nglob * wavetype;

    // allocate kernels
    int size = model.ibool.size();
    int nker = frekl.size() / (model.ibool.size());
    std::vector<double> kernels(nker * nlayer * nt);

    // z coordinates
    std::vector<float> zlist(nlayer);
    zlist[0] = 0.;
    for(int i = 0; i < nlayer - 1; i ++) {
        zlist[i + 1] = zlist[i] + thkm[i];
    }

    // loop every period 
    for(int it = 0; it < nt; it ++) {
        double freq = 1. / t[it];
        model.create_database(freq,nlayer,vphm,vpvm,
                                vshm,vsvm,etam,rhom,thkm,
                                true);
        model.prepare_matrices(wavetype);
        switch (wavetype)
        {
        case 1:
            model.compute_slegn(freq,c,displ);
            break;
        
        default:
            model.compute_sregn(freq,c,displ);
            break;
        }

        // cache the one you want
        cp[it] = c[mode];

        // compute group velocity
        if(wavetype == 1) {
            cg[it] = model.compute_love_kl(freq,c[mode],&displ[mode * eig_size],frekl);
        }
        else {
            cg[it] = model.compute_rayl_kl(freq,c[mode],&displ[mode * eig_size],frekl);
        }
        model.transform_kernels(frekl);

        // project back to original model
        for(int iker = 0; iker < nker; iker ++) {
            model.project_kl(zlist.data(),&frekl[iker * size],&kernels[it * nlayer * nker + iker * nlayer]);
        }
    }

    return kernels;

}