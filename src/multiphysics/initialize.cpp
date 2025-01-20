#include "multiphysics/vti_acoustic.hpp"

#include <algorithm>
#include <iostream>

static bool 
find_media_blocks(const char *arr, int n, std::vector<int>& blockIndex) {
    int start = 0;
    bool inOneSubarray = false, inZeroSubarray = false;
    std::vector<char> flags; flags.reserve(n);
    blockIndex.reserve(n * 2);

    for (int i = 0; i < n; ++i) {
        if (arr[i] == 1) {
            if (!inOneSubarray) {
                start = i;
                inOneSubarray = true;
            }
            inZeroSubarray = false;
        } else if (arr[i] == 0) {
            if (!inZeroSubarray) {
                start = i;
                inZeroSubarray = true;
            }
            inOneSubarray = false;
        }

        if (i == n - 1 || arr[i] != arr[i + 1]) {
            if (inOneSubarray) {
                blockIndex.push_back(start);
                blockIndex.push_back(i);
                flags.push_back(true);  // True for 1-subarray
                inOneSubarray = false;
            } else if (inZeroSubarray) {
                blockIndex.push_back(start);
                blockIndex.push_back(i);
                flags.push_back(false);  // False for 0-subarray
                inZeroSubarray = false;
            }
        }
    }

    return flags[0];
}



/**
 * @brief Create layered model for multiphysics: acoustic + vti elastic 
 * 
 * @param freq frequency used
 * @param nlayer # of layers used
 * @param rho,vpv,vph,vsv,vsh,eta model parameters. shape(nlayer) 
 * @param thk thickness of each layer, shape(nlayer)
 * @param is_layer the input model is a layered model
 */
void LayerModelMultiPhyVTI:: 
create_database(double freq,int nlayer, const float *rho,
                const float *vpv, const float* vph,
                const float *vsv, const float *vsh, const float *eta,
                const float *thk,bool is_layer)
{
    // get zlist
    std::vector<float> zlist(nlayer);
    zlist[0] = 0.;
    for(int i = 0; i < nlayer - 1; i ++) {
        zlist[i + 1] = zlist[i] + thk[i];
    }

    // find # of layers in the model, and detect which one is acoustic
    float tol = 1.0e-4;
    is_el_layer.resize(nlayer);
    int nlayer_el = 0, nlayer_ac = 0;
    for(int i = 0; i < nlayer; i ++) {
        if(vsh[i] * vsv[i] < tol) {
            is_el_layer[i] = 0;
            nlayer_ac += 1;

            // check vpv == vph
            if(std::abs(vpv[i] - vph[i]) > tol) {
                printf("You should make sure vpv == vph and vsv=vsh=0 in acoustic layers!\n");
                printf("vpv,vph = %f %f\n",vpv[i],vph[i]);
                exit(1);
            }
        }
        else {
            is_el_layer[i] = 1;
            nlayer_el += 1;
        }
    }

    // find min/max vs
    PHASE_VELOC_MAX = -1.;
    PHASE_VELOC_MIN = 1.0e20;
    for(int i = 0; i < nlayer; i ++) {
        double vmin = std::min(vsv[i],vsh[i]);
        double vmax = std::max(vsv[i],vsh[i]);
        if(!is_el_layer[i]) {
            vmin = vpv[i];
            vmax = vmin;
        }
        PHASE_VELOC_MAX = std::max(vmax,PHASE_VELOC_MAX);
        PHASE_VELOC_MIN = std::min(vmin,PHASE_VELOC_MIN);
    }
    PHASE_VELOC_MIN *= 0.85;

    // copy is_layer
    this -> IS_DICON_MODEL = is_layer;

    // half space
    nspec_grl = 1;
    nspec_ac_grl = 0; nspec_el_grl = 0;
    if(!is_el_layer[nlayer-1]) {
        nspec_ac_grl = 1;
    }
    else {
        nspec_el_grl = 1;
    }

    // determine no. of elements in each layers
    std::vector<int> nel;
    std::vector<int> mediablock;
    if(this -> IS_DICON_MODEL) {
        nspec = 0;
        nspec_ac = 0; 
        nspec_el = 0;
        nel.resize(nlayer - 1);
        for(int i = 0; i < nlayer - 1; i ++) {
            if(is_el_layer[i]) {
                float v0 = std::min(vsv[i],vsh[i]) * 0.85;
                nel[i] = thk[i] * freq / v0 + 1;
                nspec_el += nel[i];
            }
            else {
                float v0 = vpv[i] * 0.85;
                nel[i] = thk[i] * freq / v0 + 1;
                nspec_ac += nel[i];
            }
            nspec += nel[i];
        }

        // resize
        is_elastic.resize(0);
        is_elastic.reserve(nspec + nspec_grl);
        ilayer_flag.resize(nspec + nspec_grl);
        for(int i = 0; i < nlayer - 1; i ++) {
        for(int j = 0; j < nel[i]; j ++) {
            is_elastic.push_back(is_el_layer[i]);
        }}
        is_elastic.push_back(is_el_layer[nlayer-1]);

        // compute skeleton coordinates
        skel.resize(nspec * 2 + 2);
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
    }
    else {
        bool flag = find_media_blocks(is_el_layer.data(),nlayer - 1,mediablock);
        int nb = mediablock.size() / 2;
        nel.resize(nb);
        nspec = 0 ;
        for(int ib = 0; ib < nb; ib ++) {
            int i1 = mediablock[ib * 2 + 0];
            int i2 = mediablock[ib * 2 + 1];
            double maxd = zlist[i2 + 1] - zlist[i1];
            float vmin = 1.0e20;
            for(int j = i1; j <= i2; j++) {
                if(is_el_layer[j]) {
                    vmin = std::min(vmin,std::min(vsv[j],vsh[j]));
                }
                else {
                    vmin = std::min(vmin,vpv[j]);
                }
            }
            nel[ib] = 1.2 * (maxd * freq / vmin + 1);
            nspec += nel[ib];
        }
        is_elastic.resize(0);
        is_elastic.reserve(nspec + nspec_grl);
        ilayer_flag.resize(nspec + nspec_grl);
        for(int ib = 0; ib < nb; ib ++) {
            for(int i = 0; i < nel[ib]; i ++) {
                is_elastic.push_back(flag);
            }
            flag = ! flag;
        }
        is_elastic.push_back(is_el_layer[nlayer-1]);

        // compute skeleton coordinates
        skel.resize(nspec * 2 + 2);
        int id = 0;
        nspec_ac = 0; nspec_el = 0;
        for(int ib = 0; ib < nb; ib ++) {
            int i1 = mediablock[ib * 2 + 0];
            int i2 = mediablock[ib * 2 + 1];
            double h = (zlist[i2 + 1] - zlist[i1]) / nel[ib];
            for(int j = 0; j < nel[ib]; j ++ ) {
                skel[id * 2 + 0] = zlist[i1] + h  * j;
                skel[id * 2 + 1] = zlist[i1] + h  * (j+1);
                ilayer_flag[id] = i2;
                if(is_elastic[id]) {
                    nspec_el += 1;
                }
                else {
                    nspec_ac += 1;
                }
                id += 1;
            }
        }
        ilayer_flag[nspec] = nlayer - 1;
    }

    //scale factor for last layer
    double scale = PHASE_VELOC_MAX / freq / xgrl[NGRL-1] * 5;  // up to 5 wavelength

    // allocate space
    size_t size = nspec * NGLL + NGRL;
    ibool.resize(size); znodes.resize(size);
    jaco.resize(nspec + 1);
    ilayer_flag.resize(nspec + 1);

    // allocate space for parameters
    size_t size_e = nspec_el * NGLL + nspec_el_grl * NGRL;
    size_t size_a = nspec_ac * NGLL + nspec_ac_grl * NGRL;
    xrho_ac.resize(size_a); xkappa_ac.resize(size_a);
    xrho_el.resize(size_e); xA.resize(size_e);
    xC.resize(size_e); xL.resize(size_e);
    xF.resize(size_e);

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

    // min/max z for GRL layer
    float zmax = zlist[nlayer - 1];
    skel[nspec * 2 + 0] = zmax;
    skel[nspec * 2 + 1] = zmax + xgrl[NGRL-1] * scale;

    // compute coordinates and jaco in GLL layer
    for(int ispec = 0; ispec < nspec; ispec ++) {
        double h = skel[ispec * 2 + 1] - skel[ispec * 2 + 0];
        jaco[ispec] = h / 2.;
        for(int i = 0; i < NGLL; i ++) {
            double xi = xgll[i];
            znodes[ispec * NGLL + i] = skel[ispec * 2] + h * 0.5 * (xi + 1);
        }
    }

    // compute coordinates and jaco in GRL layer
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

    // create medium info ibool_*, 
    this -> create_medium_info();

    // interpolate
    std::vector<double> temp_eta(xrho_el.size());
    this -> interp_model(zlist.data(),rho,true,xrho_el);
    this -> interp_model(zlist.data(),vph,true,xA);
    this -> interp_model(zlist.data(),vpv,true,xC);
    this -> interp_model(zlist.data(),vsv,true,xL);
    this -> interp_model(zlist.data(),eta,true,temp_eta);
    this -> interp_model(zlist.data(),rho,false,xrho_ac);
    this -> interp_model(zlist.data(),vph,false,xkappa_ac);

    // check model
#ifdef SPECSWD_DEBUG
    FILE *fp = fopen("elastic.txt","w");
    for(int ispec = 0; ispec < nspec_el + nspec_el_grl; ispec ++ ) {
        int iel = el_elmnts[ispec];
        int NGL = NGLL;
        if(ispec == nspec_el) NGL = NGRL;
        for(int i = 0; i < NGL; i ++) {
            fprintf(fp,"ispec z vph = %d %f %f\n",ispec,znodes[iel * NGLL + i], xA[ispec * NGLL + i]);
        }
    }
    fclose(fp);
    fp = fopen("acoustic.txt","w");
    for(int ispec = 0; ispec < nspec_ac + nspec_ac_grl; ispec ++ ) {
        int iel = ac_elmnts[ispec];
        int NGL = NGLL;
        if(ispec == nspec_ac) NGL = NGRL;
        for(int i = 0; i < NGL; i ++) {
            fprintf(fp,"ispec z vph = %d %f %f\n",ispec,znodes[iel * NGLL + i], xkappa_ac[ispec * NGLL + i]);
        }
    }
    fclose(fp);
#endif
    
    // convert to elastic moduli
    for(int i = 0; i < (int)xrho_el.size(); i ++) {
        xA[i] = xA[i] * xA[i] * xrho_el[i];
        xC[i] = xC[i] * xC[i] * xrho_el[i];
        xL[i] = xL[i] * xL[i] * xrho_el[i];
        xF[i] = temp_eta[i] * (xA[i] - 2. * xL[i]);
    }

    // convert to acoustic moduli
    for(int i = 0; i < (int)xrho_ac.size(); i ++) {
        xkappa_ac[i] = xkappa_ac[i] * xkappa_ac[i] * xrho_ac[i];
    }

#ifdef SPECSWD_DEBUG
    // interpolate
    printf("# of elements GLL = %d, nspec_ac = %d nspec_el = %d\n",nspec,nspec_ac,nspec_el);
    printf("# of elements GRL = %d, nspec_ac_grl = %d nspec_el_grl = %d\n",nspec_grl,nspec_ac_grl,nspec_el_grl);
    printf("nglob_el nglob_ac = %d %d\n",nglob_el,nglob_ac);
#endif
}

/**
 * @brief find location of z0 in ascending list z
 * 
 * @param z depth list, shape(nlayer)
 * @param z0 current loc, must be inside z 
 * @param nlayer 
 * @return int location of z0 in z, satisfies  z0 >= z[i] && z0 < z[i + 1]
 */
static int 
find_loc(const float *z,float z0,int nz) 
{

    int i = 0;
    while(i < nz) {
        if(z0 >= z[i] && z0 < z[i + 1]) {
            break;
        }
        i += 1;
    }

    return i;
}


/**
 * @brief interpolate a model by using coordinates
 * 
 * @param z input model z coordinates, shape(nlayer)
 * @param param input model parameter, shape(nlayer)
 * @param md model required to interpolate, shape(nspec*NGLL + NGRL)
 */
void LayerModelMultiPhyVTI::
interp_model(const float *z,const float *param,bool elastic,std::vector<double> &md) const
{
    // get ilayer
    int nlay = ilayer_flag[nspec] + 1;
    int nspec_m = nspec_el;
    int nspec_m_grl = nspec_el_grl;
    if(!elastic) {
        nspec_m =  nspec_ac;
        nspec_m_grl = nspec_ac_grl;
    }
    for(int iel = 0; iel < nspec_m + nspec_m_grl ; iel ++) {
        int NGL = NGLL;
        if(iel == nspec_m) {
            NGL = NGRL;
        }

        // get layer index of this element
        int ispec;
        if(elastic) {
            ispec = el_elmnts[iel];
        }
        else{
            ispec = ac_elmnts[iel];
        }
        int ilay0 = ilayer_flag[ispec];

        // loop gll/grl nodes to interpolate
        for(int i = 0; i < NGL; i ++) {
            int id = iel * NGLL + i;
            double z0 = znodes[id];

            // set model
            if(IS_DICON_MODEL) {
                md[id] = param[ilay0];
            }
            else {
                if(z0 >= z[nlay - 1]) {
                    md[id] = param[nlay - 1];
                }
                else if (z0 <= z[0]) {
                     md[id] = param[0];
                }
                else {
                    // find location
                    int ilay = 0;
                    while(ilay < ilay0) {
                        if(is_el_layer[ilay] != elastic ) {
                            ilay += 1;
                            continue;
                        }
                        if(z0 >= z[ilay] && z0 < z[ilay + 1]) {
                            break;
                        }
                        ilay += 1;
                    }
                    float dz = z[ilay + 1] - z[ilay];
                    md[id] = param[ilay] + (param[ilay + 1] - param[ilay]) / dz * (z0 - z[ilay]);
                }
            }
        }
    }
}

/**
 * @brief Create medium information, like ibool_*, ispec_*, bdry_ac_el
 * 
 */
void LayerModelMultiPhyVTI:: 
create_medium_info()
{   
    ac_elmnts.resize(0); el_elmnts.resize(0);
    ac_elmnts.reserve(nspec_ac + nspec_ac_grl);
    el_elmnts.reserve(nspec_el + nspec_el_grl);
    for(int ispec = 0; ispec < nspec + nspec_grl; ispec ++) {
        if(is_elastic[ispec]) {
            el_elmnts.push_back(ispec);
        }
        else {
            ac_elmnts.push_back(ispec);
        }
    }

    // get nglob_el for elastic 
    ibool_el.resize(nspec_el * NGLL + nspec_el_grl * NGRL);
    nglob_el = 0;
    int idx = -1;
    for(int i = 0; i < nspec_el + nspec_el_grl; i += 1) {
        int ispec = el_elmnts[i];
        if(idx == ibool[ispec * NGLL]) nglob_el -= 1;

        int NGL = NGLL;
        if(i == nspec_el) NGL = NGRL;
        for(int igll = 0; igll < NGL; igll ++) {
            ibool_el[i * NGLL + igll] = nglob_el;
            nglob_el += 1;
        }
        idx = ibool[ispec * NGLL + NGLL-1];
    }

    // get nglob_el for  acoustic
    ibool_ac.resize(nspec_ac * NGLL + nspec_ac_grl * NGRL);
    idx = -10;
    nglob_ac = 0;
    if(!is_el_layer[0]) nglob_ac = -1; // the top point of acoustic wave is 0
    for(int i = 0; i < nspec_ac + nspec_ac_grl; i += 1) {
        int ispec = ac_elmnts[i];
        if(idx == ibool[ispec * NGLL]) nglob_ac -= 1;

        int NGL = NGLL;
        if(i == nspec_ac) NGL = NGRL;
        for(int igll = 0; igll < NGL; igll ++) {
            ibool_ac[i * NGLL + igll] = nglob_ac;
            nglob_ac += 1;
        }
        idx = ibool[ispec * NGLL + NGLL-1];
    }

    // find elastic/acoustic bdry
    nfaces_bdry = 0;
    for(int i = 0; i < nspec; i ++) {
        if(is_elastic[i] != is_elastic[i+1]) {
            nfaces_bdry += 1;
        }
    }
    ispec_bdry_loc.resize(nfaces_bdry,2);
    is_top_ac_bdry.resize(nfaces_bdry);
    idx = 0;
    for(int i = 0; i < nspec; i ++) {
        if(is_elastic[i] != is_elastic[i+1]) {
            int iloc_el = i + 1, iloc_ac = i;
            if(is_elastic[i]) {
                is_top_ac_bdry[idx] = 0;
                iloc_el = i;
                iloc_ac = i + 1;
            }
            else {
                is_top_ac_bdry[idx] = 1;
            }
            auto it = std::find(el_elmnts.begin(),el_elmnts.end(),iloc_el);
            int ispec_el = it - el_elmnts.begin();
            it = std::find(ac_elmnts.begin(),ac_elmnts.end(),iloc_ac);
            int ispec_ac = it - ac_elmnts.begin();
            ispec_bdry_loc[idx * 2 + 0] = ispec_ac;
            ispec_bdry_loc[idx * 2 + 1] = ispec_el;
            idx += 1;
        }
    }

#ifdef SPECSWD_DEBUG
    // debug
    for(int iface = 0; iface < nfaces_bdry; iface ++) {
        printf("\nface %d, ispec_ac,ispec_el = %d %d\n",iface,ispec_bdry_loc[iface*2],ispec_bdry_loc[iface*2+1]);
        printf("top is acoustic = %d\n",is_top_ac_bdry[iface]);
    }
#endif
}