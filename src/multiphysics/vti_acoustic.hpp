#ifndef SWD_LAYER_MULPHY_VTI_MODEL
#define SWD_LAYER_MULPHY_VTI_MODEL

#include <complex>
#include <vector>
#include <array>

#include "swdlayer.hpp"

class LayerModelMultiPhyVTI: public LayerModel {

public:

   LayerModelMultiPhyVTI(){};

private:
    std::vector<double> Mmat,Emat,Kmat; // matrices for SEM, om^2 M = k^2 K + E

public:
    // element type
    int nspec_ac,nspec_el;
    int nspec_ac_grl,nspec_el_grl;
    std::vector<char> is_elastic;
    std::vector<int> el_elmnts,ac_elmnts; // elements for each media, shape(nspec_? + nspec_?_grl)

    // unique array for acoustic/elastic
    int nglob_ac, nglob_el;
    std::vector<int> ibool_el, ibool_ac; // connectivity matrix, shape shape(nspec_? + nspec_?_grl)

    // density and elastic parameters
    std::vector<double> xrho_ac,xkappa_ac; // shape(nspec_ac * NGLL + nspec_ac_grl * NGRL)
    std::vector<double> xrho_el; // shape (nsepc_el * NGLL + nspec_el_grl * NGRL)
    std::vector<double> xA,xC,xL,xF; // shape(nspec_el * NGLL+ nspec_el_grl * NGRL)

    // acoustic-elastic interface
    int nfaces_bdry;
    std::vector<int> ispec_bdry_loc; // shape(nfaces_bdry,2) (i,:) = [ispec_ac,ispec_el]
    std::vector<char> is_top_ac_bdry; //if the   shape(nfaces_bdry)

private:
    std::vector<char> is_el_layer; // shape(nlayer)

public:

    // VTI model
    void create_database(double freq,int nlayer, const float *rho,
                        const float *vpv, const float* vph,
                        const float *vsv, const float *vsh, const float *eta,
                        const float *thk,bool is_layer);

    void prepare_matrices(double freq);

    void compute_egnfun(double freq, std::vector<double> &c, std::vector<double> &egn) const;
    double compute_kernels(double freq, double c,const double *egn,
                        std::vector<double> &frekl_el,std::vector<double> &frekl_ac) const;
    
    void transform_kernels(std::vector<double> &frekl) const;
    void transform_ac_egnfun();
    void interp_model(const float *z,const float *param,bool elastic,std::vector<double> &md) const;

private:
    void create_medium_info();
};

#endif