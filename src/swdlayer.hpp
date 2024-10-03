#ifndef SWDLAYER_MODEL
#define SWDLAYER_MODEL

#include <complex>
#include <vector>
#include <array>

class LayerModel {

public:
    // GLL/GRL nodes and weights
    static const int NGLL = 7, NGRL = 20;
    std::array<double,NGLL> xgll,wgll;
    std::array<double,NGRL> xgrl,wgrl;
    std::array<double,NGLL*NGLL> hprimeT,hprime; // hprimeT(i,j) = l'_i(xi_j)
    std::array<double,NGRL*NGRL> hprimeT_grl,hprime_grl;

private:
    void initialize_nodes();

public:
    // SEM Mesh
    int nspec,nspec_grl; // # of elements for gll/grl layer
    int nglob; // # of unique points
    std::vector<int> ibool; // connectivity matrix, shape(nspec * NGLL + NGRL)
    std::vector<float> skel;  // skeleton, shape(nspec * 2 + 2)
    std::vector<double> znodes; // shape(nspec * NGLL + NGRL)
    std::vector<double> jaco; // jacobian for GLL, shape(nspec + 1) dz / dxi
    std::vector<double> zstore; // shape(nglob)  

public:
    bool IS_DICON_MODEL;
    std::vector<int> ilayer_flag; // shape(nspec + 1), return layer flag 
    double PHASE_VELOC_MIN,PHASE_VELOC_MAX;

//functions
public:
    LayerModel(){};
    void initialize();
    void interp_model(const float *z,const float *param,std::vector<double> &md) const;
    void create_mesh(const int *nel, const float *thk,const float *zlist,int nlayer,double scale);
    void project_kl(const float *z, const double *param_kl, double *kl_out) const;
};

#endif