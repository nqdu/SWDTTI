#ifndef SWDLAYER_VTI_MODEL
#define SWDLAYER_VTI_MODEL

#include "swdlayer.hpp"

class LayerModelVTI: public LayerModel  {

public:
    LayerModelVTI(){};

private:
    std::vector<double> Mmat,Emat,Kmat; // matrices for SEM, om^2 M = k^2 K + E

public:

    // density
    std::vector<double> xrho;

    // vti Love parameters
    std::vector<double> xA,xC,xL,xF,xN; // shape(nspec * NGLL + NGRL)

    // VTI model
    void create_database(double freq,int nlayer, const float *rho,
                        const float *vpv, const float* vph,
                        const float *vsv, const float *vsh, const float *eta,
                        const float *thk, bool is_layer);

    void prepare_matrices(int wavetype);

    void compute_slegn(double freq,std::vector<double> &c,
                         std::vector<double> &displ) const;
    void compute_sregn(double freq,std::vector<double> &c,
                        std::vector<double> &displ) const;

    double compute_love_kl(double freq,double c,const double *displ, std::vector<double> &frekl) const;
    double compute_rayl_kl(double freq,double c,const double *displ, std::vector<double> &frekl) const;
    
    void transform_kernels(std::vector<double> &frekl) const;

private:
    void prepare_matrices_love();
    void prepare_matrices_rayl();
};

#endif
