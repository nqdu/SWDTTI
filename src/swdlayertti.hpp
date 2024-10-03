#ifndef SWD_LAYER_TTI_MODEL
#define SWD_LAYER_TTI_MODEL

#include <complex>
#include <vector>
#include <array>

#include "swdlayer.hpp"

class LayerModelTTI : public LayerModel{

typedef std::complex<double> dcmplx;
public:

    LayerModelTTI(){};

private:
    std::vector<dcmplx> Mmat,Emat,K1mat,K2mat; // matrices for SEM,shape(3*nglob,3*nglob) om^2 M = k^2 K_2 + k K_1 + E

public:

    // density
    std::vector<double> xrho;

    // tti Love parameters A,C,L,F,N, theta,phi
    std::vector<double> xA,xC,xL,xF,xN; // shape(nspec * NGLL + NGRL)
    std::vector<double> xT,xP; // theta/phi, shape(nspec *NGLL + NGRL), in rad


    // VTI model
    void create_database(double freq,int nlayer, const float *vph, const float* vpv,
                        const float *vsh, const float *vsv, const float *eta,
                        const float *theta0, const float *phi0,
                        const float *rho,const float *thk,bool is_layer);

    void prepare_matrices(double phi);

    void compute_egnfun(double freq, double phi, std::vector<double> &c, std::vector<dcmplx> &displ) const;
    std::array<double,2>
    compute_kernels(double freq, double c,double phi,
                    const std::vector<dcmplx> &displ,
                    std::vector<double> &frekl) const;
    
    void transform_kernels(std::vector<double> &frekl) const;
};

#endif