#ifndef SWD_LAYER_TTI_MODEL
#define SWD_LAYER_TTI_MODEL

#include <complex>
#include <vector>
#include <array>

class LayerModelTTI{

typedef std::complex<double> dcmplx;

private:
    // GLL/GRL nodes and weights
    static const int NGLL = 7, NGRL = 20;
    std::array<double,NGLL> xgll,wgll;
    std::array<double,NGRL> xgrl,wgrl;
    std::array<double,NGLL*NGLL> hprimeT,hprime; // hprimeT(i,j) = l'_i(xi_j)
    std::array<double,NGRL*NGRL> hprimeT_grl,hprime_grl;

    void initialize_nodes();

public:
    // SEM Mesh
    int nspec,nspec_grl; // # of elements for gll/grl layer
    int nglob; // # of unique points
    std::vector<int> ibool; // connectivity matrix, shape(nspec * NGLL + NGRL)
    std::vector<float> skel;  // skeleton, shape(nspec * 2 + 2)
    std::vector<double> znodes; // shape(nspec * NGLL + NGRL)
    std::vector<double> jaco; // jacobian for GLL, shape(nspec + 1) dz / dxi
    std::vector<double> z; // shape(nglob)

    LayerModelTTI(){};
    void initialize();

private:
    std::vector<int> ilayer_flag; // shape(nspec + 1), return layer flag 
    std::vector<dcmplx> Mmat,Emat,K1mat,K2mat; // matrices for SEM,shape(3*nglob,3*nglob) om^2 M = k^2 K_2 + k K_1 + E
    double PHASE_VELOC_MIN,PHASE_VELOC_MAX;

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
                        const float *rho,const float *thk);

    void prepare_matrices(double phi);

    void compute_egnfun(double freq, double phi, std::vector<double> &c, std::vector<dcmplx> &displ) const;
    std::array<double,2>
    compute_kernels(double freq, double c,double phi,
                    const std::vector<dcmplx> &displ,
                    std::vector<double> &frekl) const;
    
    void transform_kernels(std::vector<double> &frekl) const;
};

#endif