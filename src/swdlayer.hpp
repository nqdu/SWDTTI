#ifndef SWD_LAYER_MODEL
#define SWD_LAYER_MODEL

#include <complex>
#include <vector>
#include <array>

class LayerModel  {

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

    LayerModel(){};
    void initialize();

private:
    std::vector<int> ilayer_flag; // shape(nspec + 1), return layer flag 
    std::vector<double> Mmat,Emat,Kmat; // matrices for SEM, om^2 M = k^2 K + E
    double PHASE_VELOC_MIN,PHASE_VELOC_MAX;

public:

    // density
    std::vector<double> xrho;

    // vti Love parameters
    std::vector<double> xA,xC,xL,xF,xN; // shape(nspec * NGLL + NGRL)


    // VTI model
    void create_database(double freq,int nlayer, const float *vph, const float* vpv,
                        const float *vsh, const float *vsv, const float *eta,
                        const float *rho,const float *thk);

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