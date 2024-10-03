#include <vector>

// interface
void
__surfvti(const float *thkm,const float *vphm,const float *vpvm,
            const float *vshm,const float *vsvm,const float *etam,
            const float *rhom,int nlayer,int nt,int mode, int wavetype,
            bool only_phase,const double *t,double *cp,double *cg);
std::vector<double> 
__surfvti_kl(const float *thkm,const float *vphm,const float *vpvm,
            const float *vshm,const float *vsvm,const float *etam,
            const float *rhom,int nlayer,int nt,int mode, int wavetype,
            const double *t,double *cp,double *cg);