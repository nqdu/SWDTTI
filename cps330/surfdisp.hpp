extern "C" {
// surface wave dispersion c++ wrapper
void surfdisp96_(float *thkm,float *vpm,float *vsm,float *rhom,
                int nlayer,int iflsph,int iwave,int mode,
                int igr,int kmax,double *t,double *cg,int *ierr);
}