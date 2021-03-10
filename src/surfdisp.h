#pragma once
#include<Eigen/Core>
#include<iostream>

#define pi 3.1415926535898
class LayerModel{
    public:
    int _nlayer;
    Eigen::VectorXf _thk,_rho,_vp,_vs,_dep,_mu,_lambda;

    //==============================================
    // Initialize function
    LayerModel(){

    }

    void _initialize(float *thk,float *vp,float *vs,float *rho,int n){
        _nlayer = n;
        _vs.resize(n);
        _thk.resize(n);_vp.resize(n);
        _rho.resize(n);_mu.resize(n);
        _lambda.resize(n);_dep.resize(n);

        // copy parameters
        for(int i=0;i<n;i++){
            _thk(i) = thk[i]*1000. ; _vp(i) = vp[i] *1000.;
            _vs(i) = vs[i] *1000.;   _rho(i) = rho[i] *1000.;
            _mu(i) = _rho(i) * pow(_vs(i),2);
            _lambda(i) = _rho(i) * pow(_vp(i),2) - 2. * _mu(i);

            if(i==0)
                _dep(i) = 0.0;
            else 
                _dep(i) = _dep(i-1) + _thk(i-1);
        }
    }
    LayerModel(float *thk,float *vp,float *vs,float *rho,int n){
        _initialize(thk,vp,vs,rho,n);
    }
    // ================================================//
    // end initialize function

    //
    int read_model(std::string filename);

    private:
    void LoveMatrix(double w,double k,float mu,float rho,Eigen::Matrix2d &mat);
    void RayleighMatrix(double w,double k,float lambda,float mu,float rho,
                        Eigen::Matrix4d &mat);
    void RayleighMatrix_alt(double w,double k,float lambda,float mu,float rho,
                        Eigen::Matrix<double,5,5> &mat);

    // solve ordinary differential equations with 4-th order Runge-Kutta
    private:
    void _SolveSH(double w,double k,Eigen::VectorXd &l1,Eigen::VectorXd &l2,
                float *z,int n);
    void _SolvePSV(double w,double k,Eigen::VectorXd &r1,Eigen::VectorXd &r2,
                    Eigen::VectorXd &r3,Eigen::VectorXd &r4,float *z,int n,
                    int initial_condition = 1);
    void _SolvePSV_alt(double w,double k,Eigen::VectorXd &r1,Eigen::VectorXd &r2,
                    Eigen::VectorXd &r3,Eigen::VectorXd &r4,Eigen::VectorXd &r5,
                    float *z,int n);
    
    private:
    void _Flatten();

    // dispersion relations for SH and PSV system
    private:
    double _disper_sh(double c,double w,float *z,int n);
    double _disper_psv(double c,double w,float *z,int n);
    void _eigen_sh(double w,double k,Eigen::VectorXd &l1,Eigen::VectorXd &l2,
                double &I1,double &I2,double &I3,float *mu,float *rho,
                float *z,int n);
    void _eigen_psv(double w,double k,Eigen::VectorXd &r1,Eigen::VectorXd &r2,
                Eigen::VectorXd &r3,Eigen::VectorXd &r4,
                double &I1,double &I2,double &I3,double &I4,float *lambda,
                float *mu,float *rho,float *z,int n);

    public:
    void interpolate(float z0,float &lambda,float &mu,float &rho);
    void LovePhase(double *t,double *c,int nt,
                    float *z,int n,int mode=0);
    void RayleighPhase(double *t,double *c,int nt,
                    float *z,int n,int mode=0);
    void LoveGroup(double *t,double *c,int nt,
                float *z,int n,double *cg);
    void RayleighGroup(double *t,double *c,int nt,
                float *z,int n,double *cg);
    public:
    void kernel_sh(double *t,double *c,double *cg,int nt,float *z,int n,
                    double *sen_vs,double *sen_rho);
    void kernel_psv(double *t,double *c,double *cg,int nt,float *z,int n,
                    double *sen_vp,double *sen_vs,double *sen_rho);
    
    private: 
    void kernel_sh(double w,double k,double &I2,Eigen::VectorXd &l1,
                  Eigen::VectorXd &l2,float *mu,double *krho,double *kmu);
    void kernel_sh(double w,double k,double &I2,Eigen::VectorXd &l1,
                  Eigen::VectorXd &l2,float *mu,float *rho,
                  double *kvs,double *krho);
    void kernel_psv(double w,double k,double U,double &I1,Eigen::VectorXd &r1,
                    Eigen::VectorXd &r2,Eigen::VectorXd &r3,Eigen::VectorXd &r4,
                    float *lambda,float *mu,double *klam,double *kmu,double *krho);

    void kernel_psv(double w,double k,double U,double &I1,Eigen::VectorXd &r1,
                    Eigen::VectorXd &r2,Eigen::VectorXd &r3,Eigen::VectorXd &r4,
                    float *lambda,float *mu,float *rho,double *kvp,double *kvs,
                    double *krho);
};

void surfdisp96(float *thk,float *vp,float *vs,float *rho,
                int iwave,int igr,int mode);