#include"surfdisp.h"
#include<iostream>
#include<fstream>
using Eigen::Matrix4d;
using Eigen::Matrix2d;

/**
 * Earth Flattening Transform
 * TODO LIST
*/
void LayerModel:: _Flatten()
{

}

/**
 * read external velocity model
 * \param filename filename of external model
 */
int LayerModel:: read_model(std::string filename){
    std::ifstream infile;
    int n;
    infile.open(filename);
    infile >> n; // read no. of layers
    float thk[n],vp[n],vs[n],dep[n],rho[n];

    // read parameters
    int dummy,i;
    float d = 0.0;
    dep[0] = d;
    for(i=0;i<n-1;i++){
        infile >> dummy >> dep[i+1] >> rho[i] >> vp[i] >> vs[i];
        thk[i]= dep[i+1] - dep[i];
        d += thk[i];
    }
    i = n-1;
    infile >> dummy >> thk[i] >> rho[i] >> vp[i] >> vs[i];

    infile.close();

    _initialize(thk,vp,vs,rho,n);

    for(int i=0;i<n;i++){
        //std::cout << thk[i] << " " << dep[i] << " "<< vp[i]<<std::endl;
    }
    return 0;
};

/**
 * Interpolate function
 * \param z0  current depth
 * \param lambda,mu,rho  elastic param at z0
 */
void LayerModel:: 
interpolate(float z0,float &lambda,float &mu,float &rho)
{   
    // find index
    int i;
    if(z0 == _dep(0)){
        i = 0;
    }
    else if (z0 >= _dep(_nlayer-1)){
        i = _nlayer - 1;
    }
    else{
        i=1;
        for(;i<_nlayer-1;i++){
            if(z0 < _dep(i)) break;
        }
        i = i - 1;
    }

    // get value
    lambda = _lambda(i);
    mu = _mu(i);
    rho = _rho(i);
}

/**
 * compute Matrix form for Love wave (Aki & Richards 2002, 7.24)
 * \param w  angular frequency
 * \param k  wavenumber
 * \param mu  shear modulus
 * \param rho density g/cm^3
 * \param mat(2,2)  Love matrix
 */
void LayerModel::
LoveMatrix(double w,double k,float mu,float rho,Eigen::Matrix2d &mat)
{
    mat(0,0) = 0.0;
    mat(0,1) = 1. / mu;
    mat(1,0) = k * k * mu - w * w * rho;
    mat(1,1) = 0.;
}

/**
 * compute Matrix form for Rayleigh wave (Aki & Richards 2002, 7.28)
 * \param w  angular frequency
 * \param k  wavenumber
 * \param mu  shear modulus
 * \param lambda first Lame parameter
 * \param rho density g/cm^3
 * \param mat(4,4)  Rayleigh matrix
 */
void LayerModel:: 
RayleighMatrix(double w,double k,float lambda,float mu,float rho,
              Eigen::Matrix4d &mat)
{
    double zeta = 4. * mu * (lambda + mu) / (lambda + 2 * mu);
    mat.setZero();
    mat(0,1) = k; mat(0,2) = 1. / mu;
    mat(1,0) = -k * lambda / (lambda + 2 * mu); 
    mat(1,3) = 1. / (lambda + 2 * mu);
    mat(2,0) = k * k * zeta - w * w * rho; 
    mat(2,3) = k * lambda / (lambda + 2 * mu);
    mat(3,1) = -w * w  * rho; mat(3,2) = -k;
}

/**
 * compute Matrix form for alternative Rayleigh wave (Takeuchi & Saito 1973)
 * \param w  angular frequency
 * \param k  wavenumber
 * \param mu  shear modulus
 * \param lambda first Lame parameter
 * \param rho density g/cm^3
 * \param mat(5,5)  Rayleigh matrix
 */
void LayerModel::
RayleighMatrix_alt(double w,double k,float lambda,float mu,float rho,
                    Eigen::Matrix<double,5,5> &mat)
{
    mat.setZero();
    double lam_2mu = lambda + 2 * mu;
    mat(0,3) = 1. / mu ; mat(0,4) = -1. / lam_2mu;
    mat(1,3) = -w * w * rho; 
    mat(1,4) = w * w * rho - k * k * (lam_2mu - lambda * lambda / lam_2mu);
    mat(2,3) = k; mat(2,4) = k * lambda / lam_2mu;
    mat(3,0) = -mat(1,4); mat(3,1) = 1. / lam_2mu; mat(3,2) = -2 * k * lambda / lam_2mu;
    mat(4,0) = w * w * rho; mat(4,1) = -1. / mu; 
    mat(4,2) = -2 * k; 
}

/**
 * Solve SH ordinary equations system and find eigen functions for fixed w,k
 * \param w angular frequency
 * \param k wavenumber
 * \param l1,l2 eigen function
 * \param z depth points
 * \param n len(z)
 */
void LayerModel:: 
_SolveSH(double w,double k,Eigen::VectorXd &l1,Eigen::VectorXd &l2,
        float *z,int n)
{
    // set initial value
    l1.setZero();l2.setZero();
    l1(n-1) = 0.0; // no disp
    l2(n-1) = 1.0; // stress
    float dz = z[1] - z[0];

    // intialize runge-kutta coef
    Eigen:: Matrix2d mat;
    Eigen::Vector2d y,k1,k2,k3,k4;
    y(0) = 0.0; y(1) = 1.0;

    // loop from bottom to top
    float mu,lambda,rho;
    for(int i=n-1;i>0;i--){
        // compute k1
        interpolate(z[i],lambda,mu,rho);
        LoveMatrix(w,k,mu,rho,mat);
        k1 = mat * y;

        // compute k2
        interpolate(z[i] - dz * 0.5,lambda,mu,rho);
        LoveMatrix(w,k,mu,rho,mat);
        k2 = mat * (y - 0.5 * dz * k1);

        // compute k3
        k3 = mat * (y - 0.5 * dz * k2);

        // compute k4
        interpolate(z[i-1],lambda,mu,rho);
        LoveMatrix(w,k,mu,rho,mat);
        k4 = mat * (y - dz * k3);

        // update l1 and l2
        y = y - dz / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        l1(i-1) = y(0);
        l2(i-1) = y(1);

        // rescale to maximum to prevent overflow
        double amp = l2.array().abs().maxCoeff();
        l1 /= amp; l2 /= amp;
        y /= amp;
    }
}

/**
 * Solve P-SV ordinary equations system and find eigen functions for fixed w,k
 * \param w angular frequency
 * \param k wavenumber
 * \param r1,r2,r3,r4 eigen functions
 * \param z depth points
 * \param n len(z)
 * \param initial conditions, r3(inf) = 0 or r4(inf) = 0
 */
void LayerModel::
_SolvePSV(double w,double k,Eigen::VectorXd &r1,Eigen::VectorXd &r2,
        Eigen::VectorXd &r3,Eigen::VectorXd &r4,float *z,int n,
        int initial_condition)
{
   // set initial value
    r1.setZero();r2.setZero();
    r3.setZero();r4.setZero();
    float dz = z[1] - z[0];

    // intialize runge-kutta coef
    Eigen:: Matrix4d mat;
    Eigen::Vector4d y,k1,k2,k3,k4;
    y.setZero();
    if(initial_condition == 1){
        r3(n-1) = 1.0;
        y(2) = 1.0;
    }
    else{
        r4(n-1) = 1.0;
        y(3) = 1.0;
    }

    // integrate from bottom to top
    float mu,lambda,rho;
    for(int i=n-1;i>0;i--){
        // compute k1
        interpolate(z[i],lambda,mu,rho);
        RayleighMatrix(w,k,lambda,mu,rho,mat);
        k1 = mat * y;

        // compute k2
        interpolate(z[i] - dz * 0.5,lambda,mu,rho);
        RayleighMatrix(w,k,lambda,mu,rho,mat);
        k2 = mat * (y - 0.5 * dz * k1);

        // compute k3
        k3 = mat * (y - 0.5 * dz * k2);

        // compute k4
        interpolate(z[i-1],lambda,mu,rho);
        RayleighMatrix(w,k,lambda,mu,rho,mat);
        k4 = mat * (y - dz * k3);

        // update l1 and l2
        y = y - dz / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        r1(i-1) = y(0); r2(i-1) = y(1);
        r3(i-1) = y(2); r4(i-1) = y(3);

        // rescale to maximum to prevent overflow
        double amp;
        if(initial_condition == 1){
            amp = r3.array().abs().maxCoeff();
        }
        else{
            amp = r4.array().abs().maxCoeff();
        }

        r1 = r1 / amp; r2 = r2 / amp;
        r3 = r3 / amp; r4 = r4 / amp;
        y = y / amp;
    }
}

/**
 * Solve alternative P-SV systems, see Takeuchi & Saito 1972
 * \param w angular frequency
 * \param k wavenumber
 * \param r1,r2,r3,r4,r5 eigen functions
 * \param z depth points
 * \param n len(z)
 */
void LayerModel::
_SolvePSV_alt(double w,double k,Eigen::VectorXd &r1,Eigen::VectorXd &r2,
            Eigen::VectorXd &r3,Eigen::VectorXd &r4,Eigen::VectorXd &r5,
            float *z,int n)
{
   // set initial value
    r1.setZero();r2.setZero();
    r3.setZero();r4.setZero();
    r5.setZero();
    float dz = z[1] - z[0];

    // intialize runge-kutta coef
    Eigen:: Matrix<double,5,5> mat;
    Eigen::Matrix<double,5,1> y,k1,k2,k3,k4;
    y.setZero();
    r2(n-1) = 1.0;
    y(1) = 1.0;

    // loop from bottom to top
    float mu,lambda,rho;
    for(int i=n-1;i>0;i--){
        // compute k1
        interpolate(z[i],lambda,mu,rho);
        RayleighMatrix_alt(w,k,lambda,mu,rho,mat);
        k1 = mat * y;

        // compute k2
        interpolate(z[i] - dz * 0.5,lambda,mu,rho);
        RayleighMatrix_alt(w,k,lambda,mu,rho,mat);
        k2 = mat * (y - 0.5 * dz * k1);

        // compute k3
        k3 = mat * (y - 0.5 * dz * k2);

        // compute k4
        interpolate(z[i-1],lambda,mu,rho);
        RayleighMatrix_alt(w,k,lambda,mu,rho,mat);
        k4 = mat * (y - dz * k3);

        // update l1 and l2
        y = y - dz / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
        r1(i-1) = y(0); r2(i-1) = y(1);
        r3(i-1) = y(2); r4(i-1) = y(3);
        r5(i-1) = y(4);

        // rescale to maximum to prevent overflow
        double amp = r2.array().abs().maxCoeff();
        r1 = r1 / amp; r2 = r2 / amp;
        r3 = r3 / amp; r4 = r4/ amp;
        r5 = r5 / amp;
        y = y / amp;
    }
}

/**
 * boundary conditions for SH systems
 * \param c phase velocity m/s
 * \param w angular frequncy
 * \param z depth points
 * \param n len(z)
 */
double LayerModel:: 
_disper_sh(double c,double w,float *z,int n)
{   
    // define eigen function and 
    Eigen::VectorXd l1(n),l2(n);
    double k = w /c;
    _SolveSH(w,k,l1,l2,z,n);

    return l2(0);
}

/**
 * boundary conditions for P-SV systems
 * \param c phase velocity m/s
 * \param w angular frequncy
 * \param z depth points
 * \param n len(z)
 */
double LayerModel:: 
_disper_psv(double c,double w,float *z,int n)
{
    // define eigen function 
    double k = w / c;
    Eigen::VectorXd r1(n),r2(n),r3(n),r4(n),r5(n);
    _SolvePSV_alt(w,k,r1,r2,r3,r4,r5,z,n);

    return r2(0);
}

/**
 * Find Phase velocity for Love wave
 * \param t periods, in s 
 * \param cp phase velocity m/s
 * \param nt nt = len(t)
 * \param z depth points
 * \param n len(z)
 * \param mode mode for surface wave, =0 for fundamental
 */
void LayerModel:: 
LovePhase(double *t,double *cp,int nt,
         float *z,int n,int mode)
{
    double c,clim = _vs.maxCoeff(); 
    double dc = 1, threshold = 1.0e-3; // search length and stop threshold

    for(int i=0;i<nt;i++){
        int modecnt = 0; // mode count
        double w = 2 * pi / t[i];
        if(i == 0) {
            c = _vs.minCoeff() * 0.95; // initial trial value
        }
        else{
            c = cp[i-1] + threshold;
        }

        // bracket root
        double fun1 = _disper_sh(c,w,z,n),fun2;
        for(c=c+dc;c<=clim;c+=dc){
            fun2 = _disper_sh(c,w,z,n);
            if(fabs(fun1+fun2)==fabs(fun1)+fabs(fun2)){
                fun1=fun2;
                continue;
            }
            else{
                modecnt += 1;
                if(modecnt==mode+1){
                    break;
                }
                fun1=fun2;
            }
        }

        // find root by bisecting method
        double c1 = c- dc, c2 = c;
        double fmid;
        double cmid = 0.5 * (c1 + c2);
        int count = 0;
        while(count < 5 ){
            fmid = _disper_sh(c,w,z,n);
            if(fabs(fmid)+fabs(fun1)==fabs(fmid+fun1)){
                c1 = cmid;
                fun1 = fmid;
            }
            else{
                c2 = cmid;
                fun2 = fmid;
            }
            count += 1;
            if(fabs(cmid - 0.5 * c1 -0.5 * c2) < threshold) break;
            cmid = 0.5 * (c1 + c2);
        }

        // here we find root
        cp[i] = cmid;
    }
}

/**
 * solve Rayleigh equation in half space
 * \param a,b vp and vs for half space
 * \param c initial guess
 */
double RayleighHalf(float a,float b){
    double c = 0.95 * b;
    float kappa,gamma,k2,gk2;
    for(int j=0;j<5;j++){
        gamma = b / a,kappa = c / b;
        k2 = kappa * kappa;
        gk2 = pow(gamma * kappa,2);
        double fac1 = sqrt(1 - gk2), fac2 = sqrt(1 - k2);
        double fr = pow(2 - k2,2) -4.0 * fac1 * fac2;
        double frp = -4. * (2 - k2) * kappa 
                    + 4 * fac2 * gamma * gamma * kappa / fac1
                    + 4. * fac1 * kappa /fac2;
        frp /= b;
        c = c -fr / frp; 
    }

    return c;
}

/**
 * Find Phase velocity for Rayleigh wave
 * \param t periods, in s 
 * \param cp phase velocity m/s
 * \param nt nt = len(t)
 * \param z depth points
 * \param n len(z)
 * \param mode mode for surface wave, =0 for fundamental
 */
void LayerModel:: 
RayleighPhase(double *t,double *cp,int nt,float *z,int n,int mode)
{
    double c,clim = _vs.maxCoeff(); 
    double dc = 1, threshold = 1.0e-6; // search length and stop threshold

    for(int i=0;i<nt;i++){
        int modecnt = 0; // mode count
        double w = 2 * pi / t[i];
        if(i == 0) {
            // use half space phase velo for initial guess
            c = RayleighHalf(_vp(0),_vs(0));
        }
        else{
            c = cp[i-1] + threshold;
        }

        // bracket root
        double fun1 = _disper_psv(c,w,z,n),fun2;
        for(c=c+dc;c<=clim;c+=dc){
            fun2 = _disper_psv(c,w,z,n);
            if(fabs(fun1+fun2)==fabs(fun1)+fabs(fun2)){
                fun1=fun2;
                continue;
            }
            else{
                modecnt += 1;
                if(modecnt==mode+1){
                    break;
                }
                fun1=fun2;
            }
        }

        // find root by bisecting method
        double c1 = c- dc, c2 = c;
        double fmid;
        double cmid = (c1 + c2) * 0.5;
        int count = 0;
        while(count < 5){
            fmid = _disper_psv(c,w,z,n);
            if(fabs(fmid)+fabs(fun1)==fabs(fmid+fun1)){
                c1 = cmid;
                fun1 = fmid;
            }
            else{
                c2 = cmid;
                fun2 = fmid;
            }
            count += 1;
            if(fabs(cmid - 0.5 * (c1 + c2) < threshold ))break;
            cmid = (c1 + c2 ) * 0.5;
        }

        // here we find root
        cp[i] = cmid;
    }
}

/**
 * Find Group velocity for Love wave
 * \param t periods, in s 
 * \param cp phase velocity m/s
 * \param nt nt = len(t)
 * \param z depth points
 * \param n len(z)
 * \param cg group velocity m/s
 */
void LayerModel::
LoveGroup(double *t,double *cp,int nt,float *z,int n,double *cg)
{
    // interp at z
    float mu[n],lambda[n],rho[n];
    for(int i=0;i<n;i++){
        interpolate(z[i],lambda[i],mu[i],rho[i]);
    }

    for(int i=0;i<nt;i++){
        double c = cp[i];
        double w = 2. * pi / t[i];
        double k = w /c ;

        // compute eigen function
        Eigen::VectorXd l1(n),l2(n);
        double I1,I2,I3;
        _eigen_sh(w,k,l1,l2,I1,I2,I3,mu,rho,z,n);

        
        FILE *fp;
        std::string filename = "disp_love/";
        filename = filename + std::to_string((int)(t[i])) + ".dat";
        fp = fopen(filename.data(),"w");
        for(int kk=0;kk<n;kk++){
            fprintf(fp,"%f %g %g\n",z[kk],l1(kk),l2(kk));
        }
        fclose(fp);

        // compute group
        cg[i] = I2 / (c * I1);
        printf("Love wave, period = %f cp = %f cg = %f\n",t[i],cp[i],cg[i]);
    }
}

/**
 * Find Group velocity for Rayleigh wave
 * \param t periods, in s 
 * \param cp phase velocity m/s
 * \param nt nt = len(t)
 * \param z depth points
 * \param n len(z)
 * \param cg group velocity m/s
 */
void LayerModel::
RayleighGroup(double *t,double *cp,int nt,
            float *z,int n,double *cg)
{   
    // interp at z
    float mu[n],lambda[n],rho[n];
    for(int i=0;i<n;i++){
        interpolate(z[i],lambda[i],mu[i],rho[i]);
    }

    for(int i=0;i<nt;i++){
        double c = cp[i];
        double w = 2. * pi / t[i];
        double k = w /c ;

        // compute eigen function and energy integral
        Eigen::VectorXd r1(n),r2(n),r3(n),r4(n);
        double I1,I2,I3,I4;
        _eigen_psv(w,k,r1,r2,r3,r4,I1,I2,I3,I4,lambda,mu,rho,z,n);

        // save eigen function
        FILE *fp;
        std::string filename = "disp_ray/";
        filename = filename + std::to_string((int)(t[i])) + ".dat";
        fp = fopen(filename.data(),"w");
        for(int kk=0;kk<n;kk++){
            fprintf(fp,"%f %g %g %g %g\n",z[kk],r1(kk),r2(kk),r3(kk),r4(kk));
        }
        fclose(fp);

        // compute group
        cg[i] = (I2 + I3 / k) / (c * I1);
        printf("Rayleigh wave, period = %f cp = %f cg = %f\n",t[i],cp[i],cg[i]);
    }
}

/**
 * solve Eigen function and energy integrals for Love Wave
 * \param w  angular frequency
 * \param k  wavenumber
 * \param l1,l2 eigen function for Love wave
 * \param I1,I2,I3 energy integrals
 * \param mu,rho elastic param
 * \param z  depths to locate eigen function
 * \param n  no. of depths used
 */ 
void LayerModel::
_eigen_sh(double w,double k,Eigen::VectorXd &l1,Eigen::VectorXd &l2,
        double &I1,double &I2,double &I3,float *mu,float *rho,
        float *z,int n)
{
    assert(l1.size() == n);

    // solve eigen function
    _SolveSH(w,k,l1,l2,z,n);

    // obtain energy integrals
    I1 = 0.0; I2 = 0.0; I3 = 0.0;
    double dz = z[1] - z[0];
    for(int i=0;i<n;i++){
        double c = dz * 0.5;
        if(i==0 || i == n-1) c = c * 0.5; // trapzoid intergrate
        I1 += c * rho[i] * pow(l1(i),2);
        I2 += c * mu[i] * pow(l1(i),2);
        I3 += c * pow(l2(i),2) / mu[i];
    }
}

/**
 * solve Eigen function and energy integrals for Rayleigh Wave
 * \param w  angular frequency
 * \param k  wavenumber
 * \param l1,l2 eigen function for Love wave
 * \param I1,I2,I3 energy integrals
 * \param lambda,mu,rho elastic param
 * \param z  depths to locate eigen function
 * \param n  no. of depths used
 */ 
void LayerModel::
_eigen_psv(double w,double k,Eigen::VectorXd &r1,Eigen::VectorXd &r2,
                Eigen::VectorXd &r3,Eigen::VectorXd &r4,
                double &I1,double &I2,double &I3,double &I4,float *lambda,
                float *mu,float *rho,float *z,int n)
{
    assert(r1.size() == n);

    // solve eigen function for 2 independent systems
    Eigen::VectorXd r11(n),r21(n),r31(n),r41(n);
    Eigen::VectorXd r12(n),r22(n),r32(n),r42(n);
    _SolvePSV(w,k,r11,r21,r31,r41,z,n,1); // condi =1
    _SolvePSV(w,k,r12,r22,r32,r42,z,n,2); // condi !=1 

    // obtain the linear combination coef
    double c1 = 1., c2 = - r31(0) / r32(0);
    r1 = r11 * c1 + r12 * c2;
    r2 = r21 * c1 + c2 * r22;
    r3 = r31 * c1 + c2 * r32;
    r4 = r41 * c1 + c2 * r42;
    
    // normalize to vertical displacement at the surface
    double mm = r2(0);
    r1 = r1 / mm;
    r2 = r2 / mm;
    r3 = r3 / mm;
    r4 = r4 / mm;

    // compute energy intergrals
    double dz = z[1] - z[0];
    I1 = 0.0,I2 = 0.0,I3=0.0,I4 = 0.0;

    for(int j=0;j<n;j++){
        double c = dz * 0.5 ;
        if(j==0 || j == n-1) c *= 0.5;
        double lam_2mu = lambda[j] + 2 * mu[j];
        I1 += c * rho[j] * (pow(r1(j),2) + pow(r2(j),2));
        I2 += c * ( lam_2mu * pow(r1(j),2) + mu[j] * pow(r2(j),2));
        double I31 = -k * pow(lambda[j],2) / lam_2mu * pow(r1(j),2) 
                    - k * mu[j] * pow(r2(j),2);
        double I32 = lambda[j] /lam_2mu * r1(j) * r4(j) - 
                    r2(j) * r3(j);
        I3 += c * (I31 + I32);
    }
}

void surfdisp96(float *thk,float *vp,float *vs,float *den,
                int nlayer,int iwave,int igr,double *t,
                double *c,int nt,int mode=0,bool flatten=0)
{
    // construct layermodel class
    LayerModel model(thk,vp,vs,den,nlayer);

    // compute param vector
    // define zmax
    float zmax = model._dep(nlayer-1)  + 100.;
    float dz = model._thk.head(nlayer-1).minCoeff() / 10.0;
    int n = (int)(zmax / dz) + 1;
    dz = zmax / (n-1);

    // define param vector
    float z[n],mu[n],rho[n],lambda[n];
    for(int i=0;i<n;i++){
        z[i] = dz  * i;
        model.interpolate(z[i],lambda[i],mu[i],rho[i]);
    }
}

/**
 * solve sensitivity kernel for Love wave, Aki and Richards 2003, 7.71
 * \param w  angular frequency
 * \param k  wavenumber
 * \param I2 energy intergral 
 * \param l1,l2 eigen function for Love wave
 * \param mu shear moduli elastic param
 * \param krho,kmu  sensitivity kernel for rho and mu
 */ 
void LayerModel:: 
kernel_sh(double w,double k,double &I2,Eigen::VectorXd &l1,
          Eigen::VectorXd &l2,float *mu,double *krho,double *kmu)
{
    int n = l1.size(); // no. of points

    double deno = 4. * pow(k,3) * I2;
    for(int i=0;i<n;i++){
        double mu1 = pow(l1(i) * k,2);
        double mu2 = pow(l2(i) / mu[i],2);
        kmu[i] = w * (mu1 + mu2) / deno;
        krho[i] = -w * pow(w*l1(i),2) / deno;
    }
}

/**
 * solve sensitivity kernel for Love wave, for rho and vs
 * \param w  angular frequency
 * \param k  wavenumber
 * \param I2 energy intergral 
 * \param l1,l2 eigen function for Love wave
 * \param mu shear moduli elastic param
 * \param krho,kvp,kvs  sensitivity kernel for rho,vp and vs
 */ 
void LayerModel:: 
kernel_sh(double w,double k,double &I2,Eigen::VectorXd &l1,
         Eigen::VectorXd &l2,float *mu,float *rho,
         double *kvs,double *krho)
{
    // compute kmu and krho
    int n = l1.size(); 
    double sen_rho[n],sen_mu[n];
    kernel_sh(w,k,I2,l1,l2,mu,sen_rho,sen_mu);

    // convert kmu/krho to krho/kvs
    for(int i=0;i<n;i++){
        double vs = sqrt(mu[i] / rho[i]);
        kvs[i] = 2. * rho[i] * vs * sen_mu[i];
        krho[i] = sen_mu[i] * vs * vs + sen_rho[i];
    }
}

/**
 * solve sensitivity kernel for Rayleigh wave, Aki and Richards 2003, 7.78
 * \param w  angular frequency
 * \param k  wavenumber
 * \param U group velocity
 * \param I2 energy intergral 
 * \param r1-r4 eigen function for Love wave
 * \param lambda,mu lame elastic param
 * \param krho,klambda,kmu  sensitivity kernel for rho,lambda,mu
 */ 
void LayerModel:: 
kernel_psv(double w,double k,double U,double &I1,Eigen::VectorXd &r1,
            Eigen::VectorXd &r2,Eigen::VectorXd &r3,Eigen::VectorXd &r4,
            float *lambda,float *mu,double *klam,double *kmu,double *krho)
{
    int n = r1.size();
    double deno = 4 * k * k * U * I1;

    for(int i=0;i<n;i++){
        double lam2mu = lambda[i] + 2. * mu[i];
        klam[i] = pow((2. * k * mu[i] * r1(i) + r4(i)) / lam2mu,2);
        klam[i] /= deno;
        kmu[i] = 2. * pow(k * r1(i),2) + pow(r3(i)/mu[i],2) 
                + 2. * pow((-k * lambda[i] * r1(i) + r4(i)) / lam2mu,2);
        kmu[i] /= deno;
        krho[i] = -w * w * (pow(r1(i),2) + pow(r2(i),2));
        krho[i] /= deno;
    }
}

/**
 * solve sensitivity kernel for Rayleigh wave, for vp/vs/rho
 * \param w  angular frequency
 * \param k  wavenumber
 * \param U group velocity
 * \param I2 energy intergral 
 * \param r1-r4 eigen function for Love wave
 * \param lambda,mu,rho lame elastic param
 * \param krho,kvp,kvs  sensitivity kernel for rho,vp and vs
 */ 
void LayerModel::
kernel_psv(double w,double k,double U,double &I1,Eigen::VectorXd &r1,
            Eigen::VectorXd &r2,Eigen::VectorXd &r3,Eigen::VectorXd &r4,
            float *lambda,float *mu,float *rho,double *kvp,double *kvs,
            double *krho)
{
    int n = r1.size();
    double sen_rho[n],sen_mu[n],sen_lambda[n];

    // compute sen_rho,sen_mu,sen_lambda
    kernel_psv(w,k,U,I1,r1,r2,r3,r4,lambda,mu,sen_lambda,sen_mu,sen_rho);

    // convert to sen_vp/vs/rho
    for(int i=0;i<n;i++){
        double vs = sqrt(mu[i] / rho[i]);
        double vp = sqrt((lambda[i] + 2 * mu[i]) / rho[i]);
        kvp[i] = 2. * rho[i] * sen_lambda[i] * vp;
        kvs[i] = 2. * rho[i] * vs * (-2 * sen_lambda[i] + sen_mu[i]);
        krho[i] = sen_rho[i]  + (vp * vp -2 * vs * vs) * sen_lambda[i] 
                  + vs * vs * sen_mu[i];
        
        // relative kernel
        kvp[i] *= vp;
        kvs[i] *= vs;
        krho[i] *= rho[i];
    }
}

/**
 * solve sensitivity kernel for Love wave, for vs/rho
 * \param t  period
 * \param c,cg  phase and group velocity
 * \param z  depths to locate eigen function
 * \param n  no. of depths used
 * \param sen_vs,sen_rho  sensitivity kernel for rho and vs, shape(nz,nt), column-major
 */ 
void LayerModel::
kernel_sh(double *t,double *cp,double *cg,int nt,float *z,int n,
        double *sen_vs,double *sen_rho)
{
    // interp at z
    float mu[n],lambda[n],rho[n];
    for(int i=0;i<n;i++){
        interpolate(z[i],lambda[i],mu[i],rho[i]);
    }

    for(int i=0;i<nt;i++){
        double c = cp[i];
        double w = 2. * pi / t[i];
        double k = w /c ;

        // compute eigen function
        Eigen::VectorXd l1(n),l2(n);
        double I1,I2,I3;
        _eigen_sh(w,k,l1,l2,I1,I2,I3,mu,rho,z,n);

        // compute sensitivity kernel
        kernel_sh(w,k,I2,l1,l2,mu,rho,sen_vs + i * n,sen_rho + i * n);
    }
}

/**
 * solve sensitivity kernel for Rayleigh wave, for vp/vs/rho
 * \param t  period
 * \param cp,cg  phase and group velocity
 * \param z  depths to locate eigen function
 * \param n  no. of depths used
 * \param sen_vs,sen_rho,sen_vp sensitivity kernel for rho and vs, shape(nz,nt), column-major
 */ 
void LayerModel::
kernel_psv(double *t,double *cp,double *cg,int nt,float *z,int n,
          double *sen_vp,double *sen_vs,double *sen_rho)
{
    // interp at z
    float mu[n],lambda[n],rho[n];
    for(int i=0;i<n;i++){
        interpolate(z[i],lambda[i],mu[i],rho[i]);
    }

    for(int i=0;i<nt;i++){
        double c = cp[i];
        double w = 2. * pi / t[i];
        double k = w /c ;

        // compute eigen function and energy integral
        Eigen::VectorXd r1(n),r2(n),r3(n),r4(n);
        double I1,I2,I3,I4;
        _eigen_psv(w,k,r1,r2,r3,r4,I1,I2,I3,I4,lambda,mu,rho,z,n);

        // compute sensitivity kernel
        kernel_psv(w,k,cg[i],I1,r1,r2,r3,r4,lambda,mu,rho,
                sen_vp + i * n,sen_vs + i * n,sen_rho + i * n);
    }
}