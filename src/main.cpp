#include"surfdisp.h"
#include<fstream>
int main(){
    // mkdir to save dispersion
    int ierr = system("mkdir -p disp_ray");
    ierr = system("mkdir -p disp_love");

    LayerModel model;
    std::string filename="input/crust.txt";

    model.read_model(filename);

    int nt = 101;
    Eigen::VectorXd T(nt),cp(nt),cg(nt),crp(nt),crg(nt);

    // define zmax
    float zmax = model._dep(model._nlayer-1) + model._thk.maxCoeff();
    float dz = model._thk.head(model._nlayer-1).minCoeff() / 10;
    int n = (int)(zmax / dz) + 1;
    
    dz = zmax / (n-1);
    float z[n];
    for(int i=0;i<n;i++){
        z[i] = dz  * i;
    }

    for(int i=0;i<nt;i++){
        T(i) = 20. + i;
    }
    model.LovePhase(T.data(),cp.data(),nt,z,n);
    model.LoveGroup(T.data(),cp.data(),nt,z,n,cg.data());
    model.RayleighPhase(T.data(),crp.data(),nt,z,n);
    model.RayleighGroup(T.data(),crp.data(),nt,z,n,crg.data());

    FILE *fp = fopen("dispersion.dat","w");
    for(int i=0;i<nt;i++){
        fprintf(fp,"%lf %lf %lf %lf %lf\n",T(i),cp(i)/1000,cg(i)/1000,crp(i)/1000,
                crg(i)/1000);
        //std:: cout << T(i) << " " << cp(i) << std::endl;
    }
    fclose(fp);

    // compute sensitivity kernel
    Eigen::MatrixXd sen_vpRc(n,nt),sen_vsRc(n,nt),sen_rhoRc(n,nt);
    model.kernel_psv(T.data(),crp.data(),crg.data(),nt,z,n,sen_vpRc.data(),
        sen_vsRc.data(),sen_rhoRc.data());
    
    // save 
    for(int i=0;i<nt;i++){
        FILE *fp;
        std::string filename = "disp_ray/";
        filename = filename + std::to_string((int)(T[i])) + ".ker.dat";
        fp = fopen(filename.data(),"w");

        for(int j=0;j<n;j++){
            fprintf(fp,"%f %g %g %g\n",z[j],sen_vsRc(j,i),sen_vpRc(j,i),sen_rhoRc(j,i));
        }
        fclose(fp);
    }
    return 0;;
}
