#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include <iostream>
#include "surfdisp.hpp"

namespace py = pybind11;
using py::arg;

typedef py::array_t<float> ftensor;
typedef py::array_t<double> dtensor;
typedef std::tuple<dtensor,dtensor,dtensor,dtensor,dtensor,bool> tupledt6;
typedef std::tuple<dtensor,bool> tuple2;

/**
 * calculates the dispersion values for any layered model, any frequency, and any mode.
 * @param nlayer no. of layers
 * @param thkm,vpm,vsm,rhom model, shape(nlayer)
 * @param kmax no. of periods used
 * @param t,cp period and phase velocity, shape(kmax)
 * @param sphere true for spherical earth, false for flat earth 
 * @param wavetype one of [Rc,Rg,Lc,Lg]
 * @param mode i-th mode of surface wave, 0 fundamental, 1 first higher, ....
 * @param keep_flat keep flat earth phase/group velocity or convert it to spherical
 */
int _surfdisp(float *thk,float *vp,float *vs,float *rho,
            int nlayer,double *t,double *cg,int kmax,const std::string &wavetype,
            int mode,bool sphere,bool keep_flat)
{
    int iwave,igr,ifsph=0;
    if(wavetype=="Rc"){
        iwave = 2;
        igr = 0;
    }
    else if(wavetype == "Rg"){
        iwave = 2;
        igr = 1;
    }
    else if(wavetype=="Lc"){
        iwave = 1;
        igr = 0;
    }
    else if(wavetype=="Lg"){
        iwave = 1;
        igr = 1;
    }
    else{
        std::cout <<"wavetype should be one of [Rc,Rg,Lc,Lg]"<<std::endl;
        exit(0);
    }

    if(sphere == true) ifsph = 1;
    int ierr;
    surfdisp96_(thk,vp,vs,rho,nlayer,ifsph,iwave,mode+1,igr,kmax,t,cg,&ierr);

    if(ierr == 1){
        return ierr;
    } 

    return ierr;
}

tuple2 forward(ftensor &thk,ftensor &vp,ftensor &vs,ftensor &rho,
                dtensor &t,std::string wavetype,
                int mode=0,bool sphere=false)
{
    // check input parameters
    bool flag = wavetype == "Rc" || wavetype == "Rg" || 
                wavetype == "Lc" || wavetype == "Lg";
    if(flag == false){
        std::cout << "cnew = _flat2sphere(double t,double c,std::string wavetp)\n";
        std::cout << "parameter wavetp should be in one of [Rc,Rg,Lc,Lg]\n ";
        exit(0);
    }

    // allocate space
    int nt = t.size(), n = thk.size();
    dtensor cg(nt);
    int ierr;
    ierr = _surfdisp(thk.mutable_data(),vp.mutable_data(),vs.mutable_data(),
                rho.mutable_data(),n,t.mutable_data(),cg.mutable_data(),
                nt,wavetype,mode,sphere,false);


    bool return_flag = true;
    if(ierr == 1) return_flag = false;
    auto tt = std::make_tuple(cg,return_flag);

    return tt;
}

PYBIND11_MODULE(cps330,m){
    m.doc() = "Surface wave dispersion and sensivity kernel\n";
    m.def("forward",&forward,arg("thk"),arg("vp"),arg("vs"),
          arg("rho"),arg("period"),arg("wavetype"),
          arg("mode")=0,arg("sphere")=false,
          "Surface wave dispersion c++ wrapper");
}