#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <iostream>

#include "libswd/utils.hpp"


namespace py = pybind11;
using py::arg;

typedef py::array_t<float> ftensor;
typedef py::array_t<double> dtensor;
typedef std::tuple<dtensor,dtensor,dtensor,dtensor,dtensor,bool> tupledt6;
typedef std::tuple<dtensor,bool> tuple2;

static int 
check_wavetype(const std::string &wavetype)
{
    // check input parameters
    bool flag = wavetype == "Rc" || wavetype == "L";
    if(flag == false){
        std::cout << "parameter wavetp should be in one of [R,L]\n ";
        exit(0);
    }
    int iwave = 1;
    if(wavetype[0] == 'R') {
        iwave = 2;
    }

    return iwave;
}

auto 
swdvti(const ftensor &thk,const ftensor &vpv,const ftensor &vph,
        const ftensor &vsv,const ftensor &vsh,const ftensor &eta,const ftensor &rho, 
        const dtensor &t,const std::string &wavetype, bool only_phase,
        int mode=0)
{
    int iwave = check_wavetype(wavetype);

    // allocate space
    int nt = t.size(), n = thk.size();
    dtensor cg(nt),cp(nt);

    __surfvti(thk.data(),vph.data(),vpv.data(),vsh.data(),vsv.data(),eta.data(),
              rho.data(),n,nt,mode,iwave,only_phase,t.data(),cp.mutable_data(),
              cg.mutable_data());
    if(only_phase) cg = cp;

    return std::make_tuple(cp,cg);
}

auto 
vti_kl(const ftensor &thk,const ftensor &vpv,const ftensor &vph,
        const ftensor &vsv,const ftensor &vsh,const ftensor &eta,
        const ftensor &rho,const dtensor &t,
        const std::string &wavetype,int mode=0)
{
    int iwave = check_wavetype(wavetype);
    int nt = t.size(), n = thk.size();
    dtensor cg(nt),cp(nt);

    auto frekl = __surfvti_kl(thk.data(),vph.data(),vpv.data(),vsh.data(),vsv.data(),
                              eta.data(),rho.data(),n,nt,mode,iwave,
                              t.data(),cp.mutable_data(),cg.mutable_data());
    int nker = frekl.size() / (n * nt);
    dtensor frekl_out({nt,nker,n});
    memcpy(frekl_out.mutable_data(),frekl.data(),frekl.size() * sizeof(double));
    
    return std::make_tuple(cp,cg,frekl_out);
}



//
auto 
swdiso(const ftensor &thk,const ftensor &vp,const ftensor &vs,const ftensor &rho,
        const dtensor &t,const std::string &wavetype, bool only_phase,
        int mode=0)
{
    ftensor eta(thk.size());
    for(size_t i = 0; i < thk.size(); i++ ) eta.mutable_data()[i] = 1.;
    return swdvti(thk,vp,vp,vs,vs,eta,rho,t,wavetype,only_phase,mode);
}

auto 
iso_kl(const ftensor &thk,const ftensor &vp,const ftensor &vs,const ftensor &rho,
        const dtensor &t,const std::string wavetype,int mode=0)
{
    dtensor cp,cg,frekl;
    ftensor eta(thk.size());
    for(size_t i = 0; i < thk.size(); i++ ) eta.mutable_data()[i] = 1.;
    std::tie(cp,cg,frekl) = vti_kl(thk,vp,vp,vs,vs,eta,rho,t,wavetype,mode);

    int iwave = check_wavetype(wavetype);
    auto kl = frekl.unchecked<3>();
    
    // convert to iso type
    auto shape = frekl.shape();
    int nker = shape[1], n = shape[2], nt = shape[0];
    dtensor frekl_out;
    if(iwave == 1) { // love wave
        frekl_out.resize({nt,2,n});
        auto data = frekl_out.mutable_unchecked<3>();

        for(int it = 0; it < nt; it ++) {
            for(int i = 0; i < n; i ++) {
                data(it,0,i) = kl(it,0,i) + kl(it,1,i);
                data(it,1,i) = kl(it,2,i);
            }
        }
    }
    else {
        frekl_out.resize({nt,3,n});
        auto data = frekl_out.mutable_unchecked<3>();
        for(int it = 0; it < nt; it ++) {
            for(int i = 0; i < n; i ++) {
                data(it,0,i) = kl(it,0,i) + kl(it,1,i); // vp
                data(it,1,i) = kl(it,2,i); // vs
                data(it,2,i) = kl(it,4,i); // rho
            }
        }
    }
    
    return std::make_tuple(cp,cg,frekl_out);
}


PYBIND11_MODULE(libswd,m){
    m.doc() = "Surface wave dispersion and sensivity kernel\n";
    m.def("swdiso",&swdiso,arg("thk"),arg("vp"),arg("vs"),
          arg("rho"),arg("period"),arg("wavetype"),arg("only_phase"),
          arg("mode")=0,
          "ISO surface wave dispersion c++ wrapper");

    m.def("iso_kl",&iso_kl,arg("thk"),arg("vp"),arg("vs"),
          arg("rho"),arg("period"),arg("wavetype"),arg("mode")=0,
          "ISO surface wave kernel c++ wrapper");

    m.def("swdvti",&swdvti,arg("thk"),arg("vpv"),arg("vph"),arg("vsv"),
          arg("vsh"),arg("eta"),arg("rho"),arg("period"),arg("wavetype"),arg("only_phase"),
          arg("mode")=0,
          "VTI surface wave dispersion c++ wrapper");

    m.def("vti_kl",&vti_kl,arg("thk"),arg("vpv"),arg("vph"),arg("vsv"),
          arg("vsh"),arg("eta"),arg("rho"),arg("period"),arg("wavetype"),arg("mode")=0,
          "VTI surface wave kernel c++ wrapper");
}