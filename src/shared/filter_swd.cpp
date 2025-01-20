#include <Eigen/Core>
#include <vector>
#include <algorithm>

typedef std::complex<double> dcmplx;

/**
 * @brief filter swd in the [ vmin,vmax] and Rek >= Imk
 * 
 * @param vmin/vmax min/max velocity in this region 
 * @param om angular frequency
 * @param displ_all eigen functions, shape(size,nc_all)
 * @param k eigenvalues, wavenumber
 * @param c filtered phase velocity, shape(nc)
 * @param displ filtered eigen function, shape(nc,size)
 */
void 
filter_swd(double vmin, double vmax,double om, const Eigen::MatrixXcd &displ_all,
            const Eigen::Array<std::complex<double>,-1,1> &k,std::vector<double> &c,
           std::vector<double> &displ)
{
    Eigen::Array<double,-1,1> c_all = (om / k).real();
    //std::cout << c_all.transpose() << "\n";

    // filter swd in [vmin * 0.85,vmax] region
    Eigen::Array<double,-1,1> c_filt;
    Eigen::MatrixXd displ_filt;
    using Eigen::all;
    //auto mask = ((c_all >= 0.85 * vmin) && (c_all <= vmax)) && (k.real().abs() >= k.imag().abs());
    auto mask = ((c_all >= vmin)&& (c_all <= vmax)) && (k.real().abs() > 5 * k.imag().abs());
    // std::cout << mask << "\n";
    std::vector<int> idx0; idx0.reserve(mask.cast<int>().sum());
    for(int i = 0; i < c_all.size(); i ++) {
        if(mask[i]) {
            idx0.push_back(i);
        }
    }

    int nc = idx0.size();
    int size = displ_all.rows();
    c_filt.resize(nc); displ_filt.resize(size,nc);
    for(int i = 0; i < nc; i ++) {
        c_filt[i] = c_all[idx0[i]];
        displ_filt(all,i) = displ_all(all,idx0[i]).real();
    }
    //printf("%g %g\n",c_filt.minCoeff(),c_filt.maxCoeff());

    // sort according to ascending order 
    std::vector<int> idx;
    idx.resize(nc);
    for(int i = 0; i < nc; i ++ ) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(),
        [&c_filt](size_t i1, size_t i2) {return c_filt[i1] < c_filt[i2];}); 

    // copy to c/displ
    c.resize(nc); displ.resize(nc * size);
    for(int ic = 0; ic < nc; ic ++) {
        c[ic] = c_filt[idx[ic]];
        for(int i = 0; i < size; i ++) {
            displ[ic * size + i] = displ_filt(i,idx[ic]);
        }
    } 
}
