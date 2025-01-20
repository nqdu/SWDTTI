#include "multiphysics/vti_acoustic.hpp"
#include "swdio.hpp"

#include <iostream>
#include <fstream>

int main (int argc, char **argv){
    // read model name
    if(argc != 5 && argc != 6) {
        printf("Usage: ./surfvti_ac modelfile f1 f2 nt [is_layered=1]\n");
        printf("freqs = logspace(log10(f1),log10(f2),nt)\n");
        exit(1);
    }

    // read flag if required
    bool is_layer = true;
    if(argc == 6) {
        int flag;
        sscanf(argv[5],"%d",&flag);
        is_layer = (flag == 1);
    }

    // read model
    if(is_layer) {
        printf("reading layered velocity model:\n");
    }
    else {
        printf("reading continuous velocity model:\n");
    }
    
    printf("layer number\t thick\t rho\t vpv\t vph\t vsv\t vsh\t eta\n");

    std::vector<float> thk,vpv,vph,rho,vsv,vsh,eta;
    int nz;
    std::ifstream infile; infile.open(argv[1]);
    infile >> nz;
    thk.resize(nz); vsv.resize(nz); rho.resize(nz);
    vpv.resize(nz); vph.resize(nz); eta.resize(nz);
    vsh.resize(nz);
    for(int i = 0; i < nz; i ++) {
        infile >> thk[i] >> rho[i] >> vpv[i] >>
                  vph[i] >> vsv[i] >> vsh[i] >> 
                  eta[i];
        printf("layer %03d %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
                i + 1,thk[i],rho[i],vpv[i],
               vph[i],vsv[i],vsh[i],eta[i]);
    }
    infile.close();

    // Period
    int nt;
    float f1,f2;
    sscanf(argv[2],"%g",&f1); sscanf(argv[3],"%g",&f2);
    sscanf(argv[4],"%d",&nt);
    f1 = std::log10(f1); f2 = std::log10(f2);
    if(f1 > f2) std::swap(f1,f2);
    std::vector<double> freq(nt);
    for(int it = 0; it < nt; it ++) {
        double coef = (nt - 1);
        if(coef == 0.) coef = 1.;
        coef = 1. / coef;
        double f = f1 + (f2 - f1) * coef * it;
        freq[it] = std::pow(10,f);
    }

    // create database
    printf("\ncomputing dispersions for Rayleigh wave, with acoustic layers ...\n");
    printf("freqmin = %g freqmax = %g\n",freq[0],freq[nt-1]);

    LayerModelMultiPhyVTI model;
    model.initialize();

    // open file to write out data
    FILE *fp = fopen("out/swd.txt","w");
    FILE *fio = fopen("out/database.bin","wb");

    // write period vector into fp
    for(int it = 0; it < nt; it ++) {
        fprintf(fp,"%g ",1. / freq[it]);
    }
    fprintf(fp,"\n");

    // write meta info into database
    int nkers = 5, ncomp = 2;
    write_binary_f(fio,&nkers,1);
    write_binary_f(fio,&ncomp,1);
    
    for(int it = 0; it < nt; it ++) {
        std::vector<double> c,egn,u;

        // phase velocity/eigenfunctions
        model.create_database(freq[it],nz,rho.data(),vpv.data(),vph.data(),
                              vsv.data(),vsh.data(),eta.data(),thk.data(),
                              is_layer);
        model.prepare_matrices(freq[it]);
        model.compute_egnfun(freq[it],c,egn);

        // alloc space for group/kernels
        std::vector<double> frekl_el,frekl_ac;
        int nc = c.size();
        u.resize(nc);

        if(it == 0) {
            for (auto a : model.ibool_ac) {
                printf("%d ",a);
            }
            std::cout << "\n" <<  model.ac_elmnts[0] << "\n";
        }

        // get some consts
        int ng = model.nglob_ac + model.nglob_el * 2;
        int npts = model.ibool.size();

        // write coordinates
        write_binary_f(fio,model.znodes.data(),npts);

        for(int ic = 0; ic < nc; ic ++) {
            u[ic] = model.compute_kernels(freq[it],c[ic],&egn[ic*ng],frekl_el,frekl_ac);

            // write swd T,c,U,mode
            fprintf(fp,"%d %lf %lf %d\n",it,c[ic],u[ic],ic);

        //     // write displ
        //     std::vector<double> temp(npts*ncomp);
        //     for(int i = 0; i < npts; i ++) {
        //         int iglob = model.ibool[i];
        //         for(int j = 0; j < ncomp; j ++) {
        //             temp[j * ncomp + i] = displ[ic * nglob * ncomp + j * nglob + iglob];
        //         }
        //     }
        //     write_binary_f(fio,temp.data(),npts*ncomp);

        //     // write kernels
        //     write_binary_f(fio,&frekl[0],npts*nkers);
        }
    }
    fclose(fp);
    fclose(fio);
}