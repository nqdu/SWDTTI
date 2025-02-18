#include "vti/swdlayervti.hpp"
#include "swdio.hpp"

#include <iostream>
#include <fstream>

int main (int argc, char **argv){
    // read model name
    if(argc != 6 && argc != 7) {
        printf("Usage: ./surfvti modelfile wavetype[1 for Love and 2 for Rayleigh] f1 f2 nt [is_layered=1]\n");
        printf("freqs = logspace(log10(f1),log10(f2),nt)\n");
        exit(1);
    }
    int wavetype;
    sscanf(argv[2],"%d",&wavetype);

    // read flag if required
    bool is_layer = true;
    if(argc == 7) {
        int flag;
        sscanf(argv[6],"%d",&flag);
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
        printf("layer %d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",
                i + 1,thk[i],rho[i],vpv[i],
               vph[i],vsv[i],vsh[i],eta[i]);
    }
    infile.close();

    // Period
    int nt;
    float f1,f2;
    sscanf(argv[3],"%g",&f1); sscanf(argv[4],"%g",&f2);
    sscanf(argv[5],"%d",&nt);
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
    if(wavetype == 1) {
        printf("\ncomputing dispersions for Love wave ...\n");
    }
    else {
        printf("\ncomputing dispersions for Rayleigh wave ...\n");
    }
    printf("freqmin = %g freqmax = %g\n",freq[0],freq[nt-1]);

    LayerModelVTI model;
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
    int nkers = 3, ncomp = 1;
    if(wavetype != 1) {
        nkers = 5;
        ncomp = 2;
    }
    write_binary_f(fio,&nkers,1);
    write_binary_f(fio,&ncomp,1);
    
    for(int it = 0; it < nt; it ++) {
        std::vector<double> c,displ,u;

        // phase velocity/eigenfunctions
        model.create_database(freq[it],nz,rho.data(),vpv.data(),vph.data(),
                              vsv.data(),vsh.data(),eta.data(),thk.data(),
                              is_layer);
        model.prepare_matrices(wavetype);
        switch (wavetype)
        {
        case 1:
            model.compute_slegn(freq[it],c,displ);
            break;
        
        default:
            model.compute_sregn(freq[it],c,displ);
            break;
        }
    
        // alloc space for group/kernels
        std::vector<double> frekl;
        int nc = c.size();
        u.resize(nc);

        // get some consts
        int nglob = model.nglob;
        int npts = model.ibool.size();

        // write coordinates
        write_binary_f(fio,model.znodes.data(),npts);

        for(int ic = 0; ic < nc; ic ++) {
            if(wavetype == 1) {
                u[ic] = model.compute_love_kl(freq[it],c[ic],&displ[ic * nglob],frekl);
            }
            else {
                u[ic] = model.compute_rayl_kl(freq[it],c[ic],&displ[ic * nglob * 2],frekl);
            }
            model.transform_kernels(frekl);

            // write swd T,c,U,mode
            fprintf(fp,"%d %lf %lf %d\n",it,c[ic],u[ic],ic);

            // write displ
            std::vector<double> temp(npts*ncomp);
            for(int i = 0; i < npts; i ++) {
                int iglob = model.ibool[i];
                for(int j = 0; j < ncomp; j ++) {
                    temp[j * npts + i] = displ[ic * nglob * ncomp + j * nglob + iglob];
                }
            }
            write_binary_f(fio,temp.data(),npts*ncomp);

            // write kernels
            write_binary_f(fio,&frekl[0],npts*nkers);
        }
    }
    fclose(fp);
    fclose(fio);
}