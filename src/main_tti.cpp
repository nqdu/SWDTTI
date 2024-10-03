#include "swdlayertti.hpp"
#include <iostream>

int main (int argc, char **argv){
    // read model name
    if(argc != 2) {
        printf("Usage: ./surf_tti modelfile ");
        exit(1);
    }

    // read model
    printf("reading velocity model %s:\n",argv[1]);
    printf("layer number\t thick\t rho\t vsv\t vsh\t vpv\t vph\t theta\t phi  \n");
    std::vector<float> thk,vpv,vph,vsv,vsh,rho,theta0,phi0,eta;
    int nz;
    FILE *fp = fopen(argv[1],"r");
    fscanf(fp,"%d",&nz);
    thk.resize(nz); vsv.resize(nz); vsh.resize(nz); rho.resize(nz);
    vpv.resize(nz); vph.resize(nz); eta.resize(nz); theta0.resize(nz);
    phi0.resize(nz);
    for(int i = 0; i < nz; i ++) {
        fscanf(fp,"%f%f%f%f%f%f%f%f",&thk[i],&rho[i],&vsv[i],
                &vsh[i],&vpv[i],&vph[i],&theta0[i],&phi0[i]);
        printf("layer %d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",
                i + 1,thk[i],rho[i],vsv[i],
               vsh[i],vpv[i],vph[i],theta0[i],phi0[i]);
        eta[i] = 1.;
        theta0[i] *= M_PI / 180.;
        phi0[i] *= M_PI / 180.;
    }
    fclose(fp);

    // Period
    int nt = 120;
    std::vector<double> freq(nt);
    for(int it = 0; it < nt; it ++) {
        double f = -2.2 + 2.2 / (nt - 1) * it;
        freq[nt-1 - it] = std::pow(10,f);
    }

    // create database
    printf("\ncomputing dispersions ...\n");
    LayerModelTTI model;
    model.initialize();

    // angles for each direction
    const int ndirec = 100;
    float phi[ndirec];
    for(int i = 0; i < ndirec; i ++) {
        phi[i] = 360. / ndirec * i;
    }

    // loop every direction
    for(int id = 0; id < ndirec; id ++) {
        printf(" direction %d of %d\n",id + 1,ndirec);
        
        // open file to save swd for each direction
        std::string filename = "out/swd_direc" + std::to_string(id) + ".txt";
        fp = fopen(filename.c_str(),"w");
        
        std::vector<double> c,u;
        std::vector<std::complex<double>> displ;
        for(int it = 100; it < 101; it ++) {
            model.create_database(
                freq[it],nz,vph.data(),vpv.data(),vsh.data(),
                vsv.data(),eta.data(),theta0.data(),phi0.data(),
                rho.data(),thk.data(),true
            );
            model.prepare_matrices(phi[id]);
            model.compute_egnfun(freq[it],phi[id],c,displ);

            // write dispersion in file
            int nc = c.size();
            for(int ic = 0; ic < nc; ic ++) {
                fprintf(fp,"%g %g %g\n",phi[id],1. / freq[it],c[ic]);
            }

            // // group
            // std::vector<double> frekl;
            // int nc = c.size();
            // u.resize(nc);
            // int nglob = model.nglob;
            // for(int ic = 0; ic < nc; ic ++) {
            //     u[ic] = model.compute_kernels(freq[it],c[ic],phi[id],displ,frekl);
            //     model.transform_kernels(frekl);

            //     // write swd
            //     fprintf(fp,"%g %g %g\n",1./freq[it],c[ic],u[ic]);

            //     if(ic == 0) {
            //         // get constants
            //         int n = model.ibool.size();
            //         int nkers = frekl.size() / n;
            //         int ncomp = displ.size() / (nglob * nc);

            //         std::string filename = "out/" + std::to_string(it) + ".txt";
            //         FILE *fp1 = fopen(filename.c_str(),"w");
                    
            //         for(int i = 0; i < n; i ++) {
            //             int iglob = model.ibool[i];
                        
            //             fprintf(fp1,"%g ",model.znodes[i]);
            //             for(int j = 0; j < nkers; j ++) {
            //                 fprintf(fp1,"%g ",frekl[j * n + i]);
            //             }
            //             fprintf(fp1,"\n");
                        
            //         }
            //         fclose(fp1);
            //     }
            // }
        }

        fclose(fp);
    }
    
}