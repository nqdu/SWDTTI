#include "swdlayer.hpp"
#include <iostream>

int main (int argc, char **argv){
    // read model name
    if(argc != 3) {
        printf("Usage: ./test_love modelfile wavetype[1 for Love and 2 for Rayleigh] ");
        exit(1);
    }
    int wavetype;
    sscanf(argv[2],"%d",&wavetype);

    // read model
    printf("reading velocity model:\n");
    printf("layer number\t thick\t rho\t vs\t vp  \n");
    std::vector<float> thk,vp,rho,vs,eta;
    int nz;
    FILE *fp = fopen(argv[1],"r");
    fscanf(fp,"%d",&nz);
    thk.resize(nz); vs.resize(nz); rho.resize(nz);
    vp.resize(nz);  eta.resize(nz);
    for(int i = 0; i < nz; i ++) {
        fscanf(fp,"%f%f%f%f",&thk[i],&rho[i],&vs[i],&vp[i]);
        printf("layer %d\t %f\t %f\t %f\t %f\n",i + 1,thk[i],rho[i],vs[i],vp[i]);
        eta[i] = 1.;
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
    LayerModel model;
    model.initialize();

    // love wave swd
    fp = fopen("out/swd.txt","w");
    std::vector<double> c,displ,u;
    for(int it = 0; it < nt; it ++) {
        model.create_database(freq[it],nz,vp.data(),vp.data(),
                                vs.data(),vs.data(),eta.data(),
                                rho.data(),thk.data());
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

        // group
        std::vector<double> frekl;
        int nc = c.size();
        u.resize(nc);
        int nglob = model.nglob;
        //printf("%d of %d %f\n",it,nt,1./freq[it]);
        for(int ic = 0; ic < nc; ic ++) {
            if(wavetype == 1) {
                u[ic] = model.compute_love_kl(freq[it],c[ic],&displ[ic * nglob],frekl);
            }
            else {
                u[ic] = model.compute_rayl_kl(freq[it],c[ic],&displ[ic * nglob * 2],frekl);
            }
            model.transform_kernels(frekl);

            // write swd
            fprintf(fp,"%g %g %g\n",1./freq[it],c[ic],u[ic]);

            if(ic == 0) {
                // get constants
                int n = model.ibool.size();
                int nkers = frekl.size() / n;
                int ncomp = displ.size() / (nglob * nc);

                std::string filename = "out/" + std::to_string(it) + ".txt";
                FILE *fp1 = fopen(filename.c_str(),"w");
                
                for(int i = 0; i < n; i ++) {
                    int iglob = model.ibool[i];
                    
                    fprintf(fp1,"%g ",model.znodes[i]);
                    for(int j = 0; j < nkers; j ++) {
                        fprintf(fp1,"%g ",frekl[j * n + i]);
                    }
                    for(int j = 0; j < ncomp; j ++) fprintf(fp1,"%g ",displ[ic * nglob * 2 + j * nglob +  iglob]);
                    fprintf(fp1,"\n");
                    
                }
                fclose(fp1);
            }
        }
    }
    fclose(fp);
}