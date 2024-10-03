import numpy as np 
from scipy.io import FortranFile
import h5py
import sys 

def main():
    if len(sys.argv) != 4:
        print("Usage: python binary2h5.py binfile swdfile outfile")
        exit(1)
    
    # get input 
    binfile = sys.argv[1]
    swdfile = sys.argv[2]
    outfile = sys.argv[3]

    # open outfile
    fout:h5py.File = h5py.File(outfile,"w")

    # open swd file to read T and swd
    T = np.loadtxt(swdfile,max_rows=1)
    data = np.loadtxt(swdfile,skiprows=1)
    max_modes = int(np.max(data[:,-1])) + 1
    fout.create_group("swd")
    for imode in range(max_modes):
        gname = f"swd/mode{imode}/"
        fout.create_group(f"{gname}")
        idx = data[:,-1] == imode 
        data1 = data[idx,:]
        nt1 = data1.shape[0]

        fout.create_dataset(f"{gname}/T",shape = (nt1),dtype='f4')
        fout.create_dataset(f"{gname}/c",shape = (nt1),dtype='f4')
        fout.create_dataset(f"{gname}/u",shape = (nt1),dtype='f4')
        fout[f'{gname}/T'][:] = T[np.int32(data[idx,0])] 
        fout[f'{gname}/c'][:] = data1[:,1]
        fout[f'{gname}/u'][:] = data1[:,2]

    # write kernels
    fin:FortranFile = FortranFile(binfile,"r")
    nkers = fin.read_ints('i4')[0]
    ncomps = fin.read_ints('i4')[0]
    assert(nkers in [3,5,8])
    assert(ncomps in [1,2,3])
    if nkers == 3:
        kl_name = ['vsh_kl','vsv_kl','rho_kl']
        comp_name = ['W']
        fout.attrs['WaveType'] = 'Love'
        fout.attrs['ModelType'] = 'VTI'
    elif nkers == 5:
        comp_name = ['U','V']
        kl_name = ['vph_kl','vpv_kl','vsv_kl','eta_kl','rho_kl']
        fout.attrs['WaveType'] = 'Rayleigh'
        fout.attrs['ModelType'] = 'VTI'
    else:
        comp_name = ['U','V','W']
        kl_name = ['vph_kl','vpv_kl','vsh_kl','vsv_kl','eta_kl','theta_kl','phi_kl','rho_kl']
        fout.attrs['WaveType'] = 'Full'
        fout.attrs['ModelType'] = 'TTI'
    fout.create_group("kernels")

    for it in range(len(T)):
        idx = data[:,0] == it 
        data1 = data[idx,:]
        max_mode = int(np.max(data1[:,-1])) + 1
        #fout.attrs[f"kernels/{it}/T"] = T[it]

        # read coordinates
        zcords = fin.read_reals('f8')
        npts = zcords.size
        fout.create_dataset(f"kernels/{it}/zcords",dtype='f4',shape =(npts))
        fout[f'kernels/{it}/zcords'][:] = zcords[:]
        for imode in range(max_mode):
            gname = f"kernels/{it}/mode{imode}"
            fout.create_group(gname)

            # read eigenfuncs
            displ = fin.read_reals('f8')
            displ = displ.reshape((ncomps,npts))
            for icomp in range(ncomps):
                fout.create_dataset(f"{gname}/{comp_name[icomp]}",dtype='f8',shape=(npts))
                fout[f"{gname}/{comp_name[icomp]}"][:] = displ[icomp,:]
            
            # read kernels
            kernel = fin.read_reals('f8').reshape((nkers,npts))
            for iker in range(nkers):
                fout.create_dataset(f"{gname}/{kl_name[iker]}",dtype='f8',shape=(npts))
                fout[f"{gname}/{kl_name[iker]}"][:] = kernel[iker,:]
            

    # close 
    fin.close()
    fout.close()

    

if __name__ == "__main__":
    main()