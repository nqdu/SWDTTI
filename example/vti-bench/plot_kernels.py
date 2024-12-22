import numpy as np 
import matplotlib.pyplot as plt 
import h5py
import sys 

def main():
    if len(sys.argv) != 4:
        print("Usage: ./plot_kernel.py filename period_id mode")
        exit(1)
    
    # open file 
    fio = h5py.File(sys.argv[1],"r")
    Tid = int(sys.argv[2])
    mode = int(sys.argv[3])

    # fetch data
    if f'kernels/{Tid}/mode{mode}' not in fio:
        print("Tid = {Tid} and mode ={mode} don't exist!")
        exit(1)
    z = fio[f'kernels/{Tid}/zcords'][:]
    kl_names = [['vsv_kl','vpv_kl','vsh_kl','vph_kl'],
               ['eta_kl','phi_kl','theta_kl','rho_kl']]
    U = fio[f'kernels/{Tid}/mode{mode}/U'][:]
    V = fio[f'kernels/{Tid}/mode{mode}/V'][:]
    W = fio[f'kernels/{Tid}/mode{mode}/W'][:]

    for i in range(2):
        fig,axes = plt.subplots(1,4,sharey=True,figsize=(12,16))
        for j in range(4):
            name = kl_names[i][j]
            kl = fio[f'kernels/{Tid}/mode{mode}/{name}'][:]

            axes[j].plot(kl,-z,label=f'{name}')
            if j== 0:
                axes[j].set_ylabel("depth,km")
            axes[j].legend()
            axes[j].set_xlabel(f"{name}")
        fig.savefig(f"kernel{i+1}.mode{mode}.jpg",dpi=300)

    plt.figure(3,figsize=(14,20))
    plt.subplot(131)
    plt.plot(U.real,-z,label='U')
    plt.ylabel("depth,km")
    plt.legend()

    plt.subplot(132)
    plt.plot(W.real,-z,label='W')
    #plt.ylabel("depth,km")
    plt.legend()

    plt.subplot(133)
    plt.plot(V.real,-z,label='V')
    #plt.ylabel("depth,km")
    plt.legend()

    plt.savefig(f"eigenfunc.mode{mode}.jpg",dpi=300)


main()