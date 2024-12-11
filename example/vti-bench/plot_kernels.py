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
    vsv_kl = fio[f'kernels/{Tid}/mode{mode}/vsv_kl'][:]
    vsh_kl = fio[f'kernels/{Tid}/mode{mode}/vsh_kl'][:]
    rho_kl = fio[f'kernels/{Tid}/mode{mode}/rho_kl'][:]
    phi_kl = fio[f'kernels/{Tid}/mode{mode}/phi_kl'][:]
    theta_kl = fio[f'kernels/{Tid}/mode{mode}/theta_kl'][:]
    eta_kl = fio[f'kernels/{Tid}/mode{mode}/eta_kl'][:]
    U = fio[f'kernels/{Tid}/mode{mode}/U'][:]
    V = fio[f'kernels/{Tid}/mode{mode}/V'][:]
    W = fio[f'kernels/{Tid}/mode{mode}/W'][:]


    plt.figure(1,figsize=(14,20))
    plt.subplot(131)
    plt.plot(vsv_kl,-z,label='vsv_kl')
    plt.ylabel("depth,km")
    plt.legend()

    plt.subplot(132)
    plt.plot(rho_kl,-z,label='rho_kl')
    #plt.ylabel("depth,km")
    plt.legend()

    plt.subplot(133)
    plt.plot(phi_kl,-z,label='phi_kl')
    #plt.ylabel("depth,km")
    plt.legend()

    plt.savefig("kernel1.jpg")

    plt.figure(2,figsize=(14,20))
    plt.subplot(131)
    plt.plot(vsh_kl,-z,label='vsh_kl')
    plt.ylabel("depth,km")
    plt.legend()

    plt.subplot(132)
    plt.plot(eta_kl,-z,label='eta_kl')
    #plt.ylabel("depth,km")
    plt.legend()

    plt.subplot(133)
    plt.plot(theta_kl,-z,label='theta_kl')
    #plt.ylabel("depth,km")
    plt.legend()

    plt.savefig("kernel2.jpg")

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

    plt.savefig("eigenfunc.jpg")


main()