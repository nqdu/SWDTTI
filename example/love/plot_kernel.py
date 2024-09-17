import numpy as np 
import matplotlib.pyplot as plt 
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: ./plot_kernel.py filename")
        exit(1)
    
    d = np.loadtxt("out/" + sys.argv[1])

    plt.figure(1,figsize=(14,20))
    plt.subplot(131)
    plt.plot(d[:,1],-d[:,0],label='vsv_kl')
    plt.ylabel("depth,km")
    plt.legend()

    plt.subplot(132)
    plt.plot(d[:,2],-d[:,0],label='rho_kl')
    #plt.ylabel("depth,km")
    plt.legend()

    plt.subplot(133)
    plt.plot(d[:,3],-d[:,0],label='W')
    #plt.ylabel("depth,km")
    plt.legend()

    plt.savefig("kernel.jpg")


main()