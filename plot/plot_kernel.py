import numpy as np 
import matplotlib.pyplot as plt 

maxdep = 400
fig = plt.figure(1,figsize=(12,6))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
for i,t in enumerate([20,30,40,50,60]):
    
    d = np.loadtxt("disp_ray/" + str(t) + ".ker.dat")
    r = d[:,0] / 1000
    idx = r < maxdep
    ax1.plot(d[idx,1],r[idx],label=str(t))
    ax1.set_xlim([-0.1,0.1])
    ax1.set_ylim([0,maxdep])
    ax1.set_xlabel(r"$\partial{c}/\partial{\beta}$")
    ax1.invert_yaxis()
    ax1.set_ylabel("Depth, km")
    ax1.legend()

    d = np.loadtxt("disp_ray/" + str(t) + ".ker.dat")
    ax2.plot(d[idx,2],r[idx],label=str(t))
    ax2.set_xlim([-0.1,0.1])
    ax2.set_ylim([0,maxdep])
    ax2.set_xlabel(r"$\partial{c}/\partial{\alpha}$")
    ax2.invert_yaxis()
    ax2.set_yticks([])
    ax2.legend()

    d = np.loadtxt("disp_ray/" + str(t) + ".ker.dat")
    ax3.plot(d[idx,3],r[idx],label=str(t))
    ax3.set_xlim([-0.1,0.1])
    ax3.set_ylim([0,maxdep])
    ax3.set_xlabel(r"$\partial{c}/\partial{\rho}$")
    ax3.invert_yaxis()
    ax3.set_yticks([])
    ax3.legend()

fig.suptitle("Sensitivity Kernel for Gutenberg Model")
plt.show()