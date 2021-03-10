import numpy as np 
import matplotlib.pyplot as plt 

maxdep=100000

fig = plt.figure(1)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
for i in [20,30,40,50,60]:
    d = np.loadtxt("disp_ray/"+str(i) +".dat")
    m = d[0,2]
    d[:,1:] = d[:,1:] / m
    idx = d[:,0] <=maxdep
    ax1.plot(d[idx,1],d[idx,0]/1000,label=str(i))
    ax1.set_xlabel("$r_1$")
    ax1.invert_yaxis()
    ax1.set_yticks([0,19,38,100])
    ax1.set_ylabel("depth,km")
    ax1.set_title("Horizontal Displacement ")
    for j in [19,38]:
        ax1.axhline(j,ls="--",color='k')
    ax1.legend()

    ax2.plot(d[idx,3]/1.0e6,d[idx,0]/1000,label=str(i))
    ax2.set_xlabel("$r_3$")
    ax2.invert_yaxis()
    ax2.set_yticks([])
    #ax2.set_xticks([-2.5,0])
    ax2.set_title("Shear Stress")
    for j in [19,38]:
        ax2.axhline(j,ls="--",color='k')
    ax2.legend()

fig = plt.figure(2)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
for i in [20,30,40,50,60]:
    d = np.loadtxt("disp_ray/"+str(i) +".dat")
    idx = d[:,0] <=maxdep
    ax1.plot(d[idx,2],d[idx,0]/1000,label=str(i))
    ax1.set_xlabel("$r_2$")
    ax1.invert_yaxis()
    ax1.set_yticks([0,19,38,100])
    ax1.set_ylabel("depth,km")
    ax1.set_title("Vertical Displacement")
    for j in [19,38]:
        ax1.axhline(j,ls="--",color='k')
    ax1.legend()

    ax2.plot(d[idx,4]/1.0e6,d[idx,0]/1000,label=str(i))
    ax2.set_yticks([])
    ax2.set_xlabel("$r_4$")
    ax2.invert_yaxis()
    ax2.set_title("Normal Stress")
    #ax2.set_xticks([-2.5,0])
    for j in [19,38]:
        ax2.axhline(j,ls="--",color='k')
    ax2.legend()

fig = plt.figure(3)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
for i in [20,30,40,50,60]:
    d = np.loadtxt("disp_love/"+str(i) +".dat")
    m = d[0,1]
    d[:,1:] = d[:,1:] / m
    idx = d[:,0] <=maxdep
    ax1.plot(d[idx,1],d[idx,0]/1000,label=str(i))
    ax1.set_xlabel("$l_1$")
    ax1.invert_yaxis()
    ax1.set_yticks([0,19,38,100])
    ax1.set_ylabel("depth,km")
    ax1.set_title("Displacement")
    for j in [19,38]:
        ax1.axhline(j,ls="--",color='k')
    ax1.legend()

    ax2.plot(-d[idx,2]/1.0e6,d[idx,0]/1000,label=str(i))
    ax2.set_xlabel("$l_2$")
    ax2.invert_yaxis()
    ax2.set_yticks([])
    ax2.set_title("Stress")
    for j in [19,38]:
        ax2.axhline(j,ls="--",color='k')
    ax2.legend()

plt.show()
