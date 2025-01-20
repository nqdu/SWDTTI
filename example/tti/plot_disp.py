import numpy as np 
import matplotlib.pyplot as plt 

def fetch_data(T,data,nc,col_id):
    nt = len(T)
    c_all = np.zeros((nt,nc)) - 1

    for it in range(len(T)):
        idx = data[:,0] == it
        c = data[idx,col_id]
        for i in range(len(c)):
            if i >= c_all.shape[1] : break 
            c_all[it,i] = c[i]

    return c_all

# load swd data
data = np.loadtxt("out/swd.txt",skiprows=1)

# find unique T
T = np.loadtxt("out/swd.txt",max_rows=1)
nt = len(T)

# plot phase
nc = 3
plt.figure(1,figsize=(14,5))

c_all = fetch_data(T,data,nc,1)
u_all = fetch_data(T,data,nc,2)

for i in range(c_all.shape[1]):
    idx = c_all[:,i] > 0
    if np.sum(idx) != 0:
        plt.plot(1./T[idx],c_all[idx,i],label=f'phase_{i}',ls='--',color='k')
        plt.plot(1./T[idx],u_all[idx,i],label=f'group_{i}',color='r')


plt.legend()
plt.xlabel("Frequency,Hz")
plt.ylabel("Phase/group Velocity, km/s")
plt.savefig("swd.jpg")
plt.clf()


    