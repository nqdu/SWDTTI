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
data1 = np.loadtxt("out/swd.cps.txt",skiprows=1)

# find unique T
T = np.loadtxt("out/swd.txt",max_rows=1)
T1 = np.loadtxt("out/swd.cps.txt",max_rows=1)
nt = len(T)

# plot phase
nc = int(np.max(data[:,-1])) + 1
plt.figure(1,figsize=(14,5))
#plt.scatter(1./data[:,0],data[:,1],s=10,color='k')
c_all = fetch_data(T,data,nc,1)
for i in range(c_all.shape[1]):
    idx = c_all[:,i] > 0
    if np.sum(idx) != 0:
        plt.plot(1./T[idx],c_all[idx,i])

c1_all = fetch_data(T1,data1,nc,1)
for i in range(c1_all.shape[1]):
    idx = c1_all[:,i] > 0
    if np.sum(idx) != 0:
        plt.plot(1./T1[idx],c1_all[idx,i],color='k',ls='--',label='cps')

plt.legend()
plt.xlabel("Frequency,Hz")
plt.ylabel("Phase Velocity, km/s")
plt.savefig("phase.jpg")
plt.clf()

# plot group
plt.figure(1,figsize=(14,5))
#plt.scatter(1./data[:,0],data[:,2],s=10,color='k')
c_all = fetch_data(T,data,nc,2)

# plot
for i in range(c_all.shape[1]):
    idx = c_all[:,i] > 0
    if np.sum(idx) != 0:
        plt.plot(1./T[idx],c_all[idx,i])

c1_all = fetch_data(T1,data1,nc,2)
for i in range(c1_all.shape[1]):
    idx = c1_all[:,i] > 0
    if np.sum(idx) != 0:
        plt.plot(1./T1[idx],c1_all[idx,i],color='k',ls='--',label='cps')
plt.legend()
plt.xlabel("Frequency,Hz")
plt.ylabel("Group Velocity, km/s")
plt.savefig("group.jpg")

    