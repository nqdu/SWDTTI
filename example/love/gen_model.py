import numpy as np

# z = np.linspace(0,50,10)
# nz = len(z)
# vs = 3.5 + 0.02 * z 
z = np.array([0,30])
nz = len(z)
vs = np.array([3.0,4.0])

# vs[10:20] = 1.4
# vs[40:50] = 2.5

vp = 1.732 * vs 
rho = 0.3601 * vp + 0.541 
thk = np.zeros_like(z)
thk[0:nz-1] = np.diff(z)

f = open("model.txt","w")
f.write("%d\n"%(nz))
for i in range(nz):
    f.write("%f %f %f %f %f %f 1.\n"%(thk[i],rho[i],vp[i],vp[i],vs[i],vs[i]))
f.close()