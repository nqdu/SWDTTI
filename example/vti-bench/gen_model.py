import numpy as np

#z = np.linspace(0,50,5)
z = np.array([0,30.])
nz = len(z)
#vs = 3.5 + 0.02 * z 
vs = np.array([3.5,4.5])

# vs[10:20] = 1.4
# vs[40:50] = 2.5

vp = 1.732 * vs 
rho = 0.3601 * vp + 0.541 
thk = np.zeros_like(z)
thk[0:nz-1] = np.diff(z)

# vpv = vp * 1.
# vph = vp * 1.2
# vsv = vs * 1.
# vsh = vs * 1.2

vpv = vp * 1.
vph = vp * 1.
vsv = vs * 1.
vsh = vs * 1.

f = open("model.txt","w")
f.write("%d\n"%(nz))
for i in range(nz):
    f.write("%f %f %f %f %f %f 1. %f %f\n"%(thk[i],rho[i],vsv[i],vsh[i],vpv[i],vph[i],0.,0))
f.close()