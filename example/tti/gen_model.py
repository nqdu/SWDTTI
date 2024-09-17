import numpy as np

z = np.linspace(0,50,5)
nz = len(z)
vs = 3.5 + 0.02 * z 

# vs[10:20] = 1.4
# vs[40:50] = 2.5

vp = 1.732 * vs 
rho = 0.3601 * vp + 0.541 
thk = np.zeros_like(z)
thk[0:nz-1] = np.diff(z)

vpv = vp * 1.
vph = vp * 1.2
vsv = vs * 1.
vsh = vs * 1.2

f = open("modelhti.txt","w")
f.write("%d\n"%(nz))
for i in range(nz):
    f.write("%f %f %f %f %f %f %f %f\n"%(thk[i],rho[i],vsv[i],vsh[i],vpv[i],vph[i],90.,0))
f.close()