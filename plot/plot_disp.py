import numpy as np 
import matplotlib.pyplot as plt 
plt.figure(1)
d = np.loadtxt("dispersion.dat")
plt.plot(d[:,0],d[:,1],d[:,0],d[:,2])
plt.legend(["Phase","Group"])
plt.title("Dispersion Curve for Love Wave")

plt.figure(2)
plt.plot(d[:,0],d[:,3],d[:,0],d[:,4])
plt.legend(["Phase","Group"])
plt.title("Dispersion Curve for Rayleigh Wave")
plt.show()