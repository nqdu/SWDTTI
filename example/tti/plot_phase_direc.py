import numpy as np
import matplotlib.pyplot as plt 
from glob import glob

filenames = glob("out/swd_direc*")
nd = len(filenames)

deg = np.zeros((nd + 1))
c = np.zeros((nd + 1))

for i in range(nd):
    data = np.loadtxt(f"out/swd_direc{i}.txt",ndmin=2)
    deg[i],c[i] = data[0,0],data[0,2]
deg[-1] = deg[0]
c[-1] = c[0]
print(deg)
# plot
fig = plt.figure(1,figsize=(8,8))
ax = fig.add_subplot(1,1,1,projection='polar')
x = np.cos(deg * np.pi / 180) * c 
y = np.sin(deg * np.pi / 180) * c 
x1 = np.cos(deg * np.pi / 180) * c[0]
y1 = np.sin(deg * np.pi / 180) * c[0]
ax.plot(deg * np.pi / 180,c,label='phase velocity')
ax.plot(deg * np.pi / 180,c*0+c[0],label='phase Ref')
#ax.set_yticks([])
#plt.scatter(x,y,label='phase velocity')
#plt.plot(x1,y1,label='phase Ref')
plt.title("HTI model")
plt.legend()
plt.gca().set_aspect('equal')
plt.savefig("test.jpg")