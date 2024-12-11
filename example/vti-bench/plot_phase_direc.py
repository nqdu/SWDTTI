import numpy as np
import matplotlib.pyplot as plt 
from glob import glob

filenames = glob("out/swd.*.txt")
nd = len(filenames)

deg = np.zeros((nd))
deg_u = deg * 1.
c = np.zeros((nd))
u = c * 1.

for i in range(nd):
    data = np.loadtxt(f"out/swd.{i}.txt",skiprows=1)
    deg[i] = data[0,3]
    c[i],u[i] = data[0,1],data[0,2]
    deg_u[i] = data[0,4]
# deg[-1] = deg[0]
# c[-1] = c[0]
# u[-1] = u[0]
# deg_u[-1] = deg_u[0]

#print(deg)
print(deg_u, deg)
# plot
fig = plt.figure(1,figsize=(8,8))
ax = fig.add_subplot(1,1,1,projection='polar')
ax.plot(deg * np.pi / 180,c,label='phase velocity')
ax.plot(deg * np.pi / 180,c*0+c[0],label='phase Ref')

plt.title("HTI model")
plt.legend()
plt.gca().set_aspect('equal')
plt.savefig("phase-direc.jpg")
plt.clf()

fig = plt.figure(2,figsize=(8,8))
ax = fig.add_subplot(1,1,1,projection='polar')
ax.plot(deg * np.pi / 180,c,label='phase velocity')
ax.plot(deg * np.pi / 180,c*0+c[0],label='phase Ref')
ax.plot(deg_u * np.pi / 180,u,label="group velocity")
#ax.plot(deg * np.pi / 180,u*0+u[0],label='group Ref')

plt.title("HTI model")
plt.legend()
plt.gca().set_aspect('equal')
plt.savefig("group-direc.jpg")

fig = plt.figure(3,figsize=(8,8))
ax = fig.add_subplot(1,1,1)
x = np.cos(deg * np.pi / 180)
y = np.sin(deg * np.pi / 180)
x1 = np.cos(deg_u * np.pi / 180)
y1 = np.sin(deg_u * np.pi / 180)
ax.quiver(x,y,x,y,color='r',label='phase velocity',)
ax.quiver(x,y,x1,y1,label='group')
# ax.plot(deg * np.pi / 180,c,)
# ax.plot(deg * np.pi / 180,c*0+c[0],label='phase Ref')
# ax.plot(deg_u * np.pi / 180,u,label="group velocity")
#ax.plot(deg * np.pi / 180,u*0+u[0],label='group Ref')

plt.title("HTI model")
plt.legend()
plt.gca().set_aspect('equal')
plt.savefig("vector.jpg")