import numpy as np 
import sys

def get_disp(thk,vp,vs,rho,T,mode,wavetype,sphere=False):
    import cps330
    thk0 = np.require(thk,dtype='f4',requirements='C')
    vp0 = np.require(vp,dtype='f4',requirements='C')
    vs0 = np.require(vs,dtype='f4',requirements='C')
    rho0 = np.require(rho,dtype='f4',requirements='C')
    t = np.require(T,dtype='f8',requirements='C')

    cg,_ = cps330.forward(thk0,vp0,vs0,rho0,t,wavetype,mode,sphere)
    return cg


def main():
    if len(sys.argv) != 2:
        print("Usage: python bench_cps.py wavetype")
        exit(1)
    wavetype = int(sys.argv[1])

    data = np.loadtxt("model.txt",skiprows=1,dtype=np.float32)
    thk = data[:,0]
    rho = data[:,1]
    vs = data[:,2]
    vp = data[:,4]

    # load frequency
    T = np.loadtxt("out/swd.txt",max_rows=1)

    # open file
    fp = open("out/swd.cps.txt","w")
    for it in range(len(T)):
        fp.write("%f " %(T[it]))
    fp.write("\n")

    if wavetype == 2:
        wtp = 'R'
    else:
        wtp = 'L'
    for i in range(0,4):
        c = get_disp(thk,vp,vs,rho,T,i,f"{wtp}c",False)
        u = get_disp(thk,vp,vs,rho,T,i,f"{wtp}g",False)

        for it in range(len(T)):
            fp.write("%d %g %g %d\n" %(it,c[it],u[it],i))
    fp.close()

main()
