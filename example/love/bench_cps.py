import numpy as np 


def get_disp(thk,vp,vs,rho,T,mode,wavetype,sphere=False):
    import libsurf
    thk0 = np.require(thk,dtype='f4',requirements='C')
    vp0 = np.require(vp,dtype='f4',requirements='C')
    vs0 = np.require(vs,dtype='f4',requirements='C')
    rho0 = np.require(rho,dtype='f4',requirements='C')
    t = np.require(T,dtype='f8',requirements='C')

    cg,ier = libsurf.forward(thk0,vp0,vs0,rho0,t,wavetype,mode,sphere)
    return cg


def main():
    # z = np.linspace(0,50,60)
    # nz = len(z)
    # vs = 3.5 + 0.02 * z 

    # # vs[10:20] = 1.4
    # # vs[40:50] = 2.5
    
    # vp = 1.732 * vs 
    # rho = 0.3601 * vp + 0.541 
    # thk = np.zeros_like(z)
    # thk[0:nz-1] = np.diff(z)
    # #T = np.linspace(4,150,120)
    # read model
    data = np.loadtxt("model.txt",skiprows=1,dtype=np.float32)
    thk = data[:,0]
    rho = data[:,1]
    vs = data[:,2]
    vp = data[:,3]

    f = np.linspace(-2.2,0,120)
    T = 1. / 10**f 

    # open file
    fp = open("out/swd.cps.txt","w")

    for i in range(4):
        c = get_disp(thk,vp,vs,rho,T,i,"Lc",False)
        u = get_disp(thk,vp,vs,rho,T,i,"Lg",False)

        for it in range(len(T)):
            fp.write("%g %g %g\n" %(T[it],c[it],u[it]))
    fp.close()

main()
