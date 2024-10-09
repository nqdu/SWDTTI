import libswd 
import numpy as np 

model = np.loadtxt("dafaf")

def swdiso(thk,vp,vs,rho,t,wavetype:str,only_phase:bool,mode:int =0):
    thk_in = np.require(thk,"f4",requirements='C_CONTIGUOUS')
    vp_in = np.require(vp,"f4",requirements='C_CONTIGUOUS')
    vs_in = np.require(vs,"f4",requirements='C_CONTIGUOUS')
    rho_in = np.require(rho,"f4",requirements='C_CONTIGUOUS')
    t_in = np.require(t,"f8",requirements='C_CONTIGUOUS')

    return libswd.swdiso(thk_in,vp_in,vs_in,rho_in,t_in,wavetype,only_phase,mode)
