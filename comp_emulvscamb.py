import subprocess
import numpy,pylab,scipy
import scipy.integrate 
from scipy.interpolate import spline
import scipy.special
import sys
import pickle
import powfunc


if __name__=="__main__":
    
    k=numpy.loadtxt('kfile.dat')
    z=numpy.loadtxt('zfile.dat')

    ksize=k.shape[0]
    nz=z.shape[0]

    Plong=numpy.loadtxt('emupow.dat')
    
    P=numpy.zeros((ksize,nz),dtype='double')

    for i in xrange(nz):
        P[:,i]=Plong[i*ksize:i*ksize+ksize]
        

    X=numpy.loadtxt('MP.linear.base')
    z_camb=numpy.loadtxt('MP.redshifts.base')
    nk=X.shape[0]
    
    

    k_camb=X[:,0]
    P_camb=X[:,1:1001]

    Cl_emul=pickle.load(open('Cll_emu.pkl'))
    Cl_camb=pickle.load(open('Cll_camb.pkl'))

    Cl1_emul=Cl_emul[0][:,0]
    Cl2_emul=Cl_emul[0][:,1]
    
    l_emul=Cl_emul[1]
    
    Cl1_camb=Cl_camb[0][:,0]
    Cl2_camb=Cl_camb[0][:,1]
    
    l_camb=Cl_camb[1]

    pylab.loglog(l_emul,Cl2_emul)
    pylab.loglog(l_camb,Cl2_camb)

    pylab.xlabel('$\ell$',size=20)
    pylab.ylabel('$C_\ell$',size=20)
    
    pylab.legend(('Emu','CAMB'),loc=0)


    pylab.savefig('Cl2_emucamb.png')

    pylab.close()
    id_camb=90
    id_emu=178

    pylab.loglog(k,P[:,90])
    pylab.loglog(k_camb*0.75,P_camb[:,178]/0.75**3)
    

    pylab.xlabel('$k (Mpc^{-1})$',size=20)
    pylab.ylabel('$P(k)$',size=20)

    pylab.legend(('Emu','CAMB'),loc=0)

    pylab.savefig('Pk_emucamb_z0.36.png')
    pylab.close()
