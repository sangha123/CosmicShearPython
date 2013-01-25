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
    
    kcamb,Pkcamb=numpy.loadtxt('pow_test_matterpower.dat',unpack=True)
    
    pylab.loglog(kcamb,Pkcamb,'o',ms=0.5)
    h=0.72
    pylab.loglog(k/h,P[:,0]*h**3)
    pylab.legend(("camb","emu"),loc=0)
    pylab.savefig('emucamb.png')
    sys.exit()
        
    C=powfunc.Cosmology(h=0.75)
    nl=100
    l=numpy.arange(nl)+1.
    
    alpha=2
    beta=2
    z0=0.374
    na=40
    nz= 1./(2*z0**3) *z**2 *numpy.exp(-z/z0)

    na=40 # no. of galaxies per stradian
    B=na/numpy.sum(nz)
    nz=nz*B

    zbin=numpy.array([0.2,0.65,0.85,1.0,1.4])
    num_bin=zbin.shape[0]
    sigmaz=0.05

    num_bin=zbin.shape[0]
    sigmaz=0.05
    zbias=0.
    nbinz=C.nbinz(z,nz,sigmaz,zbias,zbin)
    nia=numpy.zeros(num_bin-1,dtype='double')
    dz=z[1]-z[0]

    for i in xrange(num_bin-1):
        nbinz[:,i]=nbinz[:,i]/numpy.sum(nbinz[:,i]*dz)

    Cll=C.Cl_shear(z,P,k,l,nbinz)
    Cl1=Cll[0][:,0]
    Cl2=Cll[0][:,1]
    Cl3=Cll[0][:,2]
    l=Cll[1]
    W=Cll[2]
    eta=Cll[3]
    k_new=Cll[4]
    pickle.dump(Cll,open('Cll_emu.pkl',mode='w'))

