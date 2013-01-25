import subprocess
import numpy,pylab,scipy
import scipy.integrate 
from scipy.interpolate import spline
import scipy.special
import sys
import pickle
import powfunc

if __name__=="__main__":

    X=numpy.loadtxt('MP.linear.base')
    z=numpy.loadtxt('MP.redshifts.base')
    nk=X.shape[0]
    
    nl=500
    l=numpy.arange(nl)+1
    k=X[:,0]
    P=X[:,1:1001]

    C=powfunc.Cosmology(h=1)

    '''z_new=z
    
    z0=z[200]
    eta=C.ang_dis(0,z0)
    rep = scipy.interpolate.splrep(k,P[:,200])
    Pz=numpy.zeros(len(l),dtype='double')
    k0=numpy.zeros(len(l),dtype='double')
    for i in xrange(len(l)):
        k0[i]=l[i]/eta
        Pz[i]=scipy.interpolate.splev(k0[i],rep)
        
    kmin=0.1
    kmax=100.
    sigma_8=0.8

    A=(sigma_8**2/C.P0_int(kmin,kmax))
    
    P0=C.Powspec_0(A,k)
    
    G=C.Growth_F(z_new)
    size_z=len(z_new)
    size_k=len(k)
    P_k_z=numpy.zeros((size_z,size_k),dtype='double')
    for i in xrange(size_z):
        P_k_z[i,:]=G[i]**2*P0

    '''
        
    alpha=2
    beta=2
    z0=1
    z=z[0:500]
    P=P[:,0:500]
    z0=0.374
    na=40
    nz= 1./(2*z0**3) *z**2 *numpy.exp(-z/z0)
    #nz=z**alpha*numpy.exp(-(z/z0)**beta)
    
    dz=z[1]-z[0]
    na=40 # no. of galaxies per stradian
    B=na/numpy.sum(nz*dz)
    nz=nz*B
    zbin=numpy.array([0.2,0.65,0.85,1.0,1.4])
    
    num_bin=zbin.shape[0]
    sigmaz=0.05

    num_bin=zbin.shape[0]
    sigmaz=0.05
    zbias=0.
    nbinz=C.nbinz(z,nz,sigmaz,zbias,zbin)
    nia=numpy.zeros(num_bin-1,dtype='double')
   
    for i in xrange(num_bin-1):
        nbinz[:,i]=nbinz[:,i]/numpy.sum(nbinz[:,i]*dz)
    

    Cll=C.Cl_shear(z,P,k,l,nbinz)
    Cl1=Cll[0][:,0]
    Cl2=Cll[0][:,1]
    Cl3=Cll[0][:,2]
    l=Cll[1]
    W=Cll[2]
    eta=Cll[3]
    pickle.dump(Cll,open('Cll_camb.pkl',mode='w'))
    
    
