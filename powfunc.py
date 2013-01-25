import subprocess
import numpy,pylab,scipy
import scipy.integrate 
from scipy.interpolate import spline
import scipy.special
import sys
import pickle

class Cosmology:

    def __init__(self,Omega_M=0.3,Omega_L=0.7,Omega_k=1.e-10,h=0.7,Omega_R=8.4e-5,keq=0.073,ns=1,z_CMB=1100.0,scale_sm=8.,H0=100.,w=-1):
        """
        Initializing variables
        """

        self.Omega_L = Omega_L
        self.Omega_M = Omega_M
        self.Omega_k = Omega_k
        self.Omega_R = Omega_R
        self.h = h
        self.keq= keq*self.Omega_M*self.h
        self.ns=ns
        self.z_CMB=z_CMB
        self.scale_sm=scale_sm # 8 h^-1 Mpc
        self.H0=H0
        self.coverh0=3000/h
        self.w=w
        
    def friedmann(self,a):
        """ Returns (H/H_0)^2 evaluated for a definite value of a """
        return (self.Omega_M/a**3 + self.Omega_R/a**4 + self.Omega_L)

    def G_integrand(self,a):
        # Linder 2005
        
        Omega_M_z=self.Omega_M/a**3/self.friedmann(a)
        gamma=0.55
        return ((Omega_M_z**gamma-1)/a)

    def G0(self):
        return (numpy.exp(scipy.integrate.quad(self.G_integrand,1/(1+self.z_CMB),1.0))[0])

    def Growth_F(self,z):
        a=1./(1.+z) 
        gz=numpy.zeros(len(z),dtype='double')
        for i in xrange(len(z)):
            gz[i]=numpy.exp(scipy.integrate.quad(self.G_integrand,1/(1+self.z_CMB),a[i]))[0]/self.G0()
            
        return gz*a

    def Transfer(self,k):
        x=k/self.keq
        return (numpy.log(1+0.171*x)/(0.171*x) *(1+0.284*x+(1.81*x)**2+(0.399*x)**3+(0.490*x)**4)**(-0.25)) # T(k) Transfer function

    def P0W(self,k0):
        return ((k0**2/(2*numpy.pi**2))*self.Transfer(k0)**2*k0**self.ns*(3*(numpy.sin(k0*self.scale_sm)-k0*self.scale_sm*numpy.cos(k0*self.scale_sm))/(k0*self.scale_sm)**3))

    def P0_int(self,kmin,kmax):
        return(scipy.integrate.quad(self.P0W,kmin,kmax)[0])

    def Powspec_0(self,A,k):
        return (A*k**self.ns*self.Transfer(k)**2)

    def ang_dis_int(self,a):
        return self.coverh0/(a*a*numpy.sqrt(self.friedmann(a)))

    def ang_dis(self,z1,z2):
        #print z1,z2
        #assert(z1<z2)
        if ((z1>z2) or (z1==z2)): 
            return 0
        else:
            a1=1/(1+z2)
            a2=1/(1+z1)
            return (scipy.integrate.quad(self.ang_dis_int,a1,a2)[0])
    
    def RedshiftfromDist(self,D):
        #In the flat universe the qwantity used is : D_A*(1+z)=comoving distance
        z=D/self.coverh0
        eps=0.01
        del_D=1236.
        while (numpy.absolute(del_D)>eps):
            
            a=1./(1.+z)
            
            D_new=scipy.integrate.quad(self.ang_dis_int,a,1)[0]
            
            del_D=D-D_new
            delz=0.01*del_D/D
            z=z+delz

        return(z)
    
    def W_lens(self,z0,ni,z):

        size_z=len(z)
        eta_z0z=numpy.zeros(size_z,dtype='double')
        eta_z=numpy.zeros(size_z,dtype='double')
        eta_z0=self.ang_dis(0,z0)
        dz=z[1]-z[0]
        
        for i in xrange(len(z)):
            eta_z0z[i]=self.ang_dis(z0,z[i])
            eta_z[i]=self.ang_dis(0,z[i])
                    
        num_bin=ni.shape[1]
        
        int_result=numpy.zeros(num_bin,dtype='double')
        for i in xrange(num_bin):
                        
            int_result[i]=1.5*self.Omega_M*(1+z0)*eta_z0*numpy.sum(dz*ni[:,i]*eta_z0z/eta_z)/self.coverh0**2

        return int_result
    
    def nbinz(self,z,nz,sigmaz,zbias,zbin):
        size_zbin=zbin.shape[0]
        ni=numpy.zeros((len(z),size_zbin-1),dtype='double')
        for i in xrange(size_zbin-1):
            xi1=(zbin[i+1]-z+zbias)/(2**0.5*sigmaz)
            xi=(zbin[i]-z+zbias)/(2**0.5*sigmaz)
            ni[:,i]= 0.5*nz*(scipy.special.erf(xi1)-scipy.special.erf(xi))
        return ni
    
    def Cl_shear(self,z,P,k,ell,ni):
        
        dz=z[1]-z[0]
        a=1./(1+z)
        da=1./(1+z)**2 *dz
        deta=da*self.ang_dis_int(a)
        num_bin=ni.shape[1]
        size_z=len(z)

        W=numpy.zeros((size_z,num_bin),dtype='double')
        eta=numpy.zeros(size_z,dtype='double')
        
        size_ell=len(ell)
        k_new=numpy.zeros((size_z,size_ell),dtype='double')
        P_new=numpy.zeros((size_z,size_ell),dtype='double')
        #print k.shape, P.shape, size_z,size_ell
        for i in xrange(size_z):
            eta[i]=self.ang_dis(0,z[i])
            WW = self.W_lens(z[i],ni,z)
            k_new[i,:]=ell/eta[i]
            
            Pfunc_spl=scipy.interpolate.splrep(k,P[:,i])
            P_new[i,:]=scipy.interpolate.splev(k_new[i,:],Pfunc_spl)
            for j in xrange(num_bin):
                W[i,j]=WW[j]
                
        C_ll=numpy.zeros((size_ell,num_bin),dtype='double')
        #print W.shape,P_new.shape
        for j in xrange(num_bin):
            for i in xrange(size_ell):
                #print i,j
                C_ll[i,j]=numpy.sum(W[:,j]**2*P_new[:,i]*deta/eta**2)
            
        return(C_ll,ell,W,eta,k_new)

        
    
