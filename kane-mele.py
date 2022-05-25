# Kane-Mele model as in cond-mat/0506581
import sys
import os
if'../lib/' not in sys.path:
    sys.path.append('../lib/')
import numpy as np
import math
#from numpy import linalg as la
from scipy import linalg as la
#from npext import *
import npext
import cmath

one2 = np.identity(2, dtype=np.complex)
one4 = np.identity(4, dtype=np.complex)
sx = np.array([[0,1],[1,0]], dtype=np.complex)
sy = np.array([[0,-1j],[1j,0]], dtype=np.complex)
sz = np.array([[1,0],[0,-1]], dtype=np.complex)
gg1 = np.kron(sx,one2)
gg2 = np.kron(sz, one2)
gg3 = np.kron(sy,sx)
gg4 = np.kron(sy,sy)
gg5 = np.kron(sy,sz)
gg12 = -1j * np.dot(gg1,gg2)
gg15 = -1j * np.dot(gg1,gg5)
gg23 = -1j * np.dot(gg2,gg3)
gg24 = -1j * np.dot(gg2,gg4)

def n2kgrid(nx,ny):
    """ given nx,ny,nz, return the k grid """
    twopi = np.pi*2
    dkx,dky = twopi / np.array([nx,ny])
    ntot = nx*ny
    kgrid = np.mgrid[0:twopi:dkx, 0:twopi:dky].reshape(2,ntot).T
    return kgrid
def label2header(labels):
    return '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
    
def decomp(hmat, gamma):
    return np.trace(np.dot(hmat,gamma))/4
def full_decomp(hmat):
    """ fully decompose a hermitian 4x4 matrix in terms of the gamma matrices """
    print("id\t%s"%(decomp(hmat,one4)))
    for i in range(1,6):
        print("%d\t%s"%(i,decomp(hmat,eval('gg%d'%i))))
    for i in range(1,6):
        for j in range(i+1,6):
            print("%d,%d\t%s"%(i,j,decomp(hmat,gg(i,j))))

def gggg(a,b,c,d):
    gga,ggb,ggc,ggd = eval('gg%d,gg%d,gg%d,gg%d'%(a,b,c,d))
    return -np.dot(np.dot(gga,ggb),np.dot(ggc,ggd))
def ggg(a,b,c):
    gga,ggb,ggc = eval('gg%d,gg%d,gg%d'%(a,b,c))
    return -1j * np.dot(gga,np.dot(ggb,ggc))
def gg(a,b):
    gga,ggb = eval('gg%d,gg%d'%(a,b))
    return -1j * np.dot(gga,ggb)


### defaults
V=0.1
R=0.05
SO=0.06
KY_SHIFT=0.00                 # so as to avoid exact 0 and pi's, etc.
class KM:
    def __init__(self,v=V,r=R,so=SO,nx=100,ny=100):
        self.v,self.r,self.so = v,r,so
        self.nx,self.ny = nx,ny
        self.dim = 4
    def hk(self,k1,k2):
        x,y = (k1-k2)/2, (k1+k2)/2
        cx,cy = np.cos((x,y))
        sx,sy = np.sin((x,y))

        d1 = 1 + 2*cx*cy
        d2 = self.v
        d3 = self.r * (1-cx*cy)
        d4 = -np.sqrt(3)*self.r * sx*sy
        d12 = -2*cx*sy
        d15 = 4*self.so * sx * (cx - cy)
        d23 = -self.r * cx*sy
        d24 = np.sqrt(3)*self.r * sx*cy

        return d1*gg1 + d2*gg2 + d3*gg3 + d4*gg4 + d12*gg12 + d15*gg15 + d23*gg23 + d24*gg24
    def export_erg(self,fn='km-erg'):
        par = open(fn+'.par','w')
        print("""
        v = %g
        r = %g
        so = %g
        nx = %d
        ny = %d
        """%(self.v,self.r,self.so,self.nx,self.ny), file=par)
        par.close()

        dat = open(fn+'.dat','w')

        labels = 'kx ky E1 E2 E3 E4'
        header = label2header(labels)

        for kx,ky in n2kgrid(self.nx,self.ny):
            print('\rkx = %g , ky = %g                    '%(kx,ky), end='')
            hh = self.hk(kx,ky)
            eig = la.eigvalsh(hh)
            print('%g\t%g\t%s'%(kx,ky, '\t'.join(['%g'%ee for ee in eig])), file=dat)
    def memo_eig_kx(self, kx):
        """ memoize eigen expansion of the Hamiltonian at given kx """

        dim=self.dim
        u_kx = np.zeros(dim*dim*self.ny,dtype=np.complex).reshape((self.ny,dim,dim))
        eig_kx = np.zeros(dim * self.ny).reshape((self.ny,dim))
        
        dky = 2*np.pi/self.ny
        for y in np.arange(self.ny):
            ky = dky * y + KY_SHIFT
            eig_kx[y], u_kx[y] = la.eigh(self.hk(kx,ky))
        self.eig_kx = eig_kx
        self.u_kx = u_kx
    def memo_band_proj(self,nbands):
        """ assuming memo_eig_kx has been invoked already, projector onto lowest 'nbands' bands"""

        proj = np.zeros(self.u_kx.shape, dtype=np.complex)

        for y,u in enumerate(self.u_kx):
            occu = u[:,:nbands]
            proj[y] = np.dot(occu, npext.dagger(occu))
        self.proj = proj
        self.sqrt_rho = proj
    def wannier(self,nbands):
        """ assuming memo_eig_kx has been populated, compute the Wannier phase of the lower nbands """

        self.memo_band_proj(nbands)
        dim = self.proj.shape[1]
        res = np.eye(dim,dtype=np.complex)

        for p in self.proj:
            res = np.dot(res,p)
        eig = la.eigvals(res)
        eig = eig[np.argsort(np.abs(eig))]
        return eig[dim - nbands:]
    def export_wannier(self,nbands,fn='km-wannier.dat'):
        dkx = np.pi*2/self.nx

        p_old = None
        with open(fn, 'w') as dat:
            labels = 'nx kx prod_amp tot_phase tot_winding ' + ' '.join(['amp_%d phase_%d'%(i+1,i+1) for i in range(nbands)])
            headline = label2header(labels)
            dat.write('#' + headline + '\n')

            wind = 0
            for x in np.arange(self.nx):
                print('\r %d / %d       '%(x,self.nx), end='')
                kx = x*dkx
                self.memo_eig_kx(kx)

                ww = self.wannier(nbands)
                aa = abs(ww)
                pp = np.angle(ww)

                prod_amp = np.prod(aa)
                tot_phase = np.sum(pp)
                
                p_old = p_old or tot_phase
                dp = tot_phase - p_old

                if abs(dp/np.pi) > 0.95:
                    wind = wind + (dp/abs(dp))
                print('%d\t%g\t%g\t%g\t%g\t%s'%(x,kx,prod_amp,tot_phase,wind,'\t'.join(['%g\t%g'%(a,p) for a,p in zip(aa,pp)])), file=dat)
                p_old = tot_phase
        None
    def memo_sqrt_rho_spin(self,rot_ky,sub=[0,2]):
        """ assuming memo_eig_kx has been invoked already,
        compute corresponding _square root_ of the spin RDM

        rot_ky: use rot_ky instead of 0 as the starting point of ky
        """
        ky_idx_shift = int(rot_ky / (2.0 * np.pi) * self.ny)
        dim = len(sub)
        sqrt_rho = np.zeros(self.ny * dim * dim, dtype=np.complex).reshape((self.ny,dim,dim))
        
        for y,(eig,u) in enumerate(zip(np.roll(self.eig_kx, ky_idx_shift, axis=0), np.roll(self.u_kx, ky_idx_shift, axis=0))):
            occu = u[:,:2]                # 2: lowest two bands are occupied
            # take the spin-up-RDM, which is the first and third
            # col/row (as indicated by the "sub" argument
            srdm = np.dot(occu,npext.dagger(occu))[sub][:,sub]
            ee,uu = la.eigh(srdm)
            # manually set negative eigenvalues to zero
            ee[ee < 0] = 0
            sqee = np.sqrt(ee)
            sqrt_rho[y] = np.dot(sqee * uu, npext.dagger(uu))
        self.sqrt_rho = sqrt_rho
        
    def uhlmann_spin(self,rot_ky = 0, sub=[0,2], prepend_rho=True):
        """ assuming memo_eig_kx has been populated, compute the spin-up uhlmann phase AT T=0. The block indices corresponding to "spin-up" should be specified by sub """

        self.memo_sqrt_rho_spin(rot_ky, sub)
        #self.memo_band_proj()

        dim = self.sqrt_rho.shape[1]
        res = np.eye(dim,dtype=np.complex)

        sq1 = self.sqrt_rho               # a0 a1 a2 ... a(N-2) a(N-1)
        sq2 = np.roll(sq1,-1,axis=0)      # a1 a2 a3 ... a(N-1) a0
        for s1,s2 in zip(sq1,sq2):
            ss = np.dot(s1,s2)
            #print('s1 = %s\ns2 = %s\n'%(s1,s2))
            u,s,vd = la.svd(ss)
            #print('u.shape = ', u.shape, 'vd.shape = ', vd.shape, 's.shape = ', s.shape)
            uvd = np.dot(u,vd)
            #res = np.dot(uvd, res)
            res = np.dot(res,uvd)

        rho0 = np.dot(sq1[0], sq1[0])
        if prepend_rho:
            # this serves to suppress low probability states in rho
            res = np.dot(rho0, res)
        
        eigs = la.eigvals(res)
        return ( eigs[np.argsort(np.abs(eigs))], np.trace(rho0) )    # sort by magnitude
    def export_uhlmann_spin(self,prepend_rho=True, rot_ky = 0, sub=[0,2], fn='km-uhlmann-spin'):
        dkx = np.pi*2/self.nx

        with open(fn+'.par', 'w') as par:
            txt = """
            v = %g
            r = %g
            so = %g
            nx = %d
            ny = %d
            rot_ky = %g * pi
            prepend_rho = %d
            """%(self.v,self.r,self.so,self.nx,self.ny,rot_ky/np.pi,prepend_rho)
            print(txt,file=par)

        labels = 'nx kx 1.abs 1.angle 2.abs 2.angle tr.abs tr.angle, rho_kappa.tr'
        header = label2header(labels)
        with open(fn+'.dat', 'w') as dat:
            print(header, file=dat)
            for x in np.arange(self.nx):
                print('\r %d / %d       '%(x,self.nx), end='')
                kx = x*dkx
                self.memo_eig_kx(kx)

                uu, rho_tr = self.uhlmann_spin(rot_ky, sub, prepend_rho)
                aa = abs(uu)
                pp = np.angle(uu)
                tr = np.sum(uu)
                print('%g\t%g\t%s\t%g\t%g\t%g'%(x,kx,'\t'.join([ '%g\t%g'%(a,p) for a,p in zip(aa,pp) ]), abs(tr),np.angle(tr),rho_tr), file=dat)

        None
