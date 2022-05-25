# Full BHZ, i.e. 4x4 with off-diagonal inversion-breaking term
# based on bhz-full.py
# purpose: compare two things when changing the inversion-breaking delta: (1) the Z2 idx, (2) the Uhlmann winding of the spin-up (say) component
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
gg1 = np.kron(sz,sx)
gg2 = np.kron(one2, sy)
gg0 = np.kron(one2,sz)
ggb = np.kron(sx, np.diag(np.array([1.0,0.0])))
gyy = np.kron(sy,sy)

def fermi(e,mu,t):
    return 1.0/ ( np.exp((e-mu)/t) + 1 )
def n2kgrid(nx,ny):
    """ given nx,ny,nz, return the k grid """
    twopi = np.pi*2
    dkx,dky = twopi / np.array([nx,ny])
    ntot = nx*ny
    kgrid = np.mgrid[0:twopi:dkx, 0:twopi:dky].reshape(2,ntot).T
    return kgrid

# defaults
M = 1.1
BX = 0.5
KY_SHIFT=0.00                 # so as to avoid exact 0 and pi's, etc.
DEL = 0.1
class BHZ:
    def __init__(self, nx=50, ny=50, m=M, delta=DEL):
        self.nx, self.ny = nx,ny
        self.dim = 4
        self.m = m
        self.delta = delta
    def hk(self,kx,ky):
        """ full BHZ """
        mk = 2 - self.m - np.sum(np.cos((kx,ky)))
        res = np.sin(kx) * gg1 + np.sin(ky) * gg2 + mk * gg0 + self.delta * gyy
        return res
    def hkhalf(self,kx,ky):
        mk = 2 - self.m - np.sum(np.cos((kx,ky)))
        res = np.sin(kx) * sx + np.sin(ky) * sy + mk * sz
        return res
    def export_hk(self, fn='erg.dat'):
        """ export energy spectrum """
        kgrid = n2kgrid(self.nx,self.ny)

        last_kx = kgrid[0,0]
        with open(fn, 'w') as dat:
            for kx,ky in kgrid:
                print('\rkx, ky = %g , %g        '%(kx,ky), end='')
                hh = self.hk(kx,ky)
                eig = la.eigvalsh(hh)
                if kx != last_kx:
                    dat.write('\n')
                print('%g\t%g\t%s'%(kx,ky,'\t'.join([ '%g'%v for v in eig ])), file=dat)
                last_kx = kx
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
    def export_wannier(self,nbands,fn='wannier.dat'):
        dkx = np.pi*2/self.nx

        p_old = None
        with open(fn, 'w') as dat:
            labels = 'nx kx prod_amp tot_phase tot_winding ' + ' '.join(['amp_%d phase_%d'%(i+1,i+1) for i in range(nbands)])
            headline = '\t'.join([ '%d:%s'%(i+1,t) for i,t in enumerate(labels.split()) ])
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
    def memo_sqrt_rho_thermal(self, mu, t, rot_ky, normalize=None):
        """ assuming memo_eig_kx has been invoked already,
        compute corresponding _square root_ of density matrices
        with fermi energy 'mu' and temperature 't'

        rot_ky: use rot_ky instead of 0 as the starting point of ky
        """

        ky_idx_shift = int(rot_ky / (2.0 * np.pi) * self.ny)
        #print('ky_idx_shift = ', ky_idx_shift)
        
        dim = self.eig_kx.shape[1]
        sqrt_rho = np.zeros(self.ny * dim * dim, dtype=np.complex).reshape((self.ny,dim,dim))
        for y,(eig,u) in enumerate(zip(np.roll(self.eig_kx, ky_idx_shift, axis=0), np.roll(self.u_kx, ky_idx_shift, axis=0))):
            f = fermi(eig,mu,t)
            if normalize == True:
                f = f/np.sum(f)               # normalize
            sqf = np.sqrt(f)
            sqrt_rho[y] = np.dot(sqf * u, npext.dagger(u))
        self.sqrt_rho = sqrt_rho
    def uhlmann(self,mu,t,rot_ky = 0, prepend_rho=True, normalize=False):
        """ assuming memo_eig_kx has been populated,

        compute the uhlmann phase(s) """

        self.memo_sqrt_rho_thermal(mu,t, rot_ky, normalize)
        #self.memo_band_proj()

        dim = self.sqrt_rho.shape[1]
        res = np.eye(dim,dtype=np.complex)

        sq1 = self.sqrt_rho               # a0 a1 a2 ... a(N-2) a(N-1)
        sq2 = np.roll(sq1,-1,axis=0)      # a1 a2 a3 ... a(N-1) a0
        for s1,s2 in zip(sq1,sq2):
            ss = np.dot(s1,s2)
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

    def export_uhlmann(self,mu,t,prepend_rho=True, rot_ky = 0, normalize = False, fn='uhlmann'):
        dkx = np.pi*2/self.nx

        with open(fn+'.par', 'w') as par:
            txt = """
            m = %g
            nx = %d
            ny = %d
            mu = %g
            T = %g
            rot_ky = %g * pi
            prepend_rho = %d
            normalize = %d
            """%(self.m,self.nx,self.ny,mu,t,rot_ky/np.pi,prepend_rho, normalize)
            print(txt,file=par)

        labels = 'nx kx 1.abs 1.angle 2.abs 2.angle 3.abs 3.angle 4.abs 4.angle tr.abs tr.angle, rho_kappa.tr'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        with open(fn+'.dat', 'w') as dat:
            print(header, file=dat)
            for x in np.arange(self.nx):
                print('\r %d / %d       '%(x,self.nx), end='')
                kx = x*dkx
                self.memo_eig_kx(kx)

                uu, rho_tr = self.uhlmann(mu,t, rot_ky, prepend_rho, normalize)
                aa = abs(uu)
                pp = np.angle(uu)
                tr = np.sum(uu)
                print('%g\t%g\t%s\t%g\t%g\t%g'%(x,kx,'\t'.join([ '%g\t%g'%(a,p) for a,p in zip(aa,pp) ]), abs(tr),np.angle(tr),rho_tr), file=dat)

        None
    def export_uhlmann_k_m(self,kx = None, ky = 0, mu=0, tlist=np.arange(0.02,0.52,0.01),mlist = np.arange(-1,5,0.1),prepend_rho=True, normalize=True, fn='uhlmann_kc_m'):
        """ At a given k point (defaults to the edge crossing kx and ky=0),
        compute the eigenvalues of rhoU and their sum (i.e., Tr(rhoU)) with varying m and T"""


        kx = kx or 0

        dky = np.pi*2/self.ny

        with open(fn+'.par', 'w') as par:
            txt = """
            kx = %g * pi
            mu = %g
            nx = %g
            ny = %g
            prepend_rho = %d
            normalize = %d
            """%(kx/np.pi, mu,self.nx,self.ny,prepend_rho, normalize)
            print(txt,file=par)

        labels = 'm T 1.abs 1.angle 2.abs 2.angle tr.abs tr.angle rho_ky.tr'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        with open(fn+'.dat', 'w') as dat:

            print(header, file=dat)
            for mm in mlist:
                self.m = mm
                self.memo_eig_kx(kx)
                for t in tlist:
                    print('\r m = %g   ,   T = %g                                       '%(mm, t), end='')
                    uu, rho_tr = self.uhlmann(mu,t, ky, prepend_rho, normalize)
                    aa = abs(uu)
                    pp = np.angle(uu)
                    tr = np.sum(uu)
                    print('%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g'%
                          (mm,t,
                           aa[0],
                          pp[0],
                          aa[1],
                          pp[1],
                          abs(tr),np.angle(tr), rho_tr), file=dat)
                dat.write('\n')
        None

    def memo_sqrt_rho_spin(self,rot_ky,normalize=None):
        """ assuming memo_eig_kx has been invoked already,
        compute corresponding _square root_ of the spin RDM

        rot_ky: use rot_ky instead of 0 as the starting point of ky
        """
        ky_idx_shift = int(rot_ky / (2.0 * np.pi) * self.ny)
        dim = self.eig_kx.shape[1]/2
        sqrt_rho = np.zeros(self.ny * dim * dim, dtype=np.complex).reshape((self.ny,dim,dim))
        
        for y,(eig,u) in enumerate(zip(np.roll(self.eig_kx, ky_idx_shift, axis=0), np.roll(self.u_kx, ky_idx_shift, axis=0))):
            occu = u[:,:2]                # 2: lowest two bands are occupied
            # take the top-left 2x2 submat (spin-up RDM)
            srdm = np.dot(occu,npext.dagger(occu))[:2,:2]
            ee,uu = la.eigh(srdm)
            # manually set negative eigenvalues to zero
            ee[ee < 0] = 0
            sqee = np.sqrt(ee)
            sqrt_rho[y] = np.dot(sqee * uu, npext.dagger(uu))
        self.sqrt_rho = sqrt_rho
    def uhlmann_spin(self,rot_ky = 0, prepend_rho=True, normalize=False):
        """ assuming memo_eig_kx has been populated, compute the spin-up uhlmann phase AT T=0 """

        self.memo_sqrt_rho_spin(rot_ky, normalize)
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

    def export_uhlmann_spin(self,prepend_rho=True, rot_ky = 0, normalize = False, fn='uhlmann-spin'):
        dkx = np.pi*2/self.nx

        with open(fn+'.par', 'w') as par:
            txt = """
            m = %g
            nx = %d
            ny = %d
            rot_ky = %g * pi
            prepend_rho = %d
            normalize = %d
            """%(self.m,self.nx,self.ny,rot_ky/np.pi,prepend_rho, normalize)
            print(txt,file=par)

        labels = 'nx kx 1.abs 1.angle 2.abs 2.angle tr.abs tr.angle, rho_kappa.tr'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        with open(fn+'.dat', 'w') as dat:
            print(header, file=dat)
            for x in np.arange(self.nx):
                print('\r %d / %d       '%(x,self.nx), end='')
                kx = x*dkx
                self.memo_eig_kx(kx)

                uu, rho_tr = self.uhlmann_spin(rot_ky, prepend_rho, normalize)
                aa = abs(uu)
                pp = np.angle(uu)
                tr = np.sum(uu)
                print('%g\t%g\t%s\t%g\t%g\t%g'%(x,kx,'\t'.join([ '%g\t%g'%(a,p) for a,p in zip(aa,pp) ]), abs(tr),np.angle(tr),rho_tr), file=dat)

        None
