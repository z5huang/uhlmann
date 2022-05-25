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

def fermi(e,mu,t):
    return 1.0/ ( np.exp((e-mu)/t) + 1 )
def n2kgrid(nx,ny):
    """ given nx,ny,nz, return the k grid """
    twopi = np.pi*2
    dkx,dky = twopi / np.array([nx,ny])
    ntot = nx*ny
    kgrid = np.mgrid[0:twopi:dkx, 0:twopi:dky].reshape(2,ntot).T
    return kgrid

class Hof:
    def __init__(self,p=3,q=7,nx=50,ny=50):
        self.nx, self.ny = nx,ny
        self.p, self.q = p,q
        self.phi = 2*np.pi * p/q

    def hk(self,kx,ky):
        res = np.diag(np.ones(self.q-1, dtype=np.complex),1)
        res[0,-1] += np.exp(-1j * ky)
        res += npext.dagger(res)

        dd = (1 + np.arange(self.q)) * self.phi
        res += np.diag(2*np.cos(kx + dd))

        return res
    def export_hk(self, fn='erg.dat'):
        """ export energy spectrum """
        kgrid = n2kgrid(self.nx,self.ny)

        with open(fn, 'w') as dat:
            for kx,ky in kgrid:
                print('\rkx, ky = %g , %g        '%(kx,ky), end='')
                hh = self.hk(kx,ky)
                eig = la.eigvalsh(hh)
                print('%g\t%g\t%s'%(kx,ky,'\t'.join([ '%g'%v for v in eig ])), file=dat)
    def memo_eig_kx(self, kx):
        """ memoize eigen expansion of the Hamiltonian at given kx """

        dim=self.q
        u_kx = np.zeros(dim*dim*self.ny,dtype=np.complex).reshape((self.ny,dim,dim))
        eig_kx = np.zeros(dim * self.ny).reshape((self.ny,dim))
        
        dky = 2*np.pi/self.ny
        for y in np.arange(self.ny):
            ky = dky * y
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
        """ assuming memo_eig_kx has been populated, compute the Wannier phase of the lower band """

        self.memo_band_proj(nbands)
        dim = self.proj.shape[1]
        res = np.eye(dim,dtype=np.complex)

        for p in self.proj:
            res = np.dot(res,p)
        eig = la.eigvals(res)
        eig = eig[np.argsort(np.abs(eig))]
        return (eig[self.q - nbands:])[::-1]
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
    def uhlmann(self,mu,t,rot_ky = 0, prepend_rho=True, normalize=None):
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
        if prepend_rho:
            # this serves to suppress low probability states in rho
            rho0 = np.dot(sq1[0], sq1[0])
            res = np.dot(rho0, res)
        
        eigs = la.eigvals(res)
        return eigs[np.argsort(np.abs(eigs))[::-1]]   # sort by magnitude
        #return eigs[np.argsort(np.imag(eigs))[::-1]]   # sort by imaginary part
    def export_uhlmann(self,mu,t,prepend_rho=True, rot_ky = 0, normalize = None, fn='uhlmann'):
        dkx = np.pi*2/self.nx

        with open(fn+'.par', 'w') as par:
            txt = """
            p = %d
            q = %d
            nx = %d
            ny = %d
            mu = %g
            T = %g
            rot_ky = %g * pi
            prepend_rho = %d
            """%(self.p,self.q,self.nx,self.ny,mu,t,rot_ky/np.pi,prepend_rho)
            print(txt,file=par)

        labels = 'nx kx Tr.amp Tr.phase ' + ' '.join(['amp_%d phase_%d'%(i+1,i+1) for i in range(self.q)]) +' ' +  ' '.join(['phasesum_highest_%d'%(i+1) for i in range(self.q)])
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        with open(fn+'.dat', 'w') as dat:
            print(header, file=dat)
            for x in np.arange(self.nx):
                print('\r %d / %d       '%(x,self.nx), end='')
                kx = x*dkx
                self.memo_eig_kx(kx)

                uu = self.uhlmann(mu,t, rot_ky, prepend_rho, normalize)
                aa = abs(uu)
                pp = np.angle(uu)
                pp_tot = np.angle([ np.prod(uu[:i+1]) for i in range(self.q)])
                
                tr = np.sum(uu)
                print('%d\t%g\t%g\t%g\t%s\t%s'%
                      (x,kx,abs(tr),np.angle(tr),
                       '\t'.join([ '%g\t%g'%(a,p) for a,p in zip(aa,pp) ]),
                      '\t'.join([ '%g'%ptot for ptot in pp_tot])
                      ), file=dat)
        None
