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

# Defaults
M=0.5
PHI=0.3 * np.pi
T1=0.3
T2=0.3
T3=0.3

def fermi(e,mu,t):
    return 1.0/ ( np.exp((e-mu)/t) + 1 )
def n2kgrid(nx,ny):
    """ given nx,ny,nz, return the k grid """
    twopi = np.pi*2
    dkx,dky = twopi / np.array([nx,ny])
    ntot = nx*ny
    kgrid = np.mgrid[0:twopi:dkx, 0:twopi:dky].reshape(2,ntot).T
    return kgrid
def proj(k,dim):
    """ |k><k| in real space """
    kstate = np.exp(1j * k * np.arange(dim))
    res = npext.v2proj(kstate)
    return res/dim
class HD:
    def __init__(self,nx=50,ny=50,m=M,phi=PHI,t1=T1,t2=T2,t3=T3):
        self.nx,self.ny = nx,ny
        self.dim = 2
        self.m,self.phi = m,phi
        self.t1,self.t2,self.t3 = t1,t2,t3
    def kc(self):
        """ return the kx point where the two edge modes will cross """
        k1 = -np.arcsin( self.m/(2*(self.t1+self.t2+self.t3)*np.sin(self.phi)) )
        k2 = np.pi - k1
        return (k1,k2)
        
    def hk(self,kx,ky):
        """ haldane model """
        bx = -1 - np.cos(kx) - np.cos(ky)
        by = -np.sin(kx) - np.sin(ky)
        bz = self.m + 2*np.sin(self.phi)*(
            self.t1 * np.sin(kx)
            - self.t2 * np.sin(ky)
            + self.t3 * np.sin(ky - kx))
        omega = -2*np.cos(self.phi)*(
            self.t1 * np.cos(kx)
            + self.t2 * np.cos(ky)
            + self.t3 * np.cos(ky - kx))
        return omega*one2 + bx*sx + by*sy + bz*sz
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

        dim=self.dim
        u_kx = np.zeros(dim*dim*self.ny,dtype=np.complex).reshape((self.ny,dim,dim))
        eig_kx = np.zeros(dim * self.ny).reshape((self.ny,dim))
        
        dky = 2*np.pi/self.ny
        for y in np.arange(self.ny):
            ky = dky * y
            eig_kx[y], u_kx[y] = la.eigh(self.hk(kx,ky))
        self.eig_kx = eig_kx
        self.u_kx = u_kx
    def memo_pky(self,ymax=None):
        """ memoize ky projectors. Only the top-left ymax x ymax block is recorded """
        ny = self.ny
        ymax = ymax or ny
        pky = np.zeros(ny*ymax*ymax, dtype=np.complex).reshape((ny,ymax,ymax))
        dky = 2*np.pi/ny
        for y in range(ny):
            ky = dky * y
            ky_state = np.exp(1j * np.arange(ymax) * ky)   # not normalized
            pp = np.outer(ky_state, npext.dagger(ky_state)) / ny
            pky[y] = pp
        self.proj_ky = pky

    def memo_band_proj(self):
        """ assuming memo_eig_kx has been invoked already, compute corresponding _square root_ of density matrices with fermi energy 'mu' and temperature 't' """

        proj = np.zeros(self.u_kx.shape, dtype=np.complex)

        for y,u in enumerate(self.u_kx):
            proj[y] = npext.v2proj(u[:,0])
        self.proj = proj
        self.sqrt_rho = proj
    def memo_corr_thermal(self, mu, t, rot_ky, normalize=False):
        """ assuming memo_eig_kx has been invoked already,
        compute corresponding correlation matrix
        with fermi energy 'mu' and temperature 't'

        rot_ky: use rot_ky instead of 0 as the starting point of ky
        """

        ky_idx_shift = int(rot_ky / (2.0 * np.pi) * self.ny)
        #print('ky_idx_shift = ', ky_idx_shift)
        
        dim = self.eig_kx.shape[1]
        corr = np.zeros(self.ny * dim * dim, dtype=np.complex).reshape((self.ny,dim,dim))
        for y,(eig,u) in enumerate(zip(np.roll(self.eig_kx, ky_idx_shift, axis=0), np.roll(self.u_kx, ky_idx_shift, axis=0))):
            f = fermi(eig,mu,t)
            corr[y] = np.dot(f * u, npext.dagger(u))
        self.corr = corr

    def wannier(self):
        """ assuming memo_eig_kx has been populated, compute the Wannier phase of the lower band """

        self.memo_band_proj()
        dim = self.proj.shape[1]
        res = np.eye(dim,dtype=np.complex)

        for p in self.proj:
            res = np.dot(res,p)
        return np.trace(res)
    def wannier_thermal(self, mu,t,rot_ky,normalize=False):
        """ assuming memo_eig_kx is populated, compute the Wannier stuff """
        self.memo_corr_thermal(mu,t,rot_ky,normalize)
        dim = self.corr.shape[1]
        res = np.eye(dim, dtype=np.complex)
        for p in self.corr:
            res = np.dot(res,p)
        eig = la.eigvals(res)
        return eig
        
    def export_wannier(self,fn='wannier.dat'):
        dkx = np.pi*2/self.nx

        p_old = None
        with open(fn, 'w') as dat:
            labels = 'nx kx amp phase winding'
            headline = '\t'.join([ '%d:%s'%(i+1,t) for i,t in enumerate(labels.split()) ])
            dat.write('#' + headline + '\n')

            wind = 0
            for x in np.arange(self.nx):
                print('\r %d / %d       '%(x,self.nx), end='')
                kx = x*dkx
                self.memo_eig_kx(kx)

                ww = self.wannier()
                aa = abs(ww)
                pp = np.angle(ww)
                
                p_old = p_old or pp
                dp = pp - p_old

                if abs(dp/np.pi) > 0.95:
                    wind = wind + (dp/abs(dp))
                print('%g\t%g\t%g\t%g\t%g'%(x,kx,aa,pp,wind), file=dat)
                p_old = pp
        None
    def export_wannier_thermal(self,mu,t,rot_ky=0,normalize=False, fn='wannier_thermal.dat'):
        dkx = np.pi*2/self.nx

        p_old = None
        with open(fn, 'w') as dat:
            labels = 'nx kx amp phase amp1 phase1 amp2 phase2'
            headline = '\t'.join([ '%d:%s'%(i+1,t) for i,t in enumerate(labels.split()) ])
            dat.write('#' + headline + '\n')

            for x in np.arange(self.nx):
                print('\r %d / %d       '%(x,self.nx), end='')
                kx = x*dkx
                self.memo_eig_kx(kx)

                ww = self.wannier_thermal(mu,t,rot_ky,normalize)
                aa = abs(ww)
                pp = np.angle(ww)
                tr = np.sum(ww)
                
                print('%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g'%(
                    x,kx,
                    abs(tr), np.angle(tr),
                    aa[0], pp[0], aa[1], pp[1]
                    ), file=dat)
        None
    
    def sqrt_rho_k_t(self,kx,ky,mu,t, normalize=False):
        """ square root of rho(k) at temperature t with chemical potential mu"""
        eig,u = la.eigh(self.hk(kx,ky))
        f = fermi(eig,mu,t)
        if normalize:
            f = f / np.sum(f)                 # normalize
        sqf = np.sqrt(f)
        res = np.dot(u, npext.dagger(sqf*u))
        return res
    
    def memo_sqrt_rho_thermal(self, mu, t, rot_ky, normalize=False):
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
            if normalize:
                f = f/np.sum(f)               # normalize
            sqf = np.sqrt(f)
            sqrt_rho[y] = np.dot(sqf * u, npext.dagger(u))
        self.sqrt_rho = sqrt_rho
    def corr_thermal_yhalf(self,mu,t,ny_half=None):
        """ assume memo_eig_kx and memo_pky are populated,
        compute the correlation matrix in the y representation,
        PROJECTED ONTO THE HALF SYSTEM, i.e. from y=0 to ny/2-1 """
        ny = self.ny
        ny_half = ny_half or int(ny/2)               # default is half or slightly smaller
        dim = self.eig_kx.shape[1] * ny_half
        res = np.zeros(dim*dim, dtype=np.complex).reshape(dim,dim)

        fsum = 0

        for y,(eig,u,pp) in enumerate(zip(self.eig_kx, self.u_kx, self.proj_ky)):
            f = fermi(eig,mu,t)
            corr = np.dot(f * u, npext.dagger(u))
            res += np.kron(pp[:ny_half, :ny_half], corr)
        return res
    def export_ent(self,mu,t,ny_half=None,fn='ent.dat'):
        ny = self.ny
        ny_half = ny_half or int(ny/2)               # default is half or slightly smalle
        self.memo_pky(ymax=ny_half)
                
        dkx = np.pi*2/self.nx
        with open(fn, 'w') as dat:
            labels = 'nx kx n f(n)'
            headline = '\t'.join([ '%d:%s'%(i+1,t) for i,t in enumerate(labels.split()) ])
            dat.write('#' + headline + '\n')
            for x in np.arange(self.nx):
                print('\r %d / %d       '%(x,self.nx), end='')
                kx = x*dkx
                self.memo_eig_kx(kx)

                corr = self.corr_thermal_yhalf(mu,t,ny_half)
                eig = la.eigvalsh(corr)
                eig[eig > 1] = 1
                eig[eig < 0] = 0
                for nf,f in enumerate(eig):
                    print('%d\t%g\t%d\t%g'%(x,kx,nf,f), file=dat)
        None

    def sqrt_rho_thermal_yhalf(self,mu,t,ny_half=None, normalize=False):
        """ assume memo_eig_kx and memo_pky are populated,
        compute the square root of the density matrix in the y representation,
        PROJECTED ONTO THE HALF SYSTEM, i.e. from y=0 to ny/2-1 """
        ny = self.ny
        ny_half = ny_half or int(ny/2)               # half or slightly smaller
        dim = self.eig_kx.shape[1] * ny_half
        res = np.zeros(dim*dim, dtype=np.complex).reshape(dim,dim)

        fsum = 0

        for y,(eig,u,pp) in enumerate(zip(self.eig_kx, self.u_kx, self.proj_ky)):
            f = fermi(eig,mu,t)
            fsum += np.sum(f)
            sqf = np.sqrt(f)
            sqrho = np.dot(sqf * u, npext.dagger(u))
            res += np.kron(pp[:ny_half, :ny_half], sqrho)
        if normalize:
            res /= np.sqrt(fsum)
        return res
    def memo_sqrt_rho_thermal_yhalf(self,mu,t,rot_kx=0,ny_half=None,normalize=False):
        """ assume memo_pky is populated, memoize the kx series of reduced sqrt_rho in half-f system """
        nx,ny = self.nx,self.ny
        ny_half = ny_half or int(ny/2)

        dim = self.dim * ny_half
        res = np.zeros(nx * dim**2, dtype=np.complex).reshape(nx,dim,dim)

        dkx = 2*np.pi/nx
        for x in range(nx):
            kx = x*dkx
            print('\r memoizing %d / %d                      '%(x,nx), end='')
            self.memo_eig_kx(kx+rot_kx)
            res[x] = self.sqrt_rho_thermal_yhalf(mu,t,ny_half,normalize)
        self.sqrt_rho_yhalf = res
        
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
    def uhlmann_yhalf(self,mu,t,rot_kx=0, ny_half=None,prepend_rho=True, normalize=False):
        """ compute the uhlmann values for the set of _entanglement_ correlation matrices at different kx,
        the entanglement cut is in the y direction and takes y = 0 to ny/2-1
        """

        nx,ny = self.nx, self.ny
        ny_half = ny_half or int(ny_half/2)
        self.memo_sqrt_rho_thermal_yhalf(mu,t,rot_kx,ny_half,normalize)

        dim = self.sqrt_rho_yhalf.shape[1]
        res = np.eye(dim,dtype=np.complex)
        sq1 = self.sqrt_rho_yhalf         # a0 a1 a2 ... a(N-2) a(N-1)
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

    def export_uhlmann(self,mu,t,prepend_rho=True, rot_ky = 0, normalize=False,fn='uhlmann'):
        dkx = np.pi*2/self.nx

        with open(fn+'.par', 'w') as par:
            txt = """
            m = %g
            phi = %g * pi
            t1 = %g # second neighbor hopping along dir 1, etc.
            t2 = %g
            t3 = %g
            nx = %g
            ny = %g
            mu = %g
            T = %g
            prepend_rho = %d
            rot_ky = %g * pi
            normalize = %d
            """%(self.m,self.phi/np.pi,self.t1,self.t2,self.t3,self.nx,self.ny,mu,t,prepend_rho, rot_ky/np.pi, normalize)
            print(txt,file=par)

        labels = 'nx kx 1.abs 1.angle 2.abs 2.angle tr.abs tr.angle, rho_kappa.tr'
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
    def export_uhlmann_yhalf(self,mu,t,ny_half=None, prepend_rho=True, rot_kx = 0, normalize=False,fn='uhlmann_yhalf'):
        nx,ny = self.nx, self.ny
        ny_half = ny_half or int(ny/2)

        dkx = np.pi*2/self.nx

        labels = 'n abs ang'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        rho_tr = None
        with open(fn+'.dat', 'w') as dat:
            print(header, file=dat)

            self.memo_pky(ny_half)
            uu, rho_tr = self.uhlmann_yhalf(mu,t, rot_kx, ny_half, prepend_rho, normalize)
            uu_tr = np.sum(uu)
            for n,u in enumerate(uu):
                print('%d\t%g\t%g'%(n,abs(u),np.angle(u)), file=dat)
        with open(fn+'.par', 'w') as par:
            txt = """
            tr_amp = %g
            tr_angle = %g
            m = %g
            ny_half = %d
            phi = %g * pi
            t1 = %g # second neighbor hopping along dir 1, etc.
            t2 = %g
            t3 = %g
            nx = %g
            ny = %g
            mu = %g
            T = %g
            prepend_rho = %d
            rot_kx = %g * pi
            normalize = %d
            """%(abs(uu_tr),np.angle(uu_tr), self.m,ny_half, self.phi/np.pi,self.t1,self.t2,self.t3,self.nx,self.ny,mu,t,prepend_rho, rot_kx/np.pi, normalize)
            print(txt,file=par)

        None
    def export_uhlmann_kx(self,mu,kx = None, tlist=np.arange(0.02,1.02,0.02),prepend_rho=True, normalize=False, fn='uhlmann_kc'):
        """ At a given kx point (defaults to the edge crossing point),
        compute the eigenvalues of rhoU with varying kappa (i.e., head/tail point in the ky
        Wilson loop) and their sum (i.e., Tr(rhoU)) """


        kx = kx or self.kc()[1]

        dky = np.pi*2/self.ny

        with open(fn+'.par', 'w') as par:
            txt = """
            kx = %g * pi
            m = %g
            phi = %g * pi
            t1 = %g # second neighbor hopping along dir 1, etc.
            t2 = %g
            t3 = %g
            nx = %g
            ny = %g
            mu = %g
            prepend_rho = %d
            normalize = %d
            """%(kx/np.pi, self.m,self.phi/np.pi,self.t1,self.t2,self.t3,self.nx,self.ny,mu,prepend_rho, normalize)
            print(txt,file=par)

        labels = 'ny kappa_y T 1.abs 1.angle 2.abs 2.angle tr.abs tr.angle, rho_kappa.tr'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        with open(fn+'.dat', 'w') as dat:

            self.memo_eig_kx(kx)
            
            print(header, file=dat)
            for y in np.arange(self.ny):
                print('\r %d / %d       '%(y,self.ny), end='')
                kappa_y = y*dky

                for t in tlist:
                    print('\r kappa = %d / %d   ,   T = %g     '%(y, self.ny, t), end='')
                    uu, rho_tr = self.uhlmann(mu,t, kappa_y, prepend_rho, normalize)
                    aa = abs(uu)
                    pp = np.angle(uu)
                    tr = np.sum(uu)
                    print('%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g'%
                          (y,kappa_y,t,aa[0],pp[0],aa[1],pp[1],abs(tr),np.angle(tr),rho_tr), file=dat)
                dat.write('\n')
        None
    def export_uhlmann_k_mu(self,kx = None, ky = 0, tlist=np.arange(0.02,1.02,0.02),mulist = np.arange(-0.2,1.2,0.1),prepend_rho=True, normalize=False, fn='uhlmann_kc_mu'):
        """ At a given k point (defaults to the edge crossing kx and ky=0),
        compute the eigenvalues of rhoU and their sum (i.e., Tr(rhoU)) with varying mu"""


        kx = kx or self.kc()[1]

        dky = np.pi*2/self.ny

        with open(fn+'.par', 'w') as par:
            txt = """
            kx = %g * pi
            m = %g
            phi = %g * pi
            t1 = %g # second neighbor hopping along dir 1, etc.
            t2 = %g
            t3 = %g
            nx = %g
            ny = %g
            prepend_rho = %d
            normalize = %d
            """%(kx/np.pi, self.m,self.phi/np.pi,self.t1,self.t2,self.t3,self.nx,self.ny,prepend_rho, normalize)
            print(txt,file=par)

        labels = 'mu T 1.abs 1.angle 2.abs 2.angle tr.abs tr.angle rho_ky.tr'
        header = '#' + '\t'.join(['%d:%s'%(i+1,txt) for i,txt in enumerate(labels.split())])
        with open(fn+'.dat', 'w') as dat:

            self.memo_eig_kx(kx)
            
            print(header, file=dat)
            for mu in mulist:
                for t in tlist:
                    print('\r mu = %g   ,   T = %g     '%(mu, t), end='')
                    uu, rho_tr = self.uhlmann(mu,t, ky, prepend_rho, normalize)
                    aa = abs(uu)
                    pp = np.angle(uu)
                    tr = np.sum(uu)
                    print('%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g'%
                          (mu,t,aa[0],pp[0],aa[1],pp[1],abs(tr),np.angle(tr), rho_tr), file=dat)
                dat.write('\n')
        None

#    def norm_rho_prod(self,mu,t):
#        """ assuming memo_eig_kx has been populated,
#        compute the product of normalized rho (normalized by its determinant)
#        """
#
#        dim = self.eig_kx.shape[1]
#
#        res = np.zeros(self.ny * dim * dim, dtype=np.complex).reshape((self.ny,dim,dim))
#
#        for y, (eig,u) in enumerate( zip(self.eig_kx, self.u_kx) ):
#            f = fermi(eig,mu,t)
#            f = f/np.prod
    def memo_eigs(self):
        dim = self.dim
        nx,ny = self.nx,self.ny
        all_u = np.zeros(dim*dim*nx*ny, dtype=np.complex).reshape(nx,ny,dim,dim)
        all_eig = np.zeros(dim*nx*ny).reshape(nx,ny,dim)

        dkx,dky = 2*np.pi/nx, 2*np.pi/ny

        for x in np.arange(nx):
            kx = dkx * x
            for y in np.arange(ny):
                ky = dky * y
                all_eig[x,y], all_u[x,y] = la.eigh(self.hk(kx,ky))
        self.all_eig = all_eig
        self.all_u = all_u
    def memo_sqrt_rho(self,mu,t,normalize=False):
        """ populate sqrt(rho) at each point in the BZ """
        print('populating sqrt(rho) for the full BZ')
        dim = 2
        nx,ny = self.nx,self.ny
        self.all_sqrt_rho = np.zeros(dim*dim*nx*ny, dtype=np.complex).reshape(nx,ny,dim,dim)
        dkx,dky = 2*np.pi/nx, 2*np.pi/ny

        for x in np.arange(nx):
            kx = dkx * x
            for y in np.arange(ny):
                ky = dky * y
                print('\rkx = %d / %d,  ky = %d / %d        '%(x,nx,y,ny), end='')
                eig,u = la.eigh(self.hk(kx,ky))
                f = fermi(eig,mu,t)
                if normalize:
                    f = f/np.sum(f)
                sqf = np.sqrt(f)
                self.all_sqrt_rho[x,y] = np.dot(sqf * u, npext.dagger(u))
        print('\nDone')
    def rhoU_plaquette(self,x,y):
        """ with memo_sqrt_rho done,
        compute the matrix rho_TR * U for
        a plaquette in BZ whose top-right corner
        is (x,y). rho_TR is the DM at (x,y) and U
        is the Uhlmann connection around the plaquette

        returns (U, rho_TR, rho_TR * U)
         """

        r1,r2,r3,r4 = \
          self.all_sqrt_rho[x,y], \
          self.all_sqrt_rho[x-1,y], \
          self.all_sqrt_rho[x-1,y-1], \
          self.all_sqrt_rho[x,y-1]
        u,s,vd = la.svd(np.dot(r1,r2))
        u12 = np.dot(u,vd)

        u,s,vd = la.svd(np.dot(r2,r3))
        u23 = np.dot(u,vd)

        u,s,vd = la.svd(np.dot(r3,r4))
        u34 = np.dot(u,vd)

        u,s,vd = la.svd(np.dot(r4,r1))
        u41 = np.dot(u,vd)

        rho1 = np.dot(r1,r1)

        u12341 = np.dot(np.dot(np.dot(u12, u23), u34), u41)

        return (u12341, rho1, np.dot(rho1, u12341))

    def export_rhoU(self,mu,t,normalize=False,fn='rhoU'):
        dim = 2

        self.memo_sqrt_rho(mu,t,normalize)

        nx,ny = self.nx, self.ny
        dkx,dky = 2*np.pi/nx, 2*np.pi/ny

        # sum over BZ of each phase band
        tot_phase_U = np.zeros(dim)
        tot_phase_rhoU = np.zeros(dim)

        with open(fn+'.dat', 'w') as dat:
            
            for x in np.arange(nx):
                kx = x*dkx
                for y in np.arange(ny):
                    print('\rkx = %d / %d,  ky = %d / %d        '%(x,nx,y,ny), end='')
                    ky = y*dky
                    
                    U,rho,rhoU = self.rhoU_plaquette(x,y)
    
                    eig_U = la.eigvals(U)
                    eig_U = eig_U[np.argsort(np.abs(eig_U))]
                    eig_rhoU = la.eigvals(rhoU)
                    eig_rhoU = eig_rhoU[np.argsort(np.abs(eig_rhoU))]

                    tot_phase_U = tot_phase_U + np.angle(eig_U)
                    tot_phase_rhoU = tot_phase_rhoU + np.angle(eig_rhoU)

                    s_U = '\t'.join([ '%g\t%g'%(abs(v),np.angle(v)) for v in eig_U ])
                    s_rhoU = '\t'.join([ '%g\t%g'%(abs(v),np.angle(v)) for v in eig_rhoU ])
                    print('%d\t%d\t%g\t%g\t%s\t%s'%(x,y,kx,ky,s_U,s_rhoU), file=dat)
                dat.write('\n')
        with open(fn+'.par', 'w') as par:
            txt = """
            m = %g
            phi = %g * pi
            t1 = %g # second neighbor hopping along dir 1, etc.
            t2 = %g
            t3 = %g
            nx = %g
            ny = %g
            mu = %g
            T = %g
            """%(self.m,self.phi/np.pi,self.t1,self.t2,self.t3,self.nx,self.ny,mu,t)
            print(txt, file=par)
            print('\n'.join(['sum_of_U_arg%d = %g'%(i,v) for i,v in enumerate(tot_phase_U)]), file=par)
            print('\n'.join(['sum_of_rhoU_arg%d = %g'%(i,v) for i,v in enumerate(tot_phase_rhoU)]), file=par)
