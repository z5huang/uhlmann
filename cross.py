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

def hk(k):
    return (1+np.cos(k))*sx + np.sin(k)*sy
def eigk(k):
    hh = hk(k)
    eig,u = la.eigh(hh)
    return eig,u
def populate_eig(nk=100):
    all_u = [ eigk(k)[1] for k in np.linspace(0,2*np.pi,100) ]
    all_conn = [ np.dot(npext.dagger(all_u[i]), all_u[i+1]) for i in np.arange(nk - 1)  ]
    #conn = np.eye(2, dtype=np.complex)
    
    return (all_u,all_conn)

        
        
