"""
Created on Wed Jul 13 09:54:38 2016

@author: SamWitte
"""
import numpy as np
import os
import scipy.integrate as integrate
import scipy.special as special
from scipy.optimize import fminbound, minimize_scalar, minimize, brentq
from scipy.interpolate import interp1d, interpn, LinearNDInterpolator
import warnings

warnings.filterwarnings('error')

try:
    MAIN_PATH = os.environ['GC_HIDDEN_PATH']
except KeyError:
    print 'GC_HIDDEN_PATH environment variable not defined, defaulting to cwd'
    MAIN_PATH = os.getcwd() + '/'

kpctocm = 3.08568 * 10 ** 21.

def dNdE_no_norm(mx, dn, Jfac, maj=True):
    dnde = dn / (8. * np.pi * mx**2.) * Jfac
    if maj:
        return dnde
    else:
        return dnde / 2.


class Gen_NFW(object):
    def __init__(self, gamma=1.2, scale_r=20., rfix=8.5, rho_fix=0.4):

        self.gamma = gamma
        self.scale_radius = scale_r
        self.rfix = rfix
        self.rho_fix = rho_fix

        self.scale_density = rho_fix * (1. + rfix/scale_r)**(3.-gamma) / (rfix/scale_r)**(-gamma) # GeV/cm^3


    def density(self, r):
        return (self.scale_density * (r/self.scale_radius)**(-self.gamma)/
                    (1. + r/self.scale_radius)**(3.-self.gamma))
    

    def J(self, l_min=0., l_max=20., b_min=2., b_max=20.):
        """
        """
        
        dist = 8.5
        l_min = l_min * np.pi / 180.
        l_max = l_max * np.pi / 180.
        b_min = b_min * np.pi / 180.
        b_max = b_max * np.pi / 180.
        
        thetat = np.linspace(np.pi * 0.01, np.pi, 200.)
        
        jfac = np.zeros_like(thetat)
        jfac2 = np.zeros_like(thetat)
        for i, th in enumerate(thetat):
            jfac[i] = integrate.quad(self.J_integrand_theta, 0., 100., args=(th),
                                     limit=200, points=np.array([0., 8.5, 9.5, 11., 100.]),
                                     maxp1=100)[0]
            #print jfac[i], jfac2[i]
        self.jt_interp = interp1d(thetat, jfac, kind = 'cubic', fill_value = 0., bounds_error=False)
        
        jfull = integrate.dblquad(self.J_bl, l_min, l_max, lambda b: b_min, lambda b: b_max)
        
        angles = 4. * (l_max - l_min) * (np.sin(b_max) - np.sin(b_min))
        return jfull[0] / angles
        

    def J_bl(self, b, l):
        theta = np.arccos(np.cos(b) * np.cos(l))
        return 4. * (kpctocm * np.cos(b) * self.jt_interp(theta))

    def J_integrand_theta(self, x, theta):
        dist = 8.5
        r = np.sqrt(dist**2. + x**2. - 2.*x*dist*np.cos(theta))
        return self.density(r)**2.



