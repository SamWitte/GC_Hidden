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
import pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl

import matplotlib.patheffects as PathEffects
import matplotlib.gridspec as gridspec
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times','Palatino']})
rc('text', usetex=True)

mpl.rcParams['xtick.major.size']=8
mpl.rcParams['ytick.major.size']=8
mpl.rcParams['xtick.labelsize']=18
mpl.rcParams['ytick.labelsize']=18


#warnings.filterwarnings('error')

try:
    MAIN_PATH = os.environ['GC_HIDDEN_PATH']
except KeyError:
    print 'GC_HIDDEN_PATH environment variable not defined, defaulting to cwd'
    MAIN_PATH = os.getcwd() + '/'

kpctocm = 3.08568 * 10. ** 21.


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

        self.scale_density = rho_fix * (1. + rfix/scale_r)**(3.-gamma) / (rfix/scale_r)**(-gamma)  # GeV/cm^3

    def density(self, r):
        return (self.scale_density * (r/self.scale_radius)**(-self.gamma)/
                    (1. + r/self.scale_radius)**(3.-self.gamma))

    def J(self, l_min=0., l_max=20., b_min=2., b_max=20.):
        """
        """
        f_name = 'JFactor_Gamma_{:.2f}_ScaleR_{:.2f}_RFix_{:.1f}_RhoFix{:.2f}_'.format(self.gamma, self.scale_radius,
                                                                                      self.rfix, self.rho_fix)
        f_name += 'lmin_{:.1f}_lmax_{:.1f}_bmin_{:.1f}_bmax_{:.1f}.dat'.format(l_min,l_max,b_min,b_max)

        if os.path.isfile(MAIN_PATH + '/FileHolding/' + f_name):
            j = np.loadtxt(MAIN_PATH + '/FileHolding/' + f_name)
            return j

        l_min = l_min * np.pi / 180.
        l_max = l_max * np.pi / 180.
        b_min = b_min * np.pi / 180.
        b_max = b_max * np.pi / 180.
        
        thetat = np.logspace(np.log10(np.pi * 0.001), np.log10(np.pi), 300)
        #thetat = np.linspace(np.pi * 0.001, np.pi, 1000)
        
        jfac = np.zeros_like(thetat)
        for i, th in enumerate(thetat):
            jfac[i] = integrate.quad(self.J_integrand_theta, 0., 100., args=th,
                                     limit=500, points=np.array([0., 100.]),
                                     maxp1=500, limlst=500)[0]
    
        #self.jt_interp = interp1d(thetat, jfac, kind='cubic', fill_value=0., bounds_error=False)
        self.jt_interp = interp1d(np.log10(thetat), np.log10(jfac), kind='cubic', fill_value=0., bounds_error=False)
        
        #jfull = integrate.dblquad(self.J_bl, b_min, b_max, lambda l: l_min, lambda b: l_max)
        jfull = integrate.dblquad(self.J_bl, l_min, l_max, lambda b: b_min, lambda b: b_max)
        
        angles = 4. * (l_max - l_min) * (np.sin(b_max) - np.sin(b_min))
        jsv = jfull[0] / angles
        np.savetxt(MAIN_PATH + '/FileHolding/' + f_name, np.array([jsv]))
        return jsv

    def J_bl(self, b, l):
        theta = np.arccos(np.cos(b) * np.cos(l))
        return 4. * (kpctocm * np.cos(b) * 10.**self.jt_interp(np.log10(theta)))

    def J_integrand_theta(self, x, theta):
        dist = 8.5
        r = np.sqrt(dist**2. + x**2. - 2.*x*dist*np.cos(theta))
        return self.density(r)**2.

    def plot_j_theta(self, fs=18):
        thetat = np.logspace(np.log10(np.pi * 10.**-4.), np.log10(np.pi), 1000)
        jfac = np.zeros_like(thetat)
        for i, th in enumerate(thetat):
            jfac[i] = integrate.quad(self.J_integrand_theta, 0., 100., args=th,
                                     limit=500, points=np.array([0., 100.]),
                                     maxp1=500, limlst=500)[0]

        thetat *= 180. / np.pi
        jfac /= (self.rfix * self.rho_fix**2.)

        pl.figure()
        ax = pl.gca()

        ax.set_xlabel(r'$\theta$  [deg]', fontsize=fs)
        ax.set_ylabel(r'$J(\theta)$', fontsize=fs)

        pl.plot(thetat, jfac, 'b', lw=1)

        plt.xlim(xmin=thetat[0], xmax=thetat[-1])
        plt.ylim(ymin=0.3, ymax=10.**3.)
        ax.set_xscale("linear")
        ax.set_yscale("log")

        pl.savefig("/Users/SamWitte/Desktop/GC_Hidden/src/Test_Plots/J_factor.pdf")

        pl.figure()
        ax = pl.gca()

        ax.set_xlabel(r'$\theta$  [deg]', fontsize=fs)
        ax.set_ylabel(r'$J(\theta)$', fontsize=fs)

        pl.plot(thetat, jfac, 'b', lw=1)

        plt.xlim(xmin=10.**-9, xmax=10.)
        plt.ylim(ymin=1., ymax=10. ** 14.)
        ax.set_xscale("log")
        ax.set_yscale("log")

        pl.savefig("/Users/SamWitte/Desktop/GC_Hidden/src/Test_Plots/J_factor_zoom.pdf")


        return
