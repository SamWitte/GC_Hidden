###############################################




###############################################


import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import os
import matplotlib as mpl
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import rc

from profiles import *

rc('font', **{'family': 'serif', 'serif': ['Times', 'Palatino']})
rc('text', usetex=True)

mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

fs=18.

try:
    MAIN_PATH = os.environ['GC_HIDDEN_PATH']
except KeyError:
    print 'GC_HIDDEN_PATH environment variable not defined, defaulting to cwd'
    MAIN_PATH = os.getcwd() + '/'

# File Contents of GC_stats.dat:
# E_low, E_mid, E_high, Stat err


def chi_covariance(log_s, spec='BB_direct_mx_50GeV.dat', maj=True, l_min=0.,
                   l_max=20., b_min=2., b_max=20., gamma=1.2):
    gce_file = np.loadtxt(MAIN_PATH + 'GC_stat_err.dat')
    gce_dat = np.loadtxt(MAIN_PATH + 'GC_dat.dat')

    #sky = (4. * np.pi / 180. * (l_max - l_min) * (np.sin(b_max * np.pi / 180.) - np.sin(b_min * np.pi / 180.)))
    sky = 1.
    
    Jfac = Gen_NFW(gamma=gamma).J(l_min=l_min, l_max=l_max, b_min=b_min, b_max=b_max)
    
    try:
        dn = np.loadtxt(MAIN_PATH + '/Spectrum/' + spec)
    except FileNotFoundError:
        print 'Spectrum File Not Found.'
        raise

    findmx = spec.find('mx_')
    findGeV = spec.find('GeV')
    mx = float(spec[findmx+3:findGeV])

    n_obs = np.copy(gce_dat[:, 1])
    model = np.zeros_like(n_obs)

    for i in range(len(n_obs)):
        model[i] = dNdE_no_norm(mx, dn[i], Jfac, maj=maj)
        n_obs[i] /= gce_dat[i, 0]**2.
        n_obs[i] *= sky

    sigma = build_covariance_matrix(gce_file[:,1], gce_file[:,3],
                                    sky, gce_file[:,0], gce_file[:, 2])

    chi = 0.
    for i in range(len(model)):
        for j in range(len(model)):
            chi += (10**log_s*model[i] - n_obs[i]) * (10**log_s*model[j] - n_obs[j]) / sigma[i, j]
    
    return chi


def build_covariance_matrix(eng, stat, sky, bnds1, bnds2, dim=24):

    sigma = np.zeros((dim, dim))

    for i in range(dim):
        for j in range(dim):
            
            # Statistical Err
            if i == j:
                sigma[i,j] += (stat[i] / eng[i]**2. * sky)**2.

            # Residual Err
            if eng[i] <= 10. and eng[j] <= 10.:
                res_1 = res_flux(bnds1[i], bnds2[i]) * sky
                res_2 = res_flux(bnds1[i], bnds2[i]) * sky
                if i == j:
                    sigma[i, j] += 2. * res_1 * res_2
                else:
                    sigma[i, j] += res_1 * res_2
            
            # Model trunc
            sigma[i, j] += model_pc_approx(bnds1[i], bnds2[i], bnds1[j], bnds2[j]) * sky**2.

    return sigma


def res_flux(eng1, eng2):
    # Takes energy in GeV, returns differential flux in units of GeV^-1 cm^-2 s^-1 sr^-1
    integ = quad(lambda x: 6. * 10. ** -8. * x ** -3., eng1, eng2)
    return integ[0] / (eng2 - eng1)

def ICS_flux(x, log_t=False, Eref=3.):
    ics = 4.0 * 10. ** -6. * x ** -2.3
    if log_t:
        ics *= np.log(x/Eref)
    return ics

def pi0_B_flux(x, log_t=False, Eref=3.):
    pi0 = 1.0 * 10. ** -5. * x ** -2. / np.sqrt(x ** 1.2 + x ** -0.4)
    if log_t:
        pi0 *= np.log(x/Eref)
    return pi0

def model_pc_approx(E11,E12, E21, E22, Eref=3.):
    alpha_ics = 0.025
    alpha_pi0_b = 0.031
    gamma_ics = 0.0093
    gamma_pi0_b = 0.
    
    ics_1_nt = quad(ICS_flux, E11, E12, args=(False,Eref))
    ics_1_t = quad(ICS_flux, E11, E12, args=(True,Eref))
    ics_2_nt = quad(ICS_flux, E21, E22, args=(False,Eref))
    ics_2_t = quad(ICS_flux, E21, E22, args=(True,Eref))

    p0_1_nt = quad(pi0_B_flux, E11, E12, args=(False,Eref))
    p0_1_t = quad(pi0_B_flux, E11, E12, args=(True,Eref))
    p0_2_nt = quad(pi0_B_flux, E21, E22, args=(False,Eref))
    p0_2_t = quad(pi0_B_flux, E21, E22, args=(True,Eref))

    ics = (alpha_ics**2.*ics_1_nt[0]*ics_2_nt[0] + gamma_ics**2.*ics_1_t[0]*ics_2_t[0])
    pi0_b = (alpha_pi0_b**2.*p0_1_nt[0]*p0_2_nt[0] + gamma_pi0_b**2.*p0_1_t[0]*p0_2_t[0])

    return (ics + pi0_b) / ((E22 - E21) * (E12 - E11))

def plot_model_vs_data(log_s, spec='BB_direct_mx_50GeV_FLUX.dat', maj=True):
    gce_file = np.loadtxt(MAIN_PATH + 'GC_stat_err.dat')
    gce_dat = np.loadtxt(MAIN_PATH + 'GC_dat.dat')
    try:
        dn = np.loadtxt(MAIN_PATH + '/Spectrum/' + spec)
    except FileNotFoundError:
        print 'Spectrum File Not Found.'
        raise
    
    findmx = spec.find('mx_')
    findGeV = spec.find('GeV')
    mx = float(spec[findmx+3:findGeV])

    n_obs = gce_dat[:, 1]
    model = np.zeros_like(n_obs)
    
    for i in range(len(n_obs)):
        model[i] = dNdE_no_norm(mx, dn[i], maj=maj)

    fig = plt.figure()
    ax = plt.gca()

    dm = 10.**log_s * model
    print np.column_stack((gce_dat[:,0], dm))
    plt.plot(gce_dat[:,0], dm, color='blue')
    plt.plot(gce_dat[:,0], gce_dat[:,1], 'o', mfc='r')

    ax.set_xscale('log')
    ax.set_yscale('log')
    pl.xlim([0.3, 500.])
    pl.ylim([10 ** -8., 10 ** -5.])
    pl.ylabel(r'', fontsize=fs)
    pl.xlabel(r'$E_\gamma$  [$GeV$]', fontsize=fs)

    figname = MAIN_PATH + '/Test_Plots/Look_at_example_spectrum.pdf'
    fig.set_tight_layout(True)
    pl.savefig(figname)

    return




