###############################################




###############################################


import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.optimize import minimize
import os
import matplotlib as mpl
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import rc
import glob

from profiles import *

rc('font', **{'family': 'serif', 'serif': ['Times', 'Palatino']})
rc('text', usetex=True)

mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

fs = 18.

try:
    MAIN_PATH = os.environ['GC_HIDDEN_PATH']
except KeyError:
    print 'GC_HIDDEN_PATH environment variable not defined, defaulting to cwd'
    MAIN_PATH = os.getcwd() + '/'

# File Contents of GC_stats.dat:
# E_low, E_mid, E_high, Stat err


def mass_scan(filef='BB_direct_mx_', gamma=1.2, maj=True, s_low=5.e-27,
              s_high=5.e-26, n_sigs=30,
              contour_val=np.array([2.3, 6.2, 11.8]),
              contour_name=np.array(['_1Sigma', '_2Sigma', '_3Sigma']),
              scale_r=20., rfix=8.5, rho_fix=0.4):

    all_files = glob.glob(MAIN_PATH + '/Spectrum/' + filef + '*.dat')

    bf_array = np.zeros(len(all_files))
    mx_list = np.zeros(len(all_files))
    # Make All Files...
    for i,f in enumerate(all_files):
        findmx = f.find('mx_')
        findGeV = f.find('GeV')
        mx_list[i] = np.float(f[findmx + 3:findGeV])
        f_tail = f[f.find(filef):]
        bf_array[i] = sig_contour(spec=f_tail, gamma=gamma, maj=maj, s_low=s_low, s_high=s_high, n_sigs=n_sigs,
                                  contour_val=contour_val, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix,
                                  make_file=True)
    #print mx_list
    #print bf_array
    order = np.argsort(mx_list)
    mx_list = mx_list[order]
    bf_array = bf_array[order]
    print 'Mass, BF'
    print np.column_stack((mx_list, bf_array))

    goal_look = interp1d(mx_list, bf_array, kind='cubic', bounds_error=False, fill_value=10.**5.)
    print goal_look(30.)
    goal = minimize(goal_look, np.median(mx_list))
    print 'Best Fit Mass: ', goal.x[0]
    print 'Best Fit ChiSq Value: ', goal.fun
    info_hold = np.zeros((len(all_files), 6))
    for i, f in enumerate(all_files):
        f_tail = f[f.find(filef):]
        hold = sig_contour(spec=f_tail, gamma=gamma, maj=maj, s_low=s_low, s_high=s_high, n_sigs=n_sigs,
                                   contour_val=contour_val, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix,
                                   make_file=False, goal=goal.fun)
        info_hold[i] = hold[-1]
    #print info_hold
    contours = info_hold
    tot_contours = len(contour_name)
    for i, cc in enumerate(contour_name):
        c_fname = MAIN_PATH + '/FileHolding/Contours/' + filef + cc
        c_fname += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
        c_fname += '.dat'
        consv = np.stack((mx_list, contours[:, i], contours[:, i+tot_contours]), axis=-1)
        np.savetxt(c_fname, consv)

    return


def sig_contour(spec='BB_direct_mx_50GeV.dat', gamma=1.2, maj=True,
                s_low=5.e-27, s_high=5.e-26, n_sigs=10,
                contour_val=np.array([2.3, 6.2, 11.8]),
                scale_r=20., rfix=8.5, rho_fix=0.4, make_file=True,
                goal=24.):

    file_name = MAIN_PATH + '/FileHolding/Contours/ChiSq/Tabbed_ChiSq_'
    file_name += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
    file_name += '_' + spec

    findmx = spec.find('mx_')
    findGeV = spec.find('GeV')
    mx = float(spec[findmx + 3:findGeV])

    if make_file:

        bf = chi_covariance(spec=spec, maj=maj, gamma=gamma, bf=True, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix)

        if os.path.isfile(file_name):
            return bf[-1]

        s_list = np.linspace(np.log10(s_low), np.log10(s_high), n_sigs)
        hold_tab = np.zeros((len(s_list), 2))
        for i, sig in enumerate(s_list):

            hold_tab[i] = chi_covariance(spec=spec, maj=maj, gamma=gamma, bf=False, sig=sig,
                                         scale_r=scale_r, rfix=rfix, rho_fix=rho_fix)
            print 'Sigma {:.2f}, ChiSq: {:.2f}'.format(sig, hold_tab[i, 1])

        np.savetxt(file_name, hold_tab)
        return bf[-1]
    else:
        chisq = np.loadtxt(file_name)
        chi_interp = interp1d(chisq[:, 0], chisq[:, 1], kind='cubic', fill_value=10.**5., bounds_error=False)

        bf = chi_covariance(spec=spec, maj=maj, gamma=gamma, bf=True, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix)

        sl = np.zeros(len(contour_val))
        sh = np.zeros(len(contour_val))

        bndsl = [(-30., bf[0])]
        bndsh = [(bf[0], -20.)]

        for j, cc in enumerate(contour_val):
            if bf[1] > (goal + cc):
                pass
            else:
                sl[j] = minimize(lambda x: np.abs(chi_interp(x) - cc - goal), np.array([bf[0] - 0.1]),
                                 method='SLSQP', bounds=bndsl).x
                sh[j] = minimize(lambda x: np.abs(chi_interp(x) - cc - goal), np.array([bf[0] + 0.1]),
                                 method='SLSQP', bounds=bndsh).x

        # Return format: Gamma, mx, Min chi^2, BF sigma_p, 1\sig_L, 2\sig_L, 3\sig_L, 1\sig_H, 2\sig_H,3\sig_H
        return gamma, mx, bf[1], 10.**bf[0][0], np.concatenate((10.**sl, 10.**sh))


def chi_covariance(spec='BB_direct_mx_50GeV.dat', maj=True, l_min=0.,
                   l_max=20., b_min=2., b_max=20., gamma=1.2, bf=True, sig=-26.,
                   scale_r=20., rfix=8.5, rho_fix=0.4):

    gce_file = np.loadtxt(MAIN_PATH + 'GC_stat_err.dat')
    gce_dat = np.loadtxt(MAIN_PATH + 'GC_dat.dat')
    
    Jfac = Gen_NFW(gamma=gamma, scale_r=scale_r, rfix=rfix,
                   rho_fix=rho_fix).J(l_min=l_min, l_max=l_max, b_min=b_min, b_max=b_max)

    try:
        load_s = np.loadtxt(MAIN_PATH + '/Spectrum/' + spec)
        dn = interp1d(load_s[:, 0], load_s[:, 1], kind='cubic', bounds_error=False, fill_value=0.)
    except FileNotFoundError:
        print 'Spectrum File Not Found.'
        raise FileNotFoundError

    findmx = spec.find('mx_')
    findGeV = spec.find('GeV')
    mx = float(spec[findmx+3:findGeV])

    n_obs = np.copy(gce_dat[:, 1])
    n_obs /= gce_dat[:, 0]**2.

    dim = len(n_obs)

    edm = np.logspace(np.log10(0.3), np.log10(30.), 300)
    model = np.zeros((len(edm),2))
    for i in range(len(edm)):
        if edm[i] < mx:
            dne = dn(edm[i])
            model[i,0] = edm[i]
            model[i,1] = dNdE_no_norm(mx, dne, Jfac, maj=maj)
        else:
            pass
    model = model[model[:, 1] > 0.]
    model_interp = interp1d(np.log10(model[:,0]), np.log10(model[:,1]),
                            kind='cubic', fill_value=0., bounds_error=False)

    dm_ev = 10.**model_interp(np.log10(gce_dat[:,0]))
    sigma = build_covariance_matrix(gce_file[:,1], gce_file[:,3], gce_file[:, 2])

    if bf:
        bf = minimize(chi_sq, np.array([-26.]), args=(dm_ev, n_obs, sigma), tol=0.001)
        return bf.x, bf.fun
    else:
        return sig, chi_sq(sig, dm_ev, n_obs, sigma)

def chi_sq(ln_s, model, n_obs, sigma, ab=False):
    chiinv = np.linalg.inv(sigma)
    vec = 10**ln_s*model - n_obs
    hold1 = np.matmul(chiinv, vec)
    chi = np.matmul(vec, hold1)
    return chi


def build_covariance_matrix(eng, stat, bnds2, dim=24, neg_stat=False):
    stat /= eng**2.
    sigma = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            # Statistical Err
            if i == j:
                if neg_stat:
                    pass
                else:
                    sigma[i,j] += (stat[i])**2.
        
            # Residual Err
            
            if bnds2[i] <= 1. and bnds2[j] <= 1.:
                res_1 = res_flux(eng[i])  # * eng[i]**2.
                res_2 = res_flux(eng[j])  # * eng[j]**2.
                if i == j:
                    sigma[i, j] += 2. * res_1 * res_2
                else:
                    sigma[i, j] += res_1 * res_2
            
            # Model trunc  -- Units of (GeV^-1 cm^-2 s^-1 sr^-1)^2
            sigma[i, j] += pc_error(eng[i], eng[j])  # * eng[i]**2. * eng[j]**2.
    return sigma


def res_flux(x):
    # Takes energy in GeV, returns differential flux in units of GeV^-1 cm^-2 s^-1 sr^-1
    return 6. * 10. ** -8. * x ** -3.


def pc_error(x, y):
    # Takes energy x,y in GeV, returns \sigma_corr err from PCs [each in units of GeV^-1 cm^-2 s^-1 sr^-1]
    pc1 = np.loadtxt(MAIN_PATH + '/FileHolding/PC_stat1.dat')
    pc2 = np.loadtxt(MAIN_PATH + '/FileHolding/PC_stat2.dat')
    pc3 = np.loadtxt(MAIN_PATH + '/FileHolding/PC_stat3.dat')
    pc1in = interp1d(np.log10(pc1[:, 0]), pc1[:, 1], kind='linear', bounds_error=False, fill_value='extrapolate')
    pc2in = interp1d(np.log10(pc2[:, 0]), pc2[:, 1], kind='linear', bounds_error=False, fill_value='extrapolate')
    pc3in = interp1d(np.log10(pc3[:, 0]), pc3[:, 1], kind='linear', bounds_error=False, fill_value='extrapolate')
    
    return (pc1in(np.log10(x))*pc1in(np.log10(y)) + pc2in(np.log10(x))*pc2in(np.log10(y)) +
            pc3in(np.log10(x))*pc3in(np.log10(y))) / (x**2. * y**2.)


def plot_model_vs_data(log_s, spec='BB_direct_mx_50GeV.dat', maj=True, gamma=1.2,
                       l_min=0., l_max=20., b_min=2., b_max=20.):
    
    gce_file = np.loadtxt(MAIN_PATH + 'GC_stat_err.dat')
    gce_dat = np.loadtxt(MAIN_PATH + 'GC_dat.dat')
    ylims = [10 ** -8., 10 ** -5.]
    try:
        load_s = np.loadtxt(MAIN_PATH + '/Spectrum/' + spec)
        dn = interp1d(load_s[:,0], load_s[:,1], kind='cubic')
    except FileNotFoundError:
        print 'Spectrum File Not Found.'
        raise

    findmx = spec.find('mx_')
    findGeV = spec.find('GeV')
    mx = float(spec[findmx+3:findGeV])

    n_obs = gce_dat[:, 1]

    Jfac = Gen_NFW(gamma=gamma).J(l_min=l_min, l_max=l_max, b_min=b_min, b_max=b_max)

    edm = np.logspace(np.log10(0.3), np.log10(30.), 300)
    model = np.zeros((len(edm),2))
    for i in range(len(edm)):
        if edm[i] < mx:
            dne = dn(edm[i])
            model[i,0] = edm[i]
            model[i,1] = dNdE_no_norm(mx, dne, Jfac, maj=maj)
        else:
            pass
    model = model[model[:, 1] > 0.]

    fig = plt.figure()
    ax = plt.gca()

    # Why does my minimzation currently give md2...
    md2 = model[:,0]**2. * 10.**(-25.506) * model[:,1]

    model[:,1] = model[:,0]**2. * 10.**log_s * model[:,1]

    # plot data points
    plt.plot(model[:, 0], model[:,1], color='blue')
    plt.plot(model[:, 0], md2, '--', color='blue')
    plt.plot(gce_dat[:, 0], gce_dat[:, 1], 'o', mfc='r')
    
    sigma = build_covariance_matrix(gce_file[:,1], gce_file[:,3], gce_file[:, 2],
                                    dim=len(n_obs))
    for i in range(len(gce_file[:, 0])):
        eng = gce_file[i, 1]
        cor_er = np.sqrt((pc_error(gce_file[i, 1],gce_file[i, 1]) * eng**4.) + gce_file[i, 3]**2.)
        #cor_er = np.sqrt(sigma[i,i])
        #print gce_file[i, 1], cor_er, pc_error(gce_file[i, 1])* eng**2., gce_file[i, 3]
        if (gce_dat[i, 1] - gce_file[i, 3]) > 0.:
            ylow = np.log10((gce_dat[i, 1] - gce_file[i, 3]) / ylims[0]) / np.log10(ylims[1]/ylims[0])
        else:
            ylow = 0.
        yhigh = np.log10((gce_dat[i,1] + gce_file[i,3]) / ylims[0]) / np.log10(ylims[1]/ylims[0])

        if (gce_dat[i,1] - cor_er) > 0.:
            ylow2 = np.log10((gce_dat[i, 1] - cor_er) / ylims[0]) / np.log10(ylims[1]/ylims[0])
        else:
            ylow2 = 0.
        yhigh2 = np.log10((gce_dat[i, 1] + cor_er) / ylims[0]) / np.log10(ylims[1]/ylims[0])
        plt.axvline(x=gce_file[i, 1], ymin=ylow2, ymax=yhigh2, linewidth=4, color='y',alpha=0.4)
        plt.axvline(x=gce_file[i, 1], ymin=ylow, ymax=yhigh, linewidth=1, color='k')

    ax.set_xscale('log')
    ax.set_yscale('log')
    pl.xlim([0.3, 500.])
    pl.ylim([ylims[0], ylims[1]])
    pl.ylabel(r'$E^2 \frac{dN}{dE}$   [$GeV cm^{-2}s^{-1}sr^{-1}$]', fontsize=fs)
    pl.xlabel(r'$E_\gamma$  [$GeV$]', fontsize=fs)

    figname = MAIN_PATH + '/Test_Plots/Look_at_example_spectrum.pdf'
    fig.set_tight_layout(True)
    pl.savefig(figname)

    return


def plot_errs():
    gce_file = np.loadtxt(MAIN_PATH + 'GC_stat_err.dat')
    pc1 = np.loadtxt(MAIN_PATH + '/FileHolding/PC_stat1.dat')
    pc2 = np.loadtxt(MAIN_PATH + '/FileHolding/PC_stat2.dat')
    pc3 = np.loadtxt(MAIN_PATH + '/FileHolding/PC_stat3.dat')
    ylims = [-1., 4.]
    fig = plt.figure()
    ax = plt.gca()

    engs = np.logspace(np.log10(0.3), np.log10(200.), 200)
    
    pc1plot = interp1d(np.log10(pc1[:,0]), pc1[:,1], kind='linear', bounds_error=False, fill_value='extrapolate')
    pc2plot = interp1d(np.log10(pc2[:,0]), pc2[:,1], kind='linear', bounds_error=False, fill_value='extrapolate')
    pc3plot = interp1d(np.log10(pc3[:,0]), pc3[:,1], kind='linear', bounds_error=False, fill_value='extrapolate')

    plt.plot(engs, pc1plot(np.log10(engs))/10**-7., color='k', lw=2)
    plt.plot(engs, pc2plot(np.log10(engs))/10**-7., color='r', lw=2)
    plt.plot(engs, pc3plot(np.log10(engs))/10**-7., color='g', lw=2)

    plt.plot(engs, engs**2. * res_flux(engs)/10.**-7., 'purple', lw=2)
    
    ax.fill_between(gce_file[:,1], gce_file[:,3]/10.**-7., -gce_file[:,3]/10.**-7.,
                    facecolor='blue', alpha=0.3, interpolate=True, lw=0)

    ax.set_xscale('log')
    ax.set_yscale('linear')
    pl.xlim([0.3, 500.])
    pl.ylim([ylims[0], ylims[1]])
    pl.ylabel(r'$E^2 \frac{dN}{dE}$   [$GeV cm^{-2}s^{-1}sr^{-1}$]', fontsize=fs)
    pl.xlabel(r'$E_\gamma$  [$GeV$]', fontsize=fs)
    
    figname = MAIN_PATH + '/Test_Plots/Systematic_Errors.pdf'
    fig.set_tight_layout(True)
    pl.savefig(figname)

    return


