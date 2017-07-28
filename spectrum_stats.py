###############################################




###############################################


import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d, UnivariateSpline, interp2d, bisplev, bisplrep
from scipy.optimize import minimize, fsolve
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


def mx_mphi_scroll(filef='BB_cascade_mphi_', gamma=1.2, maj=True,
                   scale_r=20., rfix=8.5, rho_fix=0.4,
                   contour_val=np.array([26.66, 33.15, 36.4]),
                   contour_name=np.array(['_p32', '_p10', '_p05']), cut_tail=False):
    if cut_tail:
        contour_val[0] = 22.4
        contour_val[1] = 28.4
        contour_val[2] = 31.4
    
    all_files = glob.glob(MAIN_PATH + '/Spectrum/' + filef + '*.dat')

    bf_array = np.zeros(len(all_files))
    # Order Mphi, Mx
    mass_list = np.zeros((len(all_files), 2))
    for i, f in enumerate(all_files):
        findmx = f.find('mx_')
        findGeV = f.find('GeV')
        findmphi = f.find('mphi_')
        mass_list[i] = [np.float(f[findmphi + 5:(findmx - 1)]), np.float(f[findmx + 3:findGeV])]
        print 'Mphi, Mx, BF'
        f_tail = f[f.find(filef):]
        bf_array[i] = sig_contour(spec=f_tail, gamma=gamma, maj=maj, scale_r=scale_r, rfix=rfix,
                                  rho_fix=rho_fix, ret_bf=True, cut_tail=cut_tail)
        print mass_list[i][0], mass_list[i][1], bf_array[i]

    labmp = filef.find('mphi_')

    ffsve = np.column_stack((mass_list[:, 1], mass_list[:, 0], bf_array))
    ffsve = ffsve[np.argsort(ffsve[:, 0])]
    if cut_tail
        f_save_n = MAIN_PATH + '/FileHolding/' + filef[:labmp] + '_CUT_Tabbed_mx_mphi.dat'
    else:
        f_save_n = MAIN_PATH + '/FileHolding/' + filef[:labmp] + 'Tabbed_mx_mphi.dat'
    np.savetxt(f_save_n, ffsve)
    mxlist = np.unique(mass_list[:, 1])
    mx_bflist = np.zeros((len(mxlist), 2))

    sig_cnt = np.zeros((len(mxlist), 6))

    for i, mx in enumerate(mxlist):
        mph_u = mass_list[:, 0][mass_list[:, 1] == mx]
        bf_temp = bf_array[mass_list[:, 1] == mx]
        ordr = np.argsort(mph_u)
        mph_u = mph_u[ordr]
        bf_temp = bf_temp[ordr]

        d1interp = interp1d(mph_u, bf_temp, kind='cubic', bounds_error=False, fill_value=1e5)
        bf_fixmx = minimize(d1interp, np.array([np.median(mph_u)]), tol=1.e-3)
        print 'Best fit point at mx {:.2f} is'.format(mx)
        print bf_fixmx
        mx_bflist[i] = [bf_fixmx.x, bf_fixmx.fun]

    mx_interp = interp1d(mxlist, mx_bflist[:, 1], kind='cubic', bounds_error=False, fill_value=1e5)
    overbf = minimize(mx_interp, np.array([np.median(mxlist)]), tol=1.e-4)
    print '~~~~~~~~~~~~~~~~~~~~~~~~\n'
    print 'Overall Best Fit:'
    print overbf
    goal = overbf.fun

    tot_contours = len(contour_name)
    for i, mx in enumerate(mxlist):
        for j, cc in enumerate(contour_val):
            print 'Goal: ', cc
            if mx_bflist[i, 1] < cc:
                mph_u = mass_list[:, 0][mass_list[:, 1] == mx]
                bf_temp = bf_array[mass_list[:, 1] == mx]
                ordr = np.argsort(mph_u)
                mph_u = mph_u[ordr]
                bf_temp = bf_temp[ordr]

                #bf_temp += - cc
                #bf_temp = np.abs(bf_temp)

                d1interp = lambda x: np.abs(10.**interp1d(np.log10(mph_u), np.log10(bf_temp), kind='cubic',
                                                   bounds_error=False, fill_value=1e5)(np.log10(x)) - cc)

                print 'Contour ChiSq: ', cc

                if bf_temp[0] < cc:
                    sig_cnt[i, j] = 0.
                    print 'Low: ', 0.
                else:
                    slch = fminbound(d1interp, np.min(mph_u), mx_bflist[i, 0], full_output=True, xtol=1e-3)
                    sig_cnt[i, j] = slch[0]
                    print 'Low: ', slch

                if bf_temp[-1] < cc:
                    sig_cnt[i, j + tot_contours] = mx
                    print 'High: ', mx
                else:
                    shch = fminbound(d1interp, mx_bflist[i, 0], np.max(mph_u), full_output=True, xtol=1e-3)
                    sig_cnt[i, j + tot_contours] = shch[0]
                    print 'High: ', shch

            print 'mx: ', mx, ' mphi contours: ', sig_cnt[i]

    fnl_arr = np.column_stack((mxlist, sig_cnt))
    print fnl_arr

    print 'Saving Files...'

    for i, cc in enumerate(contour_name):
        c_fname = MAIN_PATH + '/FileHolding/Contours/' + filef + cc
        if cut_tail:
            c_fname += '_CUT_'
        c_fname += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
        c_fname += '.dat'
        consv = np.column_stack((fnl_arr[:, 0], fnl_arr[:, i+1], fnl_arr[:, i+1+tot_contours]))
        consv = consv[np.argsort(consv[:, 0])]
        consv = consv[consv[:, -1] > 0.]
        np.savetxt(c_fname, consv)

    print 'Done.'
    return


def mass_scan(filef='BB_direct_mx_', gamma=1.2, maj=True, s_low=5.e-27,
              s_high=5.e-26, n_sigs=30,
              contour_val=np.array([2.3, 6.2, 11.8]),
              contour_name=np.array(['_1Sigma', '_2Sigma', '_3Sigma']),
              scale_r=20., rfix=8.5, rho_fix=0.4, cut_tail=False):

    all_files = glob.glob(MAIN_PATH + '/Spectrum/' + filef + '*.dat')

    bf_array = np.zeros(len(all_files))
    mx_list = np.zeros(len(all_files))
    # Make All Files...
    print 'Making ChiSq Tabs...'
    for i,f in enumerate(all_files):
        findmx = f.find('mx_')
        findGeV = f.find('GeV')
        mx_list[i] = np.float(f[findmx + 3:findGeV])
        f_tail = f[f.find(filef):]
        print 'On Mass {:.2f}...'.format(mx_list[i])
        bf_array[i] = sig_contour(spec=f_tail, gamma=gamma, maj=maj, s_low=s_low, s_high=s_high, n_sigs=n_sigs,
                                  contour_val=contour_val, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix,
                                  make_file=True, cut_tail=cut_tail)
    #print mx_list
    #print bf_array
    order = np.argsort(mx_list)
    mx_list = mx_list[order]
    bf_array = bf_array[order]
    tot_arr = np.column_stack((mx_list, bf_array))
    #print 'Mass, BF'
    #print tot_arr

    goal_look = interp1d(tot_arr[:, 0], tot_arr[:, 1], kind='cubic', bounds_error=False, fill_value=1.e5)
    goal = minimize(goal_look, np.median(tot_arr[:, 0]))
    print 'Best Fit Mass: ', goal.x[0]
    print 'Best Fit ChiSq Value: ', goal.fun
    info_hold = np.zeros((len(all_files), 6))
    print 'Constructing Contours...'
    for i, f in enumerate(all_files):
        findmx = f.find('mx_')
        findGeV = f.find('GeV')
        mx_list[i] = np.float(f[findmx + 3:findGeV])
        f_tail = f[f.find(filef):]
        hold = sig_contour(spec=f_tail, gamma=gamma, maj=maj, s_low=s_low, s_high=s_high, n_sigs=n_sigs,
                                   contour_val=contour_val, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix,
                                   make_file=False, goal=goal.fun, cut_tail=cut_tail)
        info_hold[i] = hold[-1]
        print 'On Mass {:.2f}...'.format(mx_list[i])
        print 'Info Hold: ', info_hold[i]

    print 'Saving Files...'
    contours = info_hold
    tot_contours = len(contour_name)
    for i, cc in enumerate(contour_name):
        c_fname = MAIN_PATH + '/FileHolding/Contours/' + filef + cc
        if cut_tail:
            c_fname += '_CUT_'
        c_fname += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
        c_fname += '.dat'
        consv = np.stack((mx_list, contours[:, i], contours[:, i+tot_contours]), axis=-1)
        consv = consv[np.argsort(consv[:,0])]
        np.savetxt(c_fname, consv)

    print 'Done.'
    return


def sig_contour(spec='BB_direct_mx_50GeV.dat', gamma=1.2, maj=True,
                s_low=5.e-27, s_high=5.e-26, n_sigs=10,
                contour_val=np.array([2.3, 6.2, 11.8]),
                scale_r=20., rfix=8.5, rho_fix=0.4, make_file=True,
                goal=24., ret_bf=False, ret_cs=False, cut_tail=False):

    file_name = MAIN_PATH + '/FileHolding/Contours/ChiSq/Tabbed_ChiSq_'
    file_name += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
    file_name += '_' + spec

    findmx = spec.find('mx_')
    findGeV = spec.find('GeV')
    mx = float(spec[findmx + 3:findGeV])

    if ret_bf:
        bf = chi_covariance(spec=spec, maj=maj, gamma=gamma, bf=True, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix,
                            cut_tail=cut_tail)
        if ret_cs:
            return bf
        else:
            return bf[-1]

    if make_file:

        bf = chi_covariance(spec=spec, maj=maj, gamma=gamma, bf=True, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix,
                            cut_tail=cut_tail)

        if os.path.isfile(file_name):
            print 'File already made.'
            return bf[-1]

        print 'Making File for mx {:.2f}'.format(mx)
        s_list = np.linspace(np.log10(s_low), np.log10(s_high), n_sigs)
        hold_tab = np.zeros((len(s_list), 2))
        for i, sig in enumerate(s_list):

            hold_tab[i] = chi_covariance(spec=spec, maj=maj, gamma=gamma, bf=False, sig=sig,
                                         scale_r=scale_r, rfix=rfix, rho_fix=rho_fix, cut_tail=cut_tail)
            print 'Sigma {:.2f}, ChiSq: {:.2f}'.format(sig, hold_tab[i, 1])

        np.savetxt(file_name, hold_tab)
        return bf[-1]
    else:
        print 'Finding contour for mx {:.2f}'.format(mx)
        chisq = np.loadtxt(file_name)
        chi_interp = interp1d(chisq[:, 0], np.log10(chisq[:, 1]),
                              kind='cubic', fill_value=10.**5., bounds_error=False)
        #chi_interp = interp1d(chisq[:, 0], chisq[:, 1], kind='linear', fill_value=10. ** 5., bounds_error=False)

        bf = chi_covariance(spec=spec, maj=maj, gamma=gamma, bf=True, scale_r=scale_r, rfix=rfix, rho_fix=rho_fix, cut_tail=cut_tail)
        print 'ChiSq Min: {:.2f}'.format(bf[1])

        sl = np.zeros(len(contour_val))
        sh = np.zeros(len(contour_val))

        bndsl = [(np.log10(s_low), bf[0])]
        bndsh = [(bf[0], np.log10(s_high))]
        #print 'Bounds: ', bndsl, bndsh

        for j, cc in enumerate(contour_val):
            print 'BF ChiSq: {:.2f}, Target: {:.2f}'.format(bf[1],goal+cc)
            if bf[1] < (goal + cc):
                slch = fminbound(lambda x: np.abs(10.**chi_interp(x) - cc - goal),
                                 bndsl[0][0],bndsl[0][1], full_output=True)
                shch = fminbound(lambda x: np.abs(10.**chi_interp(x) - cc - goal),
                                 bndsh[0][0],bndsh[0][1], full_output=True)


                print slch, shch
                #shch = minimize(lambda x: np.abs(chi_interp(x) - cc - goal), np.array([bf[0] + 0.5]),
                #                 method='SLSQP', bounds=bndsh, tol=1.e-4)
                if slch[2] == 0:
                    sl[j] = slch[0][0]
                if shch[2] == 0:
                    sh[j] = shch[0][0]

        # Return format: Gamma, mx, Min chi^2, BF sigma_p, 1\sig_L, 2\sig_L, 3\sig_L, 1\sig_H, 2\sig_H,3\sig_H
        return gamma, mx, bf[1], 10.**bf[0][0], np.concatenate((10.**sl, 10.**sh))


def chi_covariance(spec='BB_direct_mx_50GeV.dat', maj=True, l_min=0.,
                   l_max=20., b_min=2., b_max=20., gamma=1.2, bf=True, sig=-26.,
                   scale_r=20., rfix=8.5, rho_fix=0.4, cut_tail=False):

    gce_file = np.loadtxt(MAIN_PATH + 'GC_stat_err.dat')
    gce_dat = np.loadtxt(MAIN_PATH + 'GC_dat.dat')
    if cut_tail:
        gce_dat = gce_dat[:-4]
        gce_file = gce_file[:-4]
    
    Jfac = Gen_NFW(gamma=gamma, scale_r=scale_r, rfix=rfix,
                   rho_fix=rho_fix).J(l_min=l_min, l_max=l_max, b_min=b_min, b_max=b_max)

    try:
        load_s = np.loadtxt(MAIN_PATH + '/Spectrum/' + spec)
        dn = interp1d(load_s[:, 0], load_s[:, 1], kind='cubic', bounds_error=False, fill_value=0.)
    except:
        print 'Spectrum File Not Found.'
        print MAIN_PATH + '/Spectrum/' + spec
        raise ValueError

    findmx = spec.find('mx_')
    findGeV = spec.find('GeV')
    mx = float(spec[findmx+3:findGeV])

    n_obs = np.copy(gce_dat[:, 1])
    n_obs /= gce_dat[:, 0]**2.

    dim = len(n_obs)

    edm = np.logspace(np.log10(0.3), np.log10(mx - 1.), 300)
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
    sigma = build_covariance_matrix(gce_file[:,1], gce_file[:,3], gce_file[:, 2], dim=gce_file.shape[0])

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
    except:
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

    model[:,1] = model[:,0]**2. * 10.**log_s * model[:,1]

    # plot data points
    plt.plot(model[:, 0], model[:,1], color='blue')
    plt.plot(gce_dat[:, 0], gce_dat[:, 1], 'o', mfc='r')
    
    #sigma = build_covariance_matrix(gce_file[:,1], gce_file[:,3], gce_file[:, 2],
    #                                dim=len(n_obs))
    cor_er = np.zeros(len(gce_file[:, 0]))
    for i in range(len(gce_file[:, 0])):
        eng = gce_file[i, 1]
        cor_er[i] = np.sqrt((pc_error(gce_file[i, 1], gce_file[i, 1]) * eng**4.) + gce_file[i, 3]**2.)

        if (gce_dat[i, 1] - gce_file[i, 3]) > 0.:
            ylow = np.log10((gce_dat[i, 1] - gce_file[i, 3]) / ylims[0]) / np.log10(ylims[1]/ylims[0])
        else:
            ylow = 0.
        yhigh = np.log10((gce_dat[i,1] + gce_file[i,3]) / ylims[0]) / np.log10(ylims[1]/ylims[0])

        if (gce_dat[i,1] - cor_er[i]) > 0.:
            ylow2 = np.log10((gce_dat[i, 1] - cor_er[i]) / ylims[0]) / np.log10(ylims[1]/ylims[0])
        else:
            ylow2 = 0.
        yhigh2 = np.log10((gce_dat[i, 1] + cor_er[i]) / ylims[0]) / np.log10(ylims[1]/ylims[0])
        plt.axvline(x=gce_file[i, 1], ymin=ylow2, ymax=yhigh2, linewidth=4, color='y',alpha=0.4)
        plt.axvline(x=gce_file[i, 1], ymin=ylow, ymax=yhigh, linewidth=1, color='k')

    np.savetxt(MAIN_PATH + '/FileHolding/Diag_System_Err.dat', cor_er)

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


def bf_mphi(filef='BB_cascade_', gamma=1.2, maj=True,
                   scale_r=20., rfix=8.5, rho_fix=0.4):
    num_print = 5
    load_tab = np.loadtxt(MAIN_PATH + '/FileHolding/' + filef + 'Tabbed_mx_mphi.dat')
    for i in range(num_print):
        agm = np.argmin(load_tab[:, 2])
        print load_tab[agm]
        load_tab = np.delete(load_tab, agm, 0)

    return
