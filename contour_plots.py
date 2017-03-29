import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from spectrum_stats import *

import matplotlib.patheffects as PathEffects
import matplotlib.gridspec as gridspec
import glob
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times','Palatino']})
rc('text', usetex=True)

mpl.rcParams['xtick.major.size']=8
mpl.rcParams['ytick.major.size']=8
mpl.rcParams['xtick.labelsize']=18
mpl.rcParams['ytick.labelsize']=18

fs = 18.

try:
    MAIN_PATH = os.environ['GC_HIDDEN_PATH']
except KeyError:
    print 'GC_HIDDEN_PATH environment variable not defined, defaulting to cwd'
    MAIN_PATH = os.getcwd() + '/'


def mx_sigma_contours(end_p=np.array(['BB']), trans='_direct_',
                      contour_name=np.array(['_1Sigma', '_2Sigma', '_3Sigma']),
                      gamma=1.2, scale_r=20., rfix=8.5, rho_fix=0.4,
                      colorl=np.array(['blue', 'red'])):

    filename = MAIN_PATH + '/Test_Plots/'

    pl.figure()
    ax = pl.gca()

    for j, part in enumerate(end_p):
        filename += part + '_'
        for i,cc in enumerate(contour_name):
            foi = MAIN_PATH + '/FileHolding/Contours/' + part + trans + 'mx_' + cc
            foi += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
            foi += '.dat'

            loadf = np.loadtxt(foi)
            loadf = loadf[loadf[:, 1] != 1.]

            low = interp1d(np.log10(loadf[:, 0]), np.log10(loadf[:, 1]), kind='cubic', bounds_error=False)
            high = interp1d(np.log10(loadf[:, 0]), np.log10(loadf[:, 2]), kind='cubic', bounds_error=False)
            massr = np.linspace(np.min(loadf[:, 0]), np.max(loadf[:, 0]), 200)

            ax.fill_between(massr, 10.**low(np.log10(massr)), 10.**high(np.log10(massr)),
                            where=10.**high(np.log10(massr)) >= 10.**low(np.log10(massr)),
                            facecolor=colorl[j], interpolate=True, alpha=0.2)

    ax.set_xlabel(r'$m_\chi$  [GeV]', fontsize=fs)
    ax.set_ylabel(r'$<\sigma {\rm v}>$  [$cm^{3} s^{-1}$]', fontsize=fs)

    plt.xlim(xmin=4., xmax=100.)
    plt.ylim(ymin=1.*10**-27., ymax=1.*10**-25.)
    ax.set_xscale("log")
    ax.set_yscale("log")

    filename += trans + 'mx'
    filename += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
    filename += '.pdf'

    plt.savefig(filename)

    return


def csec_v_mphi_various_gamma(spec='Zprime_cascade_mphi', mx=22,
                              gam=np.array([1., 1.1, 1.2, 1.3, 1.4]),
                              maj=True, scale_r=20., rfix=8.5, rho_fix=0.4):

    all_files = glob.glob(MAIN_PATH + '/Spectrum/' + spec + '*_mx_' + str(mx) + 'GeV.dat')
    if len(all_files) == 0:
        all_files = glob.glob(MAIN_PATH + '/Spectrum/' + spec + '*_mx_{:.0f}'.format(mx) + '.GeV.dat')
    filen = MAIN_PATH + '/FileHolding/CrossSec_v_mphi_' + spec[:-3] + 'x_' + str(mx) + 'GeV'


    for j, gamma in enumerate(gam):
        fix_g_list = np.zeros((len(all_files), 3))
        for i, f in enumerate(all_files):
            findmphi = f.find('mphi_')
            findmx = f.find('mx_')
            findGeV = f.find('GeV')
            mph = np.float(f[findmphi + 5:(findmx - 1)])
            #print 'Mphi = {:.2f}'.format(mph)
            f_tail = f[f.find(spec):]

            bf_array = sig_contour(spec=f_tail, gamma=gamma, maj=maj, scale_r=scale_r, rfix=rfix,
                                   rho_fix=rho_fix, ret_bf=True, ret_cs=True)

            fix_g_list[i] = [mph, bf_array[0], bf_array[1]]
        fix_g_list = fix_g_list[np.argsort(fix_g_list[:, 0])]
        svname = filen + '_gamma_{:.2f}_rfix_{:.2f}_rhofix_{:.2f}.dat'.format(gamma, rfix, rho_fix)
        np.savetxt(svname, fix_g_list)

    return
