import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl


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


def mx_sigma_contours(end_p='BB', trans='_direct_', contour_name=np.array(['_1Sigma', '_2Sigma', '_3Sigma']),
                      gamma=1.2, scale_r=20., rfix=8.5, rho_fix=0.4):

    filename = MAIN_PATH + '/Test_Plots/' + end_p + trans + 'mx'
    filename += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
    filename += '.pdf'

    pl.figure()
    ax = pl.gca()

    linestyle = ['solid', 'dashed', 'dotted']
    colorl = ['blue', 'red']
    for i,cc in enumerate(contour_name):
        foi = MAIN_PATH + '/FileHolding/Contours/' + end_p + trans + 'mx_' + cc
        foi += 'Gamma_{:.2f}_ScaleR_{:.2f}_Rfix_{:.2f}_RhoFix_{:.2f}'.format(gamma, scale_r, rfix, rho_fix)
        foi += '.dat'

        loadf = np.loadtxt(foi)
        loadf = loadf[loadf[:, 1] != 1.]

        low = interp1d(np.log10(loadf[:, 0]), np.log10(loadf[:, 1]), kind='cubic', bounds_error=False)
        high = interp1d(np.log10(loadf[:, 0]), np.log10(loadf[:, 2]), kind='cubic', bounds_error=False)
        massr = np.linspace(np.min(loadf[:, 0]), np.max(loadf[:, 0]), 200)

        ax.fill_between(massr, 10.**low(np.log10(massr)), 10.**high(np.log10(massr)),
                        where=10.**high(np.log10(massr))>= 10.**low(np.log10(massr)),
                        facecolor='blue', interpolate=True, alpha=0.2)
        #pl.plot(massr, 10.**low(np.log10(massr)), 'r', lw=1, ls=linestyle[i])
        #pl.plot(massr, 10.**high(np.log10(massr)), 'r', lw=1, ls=linestyle[i])

    ax.set_xlabel(r'$m_\chi$  [GeV]', fontsize=fs)
    ax.set_ylabel(r'$<\sigma {\rm v}>$  [$cm^{3} s^{-1}$]', fontsize=fs)

    plt.xlim(xmin=4., xmax=100.)
    plt.ylim(ymin=1.*10**-27., ymax=1.*10**-25.)
    ax.set_xscale("log")
    ax.set_yscale("log")

    plt.savefig(filename)

    return