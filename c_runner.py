#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import argparse
from spectrum_stats import *

path = os.getcwd()

parser = argparse.ArgumentParser()

# eg ['BB_direct_mx_', 'TauTau_direct_mx_', 'BB_cascade_mphi_', 'Zprime_cascade_mphi_','Photon_cascade_mphi_']
parser.add_argument('--filef', default='Photon_cascade_mphi_')

parser.add_argument('--s_high', type=float, default=5.e-26)
parser.add_argument('--s_low', type=float, default=5.e-27)
parser.add_argument('--n_sigs', type=int, default=150)
parser.add_argument('--gamma', type=float, default=1.2)
parser.add_argument('--maj', default='T')
parser.add_argument('--scale_r', default=20.)
parser.add_argument('--rfix', default=8.5)
parser.add_argument('--rho_fix', default=0.4)
parser.add_argument('--contour_val', nargs='+', type=float, default=np.array([2.3, 6.2, 11.8]))
parser.add_argument('--contour_name', nargs='+', default=np.array(['_1Sigma', '_2Sigma', '_3Sigma']))

cut_tail = False

args = parser.parse_args()

MASS_SCAN = True
MX_MPHI_SCAN = False

if args.maj == 'T':
    maj = True
else:
    maj = False

if MASS_SCAN:
    mass_scan(filef=args.filef, gamma=args.gamma, maj=maj,
              s_low=args.s_low, s_high=args.s_high, n_sigs=args.n_sigs,
              contour_val=args.contour_val,
              contour_name=args.contour_name,
              scale_r=args.scale_r, rfix=args.rfix, rho_fix=args.rho_fix, cut_tail=cut_tail)

if MX_MPHI_SCAN:
    mx_mphi_scroll(filef=args.filef, gamma=args.gamma, maj=maj,
                   scale_r=args.scale_r, rfix=args.rfix, rho_fix=args.rho_fix,
                   cut_tail=cut_tail)

