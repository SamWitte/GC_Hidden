#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import argparse
from spectrum_stats import *

path = os.getcwd()

parser = argparse.ArgumentParser()

parser.add_argument('--filef', default='TauTau_direct_mx_')
parser.add_argument('--s_high', type=float, default=8.e-27)
parser.add_argument('--s_low', type=float, default=1.e-27)
parser.add_argument('--n_sigs', type=int, default=150)
parser.add_argument('--gamma', type=float, default=1.2)
parser.add_argument('--maj', default='T')
parser.add_argument('--scale_r', default=20.)
parser.add_argument('--rfix', default=8.5)
parser.add_argument('--rho_fix', default=0.4)
parser.add_argument('--contour_val', nargs='+', type=float, default=np.array([2.3, 6.2, 11.8]))
parser.add_argument('--contour_name', nargs='+', default=np.array(['_1Sigma', '_2Sigma', '_3Sigma']))

args = parser.parse_args()

MASS_SCAN = False
MX_MPHI_SCAN = False

if args.maj == 'T':
    maj = True
elif args.maj == 'F':
    maj = False

if MASS_SCAN:
    mass_scan(filef=args.filef, gamma=args.gamma, maj=True,
              s_low=args.s_low, s_high=args.s_high, n_sigs=args.n_sigs,
              contour_val=args.contour_val,
              contour_name=args.contour_name,
              scale_r=args.scale_r, rfix=args.rfix, rho_fix=args.rho_fix)

if MX_MPHI_SCAN:
    mx_mphi_scroll(filef='BB_cascade_mphi_', gamma=1.2, maj=True,
                   scale_r=20., rfix=8.5, rho_fix=0.4)

    