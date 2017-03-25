#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')
import argparse
from spectrum_stats import *

path = os.getcwd()

parser = argparse.ArgumentParser()

parser.add_argument('--filef', default='BB_direct_mx_')
parser.add_argument('--s_high', type=float, default=5.e-26)
parser.add_argument('--s_low', type=float, default=5.e-27)
parser.add_argument('--n_sigs', type=int, default=30)
parser.add_argument('--gamma', type=float, default=1.2)
parser.add_argument('--maj', default='T')
parser.add_argument('--scale_r', default=20.)
parser.add_argument('--rfix', default=8.5)
parser.add_argument('--rho_fix', default=0.4)
parser.add_argument('--contour_val', nargs='+', type=float, default=np.array([2.3, 6.2, 11.8]))
parser.add_argument('--contour_name', nargs='+', default=np.array(['_1Sigma', '_2Sigma', '_3Sigma']))

args = parser.parse_args()

if args.maj == 'T':
    maj = True
elif args.maj == 'F':
    maj = False


mass_scan(filef=args.filef, gamma=args.gamma, maj=True,
          contour_val=args.contour_val,
          contour_name=args.contour_name,
          scale_r=args.scale_r, rfix=args.rfix, rho_fix=args.rho_fix)