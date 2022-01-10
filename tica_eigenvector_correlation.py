# written by Ziyuan Zhao, 01082022
# For a molecular system, the tICA eigenvectors associated with the k largest, distinct eigenvalues gives relative
# weights to order parameters that contributes to the k slowest-relaxing / highest autocorrelation dof. Here, the
# script reads in the tICA model calculated from two molecular systems and compute Pearson correlation coeff. between
# every pair of eigenvectors with eigenvalues above a threshold.

import numpy as np
from msmbuilder.io import load_trajs, load_generic
import os, sys, getopt
import seaborn as sns
import matplotlib.pyplot as plt

# parse cmdline args in a C-like manner
try:
    opts, args = getopt.getopt(sys.argv[1:], "vd:i:I:o:")
except getopt.GetoptError as err:
    print(err)
    sys.exit(1)

verbose = False
dir_name = ""
model1_name = ""
model2_name = ""
output_name = "tica_eigenvec_correlation"

for o, a in opts:
    if o == "-v":
        verbose = True # do we even need this??
    elif o == "-d":
        dir_name = a
    elif o == "-i":
        model1_name = a
    elif o == "-I":
        model2_name = a
    elif o == "-o":
        output_name = a

# load two models and select eigenvectors
tica1 = load_generic(f'{model1_name}.pickl')
tica2 = load_generic(f'{model2_name}.pickl')
n1 = tica1.timescales_.shape[0]
n2 = tica2.timescales_.shape[0]
ev1 = tica1.eigenvectors_.T[:n1,:]
ev2 = tica2.eigenvectors_.T[:n2,:]
print(ev1.shape, ev1)
print(ev2.shape, ev2)
# compute pearson correlations for each pair
corr = np.corrcoef(ev1, ev2)[:n1, n1:n1+n2] # take upper right block, see numpy doc to understand

# make heatmap and save
fig1 = plt.figure()
ax1 = fig1.add_subplot()
print(np.abs(corr))
ax1 = sns.heatmap(np.abs(corr), cmap = 'coolwarm') # we care about strength of corr only
fig1.savefig(f"{output_name}.pdf")

