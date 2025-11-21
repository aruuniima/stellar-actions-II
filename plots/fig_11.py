"""
Plot Figure 11 from Arunima et al paper II
------------------------------------------------

This script visualizes the observed pairwise action differences logarithm distribution for two representative streams with different initial sizes: Theia 323 with initial size 0.3 pc and HSC 1964 with initial size 363.4 pc.

Assumes you have stream actions (from analysis/compute_stream_actions.py). 
CHANGE PATH HERE (IN USER CONFIGURATION BLOCK) to your local paths.

You can change the streams by changing big and small stream names and sizes in the USER CONFIGURATION BLOCK.

Author: Arunima
Date: 21-11-2025
"""

# ======================================================================
# USER CONFIGURATION (EDIT THESE ONLY)
# ======================================================================
stream_h5_file = 'path_here/stream_actions.h5' #output of compute_stream_actions.py CHANGE PATH HERE

save_fig = 'path_here/fig11.pdf' #CHANGE PATH HERE

big = 'HSC_1964'
small ='Theia_323'

big_stream_size = 363.4  #pc
small_stream_size = 0.3  #pc

####### =========================
# IMPORTS
# =========================
import numpy as np
import matplotlib.pyplot as plt

#local imports
from ..analysis.utils_inference import read_stream_actions

############################################################################
# plotting control settings
############################################################################
BIGGER_SIZE=12
plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE) 


############################################################################
#loading data and processing
delJz_obs,age_array,mu,sigma = read_stream_actions(stream_h5_file,big)
delJz_obs1,age_array1,mu1,sigma1 = read_stream_actions(stream_h5_file,small)

bins = np.linspace(
    min(np.log10(delJz_obs).min(), np.log10(delJz_obs1).min()),
    max(np.log10(delJz_obs).max(), np.log10(delJz_obs1).max()),
    25
)
############################################################################
# actual plotting - Figure 11
############################################################################
plt.close('all')
plt.figure(figsize=(6,5))

plt.hist(
    np.log10(delJz_obs),
    alpha=0.4,
    bins=bins,
    density=True,
    color='grey',
    label=rf'{big} : $d_{{\text{{init50}}}} = {big_stream_size}$ pc'
)

plt.hist(
    np.log10(delJz_obs1),
    alpha=0.4,
    bins=bins,
    density=True,
    color='fuchsia',
    label=rf'{small} : $d_{{\text{{init50}}}} = {small_stream_size}$ pc'
)
plt.legend()
plt.xlabel(r'$x = \log\,\delta_{\rm rel}\Delta J_z$')
plt.ylabel(r'Probability Density Function of $\{x_i\}$')
plt.tight_layout()
plt.savefig(save_fig)
