"""
Plot Figure 7 from Arunima et al paper II
------------------------------------------------

This script visualizes the cumulative distribution corresponding to the marginal posterior probability distribution given by equation 16 in Arunima et al paper II. 
Assumes you have streams' posterior probability distributions (from analysis/run_stream_analysis.py).
This one plots these for 'HSC_2282' as an example but you can do this for any stream.

CHANGE PATH HERE (IN USER CONFIGURATION BLOCK)
Author: Arunima
Date: 21-11-2025
"""

# ======================================================================
# USER CONFIGURATION (EDIT THESE ONLY)
# ======================================================================
h5_file = "path_here/cluster_posteriors.h5" #output of run_stream_analysis.py CHANGE PATH HERE

name = 'HSC_2282' #name of the stream
save_fig = 'path_here/fig7.pdf' #CHANGE PATH HERE

####### =========================
# IMPORTS
# =========================
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy
from scipy.interpolate import interp1d

############################################################################
# plotting control settings
############################################################################
BIGGER_SIZE=14
plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE) 

############################################################################
#loading data and processing
with h5py.File(h5_file) as f:
    P_r=np.array(f[name]['P_r'])
    r_grid = np.array(f[name]['r_eval'])

 #getting CDFs, median and credible intervals
finalcdf = scipy.integrate.cumulative_trapezoid(P_r, r_grid, initial=0)  # integrates properly
finalcdf /= finalcdf[-1]  # force normalization
finalcdf_interp=interp1d(finalcdf,r_grid)
high_1s = finalcdf_interp(0.84)
high_2s = finalcdf_interp(0.975)
med_r_init = finalcdf_interp(0.5)

############################################################################
# actual plotting
############################################################################

plt.close('all')
plt.figure(figsize=(6,5))
plt.plot(r_grid,finalcdf)
plt.axvline(med_r_init,c='g')
plt.axvline(high_1s,c='g',linestyle='--')
plt.axvline(high_2s,c='g',linestyle='--')
plt.xlim(med_r_init-0.5,high_2s+0.5)
plt.xlabel(r'$a = \log d_{\text{init}} $')
plt.ylabel(r'CDF of $\mathcal{P}(a \mid \{x_i\},\mu,\sigma)$')
plt.tight_layout()
plt.savefig(save_fig)




