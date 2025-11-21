"""
Plot Figures 6 from Arunima et al paper II
------------------------------------------------

This script visualizes the cumulative distributions of 1D conditional probability distribution of log vertical action difference given age and initial separation given by eq 13 in Arunima et al paper II. Assumes you have stream actions (from analysis/compute_stream_actions.py) and kde_dt_XXX.h5 result from analysis/kde_for_dt.py. 
This one plots these for 'HSC_2282' as an example (and age = 148 Myr so we only need kde_dt_148.h5) but you can do this for any stream.

CHANGE PATH HERE (IN USER CONFIGURATION BLOCK)
Author: Arunima
Date: 21-11-2025
"""

# ======================================================================
# USER CONFIGURATION (EDIT THESE ONLY)
# ======================================================================
stream_h5_file = 'path_here/stream_actions.h5' #output of compute_stream_actions.py CHANGE PATH HERE
kde_dir = 'path_here' #CHANGE PATH HERE 

name = 'HSC_2282' #name of the stream
save_fig = 'path_here/fig6.pdf' #CHANGE PATH HERE

####### =========================
# IMPORTS
# =========================
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#local imports
from ..analysis.utils_inference import read_stream_actions, load_kde_dt



# ======================================================================
# helper function: evaluate conditional CDF at fixed r_init
# ======================================================================
def cdfP_ifr_init(r_init,interp):
    A=np.cumsum(interp((np.log10(r_init),dj_grid)))
    return A/A[-1]

#plot comparing different r_init simulation action difference distribution to observations action diff distribution
delJz_obs,age_array,mu,sigma = read_stream_actions(stream_h5_file,name)
r_grid, dj_grid, P_cond, interp = load_kde_dt(kde_dir, int(mu))
 
r1,r2,r3,r4 = 0.1,0.3,9,500  #different initial distance plotted
colormap = cm.get_cmap('plasma', len(range(1, 5, 1)))

plt.close('all')
fig, axs = plt.subplots(2, 1, figsize=(5,7),sharex='all')
axs[0].set_ylabel(r' CDF of $ \mathcal{P}_\text{1D}(x \mid a,\tau)$')
axs[0].plot(dj_grid,cdfP_ifr_init(r1,interp),c=colormap(0),label=f'{r1} pc')
axs[0].plot(dj_grid,cdfP_ifr_init(r2,interp),c=colormap(1),label=f'{r2} pc')
axs[0].plot(dj_grid,cdfP_ifr_init(r3,interp),c=colormap(2),label=f'{r3} pc')
axs[0].plot(dj_grid,cdfP_ifr_init(r4,interp),c=colormap(3),label=f'{r4} pc')
axs[0].plot(np.sort(np.log10(delJz_obs)),np.linspace(0,1,delJz_obs.shape[0]),label=name,c='k',linewidth=3)
axs[0].fill_between(np.sort(np.log10(delJz_obs)),np.linspace(0,1,delJz_obs.shape[0]),color='grey',alpha=0.5)
axs[0].legend()
axs[0].set_xlim(-3.7,1.1)

delJz_obs,age_array,mu,sigma = read_stream_actions(stream_h5_file,name)
r_grid, dj_grid, P_cond, interp = load_kde_dt(kde_dir, int(mu),bias=True,max_J_obs=np.log10(delJz_obs.max()))

axs[1].set_ylabel(r' CDF of $ \mathcal{P}_\text{1D}(x \mid a,\tau)$')
axs[1].set_xlabel(r"$x = \log\,\delta_{\rm rel}\Delta J_z$")
axs[1].plot(dj_grid,cdfP_ifr_init(r1,interp),c=colormap(0),label=f'{r1} pc')
axs[1].plot(dj_grid,cdfP_ifr_init(r2,interp),c=colormap(1),label=f'{r2} pc')
axs[1].plot(dj_grid,cdfP_ifr_init(r3,interp),c=colormap(2),label=f'{r3} pc')
axs[1].plot(dj_grid,cdfP_ifr_init(r4,interp),c=colormap(3),label=f'{r4} pc')

axs[1].plot(np.sort(np.log10(delJz_obs)),np.linspace(0,1,delJz_obs.shape[0]),label=name,c='k',linewidth=3)
axs[1].fill_between(np.sort(np.log10(delJz_obs)),np.linspace(0,1,delJz_obs.shape[0]),color='grey',alpha=0.5)
axs[1].legend()
axs[1].set_xlim(-3.7,1.3)
plt.tight_layout()
plt.savefig(save_fig)
