"""
Plot Figures 6from Arunima et al paper II
------------------------------------------------

CHANGEEEEEEE EVERRTYTHINGGGGG
This script visualizes the time evolution of action differences between
stellar pairs born within various initial separations, using data produced
by `action_clustering.py`.

Data required:
  - HDF5 output files from action_clustering.py for JR, Jz, Jphi
  - Single-star reference .npz files (upload here or stellar-actions-I repo on github)

CHANGE PATH HERE (IN USER CONFIGURATION BLOCK)
Author: Arunima
Date: 13-11-2025
"""

# ======================================================================
# USER CONFIGURATION (EDIT THESE ONLY)
# ======================================================================
#plot comparing different r_init simulation action difference distribution to observations action diff distribution
name = 'HSC_2282'
delJz_obs,age_array,mu,sigma = prob_func.read_stream_actions(stream_h5_file,name)
r_grid, dj_grid, P_cond, prior_r, interp = load_kde_dt(int(mu))
def cdfP_ifr_init(r_init):
    A=np.cumsum(interp((np.log10(r_init),dj_grid)))
    return A/A[-1]
r1,r2,r3,r4 = 0.1,0.3,9,500
# tau=148 Myr ~median age of HSC2282
#
colormap = cm.get_cmap('plasma', len(range(1, 5, 1)))

name = 'HSC_2282'
delJz_obs,age_array,mu,sigma = prob_func.read_stream_actions(stream_h5_file,name)
r_grid, dj_grid, P_cond, prior_r, interp = load_kde_dt(int(mu))
plt.close('all')
fig, axs = plt.subplots(2, 1, figsize=(5,7),sharex='all')
axs[0].set_ylabel(r' CDF of $ \mathcal{P}_\text{1D}(x \mid a,\tau)$')
axs[0].plot(dj_grid,cdfP_ifr_init(r1),c=colormap(0),label=f'{r1} pc')
axs[0].plot(dj_grid,cdfP_ifr_init(r2),c=colormap(1),label=f'{r2} pc')
axs[0].plot(dj_grid,cdfP_ifr_init(r3),c=colormap(2),label=f'{r3} pc')
axs[0].plot(dj_grid,cdfP_ifr_init(r4),c=colormap(3),label=f'{r4} pc')

axs[0].plot(np.sort(np.log10(delJz_obs)),np.linspace(0,1,delJz_obs.shape[0]),label=name,c='k',linewidth=3)
axs[0].fill_between(np.sort(np.log10(delJz_obs)),np.linspace(0,1,delJz_obs.shape[0]),color='grey',alpha=0.5)
# plt.axvline(np.log10(delJz_obs).max())
# plt.hist(np.log10(delJz_obs),density=True,cumulative=True,alpha=0.6,label=name)
axs[0].legend()
axs[0].set_xlim(-3.7,1.1)

delJz_obs,age_array,mu,sigma = prob_func.read_stream_actions(stream_h5_file,name)
r_grid, dj_grid, P_cond, prior_r, interp = load_kde_dt(int(mu),bias=True,max_J_obs=np.log10(delJz_obs.max()))

axs[1].set_ylabel(r' CDF of $ \mathcal{P}_\text{1D}(x \mid a,\tau)$')
axs[1].set_xlabel(r"$x = \log\,\delta_{\rm rel}\Delta J_z$")
axs[1].plot(dj_grid,cdfP_ifr_init(r1),c=colormap(0),label=f'{r1} pc')
axs[1].plot(dj_grid,cdfP_ifr_init(r2),c=colormap(1),label=f'{r2} pc')
axs[1].plot(dj_grid,cdfP_ifr_init(r3),c=colormap(2),label=f'{r3} pc')
axs[1].plot(dj_grid,cdfP_ifr_init(r4),c=colormap(3),label=f'{r4} pc')

axs[1].plot(np.sort(np.log10(delJz_obs)),np.linspace(0,1,delJz_obs.shape[0]),label=name,c='k',linewidth=3)
axs[1].fill_between(np.sort(np.log10(delJz_obs)),np.linspace(0,1,delJz_obs.shape[0]),color='grey',alpha=0.5)
# plt.axvline(np.log10(delJz_obs).max())
# plt.hist(np.log10(delJz_obs),density=True,cumulative=True,alpha=0.6,label=name)
axs[1].legend()
axs[1].set_xlim(-3.7,1.3)
plt.tight_layout()
plt.savefig('/g/data/jh2/ax8338/project2/paper_plots/p1d_2.pdf')
