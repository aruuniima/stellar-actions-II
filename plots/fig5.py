"""
Plot Figure 5 from Arunima et al. Paper II: 
Streams' action-difference distributions overlaid on the theoretical
grid of action-space evolution (relative action difference change).

This script assumes you already generated HDF5 outputs from
stellar-actions-II/analysis/action_clustering.py. and optionally from compute_stream_actions.py

CHANGE PATH HERE (IN USER CONFIGURATION BLOCK)
Author: Arunima 
Date: 20-11-2025
"""
# ======================================================================
# USER CONFIGURATION (EDIT THESE ONLY)
# ======================================================================

# change these paths. You can get these outputs from stellar-actions-II/analysis/action_clustering.py
all_bins_R = 'path_here/JR.h5'
all_bins_z = 'path_here/Jz.h5'
all_bins_phi = 'path_here/Jphi.h5'

# this is the result from compute_stream_actions.py using all_actions=False
#it's also has been uploaded to stellar-actions-II/data/
clusters_result = 'path_here/final_clusters_result.csv'

save_dir = 'figures/'  # directory to save figures


# ======================================================================
# END USER CONFIGURATION 
# ======================================================================

#-----------------------------------------------------------------------------------------------------------------------
#IMPORTS
#-----------------------------------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import os


# ======================================================================
# Helper function to read data from HDF5 files
# ======================================================================
def load_results(filepath):
    """Load absolute and relative ΔJ datasets from an HDF5 file."""
    with h5py.File(filepath, "r") as f:
        abs_data = f["birth_bin_abschange_gt30pc"][:]
        rel_data = f["birth_bin_relchange_gt30pc"][:]
    return abs_data, rel_data



 # ======================================================================
# loading data and prep etc
# ======================================================================
#for data from all galaxy 
JR_abs_30, JR_rel_30 = load_results(all_bins_R)
Jz_abs_30, Jz_rel_30 = load_results(all_bins_z)
Jphi_abs_30, Jphi_rel_30 = load_results(all_bins_phi)

birth_bins = np.array((0.15,0.5,2,5,10,50,100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,5000))
n_bins = len(birth_bins) -1
delta_t_range = np.arange(1, 457, 1)

colormap = cm.get_cmap('viridis', len(range(1, n_bins+1, 1)))


obs_clusters=pd.read_csv(clusters_result)
obs_clusters=obs_clusters[obs_clusters['n_with_rv']>10]
mask = (obs_clusters['jr_84p']<0.8)&(obs_clusters['jz_84p']<0.8)&(obs_clusters['age_84p']<459)

plasma = cm.get_cmap("plasma")
cluster_colors = [plasma(x) for x in np.linspace(0.2, 0.9, 5)]


 # ======================================================================
# actual plotting
# ======================================================================

plt.close('all')
fig, axs = plt.subplots(3, 1, figsize=(5,11),sharex='all')
#plotting the simulation lines
for i in range(0,len(birth_bins) - 1,4):
    if birth_bins[i]<1:
        label = f"{birth_bins[i]:.1f}–{birth_bins[i+1]:.1f} pc"
    else:
        label = f"{birth_bins[i]:.0f}–{birth_bins[i+1]:.0f} pc"
    axs[0].plot(delta_t_range, JR_rel_30[:,i],c=colormap(i),linewidth=1.5,label=label)
    if i<np.shape(Jz_abs_30)[1]:
        axs[1].plot(delta_t_range, Jz_rel_30[:,i],c=colormap(i),linewidth=1.5,label=label)
    if i<np.shape(Jphi_abs_30)[1]:
        axs[2].plot(delta_t_range, Jphi_rel_30[:,i],c=colormap(i),linewidth=1.5)
#plotting the observed streams as points
for n in range(5):
    color = cluster_colors[n]
    [cluster_name,jr16,jr50,jr84,jphi16,jphi50,
 jphi84,jz16,jz50,jz84,age16,age50,age84,n_mem,
 n_used,size]=obs_clusters[mask].iloc[n].to_list()
    axs[0].errorbar(age50,jr50,xerr=[[age50-age16],[age84-age50]],
                    yerr=[[jr50-jr16],[jr84-jr50]],fmt='.',color=color,capsize=2,linewidth=1.5,markersize=10)
    axs[1].errorbar(age50,jz50,xerr=[[age50-age16],[age84-age50]],
                    yerr=[[jz50-jz16],[jz84-jz50]],fmt='.',color=color,capsize=2,linewidth=1.5,markersize=10)
    axs[2].errorbar(age50,jphi50,xerr=[[age50-age16],[age84-age50]],
                yerr=[[jphi50-jphi16],[jphi84-jphi50]],fmt='.',color=color,capsize=2,linewidth=1.5,markersize=10,label=cluster_name)

axs[0].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{R}\rangle$")
axs[1].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{z}\rangle$")
axs[2].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{\phi}\rangle$")
axs[2].set_xlabel(r"$\Delta t$  (Time interval; Myr)")
axs[0].legend(ncols=2)
axs[2].legend(ncols=2)
plt.xlim(0,350)

plt.tight_layout()
plt.savefig(f"{save_dir}/all_clusters.pdf")
plt.show()
