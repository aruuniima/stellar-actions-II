"""
Plot Figures 2 and 3 from Arunima et al paper II
------------------------------------------------
This script visualizes the time evolution of action differences between
stellar pairs born within various initial separations, using data produced
by `action_clustering.py`.

Data required:
  - HDF5 output files from action_clustering.py for JR, Jz, Jphi
  - Single-star reference .npz files (upload here or stellar-actions-I repo on github)

Author: Arunima
Date: 13-11-2025
"""

#imports-----------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib.cm as cm

##plotting controls------------------------------------------------------------------
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 12,
    'axes.labelsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 10,
    'figure.titlesize': 12,
})

# -----------------------------------------------------------------------------
# User paths (EDIT THESE)
# -----------------------------------------------------------------------------
output_file_R = 'results/JR.h5'
output_file_z = 'results/Jz.h5'
output_file_phi = 'results/Jphi.h5'

# Reference data 
abs_JR = np.load('data/abs_JR.npz')
abs_Jz = np.load('data/abs_Jz.npz')
abs_Jphi = np.load('data/abs_Jphi.npz')
rel_JR = np.load('data/rel_JR.npz')
rel_Jz = np.load('data/rel_Jz.npz')
rel_Jphi = np.load('data/rel_Jphi.npz')

save_dir = 'figures/'  # directory to save figures

#-------------------------------------------------------------------------------------
birth_bins = np.array((0.15,0.5,2,5,10,50,100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,5000))

bin_centres = 0.5 * (birth_bins[:-1] + birth_bins[1:])
n_bins = len(birth_bins) -1
delta_t_range = np.arange(1, 457, 1)

colormap = cm.get_cmap('viridis', len(range(1, n_bins+1, 1)))
colors = colormap(np.linspace(0, 1, len(range(1, n_bins+1, 1))))


bounds = birth_bins[1:n_bins+1]   
n_colors = len(bounds)

colormap = cm.get_cmap('viridis', n_colors)
norm = cm.colors.BoundaryNorm(boundaries=birth_bins[:n_colors+1], ncolors=n_colors)

# -----------------------------------------------------------------------------
# Load results
# -----------------------------------------------------------------------------
def load_results(filepath):
    """Load absolute and relative ΔJ datasets from an HDF5 file."""
    with h5py.File(filepath, "r") as f:
        abs_data = f["birth_bin_abschange_gt30pc"][:]
        rel_data = f["birth_bin_relchange_gt30pc"][:]
    return abs_data, rel_data

JR_abs_30, JR_rel_30 = load_results(output_file_R)
Jz_abs_30, Jz_rel_30 = load_results(output_file_z)
Jphi_abs_30, Jphi_rel_30 = load_results(output_file_phi)

# Single-star reference curves
singlestar_JR = abs_JR['Jphi']
singlestar_Jz =abs_Jz['Jphi']
singlestar_Jphi = abs_Jphi['Jphi']


#------------------------------------------------------------------------------------
#FIGURE 2 absolute action difference evolution
plt.close('all')
fig, axs = plt.subplots(3, 1, figsize=(5,11),sharex='all')

for i in range(0,len(birth_bins) - 1,2):
    if birth_bins[i]<1:
        label = f"{birth_bins[i]:.1f}–{birth_bins[i+1]:.1f} pc"
    else:
        label = f"{birth_bins[i]:.0f}–{birth_bins[i+1]:.0f} pc"
    axs[0].plot(delta_t_range, JR_abs_30[:,i],c=colors[i],label=label)
    if i<np.shape(Jz_abs_30)[1]:
        axs[1].plot(delta_t_range, Jz_abs_30[:,i],c=colors[i],label=label)
    if i<np.shape(Jphi_abs_30)[1]:
        axs[2].plot(delta_t_range, Jphi_abs_30[:,i],c=colors[i],label=label)
    
axs[0].plot(abs_JR['dt'], singlestar_JR,'--b',linewidth=2,label='single star')
axs[0].set_ylabel(r"$\langle \delta_{\text{abs}}\Delta J_{R}\rangle$ (kpc km/s)")
axs[1].plot(abs_Jz['dt'], singlestar_Jz,'--b',linewidth=2,label='single star')
axs[1].set_ylabel(r"$\langle \delta_{\text{abs}}\Delta J_{z}\rangle$ (kpc km/s)")
axs[2].plot(abs_Jphi['dt'], singlestar_Jphi,'--b',linewidth=2,label='single star')
axs[2].set_ylabel(r"$\langle \delta_{\text{abs}}\Delta J_{\phi}\rangle$ (kpc km/s)")
axs[2].set_xlabel(r"$\Delta t$ (Myr)")
axs[0].legend(ncols=2)

axs[0].set_xlim((-5,450))
axs[0].set_ylim((-1,20))
axs[1].set_ylim((-0.05,1))

plt.tight_layout()
##### CHANGE PATH HERE
plt.savefig(f"{save_dir}/fig2_abs_dJ.pdf", bbox_inches='tight')
print(f"Saved Figure 2 to {save_dir}/fig2_abs_dJ.pdf")

#------------------------------------------------------------------------------------
#FIGURE 3 relative action difference evolution
plt.close('all')
fig, axs = plt.subplots(3, 1, figsize=(5,11),sharex='all')

for i in range(0,len(birth_bins) - 1,2):
    if birth_bins[i]<1:
        label = f"{birth_bins[i]:.1f}–{birth_bins[i+1]:.1f} pc"
    else:
        label = f"{birth_bins[i]:.0f}–{birth_bins[i+1]:.0f} pc"
    axs[0].plot(delta_t_range, JR_rel_30[:,i],c=colors[i],label=label)
    if i<np.shape(Jz_abs_30)[1]:
        axs[1].plot(delta_t_range, Jz_rel_30[:,i],c=colors[i],label=label)
    if i<np.shape(Jphi_abs_30)[1]:
        axs[2].plot(delta_t_range, Jphi_rel_30[:,i],c=colors[i],label=label)

axs[0].plot(rel_JR['dt'],rel_JR['JR'] ,'--b',linewidth=2,label='single star')
axs[1].plot(rel_Jz['dt'],rel_Jz['Jz'] ,'--b',linewidth=2,label='single star')

axs[2].plot(rel_Jphi['dt'],rel_Jphi['Jphi'] ,'--b',linewidth=2,label='single star')

axs[0].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{R}\rangle$")
axs[1].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{z}\rangle$")
axs[2].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{\phi}\rangle$")
axs[2].set_xlabel(r"$\Delta t$ (Myr)")
axs[2].legend(ncols=2)
axs[0].set_xlim((-5,450))


plt.tight_layout()
##### CHANGE PATH HERE
plt.savefig(f"{save_dir}/fig3_rel_dJ.pdf", bbox_inches='tight')
print(f"Saved Figure 3 to {save_dir}/fig3_rel_dJ.pdf")
