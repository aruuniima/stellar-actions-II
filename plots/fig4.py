"""
Plot Figure 4 from Arunima et al. Paper II: 
Radial dependence of the time evolution of relative action differences

This script assumes you have already generated the HDF5 files
`gal_env_JR_allbins.h5`, `gal_env_Jz_allbins.h5`, and `gal_env_Jphi_allbins.h5`
using the script stellar-actions-II/analysis/gal_env_effect.py

Might need to change paths, marked by the comment #CHANGE PATH HERE in the following code.

Author: Arunima 
Date: 13-11-2025
"""

# ======================================================================
# Imports
# ======================================================================
import numpy as np
import h5py
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

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

#=============================================================================
#file paths #CHANGE PATH HERE

output_file_R   = "stellar-actions-II/data/gal_env_JR_allbins.h5"
output_file_z   = "stellar-actions-II/data/gal_env_Jz_allbins.h5"
output_file_phi = "stellar-actions-II/data/gal_env_Jphi_allbins.h5"


# ======================================================================
# Helper function to read data from HDF5 files
# ======================================================================
def read_relchange_by_bin(filename, dataset="birth_bin_relchange_gt30pc"):
    """Read rel action diff change arrays from each radial birth bin in a saved HDF5 file."""
    data = []
    with h5py.File(filename, "r") as f:
        for key in f.keys():
            data.append(np.array(f[key][dataset]))
    return np.array(data)


# ======================================================================
# Load data
# ======================================================================
print("Loading rel action diff change data by birth-radius bin...")
list_JR   = read_relchange_by_bin(output_file_R)
list_Jz   = read_relchange_by_bin(output_file_z)
list_Jphi = read_relchange_by_bin(output_file_phi)


# ======================================================================
# Plotting
# ======================================================================
plt.close('all')
fig, axs = plt.subplots(3, 1, figsize=(5,11),sharex='all')

#which birth bin are we plotting: stars born within 0.5-2 pc of each other
k = 1 

# JR panel -------------------------------------------------------------
x = np.arange(list_JR.shape[1])
for i in range(2,list_JR.shape[0],2):
    alpha = 0.1 + 0.9 * i / (list_JR.shape[0] - 1)
    axs[0].plot(x, list_JR[i,:,k],color='blue', alpha=alpha, label=f'{i}-{i+1} kpc')
axs[0].plot(delta_t_range, JR_rel_30[:,k],c='r',label='all data')
axs[0].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{R}\rangle$")

# Jz panel -------------------------------------------------------------
x = np.arange(list_Jz.shape[1])
for i in range(2,list_Jz.shape[0],2):
    alpha = 0.1 + 0.9 * i / (list_Jz.shape[0] - 1)
    axs[1].plot(x, list_Jz[i,:,k], color='blue', alpha=alpha, label=f'{i}-{i+1} kpc')
axs[1].plot(delta_t_range, Jz_rel_30[:,k],c='r',label='all data')
axs[1].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{z}\rangle$")

# Jphi panel -----------------------------------------------------------
x = np.arange(list_Jphi.shape[1])
for i in range(2,list_Jphi.shape[0],2):
    alpha = 0.1 + 0.9 * i / (list_Jphi.shape[0] - 1)
    axs[2].plot(x, list_Jphi[i,:,k], color='blue', alpha=alpha, label=f'{i}-{i+1} kpc')
axs[2].plot(delta_t_range, Jphi_rel_30[:,k],c='r',label='all data')
axs[2].set_ylabel(r"$\langle \delta_{\text{rel}}\Delta J_{\phi}\rangle$")

axs[0].legend(ncols=2)
axs[2].set_xlabel(r"$\Delta t$ (Myr)")
axs[0].set_xlim((-5,450))
plt.tight_layout()
#CHANGE PATH HERE
plt.savefig("stellar-actions-II/plots/fig4_radial_dependence.pdf")
plt.show()
