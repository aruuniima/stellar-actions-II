"""
 Figure 1 – 5th nearest-neighbour distance distribution for coeval stars
 can get action.h5 file from the following website and change the file paths here accordingly ('stellar-actions-I/data/actions.h5' in following code)
 https://www.mso.anu.edu.au/~arunima/stellar-actions-I-data/
Change paths etc at the user configuration block.

Author: Arunima
Date: 20–11–2025
"""

# ======================================================================
# User CONFIGURATION - CHANGE PATH HERE
# ======================================================================
# ===================================================================================
action_file = 'stellar-actions-I/data/actions.h5'
output_fig_path = 'PATH_TO/fig1.pdf'
ages_to_include = [2,5, 10, 100,200,300, 400] #Myr

####### =========================
# IMPORTS
# =========================
import matplotlib.pyplot as plt
import numpy as np
import h5py
from sklearn.neighbors import KDTree
import seaborn as sns
from tqdm import tqdm

# ======================================================================
# HELPER FUNCTIONS
# ======================================================================
# cylindrical to Cartesian coordinates
def cylind_to_cart(c_array):
      """Convert cylindrical (phi,R, z) to Cartesian (x, y, z)."""
    x = c_array[:,1] * np.cos(c_array[:,0])
    y = c_array[:,1] * np.sin(c_array[:,0])
    z = c_array[:,2]
    return np.array((x, y, z)).T

# ------------------------------------------------------------
# Plotting style
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 18

plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=15)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# -------------------------------------------------------------------------------------------------------------------------
# Read HDF5
with h5py.File(action_file, 'r') as f:
    C = f['coordinates'][:] # Cylindrical coordinates (phi,R,z)
    A = f['age'][:]  # Stellar ages (Myr)

# ------------------------------------------------------------
# getting 5th nearest neighbour distance for coeval stars at different ages
age_distance_dict = {t: [] for t in ages_to_include}  # Dictionary to store distances per age group

# Loop over all snapshots
for snap in tqdm(range(A.shape[0]), desc='Snapshots'):  
    for t in ages_to_include:
        mask = (A[snap, :] < t + 2) & (A[snap, :] > t)  # Select stars of age ~ t at this snapshot
        if np.any(mask):
            C_cart = cylind_to_cart(C[snap, mask])
            kdt = KDTree(C_cart, metric='euclidean')  # KDTree in current snapshot
            dist = kdt.query(C_cart, k=5)[0][:, 4]  # 5th nearest neighbour distance
            age_distance_dict[t].append(dist)

# ------------------------------------------------------------
# Plot KDE for each age group
plt.figure(figsize=(8, 6))
for t in ages_to_include:
    if age_distance_dict[t]:  # Check if there is data
        all_distances = np.concatenate(age_distance_dict[t])  # Combine across all snapshots
        print(f'{t} Myr age- {all_distances.shape} stars')
        sns.histplot(np.log10(all_distances), bins='auto', kde=True, 
                 label=f'{t} Myr', stat="density", element="step",alpha=0.1)

plt.xlabel(r"$\log(d_5/\mathrm{pc})$")
plt.ylabel(r"$dp/d\log d_5$")
plt.legend()
plt.xscale('linear')
plt.axvline(x=1.5,c='r',linestyle='--',alpha=0.8)
plt.tight_layout()
plt.minorticks_on()
plt.savefig(output_fig_path)
plt.show()




