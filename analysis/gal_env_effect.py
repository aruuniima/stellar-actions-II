"""
Compute action-space evolution for coeval stellar pairs
-------------------------------------------------------

This script calculates the absolute and relative changes in orbital actions (J_R, J_φ, J_z)
for pairs of stars born within given spatial (birth) separation bins, separated by >30 pc
after 100 Myr, and within the same 1 kpc Galactic radial bin.

It saves the median action difference relative and absolute change values over time (Δt = 1–456 Myr) into an HDF5 file, with datasets per radial bin. The results can be used to generate Figures 4 from the paper.

paths to be changed marked with comment #CHANGE PATH HERE and are in the USER CONFIGURATION BLOCK

Author: Arunima
Date: 13-11-2025
Data source: https://www.mso.anu.edu.au/~arunima/stellar-actions-I-data/
"""
# ======================================================================
# USER CONFIGURATION - edit these/ CHANGE PATH HERE
# ======================================================================

which_j = 1 # Select which action to analyze: 0=J_R, 1=J_PHI, 2=J_z

output_dir = "/PATH/TO_SAVE/"              # CHANGE PATH HERE
log_dir    = "/PATH/TO_SAVE/logs/"         # CHANGE PATH HERE
action_file = "stellar-actions-I/data/actions.h5"  # CHANGE PATH HERE if needed
#GET DATA FROM  https://www.mso.anu.edu.au/~arunima/stellar-actions-I-data/

# --- Output HDF5 file ---
output_file = output_dir + "gal_env_Jphi_allbins.h5"  # auto-selected name

birth_bins = np.array((0.15,0.5,2,5,10,50,100,200))
# --- Window for identifying co-born pairs ---
birth_window = 1 #Myr
delta_t_range = np.arange(1, 457, 1)
n_bins = len(birth_bins)-1
#birth bin 1 - 0.1-0.5 pc i=0
#2- 2-5 pc i=2
#3- 10-50 pc i=4
#4 - 100-200 pc i=6
# ======================================================================
# END USER CONFIGURATION
# ======================================================================

# ======================================================================
# Imports
# ======================================================================
import numpy as np
from tqdm import tqdm
import h5py
from collections import defaultdict
from scipy.spatial.distance import pdist, squareform
import pickle
import os

# ======================================================================
# Ensure output directories exist
# ======================================================================
os.makedirs(output_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

# ======================================================================
 # Helper Functions
# ======================================================================
def convert_to_cartesian(phi_R_z):
    """
    Convert (phi, R, z) to Cartesian coordinates (x, y, z).
    
    Parameters:
    - phi_R_z: ndarray, shape (N, 3), array of (phi, R, z) values for N stars
    
    Returns:
    - cartesian_coords: ndarray, shape (N, 3), array of (x, y, z) coordinates
    """
    phi = phi_R_z[:, 0]
    R = phi_R_z[:, 1]
    z = phi_R_z[:, 2]
    
    x = R * np.cos(phi)
    y = R * np.sin(phi)
    
    return np.column_stack((x, y, z))



# ======================================================================
# Load data
# ======================================================================

print('loading data')

with h5py.File(action_file, 'r') as f:
    J = f['actions'][:] #filtered data saved as JR,Jz,Jphi
    C = f['coordinates'][:] #in cylindrical coordinates
    A = f['age'][:]    #STELLAR AGES (MYR)
    ID = f['ID'][:]   #star IDs


#select which action to process
if which_j==0:
  filt_J = J['JR']
  j_label = "JR"  
elif which_j==1:
  filt_J = J['Jphi']
  j_label = "Jphi"
elif which_j==2:
  filt_J = J['Jz']
  j_label='Jz'
else:
    raise ValueError("Invalid `which_j` value. Must be 0, 1, or 2.")

print(f"Processing action: {j_label}")

#=====================================================================================
#defining radial bins (1kpc wide)
r_bins = np.arange(0,19000,1000)

IDs_in_bins = {i: [] for i in range(len(r_bins))}  # Dictionary to store IDs by bin
for t in tqdm(np.arange(A.shape[0]), desc="Looping over time snapshots"):
    # Identify stars born in the last 1 Myr
    valid = ~np.isnan(A[t, :])  # Valid stars
    recently_born = (A[t, valid] < 1) & (A[t, valid] >= 0)
    # Get the indices of these stars
    star_indices = np.where(valid)[0][recently_born]
    R_indices = np.digitize(C[t, valid, 1][recently_born],r_bins)-1
    # Assign star IDs to their respective bins
    for i, r_bin in enumerate(R_indices):
        if 0 <= r_bin < len(r_bins) - 1:  # Ensure index is within valid bin range
            IDs_in_bins[r_bin].append(star_indices[i])


# -----------------------------------------------------------------------------
# Main computation loop
# -----------------------------------------------------------------------------
for l_R in tqdm(np.arange(2, 19, 1), desc="Processing radial bins"):
    logfile = f"{log_dir}/completed_{j_label}_Rbin{l_R}.txt"
    J_all = filt_J[:,IDs_in_bins[l_R]]
    A_all = A[:, IDs_in_bins[l_R]]
    C_all = C[:,IDs_in_bins[l_R]]
    num_stars_total = A_all.shape[1]

    print(f'Total number of stars in bin {l_R} = {num_stars_total}')
    # Reshape to (T*N, 3) to treat all time+star entries as one list
    C_flat = C_all.reshape(-1, 3)
    # Apply the vectorised conversion
    C_cart_flat = convert_to_cartesian(C_flat)  # shape (T*N, 3)
    # Reshape back to (T, N, 3)
    C_cart = C_cart_flat.reshape(C_all.shape)

        
    #####Getting pairs in each birth_distance bin
    bin_centres = 0.5 * (birth_bins[:-1] + birth_bins[1:])
    # Storage
    bin_pair_indices = defaultdict(list)  # {bin_id: [(t0, i1, i2), ...]}
    #loop over all snapshots to get all pairwise distances
    for t0 in tqdm(range(C_cart.shape[0]),desc='looping to get pairwise birth distances in each bin'):
        born_now = np.where((A_all[t0, :] >= 0) & (A_all[t0, :] < birth_window))[0]
        if len(born_now) < 2: continue
        coords_now = C_cart[t0, born_now, :]  # shape (N, 3)
        dist_matrix = squareform(pdist(coords_now))  # (N, N)
        triu_inds = np.triu_indices(len(born_now), k=1)
        dists = dist_matrix[triu_inds]
        bin_indices = np.digitize(dists, birth_bins) - 1
        for idx, b in enumerate(bin_indices):
            if 0 <= b < len(birth_bins) - 1:
                i1 = born_now[triu_inds[0][idx]]
                i2 = born_now[triu_inds[1][idx]]
                bin_pair_indices[b].append((t0, i1, i2))
        
            # Initialize HDF5 output
    if not os.path.exists(output_file):
        with h5py.File(output_file, "w") as f:
            for r in range(len(r_bins)):
                group = f.create_group(f"bin_{r}")

            
                group.create_dataset(
                    "birth_bin_abschange_gt30pc",
                    shape=(len(delta_t_range),n_bins),
                    dtype='f8',
                    compression="gzip",
                    fillvalue=np.nan
                )
                group.create_dataset(
                    "birth_bin_relchange_gt30pc",
                    shape=(len(delta_t_range),n_bins),
                    dtype='f8',
                    compression="gzip",
                    fillvalue=np.nan
                )

    # Load completed timesteps if any
    if os.path.exists(logfile):
        with open(logfile, "r") as f:
            completed = set(map(int, f.read().split()))
    else:
        completed = set()

    with h5py.File(output_file, "a") as f:
        gt30pc_abs_ds = f[f"bin_{l_R}"]["birth_bin_abschange_gt30pc"]
        gt30pc_rel_ds = f[f"bin_{l_R}"]["birth_bin_relchange_gt30pc"]

        for i, delta_t in enumerate(tqdm(delta_t_range, desc="Looping over Δt")):
            if i in completed:
                print(f"Skipping Δt={delta_t} (index {i}) — already done.")
                continue  # skip if already done!
            
            print(f"\n>>> Starting Δt={delta_t} (index {i})")

            gt30pc_abs_row=[]
            gt30pc_rel_row=[]

            for b in range(n_bins):
                pairs = np.array(bin_pair_indices[b])
                if len(pairs) == 0:
                    print(f'no pairs found in birth bin {b}')
                    gt30pc_abs_row.append(np.nan)
                    gt30pc_rel_row.append(np.nan)
                    continue
            
                t0s, i1s, i2s = pairs[:, 0], pairs[:, 1], pairs[:, 2]
                valid = t0s + delta_t < filt_J.shape[0]
                t0s = t0s[valid]
                i1s = i1s[valid]
                i2s = i2s[valid]
            
                J0_1 = J_all[t0s, i1s]
                J0_2 = J_all[t0s, i2s]
                Jt_1 = J_all[t0s + delta_t, i1s]
                Jt_2 = J_all[t0s + delta_t, i2s]
            
                deltaJ0 = np.abs(J0_1 - J0_2)
                deltaJt = np.abs(Jt_1 - Jt_2)
                abs_change = np.abs(deltaJt - deltaJ0)
                rel_change = abs_change / np.sqrt(Jt_1 * Jt_2)
            
                # Filter for >30 pc separation after 100 Myr
                new_valid_mask = (t0s + 100) < filt_J.shape[0]
                t0s = t0s[new_valid_mask]
                i1s = i1s[new_valid_mask]
                i2s = i2s[new_valid_mask]
                deltaJt = deltaJt[new_valid_mask]
                abs_change = abs_change[new_valid_mask]
                rel_change = rel_change[new_valid_mask]
            
                pos1 = C_cart[t0s + 100, i1s]
                pos2 = C_cart[t0s + 100, i2s]
                sep = np.linalg.norm(pos1 - pos2, axis=1)
                mask_gt30 = sep > 30
            
                if np.sum(mask_gt30) == 0:
                    gt30pc_abs_row.append(np.nan)
                    gt30pc_rel_row.append(np.nan)
                else:
                    gt30pc_abs_row.append(np.nanmedian(abs_change[mask_gt30]))
                    gt30pc_rel_row.append(np.nanmedian(rel_change[mask_gt30]))
                    print(f'{np.sum(mask_gt30)} >30pc sample in birth bin {b}')

            
            gt30pc_abs_ds[i,:] = gt30pc_abs_row
            gt30pc_rel_ds[i,:] = gt30pc_rel_row
            
            f.flush()

            # Save the completed delta_t index
            with open(logfile, "a") as lf:
                lf.write(f"{i}\n")
