import numpy as np
import random
from tqdm import tqdm
import h5py
from scipy.signal import butter, filtfilt
from scipy.signal import freqz
import matplotlib
from collections import defaultdict
from scipy.spatial.distance import pdist, squareform
import pickle
import os

which_j = 0
output_file = '/g/data/jh2/ax8338/project2/save_data/JR_60pc_allbins.h5'
logfile = "/g/data/jh2/ax8338/project2/save_data/completed_timesteps_JR.txt"

sep_threshold=60
# birth_bins = np.array((1000,1200,1400,1600,1800,2000,5000))
birth_bins = np.array((0.15,0.5,2,5,10,50,100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,5000))
# np.logspace(-0.8,3,8)

birth_window = 1 #Myr
delta_t_range = np.arange(1, 457, 1)
n_bins = len(birth_bins)-1

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

print('loading data')
#load data # can try to directly load filtered data - filt_C, filt_J, A
A=np.load('/g/data/jh2/ax8338/project2/A.npy')
J=np.load('/g/data/jh2/ax8338/project2/filt_J.npz')
#filtered data saved as JR,Jz,Jphi
C=np.load('/g/data/jh2/ax8338/project2/C.npy')   #in cylindrical coordinates

if which_j==0:
    filt_J = J['JR']
elif which_j==1:
    filt_J = J['Jphi']
elif which_j==2:
    filt_J = J['Jz']

# Reshape to (T*N, 3) to treat all time+star entries as one list
C_flat = C.reshape(-1, 3)
# Apply the vectorised conversion
C_cart_flat = convert_to_cartesian(C_flat)  # shape (T*N, 3)
# Reshape back to (T, N, 3)
C_cart = C_cart_flat.reshape(C.shape)

#####Getting pairs in each birth_distance bin
bin_centres = 0.5 * (birth_bins[:-1] + birth_bins[1:])
# Storage
bin_pair_indices = defaultdict(list)  # {bin_id: [(t0, i1, i2), ...]}

# Loop over all snapshots to get all pairwise birth distances
for t0 in tqdm(range(C_cart.shape[0]),desc='looping to get pairwise birth distances in each bin'):
    born_now = np.where((A[t0, :] >= 0) & (A[t0, :] < birth_window))[0]
    if len(born_now) < 2: continue
    if (C[t0,born_now,1]>18000).any():
        print(f'{(C[t0,born_now,1]>18000).sum()} stars outside range')
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

##creating hdf5 file
if not os.path.exists(output_file):
    with h5py.File(output_file, "w") as f:
        f.create_dataset("bin_centres", data=0.5 * (birth_bins[:-1] + birth_bins[1:]))
        f.create_dataset("birth_bin_final_diff", shape=(len(delta_t_range), n_bins),dtype='f8', compression="gzip", fillvalue=np.nan)
        f.create_dataset("birth_bin_relchange", shape=(len(delta_t_range), n_bins),
                         dtype='f8', compression="gzip", fillvalue=np.nan)
        f.create_dataset("birth_bin_abschange", shape=(len(delta_t_range), n_bins),
                         dtype='f8', compression="gzip", fillvalue=np.nan)
        # NEW: For pairs that drifted to >30pc separation at 100 Myr
        f.create_dataset("birth_bin_final_diff_gt30pc", shape=(len(delta_t_range), n_bins), dtype='f8', compression="gzip", fillvalue=np.nan)
        f.create_dataset("birth_bin_abschange_gt30pc", shape=(len(delta_t_range), n_bins), dtype='f8', compression="gzip", fillvalue=np.nan)
        f.create_dataset("birth_bin_relchange_gt30pc", shape=(len(delta_t_range), n_bins), dtype='f8', compression="gzip", fillvalue=np.nan)

# Load completed timesteps if any
if os.path.exists(logfile):
    with open(logfile, "r") as f:
        completed = set(map(int, f.read().split()))
else:
    completed = set()

with h5py.File(output_file, "a") as f:
    final_diff_ds = f["birth_bin_final_diff"]
    abs_change_ds = f["birth_bin_abschange"]
    rel_change_ds = f["birth_bin_relchange"]

    if "birth_bin_final_diff_gt30pc" in f:
        gt30pc_final_ds = f["birth_bin_final_diff_gt30pc"]
        gt30pc_abs_ds = f["birth_bin_abschange_gt30pc"]
        gt30pc_rel_ds = f["birth_bin_relchange_gt30pc"]

    for i, delta_t in enumerate(tqdm(delta_t_range, desc="Looping over Δt")):
        if i in completed:
            print(f"Skipping Δt={delta_t} (index {i}) — already done.")
            continue  # skip if already done!
        
        print(f"\n>>> Starting Δt={delta_t} (index {i})")

            
        final_row = []
        abs_row = []
        rel_row = []
        

        gt30pc_final_row=[]
        gt30pc_abs_row=[]
        gt30pc_rel_row=[]
            
        for b in range(n_bins):
            pairs = np.array(bin_pair_indices[b])
            if len(pairs) == 0:
                final_row.append(np.nan)
                abs_row.append(np.nan)
                rel_row.append(np.nan)
                continue

            t0s, i1s, i2s = pairs[:, 0], pairs[:, 1], pairs[:, 2]
            valid = t0s + delta_t < filt_J.shape[0]
            t0s = t0s[valid]
            i1s = i1s[valid]
            i2s = i2s[valid]

            J0_1 = filt_J[t0s, i1s]
            J0_2 = filt_J[t0s, i2s]
            Jt_1 = filt_J[t0s + delta_t, i1s]
            Jt_2 = filt_J[t0s + delta_t, i2s]

            deltaJ0 = np.abs(J0_1 - J0_2)
            deltaJt = np.abs(Jt_1 - Jt_2)
            abs_change = np.abs(deltaJt - deltaJ0)
            rel_change = abs_change / np.sqrt(Jt_1 * Jt_2)

            print(f'{len(rel_change)}- total sample in this bin {b}')
            #saving stars that are >30pc apart after 100 Myr
            new_valid_mask = (t0s+100)<filt_J.shape[0]
            t0s=t0s[new_valid_mask]
            i1s = i1s[new_valid_mask]
            i2s = i2s[new_valid_mask]
            deltaJt = deltaJt[new_valid_mask]
            abs_change=abs_change[new_valid_mask]
            rel_change=rel_change[new_valid_mask]
    
            pos1 = C_cart[t0s+100, i1s]
            pos2 = C_cart[t0s+100, i2s]
            sep = np.linalg.norm(pos1 - pos2, axis=1)  # just the distance between them at the later time so for later bins, might include everything
            mask_gt30 = sep >sep_threshold
            if np.sum(mask_gt30) ==0:
                gt30pc_final_row.append(np.nan)
                gt30pc_abs_row.append(np.nan)
                gt30pc_rel_row.append(np.nan)
            else:
                gt30pc_final_row.append(np.nanmedian(deltaJt[mask_gt30]))
                gt30pc_abs_row.append(np.nanmedian(abs_change[mask_gt30]))
                gt30pc_rel_row.append(np.nanmedian(rel_change[mask_gt30]))
                print(f'{np.sum(mask_gt30)}- >30pc sample in this bin {b}')


        
        gt30pc_final_ds[i,:] = gt30pc_final_row
        gt30pc_abs_ds[i,:] = gt30pc_abs_row
        gt30pc_rel_ds[i,:] = gt30pc_rel_row
        
        f.flush()

        # Save the completed delta_t index
        with open(logfile, "a") as lf:
            lf.write(f"{i}\n")
