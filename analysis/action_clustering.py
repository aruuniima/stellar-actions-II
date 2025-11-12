#imports
import numpy as np
from tqdm import tqdm
import h5py
from collections import defaultdict
from scipy.spatial.distance import pdist, squareform
import os

#------------------------------------------------------------------------------------
#function 
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


#===================================================================================
# which action component, 0 is radial, 1 is azimuthal and 2 is vertical action
which_j = 0
# CHANGE PATH HERE
output_file = 'insert_output_file_path_here.h5'
logfile = "insert_log_file_path_here.txt"    #to keep track of which timesteps have been done

sep_threshold=30 # the threshold distance to separate unbound stars - we only want to analyse these. 30 pc comes from the Fig 1 of the paper

birth_bins = np.array((0.15,0.5,2,5,10,50,100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,5000))   #bins of birth separation between two stars - these are the values used in Fig.2,3

birth_window = 1 #Myr
delta_t_range = np.arange(1, 457, 1)
n_bins = len(birth_bins)-1

print('loading data')
#GET DATA FROM  https://www.mso.anu.edu.au/~arunima/stellar-actions-I-data/
#CHANGE PATH HERE if needed
with h5py.File('stellar-actions-I/data/actions.h5', 'r') as f:
    J = f['actions'][:] #filtered data saved as JR,Jz,Jphi
    C = f['coordinates'][:] #in cylindrical coordinates
    A = f['age'][:]    #STELLAR AGES (MYR)


    
#================================================================================
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
    born_now = np.where((A[t0, :] >= 0) & (A[t0, :] < birth_window))[0]   #get stars born in the last [birth_window] Myr
    #if not even one pair found, skip.
    if len(born_now) < 2: continue
    coords_now = C_cart[t0, born_now, :]  # shape (N, 3)
    dist_matrix = squareform(pdist(coords_now))  # (N, N)
    triu_inds = np.triu_indices(len(born_now), k=1)
    dists = dist_matrix[triu_inds] #getting distances between all pairs
    bin_indices = np.digitize(dists, birth_bins) - 1        #putting distances in the birth_bins defined earlier
    for idx, b in enumerate(bin_indices):
        if 0 <= b < len(birth_bins) - 1:
            i1 = born_now[triu_inds[0][idx]]
            i2 = born_now[triu_inds[1][idx]]
            bin_pair_indices[b].append((t0, i1, i2))        #stores snapshot number and the indices of the pair of star that fall in each birth_bin

if which_j==0:
    filt_J = J['JR']
elif which_j==1:
    filt_J = J['Jphi']
elif which_j==2:
    filt_J = J['Jz']
    
##creating hdf5 file
if not os.path.exists(output_file):
    with h5py.File(output_file, "w") as f:
          # For pairs that drifted to >30pc separation at 100 Myr, absolute and relative change in action differences
        f.create_dataset("birth_bin_abschange_gt30pc", shape=(len(delta_t_range), n_bins), dtype='f8', compression="gzip", fillvalue=np.nan)
        f.create_dataset("birth_bin_relchange_gt30pc", shape=(len(delta_t_range), n_bins), dtype='f8', compression="gzip", fillvalue=np.nan)

# Load completed timesteps if any
if os.path.exists(logfile):
    with open(logfile, "r") as f:
        completed = set(map(int, f.read().split()))
else:
    completed = set()

with h5py.File(output_file, "a") as f:

    gt30pc_abs_ds = f["birth_bin_abschange_gt30pc"]
    gt30pc_rel_ds = f["birth_bin_relchange_gt30pc"]

    for i, delta_t in enumerate(tqdm(delta_t_range, desc="Looping over Δt")):
        if i in completed:
            print(f"Skipping Δt={delta_t} (index {i}) — already done.")
            continue  # skip if already done!
        
        print(f"\n>>> Starting Δt={delta_t} (index {i})")


        gt30pc_abs_row=[]
        gt30pc_rel_row=[]
            
        for b in range(n_bins):
            pairs = np.array(bin_pair_indices[b])

            #getting the indices and snapshots of the pair of stars
            t0s, i1s, i2s = pairs[:, 0], pairs[:, 1], pairs[:, 2]
            valid = t0s + delta_t < filt_J.shape[0]  #only choosing where t0s+delta_t is within our simulation time
            t0s = t0s[valid]
            i1s = i1s[valid]
            i2s = i2s[valid]

            J0_1 = filt_J[t0s, i1s] #actiion of star 1 at earlier snapshot
            J0_2 = filt_J[t0s, i2s] #action of star 2 at earlier snapshot
            Jt_1 = filt_J[t0s + delta_t, i1s] #action of star 1 at later snapshot 
            Jt_2 = filt_J[t0s + delta_t, i2s]    #action of star 2 at later snapshot

            deltaJ0 = np.abs(J0_1 - J0_2)    #absolute difference in action of stars at previous snapshot
            deltaJt = np.abs(Jt_1 - Jt_2)    #absolute difference in actions of stars at later snapshot
            abs_change = np.abs(deltaJt - deltaJ0)    #absolute change in action differences
            rel_change = abs_change / np.sqrt(Jt_1 * Jt_2)    #relative change in action difference

            #saving stars that are >30pc apart after 100 Myr (to select unbound stars)
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
                gt30pc_abs_row.append(np.nan)
                gt30pc_rel_row.append(np.nan)
            else:
                gt30pc_abs_row.append(np.nanmedian(abs_change[mask_gt30]))
                gt30pc_rel_row.append(np.nanmedian(rel_change[mask_gt30]))
                print(f'{np.sum(mask_gt30)}- >30pc sample in this bin {b}')


        #saving everything        
        gt30pc_abs_ds[i,:] = gt30pc_abs_row
        gt30pc_rel_ds[i,:] = gt30pc_rel_row
        
        f.flush()

        # Save the completed delta_t index
        with open(logfile, "a") as lf:
            lf.write(f"{i}\n")
