import numpy as np
import random
from tqdm import tqdm
from scipy.spatial.distance import pdist, squareform
import h5py
import os
import sys


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

which_j =int(sys.argv[1]) #enter which component to process

print('loading data')
#load data # can try to directly load filtered data - filt_C, filt_J, A
A=np.load('/g/data/jh2/ax8338/project2/A.npy')
J=np.load('/g/data/jh2/ax8338/project2/filt_J.npz')
#filtered data saved as JR,Jz,Jphi
C=np.load('/g/data/jh2/ax8338/project2/C.npy')   #in cylindrical coordinates

if which_j==0:
    filt_J = J['JR']
    comp = 'R'
elif which_j==1:
    filt_J = J['Jphi']
    comp ='phi'
elif which_j==2:
    filt_J = J['Jz']
    comp = 'z'
output_file = f'/g/data/jh2/ax8338/project2/save_data/all_data_interp_{comp}.h5'
logfile = f"/g/data/jh2/ax8338/project2/save_data/completed_timesteps_J{comp}.txt"

sep_threshold = 30
birth_window = 1 #Myr
delta_t_range = np.arange(1, 457, 1)

# Reshape to (T*N, 3) to treat all time+star entries as one list
C_flat = C.reshape(-1, 3)
# Apply the vectorised conversion
C_cart_flat = convert_to_cartesian(C_flat)  # shape (T*N, 3)
# Reshape back to (T, N, 3)
C_cart = C_cart_flat.reshape(C.shape)


if os.path.exists(logfile):
    with open(logfile, "r") as lf:
        completed_pairs = set(tuple(map(int, line.split())) for line in lf)
else:
    completed_pairs = set()

#creating the hdf5 file. there are arrays with shape (n_dt,number_of_pairs(they grow))
def init_h5(fname, n_dt):
    with h5py.File(fname, "w") as f:
        # dt axis is fixed, pair axis is unlimited
        maxshape = (n_dt, None)
        kwargs = dict(maxshape=maxshape, chunks=(n_dt, 1000), compression="gzip")

        f.create_dataset("r0", shape=(0,), dtype="f4", maxshape=(None,),
                        chunks=(1000,), compression="gzip", fillvalue=np.nan)
        f.create_dataset("relJ", shape=(n_dt, 0), fillvalue=np.nan, dtype="f4", **kwargs)
        f.create_dataset("absJpn", shape=(n_dt, 0), fillvalue=np.nan, dtype="f4", **kwargs)

# --- append results from one t0 block ---
def append_block(fname, dist_array, relJ_array, absJ_array):
    M = relJ_array.shape[1]  # number of pairs in this block
    if M == 0:
        return
    with h5py.File(fname, "a") as f:
        for name, data in [
            ("relJ", relJ_array),
            ("absJpn", absJ_array),
        ]:
            ds = f[name]
            old_n = ds.shape[1]
            ds.resize(old_n + M, axis=1)
            ds[:, old_n:old_n+M] = data
        #store initial distance r0 separately:
        ds_r0 = f["r0"]
        old_m = ds_r0.shape[0]
        ds_r0.resize(old_m + M, axis=0)
        ds_r0[old_m:old_m+M] = dist_array

#create file if not present
if not os.path.exists(output_file):
    init_h5(output_file,len(delta_t_range))

#open and append
with h5py.File(output_file, "a") as f:
# C_cart.shape[0]
    for t0 in tqdm(range(2), desc="Processing births"):
        # stars born in this snapshot
        born_now = np.where((A[t0, :] >= 0) & (A[t0, :] < birth_window))[0]
        if len(born_now) < 2:
            continue

        coords_now = C_cart[t0, born_now, :]
        dist_matrix = squareform(pdist(coords_now))
        triu_inds = np.triu_indices(len(born_now), k=1)

        global_i1 = born_now[triu_inds[0]]
        global_i2 = born_now[triu_inds[1]]

        #filtering: not considering pairs with birth distances more than 2000 pc
        dist_birth = np.linalg.norm(C_cart[t0, global_i1] - C_cart[t0, global_i2], axis=1)
        birth_mask = dist_birth <= 2000  # only pairs within 2 kpc at birth
        global_i1 = global_i1[birth_mask]
        global_i2 = global_i2[birth_mask]
        dist_birth = dist_birth[birth_mask]
        
        # --- another filter step: only keep pairs that drift >30pc by dt=100 ---
        t_check = t0 + 100
        if t_check >= C_cart.shape[0]:
            print('skipping here')
            continue
        dist_future_check = np.linalg.norm(
            C_cart[t_check, global_i1] - C_cart[t_check, global_i2], axis=1
        )
        keep_mask = dist_future_check > sep_threshold
        if not np.any(keep_mask):
            continue

        global_i1 = global_i1[keep_mask]
        global_i2 = global_i2[keep_mask]
        dist_birth = dist_birth[keep_mask]

        #got the pairs
        M = len(global_i1) #number of pairs
        if M ==0:
            continue
        #allocate nan array
        relJ_block = np.full((len(delta_t_range), M), np.nan, dtype=np.float32)
        absJ_block = np.full((len(delta_t_range), M), np.nan, dtype=np.float32)


        for j,dt in enumerate(tqdm(delta_t_range)):
            #checking if t0 and delta_t already here
            if (t0,dt) in completed_pairs:
                continue
            t1 = t0 + dt
            if t1 >= C_cart.shape[0]:
                break

            dJ_t1 = abs(filt_J[t1, global_i1] - filt_J[t1, global_i2])
            dJ_t0 = abs(filt_J[t0, global_i1] - filt_J[t0,global_i2])
            # modulus of difference b/w the actions of two stars ar t1 and at t2
            # then take difference of those - if negative then difference decreased at a later time
            abs_change = dJ_t1 - dJ_t0
            rel_change = abs_change/(np.sqrt(filt_J[t1, global_i1]*filt_J[t1,global_i2]))

            relJ_block[j,:] = rel_change
            absJ_block[j,:] = abs_change
            
            with open(logfile, "a") as lf:
                    lf.write(f"{t0} {dt}\n")
                
        append_block(output_file,dist_birth,relJ_block,absJ_block)
        print(f"Appended t0={t0}, {M} pairs")

                
     
            

            
            
            




