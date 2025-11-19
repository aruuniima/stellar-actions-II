"""
Construct conditional PDFs P_1D(\Delta J | r0, \tau) using Gaussian KDEs.

This script reads the interpolated stellar-pair data from the output of the
data_for_model.py pipeline and computes, for each time separation Δt:

    P(log10 \Delta J | log10 r0, \tau)

Following Eq. (13) in Arunima et al

Outputs for each Δt:
    kde_dt_XXX.h5
which contains:
    • r_grid       – log10 r0 grid
    • dj_grid      – log10 \Delta J grid
    • P_cond       – conditional PDF, shape (Nr, NdJ)
    • metadata     – stored in HDF5 attributes

The KDE is computed in log-space for both \Delta J and r0.

change paths in USER SETTINGS BLOCK

Author: Arunima
Date: 19-11-2025
"""
###############################################################################################################
#imports
###############################################################################################################
import os
import h5py
import numpy as np
from tqdm import tqdm
from scipy.stats import gaussian_kde
from multiprocessing import Pool, cpu_count

# ==============================================================================
# USER SETTINGS - CHANGE PATH HERE
# ==============================================================================

interp_file = "INSERT_PATH/all_data_interp_z.h5" #output from data_for_model.py
out_kde_dir = "INSERT_OUTPUT_DIR/kde_by_dt"
os.makedirs(out_kde_dir, exist_ok=True)

n_r = 200    # number of grid points in r
n_j = 100    # number of grid points in \Delta J
r_grid = np.linspace(-1, 3, n_r)   # log10 r0 grid  

# ==============================================================================
# Load interpolated data
# ==============================================================================

print("Loading interpolated data:", interp_file)
with h5py.File(interp_file, "r") as f:
    relJR = f["relJ"][:]   # shape (\Delta t, Npairs)
    r_init = f["r0"][:]    # shape (Npairs,)

# ==============================================================================
# KDE processing for a single \Delta t 
# ==============================================================================


def process_dt(dt):
      """Compute conditional KDE for one snapshot"""
    out_fname = os.path.join(out_kde_dir, f"kde_dt_{int(dt):03d}.h5")
    if os.path.exists(out_fname):
        return (dt, "exists")
    X = np.vstack((r_init,relJR[dt,:]))
    X=X[:,~np.isnan(X[1,:])]  #drop nan
    X=abs(X) #absolute value
    if X.shape[1]>2e6:
      idx = np.random.choice(X.shape[1], size=2000000, replace=False) #random subsampling for speed
      X_fit = X[:,idx]
    else:
      X_fit = X
    #taking log in both delJ 
    X_fit=np.log10(X_fit)

    #fit gaussian kde with silverman bw
    kde = gaussian_kde(X_fit, bw_method='silverman')
    print('kde fit')

    # Construct 2D evaluation grid
    dj_min, dj_max = X_fit[1,:].min(), X_fit[1,:].max()   # use full data range
    r_min, r_max = r_grid.min(), r_grid.max()
    pos_X,pos_Y=np.mgrid[r_min:r_max:n_r*1j,dj_min:dj_max:n_j*1j]
    positions = np.vstack([pos_X.ravel(), pos_Y.ravel()])
    x,y=np.unique(pos_X[:,0]), np.unique(pos_Y[0,:])
    Z = np.reshape(kde(positions).T, pos_X.shape)
    #find integral over vertical slices (all log delJ values) at each log r_init value and divide the kde at each r_init by that to normalise
    #apply a floor to min Z so division doesn't give numerical errors
    eps = 1e-14
    Z+=eps
    Nj = np.trapezoid(Z,x=y,axis=1)
    P_cond = Z/Nj[:,None]
    print('kde normalised')
   
    # Save to HDF5
    with h5py.File(out_fname, "w") as hf:
        grp = hf.create_group("kde")
        grp.create_dataset("r_grid", data=x)
        grp.create_dataset("dj_grid", data=y)
        grp.create_dataset("P_cond", data=P_cond, compression="gzip")
        meta = dict(n_samples=X_fit.shape[1], dt=int(dt))
        for k, v in meta.items():
            grp.attrs[k] = v

    return (dt, "ok", X_fit.shape[1])

# ==============================================================================
# Parallel driver
# ==============================================================================

if __name__ == "__main__":
    dt = np.arange(0,relJR.shape[0])
    # choose number of processes
    nproc = max(1, min(cpu_count()-1, 100))
    print(f"Running with {nproc} processes")
    with Pool(nproc) as P:
        results = list(tqdm(P.imap_unordered(process_dt, dt_list), total=len(dt_list)))
    print("Done.")
