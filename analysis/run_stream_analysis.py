"""
Run posterior reconstruction of initial cluster sizes for all observed streams.

Uses:
    - utils_inference.py (the functions script)
    - KDE files from kde_for_dt.py
    - Observational \Delta J from compute_stream_actions.py

Outputs:
    - CSV file containing summary statistics per stream
    - HDF5 file containing full posterior PDFs per stream

Author: Arunima
Date: 20–11–2025
"""
# =====================================================================
# USER CONFIGURATION BLOCK - CHANGE PATH HERE
# =====================================================================
#input paths
stream_h5_file = 'path_here/streams_actions.hdf5'  #result of compute_stream_actions.py
kde_dir = 'path_here' #the directory where the results of kde_for_dt.py are saved
#output paths
outfile = "path_here/d_init.csv"
h5_file = "path_here/cluster_posteriors.h5"


fieldnames = ["cluster_name", "med_r_init", "high_1s", "high_2s",'l1','med_r_init_bias','1s_bias','2s_bias','l1_bias']
############################################################################################
# IMPORTS
############################################################################################
import numpy as np
import h5py
import os
import csv
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
import scipy.integrate

# === local import ===
import utils_inference as u


# =====================================================================
# WRAPPER FOR PARALLEL AGE PROCESSING
# =====================================================================

# --- Define wrapper for one iteration - both simple and biased  ---
def process_age_dual(kde_dir,age, delJz_obs):
    # Non-biased
    interp, r_grid, dj_grid, post, r_ml, r_median, high_1s, high_2s = u.process_p_post(kde_dir,age, delJz_obs)

    # Biased (truncated)
    interp_b, r_grid_b, dj_grid_b, post_b, r_ml_b, r_median_b, high_1s_b, high_2s_b = u.process_p_post(kde_dir, age, delJz_obs, bias=True)

    # L1 norms
    kde_obs = gaussian_kde(delJz_obs, bw_method='silverman')
    dj_grid_new = np.linspace(dj_grid.min(), dj_grid.max(), 1000)
    l1_norm = u.compute_L1(r_ml, dj_grid_new, kde_obs, interp)
    l1_norm_b = u.compute_L1(r_ml_b, dj_grid_new, kde_obs, interp_b)

    return dict(
        post=post,
        post_b=post_b,
        rest=[r_ml, r_median, high_1s, high_2s],
        rest_b=[r_ml_b, r_median_b, high_1s_b, high_2s_b],
        L1=l1_norm,
        L1_b=l1_norm_b,
        r_grid=r_grid,
        r_grid_b=r_grid_b
    )


# =====================================================================
# MAIN SCRIPT
# =======================================================##############

def main():
    
    #read the clusters' names
    with h5py.File(stream_h5_file,'r') as f:
        clusters_h5=list(f.keys())
    print(f"Found {len(cluster_names)} clusters in stream file.")
    # -----------------------------------------------------
    # Prepare output CSV
    # -----------------------------------------------------
    write_header = not os.path.exists(outfile)

    with h5py.File(h5_file, "a") as h5out:
        with open(outfile, "a", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=fieldnames)
            if write_header:
                writer.writeheader()
            # ============================================================
            # Main loop over clusters
            # ============================================================
            for cluster_name in tqdm(clusters_h5):
                delJz_obs,age_array,mu,sigma = u.read_stream_actions(stream_h5_file,cluster_name)
                delJz_obs[delJz_obs==0]=float(1e-10) #just putting a really small number so doesn't cause an error
                delJz_obs = np.log10(delJz_obs)
                
                # --------------------------------------------
                # Age sampling
                # --------------------------------------------
                ages = np.random.normal(loc=mu,scale=sigma,size=int(1.5e4))
                ages = ages[(ages>0)&(ages<399)]
                ages = np.array([int(age) for age in ages])
                # --- unique ages for safe HDF5 access ---
                unique_ages, inv_idx = np.unique(ages, return_inverse=True)
    
                # --- Parallel execution across unique ages ---
                results = Parallel(n_jobs=8)(
                    delayed(process_age_dual)(kde_dir,age, delJz_obs) for age in tqdm(unique_ages, desc="Computing unique ages")
                )
    
                # constructing full results from the unique age results
                posterior = [results[i]['post'] for i in inv_idx]
                posterior_b = [results[i]['post_b'] for i in inv_idx]
                L1 = [results[i]['L1'] for i in inv_idx]
                L1_b = [results[i]['L1_b'] for i in inv_idx]
    
                r_grid = results[0]['r_grid']
                r_grid_b = results[0]['r_grid_b']
                  
                # --------------------------------------------
                # Average posterior over ages
                # --------------------------------------------
                #combining all results
                #posterior normal
                P_r = np.average(posterior,axis=0)
                P_r/=np.trapezoid(P_r,r_grid) #normalise
                #posterior biased/truncated
                P_r_b = np.average(posterior_b,axis=0)
                P_r_b/=np.trapezoid(P_r_b,r_grid_b) #normalise
                
                # --------------------------------------------
                # L1 norms
                # --------------------------------------------
                L1_final = np.mean(L1)
                L1_bias_final = np.mean(L1_b)
                
                # --------------------------------------------
                # Compute credible intervals for unbiased posterior
                # --------------------------------------------
                #getting CDFs, median and credible intervals
                finalcdf = scipy.integrate.cumulative_trapezoid(P_r, r_grid, initial=0)  # integrates properly
                finalcdf /= finalcdf[-1]  # force normalization
                finalcdf_interp=interp1d(finalcdf,r_grid)
                high_1s = finalcdf_interp(0.84)
                high_2s = finalcdf_interp(0.975)
                med_r_init = finalcdf_interp(0.5)
                
                # --------------------------------------------
                # Compute credible intervals for biased posterior
                # --------------------------------------------
                #truncated/biased version
                finalcdf_bias = scipy.integrate.cumulative_trapezoid(P_r_b, r_grid_b, initial=0)  # integrates properly
                finalcdf_bias /= finalcdf_bias[-1]  # force normalization
                finalcdf_interp_bias=interp1d(finalcdf_bias,r_grid_b)
                high_1s_bias = finalcdf_interp_bias(0.84)
                high_2s_bias = finalcdf_interp_bias(0.975)
                med_r_init_bias = finalcdf_interp_bias(0.5)
                
                 # --------------------------------------------
                # Write CSV entry
                # --------------------------------------------                
                writer.writerow({
                    "cluster_name": cluster_name,
                    "med_r_init": float(med_r_init),
                    "high_1s": float(high_1s),
                    "high_2s": float(high_2s),
                    "l1":float(L1_final),
                    'med_r_init_bias':float(med_r_init_bias),
                    '1s_bias':float(high_1s_bias),
                    '2s_bias':float(high_2s_bias),
                    "l1_bias":float(L1_bias_final)
                })
                # --------------------------------------------
                # Save full posterior to HDF5
                # --------------------------------------------
                if cluster_name in h5f:
                    del h5f[cluster_name]  # overwrite if exists
                grp = h5f.create_group(cluster_name)
                grp.create_dataset("r_eval", data=r_grid)
                grp.create_dataset("P_r", data=P_r)
                grp.attrs["med_r_init"] = float(med_r_init)
                grp.attrs["high_1s"] = float(high_1s)
                grp.attrs["high_2s"] = float(high_2s)
                grp.attrs['l1_norm'] = float(L1_final)
                
                grp.create_dataset("P_r_bias", data=P_r_b)
                grp.attrs["med_r_init_bias"] = float(med_r_init_bias)
                grp.attrs["high_1s_bias"] = float(high_1s_bias)
                grp.attrs["high_2s_bias"] = float(high_2s_bias)
                grp.attrs['l1_norm_bias'] = float(L1_bias_final)


# =====================================================================
# RUN
# =====================================================================
if __name__ == "__main__":
    main()

           
