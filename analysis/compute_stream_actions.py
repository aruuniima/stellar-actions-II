"""
Compute action differences for observed moving groups and save results to HDF5.

This script:
- Reads the Hunt & Reffert (2024) moving group catalog
- Loads cleaned membership data with radial velocities
- Computes action differences for each group using functions_obs_actions.process_cluster()
- Saves ΔJR, ΔJz, ΔJphi + age metadata to an HDF5 output file
- (Optional) Saves summary percentiles to CSV if all_actions=False

Author: Arunima
Date: 14-11-2025
"""

# ======================================================================
# USER CONFIGURATION (EDIT THESE ONLY)
# ======================================================================

all_actions=True #are we saving the entire distribution of actions (along with age mu and sigma) to a hdf5 file
#False would just save summary percentiles

#define a function to write out the hdf5 file if all_actions==True:
# CHANGE PATH HERE
if all_actions==True:
    stream_h5_file = 'PATH_HERE/streams_actions.hdf5'
if all_actions==False:
    csv_file = 'PATH_HERE/final_clusters_result.csv'
clusters_file = 'stellar-actions-II/data/clusters' # CHANGE PATH HERE

mem_file =  'stellar-actions-II/data//all_stars.csv' # CHANGE PATH HERE

#-----------------------------------------------------------------------------------------------------------------------
#IMPORTS
#-----------------------------------------------------------------------------------------------------------------------

import os 
import h5py
import numpy as np
import pandas as pd
import obs_functions as func


# =============================================================================
# HDF5 Writing Utility
# =============================================================================

def write_stream_actions(filename,cluster_name, delR, delz,delphi, age_mu, age_sigma):
    """
    Save cluster's actions differences and age info to an HDF5 file.

    Parameters
    ----------
    filename : str
        Path to the HDF5 file to create/append.
    cluster_name : str
        Name of the cluster (used as group name).
    delR, delz, delphi : array-like
        Arrays of action/coordinate differences.
    age_mu, age_sigma : float
        Mean and uncertainty (1sigma) of cluster age.
    """
    with h5py.File(filename, "a") as f:  # "a" = append if exists
        # Create or overwrite group for this cluster
        if cluster_name in f:
            del f[cluster_name]   # remove old version
        g = f.create_group(cluster_name)

        # Store arrays
        g.create_dataset("delR", data=np.asarray(delR))
        g.create_dataset("delz", data=np.asarray(delz))
        g.create_dataset("delphi", data=np.asarray(delphi))

        # Store age info as attributes
        g.attrs["age_mu"] = float(age_mu)
        g.attrs["age_sigma"] = float(age_sigma)
    print(f'successfully saved stream ({cluster_name}) data to file')
    
if not os.path.exists(stream_h5_file):
    with h5py.File(stream_h5_file, "w") as f:
        pass


# =============================================================================
# Load Cluster Summary (Hunt & Reffert 2024)
# =============================================================================
cluster_summary = pd.read_fwf(
    clusters_file,  # replace with your file path
    colspecs=[
        (0, 20),    # cluster name
        (21, 25),   # internal cluster ID
        (280, 281), # type of object
        (294, 300), # number of member stars
        (460, 473), # total radius in pc
        (733, 737), # number with radial velocity
        (795, 806), # age 16th percentile
        (807, 818), # median age
        (819, 831), # age 84th percentile
        (584,598), #distance 16th percentile
        (600,614), #dist median
        (616,631), #distance 84th percentile
        (633,637), #number of stars used to get distance
    ],
    names=[
        "cluster_name", "internal_id", "object_type", "n_members",
        "radius_pc", "n_with_rv", "age_16p", "age_median", "age_84p","dist_16p","dist_median",
        "dist_84p", "n_with_dist"
    ]
)
#selecting moving groups which are younger than 457 Myr 
cluster_summary=cluster_summary[cluster_summary['object_type']=='m'] #should be a moving group (not a bound cluster)
cluster_summary['age_median']= 10**cluster_summary['age_median']/10**6
cluster_summary['age_16p']= 10**cluster_summary['age_16p']/10**6
cluster_summary['age_84p']= 10**cluster_summary['age_84p']/10**6
cluster_summary=cluster_summary[cluster_summary['age_median']<=457] #age should be less than 457 Myr because that's the simulation time we have

# =============================================================================
# Load Member Star Data
# =============================================================================

full_df = pd.read_csv(mem_file)
full_df=full_df[full_df['n_with_rv']>10]
full_df=full_df[full_df['age_84p']<456]
common_names = set(cluster_summary['cluster_name'].unique()).intersection(set(full_df['cluster_name'].unique()))
# =============================================================================
# Process All Clusters
# =============================================================================

all_clusters_summary = []
for name in common_names :
    try:
        if all_actions:
            mu,sigma,J = func.process_cluster(full_df,name,all_actions=all_actions)
            delR,delphi,delz = J
            write_stream_actions(stream_h5_file,name, delR, delz,delphi, mu, sigma)
        else:
            result = func.process_cluster(full_df,name)
            all_clusters_summary.append(result)
    except Exception as e:
        print(f"Error processing {name}:{e}")


# =============================================================================
# If requested, save percentile summary table
# =============================================================================

if all_actions==False:
    columns = [
        'cluster_name',
        'jr_16p', 'jr_50p', 'jr_84p',
        'jphi_16p', 'jphi_50p', 'jphi_84p',
        'jz_16p', 'jz_50p', 'jz_84p',
        'age_16p', 'age_50p', 'age_84p',
        'n_members', 'n_with_rv', 'radius_pc'
    ]
    
    summary_df = pd.DataFrame(all_clusters_summary, columns=columns)
    # CHANGE PATH HERE
    summary_df.to_csv(csv_file,index=False)
