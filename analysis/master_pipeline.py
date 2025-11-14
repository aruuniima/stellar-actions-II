import pandas as pd
import obs_functions as func
import os 
import h5py
import numpy as np

all_actions=True #are we saving the entire distribution of actions (along with age mu and sigma) to a hdf5 file

#define a function to write out the hdf5 file if all_actions==True:
stream_h5_file = '/g/data/jh2/ax8338/project2/obs_applications/streams_actions_pruned.hdf5'
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
        Mean and uncertainty (1σ) of cluster age.
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
        f.create_group('trial')

clusters_file = '/g/data/jh2/ax8338/project2/obs_applications/clusters'
mem_file =  '/g/data/jh2/ax8338/project2/obs_applications/members'

#reading clusters file from the Hunt and Reffert 2024 catalog
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
# cluster_summary=cluster_summary[cluster_summary['radius_pc']>=50] #radius should at least be 50 pc otherwise it's basically all in the same area and can be visually identified

# ##reading members data
# ### members data
# members = pd.read_fwf(
#     mem_file,  # replace with your file path
#     colspecs=[
#         (0, 20),    # cluster name
#         (21, 25),   # internal cluster ID
#         (26, 45),   # Gaia DR3 source ID
#         (61, 73),   # RA
#         (85, 96),   # DEC
#         (108, 120), # Galactic longitude
#         (121, 132), # Galactic latitude
#         (133, 144), # PM RA*cos(dec)
#         (156, 167), # PM Dec
#         (528, 537), # Radial velocity
#         (538, 550), # Radial velocity error
#     ],
#     names=[
#         "cluster_name", "internal_id", "source_id", "ra", "dec",
#         "glon", "glat", "pmra", "pmdec", "radial_velocity", "radial_velocity_error"
#     ]
# )
# # Create a set of valid (cluster_name, internal_id) pairs
# valid_clusters = set(zip(cluster_summary["cluster_name"], cluster_summary["internal_id"]))

# # Filter members that belong to those clusters
# filtered_members = members[
#     members.apply(lambda row: (row["cluster_name"], row["internal_id"]) in valid_clusters, axis=1)
# ]

# #combined dataframe
# full_df = pd.merge(
#     filtered_members,
#     cluster_summary,
#     on=["cluster_name", "internal_id"],
#     how="left"  # keep all member stars
# )

# full_df = pd.read_csv('/g/data/jh2/ax8338/project2/obs_applications/OCT_all_stars_moreRVs.csv')
#full df with outliers removed- radial velocities beyond mean+/-3sigma
full_df = pd.read_csv('/g/data/jh2/ax8338/project2/obs_applications/NOV_all_stars_pruned.csv')


all_clusters_summary = []
for name in cluster_summary['cluster_name'].unique():
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
    
    summary_df.to_csv('/g/data/jh2/ax8338/project2/obs_applications/final_clusters_result_OCT.csv',index=False)
