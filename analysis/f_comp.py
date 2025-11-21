"""
Getting completeness fraction by comparing the actions difference distribution from the simulation with the best fit with the assumption that observations are incomplete to that from the observations of the stream members.
See section 4.4 of the paper for details.

Assumes you have stream actions (from analysis/compute_stream_actions.py) and kde_dt_XXX.h5 result from analysis/kde_for_dt.py. 

Assumes you have d_init.csv file containing summary statistics per stream (from analysis/run_stream_analysis.py).
OR you can use stellar-actions-II/data/Table_1.csv file
This script is written assuming you are using Table_1.csv. If using the result from run_stream_analysis.py directly, you will have to change the fieldnames as follows (and all sizes are in logarithm in the result so need to change that):

fieldnames = ["cluster_name", "med_r_init", "high_1s", "high_2s",'l1','med_r_init_bias','1s_bias','2s_bias','l1_bias']

Here we have used ['Cluster Name', 'd_init50', 'd_init84', 'd_init97', 'L1_f', 'd_trunc50', 'd_trunc84', 'd_trunc97', 'L1_ftrunc'].

Output: 
- Gives a csv file with cluster names and completeness fraction.

Change paths in the user configuration section below.

Author: Arunima
Date: 21-11-2025
"""

# ======================================================================
# USER CONFIGURATION (EDIT THESE ONLY)
# ======================================================================
stream_h5_file = 'path_here/stream_actions.h5' #output of compute_stream_actions.py CHANGE PATH HERE
kde_dir = 'path_here' #CHANGE PATH HERE 

dinit_csv_file = "path_here/d_init.csv" #output of run_stream_analysis.py CHANGE PATH HERE
f_comp_file = 'path_here/completeness_fraction.csv' #output file to save completeness fractions CHANGE PATH HERE

####### =========================
# IMPORTS
# =========================
import numpy as np
import h5py
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import make_interp_spline

#local imports
from ..analysis.utils_inference import read_stream_actions, load_kde_dt

##############################################################################
def cdfP_ifr_init(r_init,interp,dj_grid):
    A=np.cumsum(interp((np.log10(r_init),dj_grid)))
    return A/A[-1]

# ======================================================================
# LOAD INPUT TABLE
# ======================================================================
d_init = pd.read_csv(dinit_csv_file)
r_bias = d_init['d_trunc50'].to_numpy()
clusters = d_init['Cluster Name'].to_numpy()

#getting completeness factor
f_comp =[]
for cluster in tqdm(clusters):
    r_fit = r_bias[clusters==cluster][0]
    delJz_obs,age_array,mu,sigma = read_stream_actions(stream_h5_file,cluster)
    r_grid, dj_grid, P_cond, interp = load_kde_dt(kde_dir,int(mu))
      # compute the CDF curve once, then spline it
    completeness=make_interp_spline(dj_grid, cdfP_ifr_init(r_fit,interp,dj_grid),k=3)
      # completeness = CDF at max observed delJz
    f_comp.append(completeness(np.log10(delJz_obs).max()))


# ======================================================================
# SAVE OUTPUT
# ======================================================================
df_out = pd.DataFrame({
    'Cluster Name': clusters,
    'f_comp': f_comp
})

df_out.to_csv(f_comp_file, index=False)
