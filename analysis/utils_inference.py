"""
This file conatins functions that can be used to get the initial sizes of the stellar streams using output from kde_for_dt.py (simulation part) and compute_stream_actions.py (observational part). The functions are used in run_stream_analysis.py.

Author: Arunima
Date: 19-11-2025
"""
####### =========================
# IMPORTS
# =========================
import numpy as np
import h5py
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import norm


# =========================
# KDE LOADING
# =========================
# using the conditional probability saved for each delta_t on the observed cluster's pairs' delta_J
# read the kde file
def load_kde_dt(kde_dir,dt, bias=False,max_J_obs=None):
    fname = f"{kde_dir}/kde_dt_{int(dt):03d}.h5" #CHANGE PATH HERE (IF NEEDED)
    with h5py.File(fname, "r") as f:
        g = f["kde"]
        r_grid = g["r_grid"][:]
        dj_grid = g["dj_grid"][:]
        P_cond = g["P_cond"][:]    # shape (n_r, n_j)
    # build interpolator that returns p(dj | r)
    interp = RegularGridInterpolator((r_grid, dj_grid), P_cond,
                                     bounds_error=False, fill_value=0.0)
    if not bias:
        return r_grid, dj_grid, P_cond, interp
    if bias:
        #truncate the P_cond and renormalise it according to max_J_obs
        if max_J_obs is None:
            raise ValueError("For biased version, max_J_obs must be specified.")

        # Find index where dj_grid <= dj_max
        mask = dj_grid <= max_J_obs
        dj_trunc = dj_grid[mask]
        P_trunc = P_cond[:, mask].copy()
    
        # Renormalize along ΔJ axis for each r
        norm = np.trapezoid(P_trunc, dj_trunc, axis=1)
        norm[norm == 0] = 1.0  # avoid divide-by-zero
        P_trunc /= norm[:, None]
    
        # Build biased interpolator
        interp_biased = RegularGridInterpolator((r_grid, dj_trunc), P_trunc,
                                                bounds_error=False, fill_value=0.0)
    
        return r_grid, dj_trunc, P_trunc, interp_biased


# =========================
# POSTERIOR & LIKELIHOOD
# =========================

#likelihood function
def likelihood_over_r(dj_obs, r_eval=None, interp_cond=None):
    """
    Compute likelihood L(r_init) given observed ΔJ values and a conditional KDE/interpolator.
    
    Parameters
    ----------
    dj_obs : array_like, shape (N_obs,)
        Observed log ΔJ values for a cluster.
    r_eval : array_like, optional
        Grid of log r_init values to evaluate the likelihood on. 
        If None, uses default from the interpolator (must be supplied otherwise).
    interp_cond : callable
        Interpolator returning p(log ΔJ | log r_init). Must accept shape (N, 2) array.
    
    Returns
    -------
    r_eval : ndarray, shape (N_r,)
        Grid of r values.
    logL : ndarray, shape (N_r,)
        Log-likelihood of observing `dj_obs` at each r_init.

    """
    dj_obs = np.atleast_1d(dj_obs)
    
    if r_eval is None:
        raise ValueError("r_eval must be provided if not embedded in interp_cond.")
    r_eval = np.atleast_1d(r_eval)
    
    # Create (r, dj) grid for querying the interpolator
    R, DJ = np.meshgrid(r_eval, dj_obs, indexing="ij")  # shape (N_r, N_obs)
    pts = np.column_stack([R.ravel(), DJ.ravel()])      # shape (N_r*N_obs, 2)
    
    # Evaluate conditional PDF
    pvals = interp_cond(pts)                            # shape (N_r*N_obs,)
    pvals = pvals.reshape(len(r_eval), len(dj_obs))    # shape (N_r, N_obs)
    
    # Avoid log(0)
    pvals = np.maximum(pvals, 1e-300)
    
    # Log-likelihood: sum over all observed ΔJ
    logL = np.sum(np.log(pvals), axis=1)               # shape (N_r,)
    
    return r_eval, logL

#defining a function that gives posterior probability,r_ml, r_med, 1sigma, 2sigma
def process_p_post(kde_dir,age,delJz_obs,bias=False):
    max_J_obs=delJz_obs.max()
    r_grid, dj_grid, P_cond,interp = load_kde_dt(kde_dir,age,bias,max_J_obs)
    r_eval, logL_r=likelihood_over_r(delJz_obs,r_eval=r_grid,interp_cond=interp)
    log_post = logL_r #assuming prior is flat in log r_init
    #getting the uncertainty for now i.e. one delta_t
    post = np.exp(log_post - log_post.max())  # subtract max to avoid overflow
    post /= np.trapezoid(post, r_eval)           # normalize so area = 1
    # cumulative distribution (CDF)
    cdf = scipy.integrate.cumulative_trapezoid(post, r_eval, initial=0)  # integrates properly
    cdf /= cdf[-1]  # force normalization
    cdf_interp=interp1d(cdf,r_eval)
    r_ml = r_eval[np.argmax(log_post)]
    r_median = cdf_interp(0.5)
    high_1s = cdf_interp(0.84)
    high_2s = cdf_interp(0.975)
    return interp,r_grid,dj_grid,post,r_ml,r_median,high_1s,high_2s

def normalise_posterior(L, prior):
    post = L * prior
    Z = np.trapz(post, x=np.arange(len(post)))
    return post / Z if Z > 0 else post

#getting L1 norm
def compute_L1(r_ml, dj_grid_new, kde_obs, interp):
    """
    Computes L1 distance.
    """
    # --- Observed PDF ---
    p_obs = kde_obs(dj_grid_new)
    p_obs = np.maximum(p_obs, 0)
    p_obs /= np.trapezoid(p_obs, dj_grid_new)

    # --- Theoretical PDF (full) ---
    p_theo = interp((r_ml, dj_grid_new))
    p_theo = np.maximum(p_theo, 0)
    p_theo /= np.trapezoid(p_theo, dj_grid_new)

    # --- L1 calculations ---
    L1_full = 0.5 * np.trapezoid(np.abs(p_obs - p_theo), dj_grid_new)

    return L1_full

# =========================
# AGE SAMPLING
# =========================

def sample_ages(mu, sigma, N=1000):
    """Sample stellar ages from Gaussian truncated at [0, 399]."""
    ages = np.random.normal(mu, sigma, N)
    ages = np.clip(ages, 0, 399)
    return ages


# =========================
# STREAM INPUT
# =========================
# function to read the hdf5 file with all the clusters' pairs' action change
def read_stream_actions(file_name,cluster_name,comp=2):
    """
    Read action-change arrays and age info for one cluster
    from an HDF5 file created by compute_stream_actions.py

    Parameters
    ----------
    file_name : str
        Path to the HDF5 file.
    cluster_name : str
        Cluster name (group name in HDF5).
    comp : int
        which component of action. 0 for R, 1 for phi, 2 for z. By default z i.e. 2

    Returns
    -------
    delJobs : array
        arrays of action changes of given component
        
    age_info : array
        array with (mu-4sigma,mu-3sigma,mu-2sigma,mu-sigma,mu,mu+sigma,mu+2sigma,mu+3sigma,mu+4sigma)
    """
    with h5py.File(file_name, "r") as f:
        if cluster_name not in f:
            raise KeyError(f"Cluster '{cluster_name}' not found in {file_name}")
        g = f[cluster_name]
        if comp==0:
            delJobs = g["delR"][:]
        elif comp==2:
            delJobs = g["delz"][:]
        elif comp==1:
            delJobs = g["delphi"][:]

        mu = g.attrs["age_mu"]
        sigma = g.attrs["age_sigma"]
        age_info = np.array((mu-4*sigma,mu-3*sigma,mu-2*sigma,mu-sigma,mu,mu+sigma,mu+2*sigma,mu+3*sigma,mu+4*sigma))

    return delJobs, age_info,mu,sigma
    








    
    
