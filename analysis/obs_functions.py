import galpy
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
from galpy.potential import evaluateR2derivs, evaluatePotentials, evaluatez2derivs
from galpy.potential import MWPotential2014, evaluateRforces, evaluatezforces, evaluateRzderivs
from galpy.orbit import Orbit
from galpy.actionAngle import actionAngle
from astropy import units
from galpy.actionAngle import estimateDeltaStaeckel
from galpy.actionAngle import actionAngleStaeckel


# functions

def select_cluster_members(clusters_df, cluster_name):
    """select members of given cluster from big datframe"""
    return clusters_df[clusters_df['cluster_name'] == cluster_name].copy()

##getting gaia bailer jones distances

from astroquery.gaia import Gaia
def get_bj_distances(source_ids):
    ids_str = ','.join(str(int(sid)) for sid in source_ids)
    query = f"""
    SELECT d.source_id, d.r_med_geo
    FROM external.gaiaedr3_distance AS d
    WHERE d.source_id IN ({ids_str})
    """
    job = Gaia.launch_job_async(query)
    res = job.get_results().to_pandas()
    return res

#merging distance data with rest of the cluster data
def merge_dist_dataframe(members_df,bj_df):
    return pd.merge(members_df,bj_df,on='source_id',how='left')

##GETTING ACTIONS AND THEIR DIFFERENCES

#change ra and dec of this catalogue to J2000 from J2016
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

# function to get action change once we have actions from galpy in galpy coordinates
def action_diff(J):
    JR=J*8*220 #kpc km/s
    # Compute outer differences and geometric means
    JR_i = JR[:, None]  # shape (N, 1)
    JR_j = JR[None, :]  # shape (1, N)
    
    # Absolute difference matrix
    delta_J = np.abs(JR_i - JR_j)
    
    # Geometric mean matrix
    gm_J = np.sqrt(JR_i * JR_j)
    
    # Relative difference matrix (with division by zero protection)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_diff_matrix = np.where(gm_J > 0, delta_J / gm_J, np.nan)
    i_inds, j_inds = np.triu_indices(len(JR), k=1)
    rel_diff_upper = rel_diff_matrix[i_inds, j_inds]
    return np.nanpercentile(rel_diff_upper,16),np.nanmedian(rel_diff_upper),np.nanpercentile(rel_diff_upper,84)

# function to get action change once we have actions from galpy in galpy coordinates - gives all pairs' action differences instead of just percentiles
def all_action_diff(J):
    JR=J*8*220 #kpc km/s
    # Compute outer differences and geometric means
    JR_i = JR[:, None]  # shape (N, 1)
    JR_j = JR[None, :]  # shape (1, N)
    
    # Absolute difference matrix
    delta_J = np.abs(JR_i - JR_j)
    
    # Geometric mean matrix
    gm_J = np.sqrt(JR_i * JR_j)
    
    # Relative difference matrix (with division by zero protection)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_diff_matrix = np.where(gm_J > 0, delta_J / gm_J, np.nan)
    i_inds, j_inds = np.triu_indices(len(JR), k=1)
    rel_diff_upper = rel_diff_matrix[i_inds, j_inds]
    return rel_diff_upper
    
#writing a function to get the action change to put on the plot given RA, dec, dist,pmra,pmdec,vlos

def action_difference_pipeline(df,j2016=True,all_actions=False):
    """
    function uses galpy to get actions of all the member stars when a dataframe of the stars is given
    the dataframe should have following columns:
    ra,dec,pmra,pmdec,dist_median,radial_velocity
    if epoch is J2000 instead of J2016, put j2016 as false, it is True by default
    Distance and radial velocity won't be available for a lot of members.
    Should have options? can take all missing rad vel and distances as mean of the ones available
    OR ignore stars that dont have radial velocities, 
    
    uses MWPotential2014 and returns 16th, median and 84th percentiles of relative action change 
    in three components. [[jr16,jr50,jr84],[jphi16,jphi50,jphi84],[jz16,jz50,jz84]]
    """
    ra = np.array(df['ra'])
    dec = np.array(df['dec'])
    dist = np.array(df['dist_median'])
    pmra = np.array(df['pmra'])
    pmdec = np.array(df['pmdec'])
    vlos = np.array(df['best_rv'])
    #removing the ones with no RVs
    finite_mask = np.isfinite(vlos)
    ra=ra[finite_mask]
    dec=dec[finite_mask]
    dist=dist[finite_mask]
    pmra=pmra[finite_mask]
    pmdec=pmdec[finite_mask]
    vlos=vlos[finite_mask]
    if j2016==True:
        # Example: convert your columns into astropy quantities
        coord_2016 = SkyCoord(
            ra=ra*u.deg,
            dec=dec*u.deg,
            distance=dist*u.pc,
            pm_ra_cosdec=pmra*u.mas/u.yr,
            pm_dec=pmdec*u.mas/u.yr,
            radial_velocity=vlos*u.km/u.s,
            obstime=Time('J2016.0')
        )
        
        # Propagate to J2000.0
        coord_2000 = coord_2016.apply_space_motion(new_obstime=Time('J2000.0'))
        
        # Get new coordinates
        ra = coord_2000.ra.deg
        dec = coord_2000.dec.deg
        pmra = coord_2000.pm_ra_cosdec.to_value(u.mas/u.yr)
        pmdec = coord_2000.pm_dec.to_value(u.mas/u.yr)
        

    o = Orbit(np.array((ra,dec,dist,pmra,pmdec,vlos)).T,radec=True,ro=8,vo=220)
    ts = np.linspace(0., 50, 10000) 
    o.integrate(ts, MWPotential2014)
    delta=estimateDeltaStaeckel(MWPotential2014,o.R(ts),o.z(ts))
    aAS= actionAngleStaeckel(pot=MWPotential2014,delta=delta,c=True)
    jr,jphi,jz = aAS(o)
    jr16,jr50,jr84 = action_diff(jr)
    jphi16,jphi50,jphi84 = action_diff(jphi)
    jz16,jz50,jz84 = action_diff(jz)
    if all_actions:
        return all_action_diff(jr), all_action_diff(jphi), all_action_diff(jz)
    return [[jr16,jr50,jr84],[jphi16,jphi50,jphi84],[jz16,jz50,jz84]]

#how can i get a cluster dataframe at the end which will have all actions percentiles, age percentiles, cluster name, number of members used for calculation, total number of members in the hunt and reffert catalog
# full_df['cluster_name','n_members','radius_pc','n_with_rv','age_16p','age_median','age_84p'] + include action percentiles from our calculation

def process_cluster(clusters_df,cluster_name,all_actions=False):
    print('selecting members of cluster')
    members = select_cluster_members(clusters_df,cluster_name)
    print('getting distances from gaia archive bailer-jones catalogue')
    bj_df = get_bj_distances(members['source_id'].tolist())
    print('got distances, merging it with datframe')
    members = merge_dist_dataframe(members,bj_df)
    #putting missing distances as mean of distances/ converting to kpc etc
    members['dist_median'] = members['r_med_geo'][members['r_med_geo'].notna()].mean()/1000 #in kpc
    members.loc[members['r_med_geo'].notna() , 'dist_median'] = members['r_med_geo']/1000
    print('getting action differences')
    [[jr16,jr50,jr84],[jphi16,jphi50,jphi84],[jz16,jz50,jz84]] = action_difference_pipeline(members)
    age16 = np.nanmedian(members['age_16p'])
    age50 = np.nanmedian(members['age_median'])
    age84=np.nanmedian(members['age_84p'])
    n_mem = np.nanmedian(members['n_members'])
    n_used = (~members['best_rv'].isna()).sum()
    size = np.nanmedian(members['radius_pc'])
    if all_actions:
        return age50, (age84-age16)/2, action_difference_pipeline(members,all_actions=True)
    print('returning everything')
    return [cluster_name,jr16,jr50,jr84,jphi16,jphi50,jphi84,jz16,jz50,jz84,age16,age50,age84,n_mem,n_used,size]
    