"""
Plot Figure 8 (or A1), 9, 10, 12 and 13 from Arunima et al paper II
------------------------------------------------

This script visualizes:
- the distribution of inferred median logarithmic initial sizes for all streams. (Fig 8 or with some work A1)
- the distribution of the L1 norm which measures the quality of our fit (Fig 9).
- the cumulative distribution of inferred logaithmic initial sizes (different percentiles, Fig 10)
- Distribution of inferred logarithmic median intitial size as a function of median age of the streams (Fig 12).
- Distribution of the completeness fraction of the streams

Assumes you have d_init.csv file containing summary statistics per stream (from analysis/run_stream_analysis.py) and completeness fraction of the streams from f_comp.py.
OR you can use stellar-actions-II/data/Table_1.csv file
This script is written assuming you are using Table_1.csv. If using the result from run_stream_analysis.py directly, you will have to change the fieldnames as follows (and all sizes are in logarithm in the result so need to change that):

fieldnames = ["cluster_name", "med_r_init", "high_1s", "high_2s",'l1','med_r_init_bias','1s_bias','2s_bias','l1_bias']

Here we have used ['Cluster Name', 'd_init50', 'd_init84', 'd_init97', 'L1_f', 'd_trunc50', 'd_trunc84', 'd_trunc97', 'L1_ftrunc','f_comp'].
Plus you will have to add the age from the all_stars.csv file.


For getting figure A1, you need to use a different d_init.csv file by running the run_stream_analysis.py and compute_stream_actions.py scripts with edited mem_file (stellar-actions-II/data/all_stars.csv) wherein you remove the members with radial velocity more than 3 sigma away from the stream mean radial velocity for each stream.

CHANGE PATH HERE (IN USER CONFIGURATION BLOCK)
Author: Arunima
Date: 21-11-2025
"""

# ======================================================================
# USER CONFIGURATION (EDIT THESE ONLY)
# ======================================================================
dinit_csv_file = "path_here/d_init.csv" #output of run_stream_analysis.py CHANGE PATH HERE

fig_dir = 'path_here/' #CHANGE PATH HERE

####### =========================
# IMPORTS
# =========================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


############################################################################
# plotting control settings
############################################################################
BIGGER_SIZE=14
plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE) 

############################################################################
#loading data and processing
d_init = pd.read_csv(dinit_csv_file)

r = d_init['d_init50'].to_numpy()
r_bias = d_init['d_trunc50'].to_numpy()

r_low = d_init['d_init84'].to_numpy()
r_high = d_init['d_init97'].to_numpy()

l1_norm = d_init['L1_f'].to_numpy()
l1_bias = d_init['L1_ftrunc'].to_numpy()

age = d_init['Median age'].to_numpy()

f_comp = d_init['f_comp'].to_numpy()


############################################################################
# actual plotting - Figure 8 (or A1) - Distribution of inferred median log initial sizes
############################################################################
mask = l1_norm<np.percentile(l1_norm,90)
mask_bias = l1_bias<np.percentile(l1_bias,90)

plt.close('all')
fig, axs = plt.subplots(2, 1, figsize=(5,7),sharex='all')
bins = np.linspace(np.sort(np.log10(r[mask]))[0],3,20)

axs[0].hist(np.log10(r),alpha=0.6,bins=bins,histtype='step',linewidth=1,color='purple',label='all')
axs[0].hist(np.log10(r[mask]),alpha=0.6,bins=bins,histtype='stepfilled',linewidth=1,color='purple',label='masked')
axs[0].set_ylabel('Number of streams')
axs[0].legend()

axs[1].hist(np.log10(r_bias),alpha=0.6,bins=bins,histtype='step',linewidth=1,color='green',label='all')
axs[1].hist(np.log10(r_bias[mask_bias]),alpha=0.6,bins=bins,histtype='stepfilled',linewidth=1,color='green',label='masked')
axs[1].set_ylabel('Number of streams')
axs[1].legend()
axs[1].set_xlabel(r'$a_{50} = \log d_\text{init50}$')
plt.tight_layout()
plt.savefig(fig_dir+ 'fig8.pdf')

############################################################################
# actual plotting - Figure 9 - Distribution of L1 norm
############################################################################
bins=np.linspace(0,0.5,30)
plt.close('all')
plt.figure()
plt.hist(l1_norm,bins=bins,alpha=0.4,color='purple',label=r'$L^1_f$')
plt.hist(l1_bias,bins=bins,alpha=0.4,color='green',label=r'$L^1_{f,\text{trunc}}$')
plt.xlabel(r'$L^1$ norm')
plt.ylabel('Number of streams')
plt.axvline(np.percentile(l1_norm,90),c='k')
plt.axvline(np.percentile(l1_bias,90),c='k',linestyle='--')
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir+'Fig_9.pdf')

############################################################################
# actual plotting - Figure 10 - Cumulative distribution of inferred logaithmic initial sizes (different percentiles)
############################################################################
## cdf of dinit50,dinit84,dinit97 plot
bins = np.linspace(np.sort(np.log10(r[mask]))[0],3,100)

plt.close('all')
plt.figure()
plt.hist(np.log10(r[mask]),weights=np.ones(np.size(r[mask]))/np.size(r[mask]),cumulative=True,bins=bins,histtype='step',label=r'$a_{50}$')
plt.hist(np.log10(r_low[mask]),weights=np.ones(np.size(r[mask]))/np.size(r[mask]),cumulative=True,bins=bins,histtype='step',label=r'$a_{84}$')
plt.hist(np.log10(r_high[mask]),weights=np.ones(np.size(r[mask]))/np.size(r[mask]),cumulative=True,bins=bins,histtype='step',label=r'$a_{97}$')
plt.ylabel('Cumulative Fraction of Streams')
plt.xlabel(r'$a$')
plt.legend()
plt.tight_layout()
plt.savefig(fig_dir+'fig10.pdf')

############################################################################
# actual plotting - Figure 12 - Distribution of inferred logarithmic median intitial size as a function of median age of the streams
############################################################################

# Create the figure
# plt.figure(figsize=(6,5))
plt.close('all')
x=np.log10(r)
y=age
plt.figure()
# Plot filled contours for dense regions
sns.kdeplot(x=x, y=y, fill=True,alpha=0.8, levels=10, cmap='Blues', thresh=0.1)

# Overlay scatter points (optionally only for low-density points)
sns.scatterplot(x=x, y=y, s=10, color='k')

plt.xlim(-0.8,3)
plt.ylim(0,250)
plt.ylabel('Median Age (Myr)')
plt.xlabel(r'$a_{50} = \log d_\text{init50}$')
plt.tight_layout()
plt.savefig(fig_dir+'fig12.pdf')

############################################################################
# actual plotting - Figure 13 - distribution of completeness fraction of the streams
############################################################################

plt.close('all')
plt.figure()
bins=np.linspace(0.4,1,20)
plt.hist(f_comp,alpha=0.3,bins=bins,label='all streams')
plt.hist(f_comp[mask],bins=bins,alpha=0.3,color='fuchsia',label=r'streams with $d_{\text{trunc}50}<10$ pc')
plt.ylabel("Number of streams")
plt.xlabel(r'$f_\text{comp}$ (completeness fraction)')
plt.tight_layout()
plt.legend()
plt.savefig(fig_dir + 'fig13_completeness.pdf')
