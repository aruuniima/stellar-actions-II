# `analysis/` directory

This directory contains all scripts used to perform the core analysis for the paper _Evolution of action-space coherence in a Milky Way-like simulation_. 

Scripts are documented in the order they appear in the full workflow.
All scripts require user edits before running. Editable sections are clearly marked inside a `USER CONFIGURATION` block with  `# CHANGE PATH HERE`. These include paths to input/output files, dataset locations, and global variables.

None of the scripts take command-line arguments. After editing the configuration block, each script can be run with:

```
python3 script_name.py
```

If used, the variable `which_j` selects the action component:  
- 0 → $J_R$  
- $1 → J_\Phi$  
- 2 → $J_z$

---

### 1. `action_clustering.py`:
Computes the evolution of **absolute and relative action differences** for stellar pairs born within spatial separation bins.  
Used for **Figures 2–3** of the paper.

Requires:
- Simulation data (actions, coordinates, ages) available at [the website](https://www.mso.anu.edu.au/~arunima/stellar-actions-I-data/).
- Editable parameters: `which_j`, output HDF5 path, logfile path, birth-separation bins, unbound-pair separation threshold.

Outputs:
- HDF5 file with median change in action differences (absolute & relative) vs. $\Delta t$ for each birth bin.
- Logfile tracking processed snapshots.

---

### 2. `gal_env_effect.py`: 
Computes action-difference evolution for pairs of stars **within 1 kpc Galactic radial bins**, for pairs:
- Born within specified spatial bins,
- Separated by >30 pc after 100 Myr.

Produces datasets used for **Figure 4**.

Inputs/outputs mirror `action_clustering.py` (same simulation data) except that the output HDF5 file has datasets per radial bin.
User-editable settings: `which_j`, radial bins, separation threshold, output HDF5 path, logfile.

---

### 3. `functions_obs_actions.py`:
Utilities used by observational-analysis scripts, including:
- Selecting cluster/stream members,
- Fetching Gaia distances (Bailer-Jones),
- Converting coordinates to J2000,
- Computing orbital actions (using `galpy`),
- Computing pairwise $\Delta J$ for stellar groups.

Used mainly by `compute_stream_actions.py`.

---

### 4. `compute_stream_actions.py`:
Computes $\Delta J$ for pairs of member stars in **observed moving groups/streams**.

Inputs (all included in the repo `data/` directory):
- `cluster.csv` (Hunt & Reffert catalog)
- `all_stars.csv` (member stars + RVs from GALAH/RAVE/APOGEE)

Functionality:
- Selects streams,
- Computes actions for all members,
- Computes pairwise $\Delta J$ distributions.

Outputs:
- If `all_actions=True`: full $\Delta J$ distributions per stream in an HDF5 file (`stream_h5_file`)
- If `all_actions=False`: summary percentiles saved to a CSV (`csv_file`)

All file paths can be edited in the configuration block.

---

### 5. `data_for_model.py`:
Processes **all stellar pairs** in the simulation to compute absolute & relative action difference change for use in the inference model used for finding the intiial sizes of real observed moving groups/stellar streams.

Notes:
- Very large output (~200 GB HDF5 file). Not included in the repository.
- Inputs: same simulation data as `action_clustering.py` / `gal_env_effect.py`.

User edits (in user configuration block):
- `DATA_PATH` - path to simulation data
- `which_j`
- Output HDF5 file path (`output_file`)
- Logfile path (`logfile`)

Used as input to `kde_for_dt.py`.

---

### 6. `kde_for_dt.py`:
Construct conditional PDFs $\mathcal{P}_{1D}(x \mid a , \tau)$ for each $\tau = \Delta t$ following equation 13 in Arunima et al paper II. 

Requires:
- Output HDF5 from `data_for_model.py` (`interp_file`)

Produces:
- One HDF5 KDE file per Δt (named `kde_dt_XXX.h5`)
- Output directory controlled via `out_kde_dir`

User-editable: number of grid points in initial separation and action difference used to calculate this probability (see paper for more details).

---

### 7. `utils_inference.py`:
Utility functions used for reconstructing **initial cluster sizes** using:
- KDE results from simulations (`kde_for_dt.py`)
- Observational pairwise $\Delta J$ from `compute_stream_actions.py`

These functions are called by `run_stream_analysis.py`.

---

### 8. `run_stream_analysis.py`:
Runs the full inference pipeline to recover initial sizes of observed streams.    

Uses:
- `utils_inference.py` 
- KDE files from `kde_for_dt.py` (path in `kde_dir`)
- Observational $\Delta J$ from `compute_stream_actions.py` (path in `stream_h5_file`)

Outputs:
- CSV file with median/percentile initial sizes per stream (path in `outfile`)
- HDF5 file with full posterior PDFs (path in `h5_file`)
  All the path variables mentioned here can be change in the user configuration block.

---

### 9. `f_comp.py`:
Computes the **completeness fraction** by comparing simulation-based non-truncated action change difference distributions found using the best-fit size ($d_\text{trunc50}$, calculated using the truncated distribution/with observational incompleteness) to the observed $\Delta J$ distributions. See section 4.4 of the paper for details.

Requires:
- Stream actions (`compute_stream_actions.py`) Path as `stream_h5_file`.
- KDE files (`kde_for_dt.py`) Path in `kde_dir`
- Initial-size catalog (`run_stream_analysis.py` or `Table_1.csv`) Path in `dinit_csv_file`

If using `run_stream_analysis.py` outputs directly, field names must be adjusted (and all sizes are in logarithm in the result so need to change that):
```
fieldnames = ["cluster_name", "med_r_init", "high_1s", "high_2s",'l1','med_r_init_bias','1s_bias','2s_bias','l1_bias']
```
This script assumes the Table 1 default:
```
['Cluster Name', 'd_init50', 'd_init84', 'd_init97', 'L1_f', 'd_trunc50', 'd_trunc84', 'd_trunc97', 'L1_ftrunc'].
```

Output: 
- CSV file with completeness fractions (`f_comp_file`)
  
All paths are set in the user configuration block as usual.
