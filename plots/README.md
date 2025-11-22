# `plots/` Directory

This directory contains all scripts used to generate the figures in  
**_Evolution of action-space coherence in a Milky Way-like simulation_**  
(Arunima et al., Paper II).

Each script corresponds to a specific figure (or group of figures) in the paper.  
All scripts contain a **`USER CONFIGURATION`** block at the top where you must update:

- Paths to input data  
- Paths to output figures  
- Stream names (where applicable)  
- Optional parameter choices  

No script takes command-line arguments.  
Running a script after updating the configuration is simply:
```
python3 fig_X.py
```
---
### `fig_1.py`
**Reproduces Figure 1**:  
Distribution of the 5th nearest-neighbour distance for coeval stars.

**Requires:**
- Simulation actions file: `actions.h5`  
  (download from https://www.mso.anu.edu.au/~arunima/stellar-actions-I-data/)

**User edits:**
- `action_file` — path to `actions.h5`  
- `output_fig_path` — save location for the PDF  
- `ages_to_include` — list of coeval ages (Myr)

---

### `fig_2_3.py`
**Reproduces Figures 2 and 3**:  
Time evolution of action differences for stellar pairs binned by initial separation.

**Requires:**
- HDF5 outputs from `analysis/action_clustering.py`:  
  - `JR.h5`, `Jz.h5`, `Jphi.h5`
- Reference single-star orbital action files (`abs_*.npz`, `rel_*.npz`)  
  (from stellar-actions-I repository)

**User edits:**
- Paths to the three action-difference HDF5 files  (`output_file_R`,`output_file_z`,`output_file_phi`)
- Paths to the six `.npz` reference files  (`abs_JR`,`abs_Jz`,`abs_Jphi`,`rel_JR`,`rel_Jz`,`rel_Jphi`)
- Output directory for saving the figures (`save_dir`)

---

### `fig_4.py`

**Reproduces Figure 4**:  
Radial dependence of action-space divergence.

**Requires:**
- HDF5 outputs from `analysis/gal_env_effect.py`:  
  - paths to be put as variables : `output_file_R`, `output_file_z`, `output_file_phi`
- HDF5 outputs from`analysis/action_clustering.py` for reference curves:
  - paths to be put as variables: `all_bins_R`, `all_bins_z`, `all_bins_phi`

**User edits:**
- Paths for both sets of HDF5 inputs  
- Output directory (`save_dir`)

---

### `fig5.py`

**Reproduces Figure 5**:  
Observed stream action-difference distributions overlaid on theoretical evolution grids.

**Requires:**
- Same HDF5 outputs from `analysis/action_clustering.py`
- `final_clusters_result.csv` (output from `compute_stream_actions.py` with `all_actions=False` or directly take from data/)

**User edits:**
- Paths to HDF5 outputs: `all_bins_R`, `all_bins_z`, `all_bins_phi`
- Path to `final_clusters_result.csv` : `clusters_result` 
- Output directory (`save_dir`)

---

### `fig6.py`
**Reproduces Figure 6**:  
Cumulative distributions of the 1D conditional PDF 

$\mathcal{P}_{1D}(x \mid a, \tau)$

for a chosen stream at a specific age.

**Requires:**
- `stream_actions.h5` (output of `compute_stream_actions.py`)  
- KDE files from `analysis/kde_for_dt.py` (`kde_dt_XXX.h5`)

**Default example:** HSC_2282 at τ = 148 Myr.

**User edits:**
- Paths to `stream_actions.h5` and `kde_dir`  : `stream_h5_file` and `kde_dir`
- `name` — stream name to plot  
- Save path for the figure

---

### `fig7.py`
**Reproduces Figure 7**:  
Posterior cumulative distribution function of the initial size  

**Requires:**
- `cluster_posteriors.h5` from `analysis/run_stream_analysis.py`

**Default example:** HSC_2282.

**User edits:**
- Path to `cluster_posteriors.h5` : `h5_file`
- Stream name (`name`)  
- Output figure path (`save_fig`)

---

### `fig8-13.py`
**Reproduces Figures 8, 9, 10, 12, 13** (and A1 with modifications)

Includes plots of:
- Distribution of inferred median initial log-sizes (Fig 8 / A1)  
- Distribution of L1 norm (Fig 9)  
- CDFs of inferred initial sizes (Fig 10)  
- Initial size vs stream age (Fig 12)  
- Completeness fractions (Fig 13)

**Requires:**
- Summary statistics (`d_init.csv` from `run_stream_analysis.py`) [see comments in the script to see the changes needed if using this]
  
  OR  
- `Table_1.csv` from the repository (default expected format)

**Optional for A1:**
- Use a modified version of `all_stars.csv` with 3 $\sigma$ radial-velocity outliers removed and rerun run_stream_analysis.py using this.

**User edits:**
- Path to `d_init.csv` : `dinit_csv_file` 
- Output directory (`fig_dir`)

---

### `fig11.py`
**Reproduces Figure 11**:  
Histogram of observed $\Delta J_z$ (vertical action differences)  
for two representative streams:
- “small” initial size (e.g. Theia 323)  
- “large” initial size (e.g. HSC 1964)

**Requires:**
- `stream_actions.h5` from `compute_stream_actions.py`

**User edits:**
- Path to `stream_actions.h5`  : `stream_h5_file`
- Two stream names (`big`, `small`)  
- Their corresponding inferred sizes  
- Output figure path (`save_fig`)

## Notes on Running the Figure Scripts

- All scripts are **standalone**: run each individually after updating paths.  
- Figures are saved as `.pdf` by default but can be changed in the config block.  
- The ordering of the scripts corresponds to the figure numbering in the paper.

