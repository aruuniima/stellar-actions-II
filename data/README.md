# `data/` Directory

This directory contains all observational data products used in the paper  
**_Evolution of action-space coherence in a Milky Way-like simulation_**,  
as well as metadata required by the analysis scripts in `analysis/`.

The data here falls into two categories:

1. **Local observational catalogs and Table 1 from the paper** (included in this repository)  
2. **External simulation datasets** (too large to store here; downloadable separately)

---

## 1. Local Files Included in This Directory

### `cluster.csv`
Source: Hunt & Reffert catalog of moving groups and stellar streams.

See the Hunt and Reffert data readme for the entire information about this file. Here I mention columns I am using in `compute_stream_actions.py` script:
- Cluster Name
- Type of object (cluster/moving group)
- Number of member stars
- Total radius in pc
- Age (16th, 50th and 84th percentiles)
- Distance median 

---

### `all_stars.csv`
Compiled catalog of stream member stars. Uses the members.csv file provided by Hunt and Reffert and adds radial velocities from GALAH, RAVE and APOGEE.

Includes:
- Gaia astrometry  
- Radial velocities (combined from **GALAH**, **RAVE**, **APOGEE**)  
- Membership assignments from Hunt & Reffert (cluster name)

---

### `Table_1.csv`
Table of inferred **initial cluster sizes** from Paper II  
(Arunima et al., stellar actions II).

Contains:
- Cluster/stream names
- Median age (Myr)
- Present size (pc)
- 50th, 84th and 97th percentiles of initial size (pc)
- L1 norm
- 50th, 84th and 97th percentiles of initial size (pc) assuming truncation of the probability  based on observational incompleteness
- L1 norm for this calculation
- Completeness fraction

Used in:
- `analysis/f_comp.py` (default input for completeness fraction)  
- Can substitute for results of `run_stream_analysis.py` if desired

Field names match those required by the script.

---

## 2. External Simulation Data (Not Included Here)

A large set of simulation outputs — orbital actions, coordinates, ages of all stars —  
is required for the analysis scripts that produce Figures 2–4 and the inference model.

These datasets are available at:

**https://www.mso.anu.edu.au/~arunima/stellar-actions-I-data/**

You will need to download the HDF5 file `actions.h5` to run:

- `analysis/action_clustering.py`  
- `analysis/gal_env_effect.py`  
- `analysis/data_for_model.py`  
- `analysis/kde_for_dt.py`

Because the full set is extremely large (hundreds of GB), it is not stored in this repository.

Each of these scripts contains a `USER CONFIGURATION` block where you can set:
- The directory where you downloaded the simulation data  
- Which action component (`which_j`) to analyze  
- Output locations for HDF5 and logs  

For more details, see the README in the `analysis/` directory.

---

## 3. Summary of Required Inputs per Script

| Script | Uses Local Data? | Uses External Simulation? | Notes |
|-------|------------------|---------------------------|-------|
| `action_clustering.py` | no | yes | used for Fig 2,3 |
| `gal_env_effect.py` | no | yes |  used for Fig 4 |
| `compute_stream_actions.py` |yes: `cluster.csv`, `all_stars.csv` | no | Produces observational $\Delta J$ distributions |
| `data_for_model.py` | no | yes | ~200 GB output |
| `kde_for_dt.py` | no | indirectly (requires output from `data_for_model.py` | Generates KDE files |
| `run_stream_analysis.py` | requires output from `compute_stream_actions.py` and `kde_for_dt.py` | no | Reconstructs initial sizes |
| `f_comp.py` | output from `compute_stream_actions.py` and `kde_for_dt.py`, output from `run_stream_analysis.py` or `Table_1.csv` | no| Computes completeness |


