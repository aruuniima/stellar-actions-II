# stellar-actions-II Repository

This repository contains the full codebase and partial datasets used in the analysis of stellar actions in a Milky Way–like hydrodynamical simulation. The work presented here is associated with the paper:

**“Evolution of action-space coherence in a Milky Way–like simulation”**  
Arunima et al. (submitted)  
[arXiv:2511.15215](https://arxiv.org/abs/2511.15215)

The repository is organised into several top-level directories, each corresponding to a major component of the project workflow. Every directory includes its own `README` with detailed descriptions of the scripts and their roles. 

The workflows in this repository are designed to be modular:  
- **`analysis/`** - Core analysis scripts. Computes action differences, clustering statistics, stream evolution measures. Contains helper utilities.
- **`data/`** -  Contains raw and processed data required for the analysis. (Some files excluded due to size.)
- **`plots/`** contains clean, lightweight scripts that reproduce the final paper figures.  

Please refer to the individual directory READMEs for specifics on running the scripts, input formats, and expected outputs.
