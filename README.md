# Network Hamiltonian Model for CSH-CSSC system

This directory is the publication-facing subset of the NHMs project.

## Structure

```text
GitHub_Publication/
├── scripts/
│   ├── csh_cssc_1node_estimation_RestFiber.R
│   ├── csh_cssc_1node_functions.R
│   └── submit_R_analysis.sh
├── plotting/
│   ├── AA_analysis_results.R
│   ├── CG_analysis_results.R
│   ├── NHM_analysis_plots.Rmd
│   ├── visualize_csh_cssc_networks.R
│   └── visualize_simulated_networks.R
├── results/
│   ├── csh_cssc_1node_fit.Rdata
│   ├── csh_cssc_1node_fit_CG.Rdata
│   ├── csh_cssc_1node_fit_RestFiber.Rdata
│   ├── final_AA_vs_CG_comparison_20250925_1857.dat
│   └── intermediate_CG_estimate_comparison_20250925_11am.dat
└── plots/
    └── Rplots.pdf
```

## Scripts

- `scripts/csh_cssc_1node_estimation_RestFiber.R`: runs the restrained-fiber NHM fitting workflow and writes result files to `results/`.
- `scripts/csh_cssc_1node_functions.R`: helper functions used by the restrained-fiber NHM fitting script.
- `scripts/submit_R_analysis.sh`: SLURM wrapper for running the restrained-fiber analysis with a user-defined main directory.
- `plotting/AA_analysis_results.R`: loads atomistic fit results and generates summary tables and GOF plots.
- `plotting/CG_analysis_results.R`: loads coarse-grained fit results and generates summary tables and GOF plots.
- `plotting/NHM_analysis_plots.Rmd`: main R Markdown driver for the bundled plotting workflow.
- `plotting/visualize_csh_cssc_networks.R`: plots molecular-dynamics network snapshots when the underlying simulation data are available locally.
- `plotting/visualize_simulated_networks.R`: plots simulated fitted networks if simulated-network result files are added separately.

## Path Setup

Set the bundle root once with `maindir` in R or `NHM_MAIN_DIR` in the shell.

```bash
export NHM_MAIN_DIR=/path/to/GitHub_Publication
export NHM_DATA_DIR=/path/to/restrained_fiber_input_data
cd /path/to/GitHub_Publication/scripts
Rscript csh_cssc_1node_estimation_RestFiber.R
```

```r
maindir <- "/path/to/GitHub_Publication"
source(file.path(maindir, "plotting", "AA_analysis_results.R"))
```

## Data Note

This bundle includes analysis/plotting code, result files used for plotting, and plot outputs. Large raw simulation data are not included.
# Network_Hamiltonian_Model_for_CSH_CSSC
