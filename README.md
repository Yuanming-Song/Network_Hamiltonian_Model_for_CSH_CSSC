# Network Hamiltonian Model for CSH-CSSC system

This directory is the publication-facing subset of the NHMs project.

## Structure

```text
GitHub_Publication/
├── scripts/
│   ├── csh_cssc_1node_estimation_RestFiber.R
│   └── csh_cssc_1node_functions.R
├── plotting/
│   ├── AA_analysis_results.R
│   ├── CG_analysis_results.R
│   ├── NHM_analysis_plots.Rmd
│   ├── visualize_csh_cssc_networks.R
│   └── visualize_simulated_networks.R
├── results/
│   ├── README.md
│   ├── csh_cssc_1node_fit.Rdata
│   ├── csh_cssc_1node_fit_CG.Rdata
│   ├── csh_cssc_1node_fit_RestFiber.Rdata
│   └── final_AA_vs_CG_comparison_20250925_1857.dat
└── plots/
    └── Rplots.pdf
```

## Scripts

- `scripts/csh_cssc_1node_estimation_RestFiber.R`: runs the restrained-fiber NHM fitting workflow and writes result files to `results/`.
- `scripts/csh_cssc_1node_functions.R`: helper functions used by the restrained-fiber NHM fitting script.
- `plotting/AA_analysis_results.R`: loads atomistic fit results and generates summary tables and GOF plots.
- `plotting/CG_analysis_results.R`: loads coarse-grained fit results and generates summary tables and GOF plots.
- `plotting/NHM_analysis_plots.Rmd`: main R Markdown driver for the bundled plotting workflow.
- `plotting/visualize_csh_cssc_networks.R`: plots molecular-dynamics network snapshots when the underlying simulation data are available locally.
- `plotting/visualize_simulated_networks.R`: plots simulated fitted networks if simulated-network result files are added separately.
- `results/README.md`: documents the structure of the bundled `.Rdata` result files and the comparison data file.

## Results

See [results/README.md](/dfs9/tw/yuanmis1/mrsec/CSH-CSSC/NHMs/GitHub_Publication/results/README.md) for the `.Rdata` directory tree and object structure.

## Data Note

This bundle includes analysis and plotting code, selected result files, and plot outputs. Large raw simulation data are not included.
