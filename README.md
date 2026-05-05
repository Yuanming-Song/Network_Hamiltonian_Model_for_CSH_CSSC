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

## Rdata Structure

- `results/csh_cssc_1node_fit.Rdata` contains `best.dev`, `best.fit`, `best.gof`, `best.terms`, `cand`, `dmax`, `maxdegs`, `seed`, and `terminc`.
- `results/csh_cssc_1node_fit_CG.Rdata` contains `best.dev`, `best.fit`, `best.gof`, `best.terms`, `cand`, `dmax`, `hist.dev`, `hist.terms`, `maxdegs`, `seed`, `term_counter`, and `terminc`.
- `results/csh_cssc_1node_fit_RestFiber.Rdata` contains `best.dev`, `best.fit`, `best.terms`, `cand`, `dmax`, `hist.dev`, `hist.terms`, `maxdegs`, `seed`, `term_counter`, and `terminc`.
- In all three fit objects, `best.fit` stores coefficient and fit summaries such as `coef`, `cov`, `se`, `ss.cov`, `ss.mean`, `ss.target`, `terms`, `dmax`, and held-out deviance fields.
- `best.gof` is present in the atomistic and coarse-grained `.Rdata` files. The atomistic file stores compact GOF summaries (`degree.q`, `degree.z`, `comp.q`, `comp.z`, `esp.q`, `esp.z`), while the coarse-grained file also stores observed, mean, standard-deviation, quantile, and type-indexed GOF outputs.
- `best.terms` is a logical vector indicating the selected model terms. Its length is 25 for the atomistic fit, 27 for the coarse-grained fit, and 32 for the restrained-fiber fit.

## Data Note

This bundle includes analysis and plotting code, selected result files, and plot outputs. Large raw simulation data are not included.
