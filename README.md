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

## ERGM Terms

The restrained-fiber fitting script evaluates the following candidate ERGM/NHM terms. A baseline `edges` term and Krivitsky size offset are added automatically by the fitting code and are not listed in `cand`.

- `kstar(2)`: two-star term; captures tendency for nodes to form multiple ties.
- `nodematch("Resname")`: overall homophily by residue type, i.e. preference for same-type ties.
- `nodematch("Resname", levels=I(resnames), diff=TRUE)`: type-specific homophily, separating same-type effects by residue class.
- `nodefactor("Resname", levels=I(resnames)[2])`: main-effect shift associated with the second residue type, `CSSC`.
- `isolates`: counts degree-0 nodes.
- `dimers`: counts isolated two-node connected components.
- `degree(0, ...)`, `degree(1, ...)`, `degree(2, ...)`, `degree(3, ...)`: residue-specific degree terms; these track how often nodes of a given residue type have degree 0, 1, 2, or 3.
- `gwdegree(0.25, fixed=TRUE)`, `gwdegree(0.5, fixed=TRUE)`, `gwdegree(1, fixed=TRUE)`, `gwdegree(2, fixed=TRUE)`, `gwdegree(3, fixed=TRUE)`: geometrically weighted degree terms; these summarize overall degree heterogeneity while downweighting very high degrees.
- `components`: counts connected components.
- `compsizesum(pow=2)`: emphasizes larger connected components by summing squared component sizes.
- `esp(0)`, `esp(1)`, `esp(2)`: edgewise shared partner terms; these count edges whose endpoints share 0, 1, or 2 common neighbors.
- `nsp(1)`, `nsp(2)`: non-edgewise shared partner terms; these count nonadjacent node pairs sharing 1 or 2 common neighbors.
- `gwesp(0.05, fixed=TRUE)`, `gwesp(0.25, fixed=TRUE)`, `gwesp(0.5, fixed=TRUE)`, `gwesp(0.75, fixed=TRUE)`, `gwesp(1, fixed=TRUE)`, `gwesp(2, fixed=TRUE)`: geometrically weighted edgewise shared partner terms; these summarize triadic closure / clustering across shared-partner counts.
- `degcor`: degree-correlation term; captures dependence between the degrees of linked nodes.

## Results

See [results/README.md](/dfs9/tw/yuanmis1/mrsec/CSH-CSSC/NHMs/GitHub_Publication/results/README.md) for the `.Rdata` directory tree and object structure.

## Data Note

This bundle includes analysis and plotting code, selected result files, and plot outputs. Large raw simulation data are not included.
