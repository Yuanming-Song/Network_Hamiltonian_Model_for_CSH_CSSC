# Results Directory

## Structure

```text
results/
├── README.md
├── csh_cssc_1node_fit.Rdata
├── csh_cssc_1node_fit_CG.Rdata
├── csh_cssc_1node_fit_RestFiber.Rdata
└── final_AA_vs_CG_comparison_20250925_1857.dat
```

## Rdata Structure

- `csh_cssc_1node_fit.Rdata`

```text
.
├── best.dev
├── best.fit
├── best.gof
├── best.terms
├── cand
├── dmax
├── maxdegs
├── seed
└── terminc
```

- `csh_cssc_1node_fit_CG.Rdata`

```text
.
├── best.dev
├── best.fit
├── best.gof
├── best.terms
├── cand
├── dmax
├── hist.dev
├── hist.terms
├── maxdegs
├── seed
├── term_counter
└── terminc
```

- `csh_cssc_1node_fit_RestFiber.Rdata`

```text
.
├── best.dev
├── best.fit
├── best.terms
├── cand
├── dmax
├── hist.dev
├── hist.terms
├── maxdegs
├── seed
├── term_counter
└── terminc
```

## Object Notes

- `best.dev` is the best held-out deviance summary used during model selection. In the fitting scripts it is stored as a two-value vector: test deviance and the corresponding standard error.
- `best.fit` is the main fitted compositionally pooled ERGM object.
  - `ERGM` means exponential random graph model.
  - `coef` stores the fitted model coefficients.
  - `cov` stores the estimated covariance matrix of the fitted coefficients.
  - `se` stores the coefficient standard errors.
  - `ss.cov` means covariance matrix of the simulated sufficient statistics used in the stochastic-approximation fit.
  - `ss.mean` means mean simulated sufficient statistics at the fitted parameter values.
  - `ss.target` means target sufficient statistics computed from the training networks.
  - `terms` stores the selected ERGM/NHM terms included in the fitted model.
  - `dmax` means the maximum allowed node degree used as the degree-bound constraint during estimation and simulation.
  - `deviance.test` is the held-out test deviance computed on the test set.
  - `deviance.test.se` is the standard error of the held-out test deviance estimate.
  - `train.sample` stores the compositionally pooled training-sample metadata, including composition weights (`w`), composition types (`types`), type indices (`type.ind`), and counts per type (`type.count`).
  - In the restrained-fiber fit, `work_msg` and `finish_msg` store progress/log text for batch candidate evaluation.
- `best.gof` is the goodness-of-fit summary object.
  - `GOF` means goodness of fit.
  - It is present in the atomistic and coarse-grained result files.
  - `degree` refers to node degree distributions.
  - `comp` refers to component statistics, i.e. connected-component size/count summaries.
  - `esp` refers to edgewise shared partners.
  - `.obs` means observed values from the held-out networks.
  - `.mean` means mean values across simulated networks from the fitted model.
  - `.sd` means standard deviation across simulated networks.
  - `.q025` and `.q975` are the 2.5% and 97.5% simulation quantiles, used as an approximate 95% simulation interval.
  - `.q` stores quantile-position summaries comparing observed statistics to the simulated reference distribution.
  - `.z` stores z-score style standardized discrepancies between observed and simulated values.
  - `types`, `type.count`, and `type.ind` index the composition classes used when pooling across systems.
- `best.terms` is a logical vector marking which candidate model terms were selected in the final model. Its length is 25 for the atomistic fit, 27 for the coarse-grained fit, and 32 for the restrained-fiber fit.
- `cand` means the full candidate term list considered during forward model selection.
- `maxdegs` stores the observed maximum degrees by system/composition before collapsing to the fitted degree bound.
- `seed` is the random-number seed used to make the fitting workflow reproducible.
- `terminc` means the current term-inclusion logical vector used during forward selection.
- `hist.dev` stores the deviance history across accepted forward-selection steps.
- `hist.terms` stores the corresponding term-inclusion history across accepted forward-selection steps.
- `term_counter` stores how many candidate terms had been tested so far in the modified batch-processing workflow.
- `final_AA_vs_CG_comparison_20250925_1857.dat` is the bundled final comparison data file comparing atomistic (`AA`) and coarse-grained (`CG`) outputs.
