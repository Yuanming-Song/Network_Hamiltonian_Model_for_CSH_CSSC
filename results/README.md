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

- `best.dev`: two-value summary of the current best held-out deviance and its standard error.
- `best.fit`: main fitted ERGM object.
- `best.gof`: goodness-of-fit summary object.
- `best.terms`: logical vector marking selected candidate terms; length 25 for atomistic, 27 for coarse-grained, 32 for restrained-fiber.
- `cand`: candidate ERGM term list used during forward selection.
- `dmax`: maximum allowed node degree used in the degree-bound constraint.
- `maxdegs`: observed maximum degrees by system/composition before collapsing to `dmax`.
- `seed`: random-number seed used for reproducibility.
- `terminc`: current term-inclusion logical vector during forward selection.
- `hist.dev`: deviance history across accepted forward-selection steps.
- `hist.terms`: term-inclusion history across accepted forward-selection steps.
- `term_counter`: number of candidate terms tested so far in the batch-processing workflow.
- `coef`: fitted model coefficients.
- `cov`: covariance matrix of the fitted coefficients.
- `se`: standard errors of the fitted coefficients.
- `ss.cov`: covariance matrix of the simulated sufficient statistics.
- `ss.mean`: mean simulated sufficient statistics at the fitted parameter values.
- `ss.target`: target sufficient statistics computed from the training networks.
- `terms`: selected ERGM/NHM terms included in the fitted model.
- `deviance.test`: held-out test deviance estimate.
- `deviance.test.se`: SE of held-out test deviance.
- `train.sample`: pooled training-sample metadata.
- `w`: composition weights used in pooled fitting.
- `types`: composition classes used in pooling.
- `type.ind`: composition-class index for each network.
- `type.count`: number of networks in each composition class.
- `work_msg`: stored progress message for a candidate-model batch job in the restrained-fiber workflow.
- `finish_msg`: stored completion message for a candidate-model batch job in the restrained-fiber workflow.
- `degree`: node degree distribution summary.
- `comp`: connected-component summary.
- `esp`: edgewise shared partner summary.
- `.obs`: observed values from held-out networks.
- `.mean`: mean values from simulated networks.
- `.sd`: standard deviation across simulated networks.
- `.q025`: 2.5% simulation quantile.
- `.q975`: 97.5% simulation quantile.
- `.q`: quantile-position summary comparing observed to simulated values.
- `.z`: standardized discrepancy summary comparing observed to simulated values.
- `final_AA_vs_CG_comparison_20250925_1857.dat`: final comparison data file for atomistic (`AA`) versus coarse-grained (`CG`) results.
