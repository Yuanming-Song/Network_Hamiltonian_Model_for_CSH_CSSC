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

- In all three files, `best.fit` contains the main fit summary, including coefficient and covariance outputs such as `coef`, `cov`, `se`, `ss.cov`, `ss.mean`, `ss.target`, `terms`, `dmax`, and held-out deviance fields.
- `best.gof` is present in the atomistic and coarse-grained result files.
- In `csh_cssc_1node_fit.Rdata`, `best.gof` contains compact GOF summaries: `degree.q`, `degree.z`, `comp.q`, `comp.z`, `esp.q`, `esp.z`.
- In `csh_cssc_1node_fit_CG.Rdata`, `best.gof` contains observed, mean, standard-deviation, quantile, and type-indexed GOF outputs.
- `best.terms` is a logical vector indicating the selected model terms. Its length is 25 for the atomistic fit, 27 for the coarse-grained fit, and 32 for the restrained-fiber fit.
- `final_AA_vs_CG_comparison_20250925_1857.dat` is the bundled final comparison data file.
