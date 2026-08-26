# Hendrickson et al. (2020) original-code comparison arm

*2026-08-26*

Supports Section 2.2.4 of `analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`.

## Provenance

`vendored-hendrickson-generateData.R` and
`vendored-hendrickson-lme_analysis.R` are downloaded verbatim from
Hendrickson et al.'s own public repository,
`https://github.com/rchendrickson/pmsimstats`, `master` branch, one
commit (`3035581`, 2026-06-28) before their paper's 2020-07-06
acceptance date. One line was patched for compatibility with a
modern `data.table` version (Hendrickson's original
`generateData.R` line 118 used a pre-2021 `data.table` indexing
idiom, `modelparam[(paste("c",c,sep="."))][,]`, that current
`data.table` rejects; the patch replaces it with the semantically
identical `modelparam[[paste("c",c,sep=".")]]`). No other lines were
changed. `lme_analysis.R` is unpatched.

The `Produce_Publication_Results_1_generate_data.Rmd` vignette in
the same repository (not vendored here) confirms the exact
`analysisparams <- expand.grid(useDE=FALSE, t_random_slope=FALSE,
full_model_out=FALSE)` call used for the published figures, and that
no `carryover_t1half` is passed to `lme_analysis()` (so its own
fallback sets it to 0, giving a Dbc treatment that does not track
the DGP's true carryover half-life). `generateData()`'s
`scalefactor` argument is left at its default (2), also matching
the vignette, which never overrides it.

## Scripts

- `01-hendrickson-orig-driver.R` — runs Hendrickson et al.'s original
  DGP and analysis (as vendored above) through the same three
  designs (CO, Hybrid, OL+BDC), carryover half-lives (0, 0.5, 1.0),
  and `c.bm` (0, 0.45) as the paper's own Section 3.1 factorial, at
  `N = 70` total per design. Trial designs, baseline/response
  parameters (`extracted_bp`, `extracted_rp`), and utility functions
  (`buildtrialdesign()`, `cumulative()`, `modgompertz()`) are reused
  from the current package, confirmed compatible by direct
  inspection (unchanged between the original repository and this
  package for the functions and data actually used). Run with
  `NREPS=<n> Rscript 01-hendrickson-orig-driver.R` (default 300;
  the paper's own production run uses 1,000). Outputs
  `../../../data/quick-sim/hendrickson-original-comparison/hendrickson-orig-summary.txt`
  and `-replicates.rds`.
- `02-pd-sweep.R` — sweeps `c.bm` from 0.1 to 0.8 (step 0.05) at the
  OL+BDC design, path A, `t1half = 1.0`, comparing whether the raw
  (pre-correction) covariance matrix is positive definite under
  Hendrickson's original compound-symmetry construction versus this
  package's AR(1) construction (`buildSigma(..., dgp_architecture =
  'mvn')`). Uses the project's own monkey-patch pattern for
  `corpcor::is.positive.definite` (see
  `analysis/scripts/figure4/05-pd-diagnostics.R` for the original
  use of this pattern). Prints its result table to stdout; not
  saved to a file (rerun to reproduce, seconds to run).

## Result summary (see report.Rmd Section 2.2.4 for full context)

- PD sweep: raw covariance matrix is positive definite up to
  `c.bm ~ 0.35` under compound symmetry, `c.bm ~ 0.50` under AR(1).
  At `c.bm = 0.45` (this paper's calibration), compound symmetry
  fails, AR(1) does not.
- Power comparison: Hendrickson-original's relative carryover loss
  (`t1half` 0 to 1.0) is 4.2% (CO), 12.5% (Hybrid), 19.0% (OL+BDC) —
  between this paper's Architecture A and Architecture B at every
  design.
