# Analysis-ladder tables and figures (paper 02)

Script: `37-ladder-tables-figures.R`. Run from the repository root.
It reads only stored data and runs no simulation (about 7 seconds).

## Outputs

CSV files in `analysis/data/carryover-sensitivity/ladder/`:

- `ladder-long.csv`: tidy long file (family, cell_group, scenario,
  level, analysis, quantity, value, mcse, n, from, reference,
  in_family_pool). `level` is `cell` or `family`. Quantities are
  `size`, `power`, `kappa`, `kappa_null`, `kappa_alt`, `size_diff`,
  `power_diff` (value is `analysis` minus `from`, mcse is the paired
  SE), and `coverage` with `coverage_lower` and `coverage_upper`.
  Figures 1 and 2 are drawn from the family-level rows.
- `ladder-table-A.csv`: one row per step, formatted cells
  `change (paired SE)` for size and power by family.
  `ladder-table-A-numeric.csv` holds the same numbers unformatted.
- `ladder-table-B.csv`: one row per analysis (families in the
  0.04 to 0.06 band, extreme cell sizes, power relative to G9, worst
  power loss against G9 by family and by cell, coverage).
- `ladder-table-C.csv`: formula, inference, assumptions and rung.

Figures in `analysis/figures/` (PDF and PNG, 6.5 inches wide):
`02-ladder-fig1-size`, `02-ladder-fig2-power`,
`02-ladder-fig3-calibration`, `02-ladder-fig4-trend`.

## Scope (`--scope full|predictive`)

`--scope full`, the default, gives the outputs above, the ladder of
paper 02 before the split. Running it after the 2026-10-31 change
reproduces the stored CSV files byte for byte.

`--scope predictive` gives the ladder of paper 02 after the split
(split plan, `revisions/split-plan-prognostic-2026-10-29.md`, Section
2.4). It writes the same four CSV files to
`analysis/data/carryover-sensitivity/ladder-predictive/` and Figures
1 to 3 with a `-predictive` suffix. It has no Figure 4. Compared with
the full scope:

- Families 5 (prognostic) and 6 (large N), the trend block and the
  block N cells with a prognostic pathway are dropped. Family 7 is the
  patient-slope cell only.
- F1, F1A0B and the fixed-half-life G8 variants are dropped.
- G9-null is labeled G10. It is rung 4 and the reference: Table B
  reports `mean_power_pct_of_g10`, and the worst-loss columns are
  against G10. G9 is a side analysis. The paired steps are `G8>G10`
  (core) and `G10>G9` (side).
- G10 replicates for families 1 to 4 and the slope block come from
  the backfill (`analysis/data/carryover-sensitivity/g10-backfill/`,
  script 41). A replicate is used only if its G9 refit reproduced the
  stored G9 (`guard.rds`). For blocks N and A, G10 is the stored
  `E7null_cr2`.
- Coverage (Figure 3, right panel) is stored for G8 and G9 only, so
  it has no G10 entry.

## Definitions

- Reject is `p_value < 0.05` among converged replicates with a
  p-value. This reproduces the stored summary rates exactly.
- Families: 1 original 30-cell grid; 2 misspecified half-life and
  decay shape (Stage 4 cells plus supplementary B1, B2); 3
  autocorrelation (B4); 4 patient-level dropout (job 33); 5 prognostic
  pathways (job 28 at strength 0.3, supplementary C1 ladder);
  6 large N (C2, null only); 7 trend pathway and patient slope
  (referee blocks T and S); S the referee Hybrid run (blocks N and A),
  the only source of G9-null and G8 at fixed half-lives.
- Size pools all null cells of a family. Power pools the reference
  alternative: family 1 c_bm 0.45 at t_half 1 (Hybrid, OL+BDC, CO,
  N 35 and 70); 2 c_bm 0.25; 3 c_bm 0.30; 4 c_bm 0.45; 5 job 28 c_bm
  0.45 (the C1 alternatives enter only the cell-level rows); 7 c_bm
  0.45; S c_bm 0.45 (t_half 0.5 and 2). Family 6 has no alternative.
  Family S size uses the cpb = 0 cells only; the cpb = 0.3 cells are
  in the long file with `in_family_pool = FALSE`.
- Pooled estimates weight cells by replicates. MCSE is
  sqrt(p(1 - p) / n) with n the total replicates, which treats cells
  as independent.
- Paired differences use the same (cell, replicate). The family value
  is the replicate-weighted mean of cell differences, and its SE is the
  replicate-weighted mean of the cell paired SEs. This is a
  conservative bound because cells within a job share random-number
  streams (blocks N to S are seeded per cell and replicate).
- Calibration ratio kappa is mean model SE over empirical SD per cell,
  averaged without weights over all cells of a family. Coverage is the
  mean over cells of the lower and upper bounds in
  `referee-reanalysis/m3-coverage.csv` (G8 and G9 only; families 2 to 5).
- The shaded band 0.04 to 0.06 in Figures 1 and 4 is a fixed
  descriptive tolerance around 0.05. It is not a confidence interval
  or a test.

## Known limits

- F1 and F1A0B are not identifiable under the CO design, so family 1
  pools them over Hybrid and OL+BDC. Table B compares them with G9 on
  the same cells.
- F1A0B is stored only as cell summaries: it has no paired
  differences and no kappa.
- The trend pathway (block T) and the slope cell (block S) stored only
  G8, G9 and F1. Figure 4 therefore shows three analyses, not every
  rung. E0, G1 and G6 were not fitted there.
- G9-null and G8 at fixed half-lives exist only in family S.
- Excluded: strength-0.6 prognostic cells, mean-moderation cells,
  occasion-level dropout (job 31), the A2 size block, `progtvgrow`.

## Isolation branch

If `analysis/data/carryover-sensitivity/e0-isolation/all.rds` exists
(methods `E0, E0n, E0nc, E1nc, E1`, same columns as the tier-1
`all.rds`), rung 1 of family 1 is split into four single-factor steps
in Figures 1 and 2 and in Table A (E0 to E0n, E0n to E0nc, E0n to
E1nc, E1nc to G1). Otherwise the combined E0 to G1 step is shown with
the note "three respects at once". Override paths with
`LADDER_ISOLATION_DIR`, `LADDER_OUT_DATA` and `LADDER_OUT_FIG`.
