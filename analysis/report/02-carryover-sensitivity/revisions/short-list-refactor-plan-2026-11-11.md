# Paper 02: refactor tables and figures to the short list
*2026-11-11 11:05 PST*

pmsimstats team

## Decisions (author, 2026-11-11)

- Short list of analyses in every table and figure: the five rungs
  only. G9 (side analysis), G3 and G4 (working-covariance
  counterparts of G8 and G9), G2, G5, G7 and the fixed-half-life G8
  variants leave the tables and figures.
- Labels in paper 02 only, by rung number. Stored codes are
  unchanged; only the display layer maps them.

| Stored code | Old display | New display | Rung |
|---|---|---|---|
| `E0` | E0 | A0 | 0, starting analysis |
| `E1` | G1 | A1 | 1, binary indicator, working-covariance SE |
| `E1cr2` | G6 | A2 | 2, binary indicator, CR2 |
| `E2cr2` | G8 | A3 | 3, exposure-weighted at the true half-life (oracle) |
| `E7null_cr2` | G10 | A4 | 4, half-life by AIC without the interaction, CR2 |

- One sentence in the paper maps A4 to G10 in the companion papers
  (01, 15).
- About 12 tables move to `supplement.Rmd`.
- "Working-covariance standard error" against "cluster-robust
  standard error (CR2)" replaces "model-based" (done 2026-11-11).

## Data availability for the short list

- Family 1 (job 27, 30 cells): G10 backfilled for all 30 cells
  (`g10-backfill/summary.rds`, source `t1`). Tables 6, 8, 9 and the
  job-27 figures can carry A4 without new simulation.
- Family 2 (job 21): G10 backfilled for 4 cells only (source `gvf`).
- Supplementary study (`g9-gap-fill`, families 2 to 4): G10
  backfilled for 60 rows (source `gf`).
- Job 26 (full decay-form factorial): no G10. Text that uses job 26
  states that A4 was not fitted there.

## Stages

1. Label layer: add the A0-A4 map to `spec-labels.R` (new vectors;
   existing G vectors untouched, other papers read them).
2. Ladder script (`37-ladder-tables-figures.R --scope predictive`):
   drop G9 and the optional rows from the figures and display tables;
   use A labels; rerun (stored data only).
3. Job-27 tables (6, 8, 9) and `40-stage-power-figure.R`,
   `12-hendrickson-heatmap.R`: add A4 from the backfill; reduce to the
   short list; rerun the figure scripts.
4. Remaining tables: reduce to the short list or move to the
   supplement (see table map below).
5. Prose: replace display codes with A labels; remove G9, G3, G4 from
   the text, keeping one Discussion sentence on why A4 selects without
   the interaction; update glossary and abstract.
6. Rebuild `summary.Rmd`; render report, summary and supplement
   through `tools/render.sh`; run `notation-lint.pl`.

## Table map

| # | Label | Action |
|---|---|---|
| 1 | `tab-families` | Keep; caption to A labels |
| 2 | `tab-params` | Keep |
| 3 | `tab-spec-summary` | Reduce to A0-A4 |
| 4 | `tab-ladder-def` | Keep; drop side and optional rows |
| 5 | `tab-ladder-inc` | Keep; core steps only |
| 6 | `stage3-null` | Short list with A4 from backfill; supplement |
| 7 | `iso-table` | Supplement |
| 8 | `tab-type1-cr2` | A1/A3 regressor rows, WC and CR2 SE columns; supplement |
| 9 | `tab-ref-cell` | Keep; short list with A4 from backfill |
| 10 | `stage2` | Supplement (working-covariance power only) |
| 11 | `stage3-pooled` | Supplement or drop (G9-centered) |
| 12 | `stage3-selection` | Supplement; A4 selections from backfill |
| 13 | `stage3-variants-null` | Supplement; A3, A4 |
| 14 | `stage3-variants-alt` | Supplement; A3, A4 |
| 15-18 | `stage4-*` | Supplement; short list |
| 19 | `stage4-dropout` | Keep; short list |
| 20 | `heavy-missing` | Supplement |
| 21 | `stage5-slope` | Supplement; A3, A4 |
| 22 | `stage3-coverage` | Keep; A3, A4 |
| 23 | `tab-ladder-score` | Keep; short list |
| 24 | `tab-size-equiv` | Keep; A4 against A3 |

## Progress

- 2026-11-11, stage 1 done: `a_code`, `a_labels`, `to_a()`,
  `a_factor()` added to `spec-labels.R`.
- 2026-11-11, stage 2 done: `37-ladder-tables-figures.R --scope
  predictive` draws only A0-A4 (no G9, no isolation analyses, no side
  gap); figure 3 shows A3 and A4 coverage (A4 from
  `m5-coverage-g10.csv`; its patient-level dropout cells are found by
  the `-pt` suffix); figure 3 family axis uses the output numbering
  (5, not 7). CSVs regenerated: values unchanged except the
  floating-point last digit and Table B's coverage columns, which now
  hold G10 instead of G9. Captions of the three ladder figures in
  `report.Rmd` updated.
- Next: stage 3 (job-27 tables and figures with A4 from the backfill).
- 2026-11-11, stage 3 done: job-27 G10 (A4) joined from the backfill in
  the report setup (`a27_g10`, `s27_g10`; 4 guard-failed replicates
  dropped; 494 to 500 replicates per cell). Tables `stage3-null` and
  `tab-ref-cell` show A0-A4; `tab-type1-cr2` has one row per
  regressor with WC and CR2 columns plus an A4 row (CR2 only).
  `40-stage-power-figure.R` rewritten for A0-A4; new
  `45-ladder-heatmaps.R` writes `02-heatmap-ladder-{a,b}.pdf` (the
  nine-specification `02xs-heatmap-hendrickson-*` files are untouched;
  the supplement, summary and bullets still use them). Purled report
  code runs without error. Finding for stage 5: at the Hybrid
  reference cell (N = 70, t1/2 = 1) A4 has kappa 0.82 (empirical SD
  1.28, mean SE 1.05).
- 2026-11-11, stage 4 table contents done (placement not yet done):
  `tab-spec-summary` is now the single A0-A4 definition table (from
  `ladder-table-C.csv`, with the companion-paper code and stored code);
  `tab-ladder-def` deleted and its references redirected (report and
  summary template). `tab-ladder-inc` shows the four core steps only.
  `stage2` and `stage3-pooled` deleted (redundant without G3, G4, G9).
  `stage3-selection` rebuilt for A4 from the backfill `h` column
  (job-27 Hybrid and supplementary-study cells, guarded). Referee-
  variant tables reduced to A3 and A4 (fixed-half-life rows dropped;
  their numbers stay available to the prose via `fxA`). Stress,
  dropout, slope, coverage, scoreboard and equivalence tables show the
  short list; `focus`/`focus_cols` now hold A0-A4. Heavy-missingness
  table shows A0-A3 (A4 never fitted to the per-occasion cells). Table
  count 24 -> 21. Purled report code runs without error.
- Prose still references deleted tables (`tab:stage2`,
  `tab:stage3-pooled`) and G codes throughout: stage 5.
- 2026-11-11, stage 4 placement and stage 5 done: shared data code in
  `shared-data.R` (sourced by report, summary and supplement); the old
  supplement archived as `archive/supplement-2026-10-21.*`; the new
  `supplement.Rmd` holds Tables S1-S12; report prose relabeled A0-A4,
  G9/G3/G4/G2/G5/G7 and the fixed-half-life finding removed (one
  sentence in Section sec-comparators keeps the selection variant's
  size excess as the reason A4 leaves the interaction out). Summary
  template rewritten for the main-text tables only.
- 2026-11-11, stage 6: rendered report (45 pp), supplement (7 pp) and
  summary through `tools/render.sh`; no unresolved references.
  OPEN: `notation-lint.pl` flags A1-A3 as `legacy-spec-label` (152
  hits; NOTATION.md retired A1-A3 on 2026-07-31 with different
  meanings). Author rendered with A0-A4 as is; letter choice pending.
- 2026-11-11: A0 dropped from every power display (author decision):
  power ladder (first change shown is A1 to A2), power against
  half-life, both heatmaps, the power columns of the dropout table and
  of Tables S8 and S10; power shown as -- for A0 in the reference-cell
  table, the A0-to-A1 power increment and the scoreboard's power
  columns. A0 stays in the size ladder, Table S1 and every Type I
  column. Captions say why. Re-rendered; no unresolved references.

## Sample-size simulations for A4 (started 2026-11-11 15:10)

Driver `46-run-sample-size.R` (commit 5819e1d): A4 and A2 power over a
grid of N, designs CO/OLBDC/Hybrid, t1/2 {0, 0.5, 1}, c_bm {0.2, 0.3,
0.4}, 500 replicates per cell, both parameter sets (author decision).
Local pilot (Hybrid, c_bm 0.4, 100 reps): A4 power at t1/2 = 1 is 0.86
at N = 100 under simcore2026 and 0.61 under published, so the
published N grid is roughly doubled. zeus was unreachable (Tailscale
timeout, home alias unresolved), so the six jobs went to:

| Host | Job | Set, design |
|---|---|---|
| athena | 21 | published OLBDC |
| joe | 14 | published Hybrid |
| zeus (sim-zeus-ts) | 66, then 67 | published CO, simcore2026 Hybrid (moved from athena 22 and joe 15 when zeus came back; commit 431f126, same driver and core as 5819e1d) |
| mercury (sim-mercury-ts) | 53, then 54 | simcore2026 CO, simcore2026 OLBDC |

athena and joe build the image first (30 to 60 min under Rosetta).
Outputs: `analysis/data/carryover-sensitivity/sample-size-<set>-<design>/`.
After `simrun fetch` / `apply`: write 47-sample-size-curves.R (probit
fit of power on sqrt(N) per design, t1/2, c_bm and analysis; invert at
0.8 with an interval), then the table and figure.

### Paper 01 G10 planning material, staged for merge (2026-11-11)

Moved from paper 01 by another session (commits c5ab2118b, 62ea70327,
c91186258): `staging/from-paper01-g10-planning.Rmd` (planning chunks
and prose, G10 relabeled A4; E9 and one G9 sentence still to remove or
cite to paper 01) and `staging/reconcile-g10-with-46.R` (probit in
c_bm*sqrt(N) at t1/2 = 0 and 1, likelihood-ratio test of a study
difference, N80 side by side). Both use the published parameter set.
Papers 01 and 15 now cite paper 02 for the A4 sample sizes, so paper 02
must carry them. Merge after joe job 14 (published Hybrid) is fetched:
run the reconciliation first; merge only if the two studies agree.
Decision needed: paper 02's other results use simcore2026, these use
published.
