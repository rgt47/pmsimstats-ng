# 06-component-decomposition

*2026-10-28 18:10 PDT* (supersedes the 2026-05-07 note)

## Working title

Three-component decomposition of treatment response in
aggregated N-of-1 trials: pharmacological, expectation-driven,
and natural-history components.

## What the paper now says

The paper formalizes the BR, PB, TV decomposition, maps its
identifiability onto the phase contrasts of the hybrid design, and
derives the omitted-variable identity for the one-component treatment
effect and biomarker slope. A corollary of the identity gives the
displacement of the one-component biomarker slope under a biomarker
correlated with the placebo response as
$-c_{bm,PB}\,\sigma_{PB}\,(\bar\eta_{on} - \bar\eta_{off})$, whose sign
is set by the design. Four exploratory simulation studies test it:

| Study | Design | What it shows | Stored as |
|---|---|---|---|
| 1 | eight-visit (N 35, 70), sixteen-visit (N 200, 300) | no bias with a biomarker correlated with BR alone; phase-augmented analysis dominated | `study-a-alt1000-null5000` |
| 2 | eight-visit N 70, sixteen-visit N 300 (100 reps) | displacement in the predicted direction on the feasible eight-visit cells; sixteen-visit cells at $c_{bm,PB} \ge 0.15$ are repair-artifact exhibits | `study-contam-alt100-null100` |
| 3 | sixteen-visit (1,000 reps) | F1 and three other analyses unbiased at $c_{bm,PB} = 0$; contaminated cells are exhibits | `study-b-recovery-pilot` |
| 4 | balanced-placebo (1,000 reps) | one-component slope displaced upward as predicted; F1 recovers the pharmacological slope at every feasible level | `study-b-balanced` |

Feasibility bounds, the repaired-process population values for the
exhibit cells, and the identity predictions are computed by
`analysis/scripts/component-decomposition/05-feasibility-bounds.R` into
`feasibility-bounds.rds`. All data directories are under
`analysis/data/derived_data/component-decomposition/`.

## Pending

`run-study-2-feasible.R` is the corrected contamination study: feasible
cells only, `stopOnNonPD = TRUE` (the default), seed 2061, 22
alternative cells at 1,000 replicates and 6 null cells at 2,000, with a
pre-flight eigenvalue check. A 3-replicate smoke test passed on
2026-10-28. Estimated cost about 45 core-hours. Not yet run; the
manuscript's section 7.7 renders a RESULTS PENDING notice until
`study-2-feasible/study-2-summary.rds` exists.

Not run and listed as such in section 7.8: the full nonlinear fit
(protocol Study A method 3 and Study B), protocol Study C
(subadditivity), the eta-mismatch, TV and dropout grids, a prognostic
$c_{bm,TV}$ cell, Study 1 at N 200 and 300 on the eight-visit design, a
distributed-lag arm.

## Revision history

- **2026-10-28.** Referee report `referee-report-2026-10-28.md`
  addressed. Finding 1: the June contamination and hybrid recovery runs
  drew their cells at $c_{bm,PB} \ge 0.15$ (sixteen-visit) and $0.45$
  (eight-visit) from silently repaired non-positive-definite matrices;
  the `+0.51` headline and the "30 MCSE", "scales monotonically" and
  "toward zero" statements were withdrawn, a feasibility table added,
  and the affected cells re-reported as exhibits against the repaired
  population value. Finding 2: the identity corollary added to section
  6.1 and used to predict each design; the N 70 cells recognized as the
  informative ones. Finding 3: convergence-rate claim for the stub
  nonlinear fit deleted; Study C, eta-mismatch, TV and dropout grids
  moved to a not-run list; Gompertz-basis recovery claim deleted;
  Reproducibility rewritten from the drivers (seeds 642/644/916/815,
  `set.seed(seed + 1000 cell + rep)`, derived_data paths, the
  2026-08-08 N migration). Finding 4: abstract rewritten under the
  30-word rule without dated notes; section 1.1 made conditional;
  `D_{it}` versus `D_{bc,it}` per NOTATION rule 1; "coupling" removed
  from prose; glossary added; "enroll". Finding 5: expectancy defined
  once in the terms table (covert on-drug 0, concealed placebo 1,
  disclosed open-label placebo not simulated); deception paragraph with
  Miller et al. 2005 and Meissner et al. 2023; Petkova 2009 and Tarpey
  2010 cited; the Study B decomposition named as the series' F1 model
  with the Cor(D, eta) difference stated. Minor comments: coverage
  0.99, 1.8 MCSE, the m_PB-row attribution, the N 70 c 0.2 cell of
  Study 4 reported, the N 300 c 0.15 phase-augmented difference
  acknowledged, worked example marked as the noiseless limit, TV
  increment quantified, `in_header: claudecode.tex` removed, bib fixed
  (nurko2022 authors, rohsenow1981 DOI). Studies numbered.
- **2026-10-08.** Dated correction notes added for: the sign of the
  d_BR interpretation; the exact aliasing of `Dbc:phase` at every N
  (the production run fitted the fallback formula); the Study B
  within-one-MCSE claim (three of four cells at 1.0 to 1.4 MCSE); the
  Cor(D, eta) figure (0.753 on the sixteen-visit design, 0.485 on the
  four-path Hybrid); the unreproduced explanation of the contamination
  sign. All absorbed into the prose on 2026-10-28.
- **2026-08-08.** Stored `N` recoded from per-path to total
  (`migrate-n-to-total.R`), no rerun.
- **2026-06-12 to 06-16.** Four referee rounds; Study A production run,
  contamination pilot, hybrid recovery run and balanced-placebo Study B
  added.
- **2026-05-07.** Manuscript, ADEMP protocol and pilot created.

## Relationship to other papers in the series

- `01-dgp-mean-moderation-vs-mvn`: shares the Gompertz parameterization
  and the correlation encoding of the interaction.
- `02-carryover-sensitivity`: the F-family analyses (F1 and its
  variants) evaluated there and in `docs/47` to `docs/49`, `docs/56`;
  `docs/65` documents the same repair artifact in that paper's drivers.
- `04-treatment-main-effect`: uses the three-component DGP for the main
  effect.

## Author

pmsimstats team
