# Migration plan: common Morris-compatible Results sections across all nine pmsimstats manuscripts

*2026-05-07 18:21 PDT*

## Objective

Bring all nine manuscripts in `analysis/report/01-` through
`09-` into a common Results-section format that meets the Morris,
White and Crowther (2019) ADEMP performance-measure
recommendations end-to-end. The format must be reproducible (a
single shared infrastructure package generates the tables and
figures rather than each paper rolling its own), comparable
(rows and columns are the same across papers, so cross-paper
synthesis is trivial), and complete (every reported estimate is
accompanied by its Monte Carlo standard error).

The plan accepts the cost of **rerunning simulations where
necessary** to satisfy the format. Some papers' existing
simulations preserved only summary outputs that lack the
per-replicate granularity needed to compute MCSE retrospectively;
those papers need their simulations rerun with the canonical
output schema.

## What 'common Morris-compatible' means

### Canonical Results-section structure

Every paper's Results section follows the same six-subsection
structure:

```
# Results

## Setup and replicate counts
   (Per-cell n_reps, total cell count, total compute time,
    parameter-grid summary; confirms what was actually run.)

## Power and type I error
   (Power for primary estimand at each cell, with binomial-
    proportion MCSE; type I error under each null cell, with
    same MCSE.)

## Bias of the primary estimand
   (Mean bias and Monte Carlo SE of the bias; reported as
    absolute and as relative to the calibrated true effect.)

## Coverage of nominal-95% CIs
   (Empirical coverage rate of nominal 95% CIs at each cell,
    with binomial-proportion MCSE; coverage of the model SE
    against the empirical SE across replicates.)

## Convergence and computational cost
   (Convergence rate per cell, mean fit time, any cells with
    convergence below 90% flagged.)

## Sensitivity and cross-cell synthesis
   (Paper-specific axis-by-axis sensitivity plots; any
    cross-paper comparisons that are made.)
```

Papers with multiple substudies (e.g., 02 with three carryover
specs, 03 with five studies, 04 with M1 plus S1-S12, 09 with
four studies) repeat the six-subsection block per substudy, with
a synthesis subsection at the end.

### Canonical results-table format

Every per-cell Results table has the same column structure:

| cell | n_reps | power | mcse(power) | bias | mcse(bias) | empSE | mcse(empSE) | modSE | coverage | mcse(coverage) | conv_rate |
|---|---|---|---|---|---|---|---|---|---|---|---|

where `cell` is a paper-specific identifier (e.g., `(N, t1/2,
design)`), `n_reps` is the per-cell replicate count, and the
remaining columns are the Morris ADEMP performance measures with
their MCSEs.

### Canonical figure conventions

- **Power-by-axis plots**: x-axis is the perturbed parameter,
  y-axis is power $\in [0, 1]$, error bars are 95% binomial
  CIs (Wilson interval), faceting by other axes.
- **Bias-by-axis plots**: x-axis is the parameter, y-axis is
  $\hat{\beta} - \beta_{\text{true}}$, error bars are
  $\pm 1.96 \cdot \text{MCSE}(\hat{\beta})$.
- **Heatmaps**: cells coloured by the primary performance
  measure (typically power); each cell annotated with the
  point estimate and its MCSE in subscript or in a tooltip.
- **Coverage plots**: x-axis is the parameter, y-axis is
  empirical coverage with a horizontal reference line at 0.95;
  error bars as binomial 95% CIs.
- **Palette**: `viridis` or `cividis` (colour-blind safe);
  identical across all nine papers.

### Per-replicate output schema

Every simulation driver writes per-replicate output to
`analysis/data/<paper-slug>/<study-slug>-replicates.rds` with
the schema:

```
data.table(
  cell_id, rep_idx,
  beta_hat, beta_se, p_value, ci_lo, ci_hi, converged, fit_time
)
```

plus paper-specific columns. The per-replicate granularity is
mandatory because all Morris performance measures (especially
MCSE of bias, empirical SE, and coverage) require recomputation
from individual replicates rather than summary-only output. Any
existing driver that writes only per-cell summaries must be
modified to write per-replicate data.

## Shared infrastructure: `tools/morris-helpers.R`

A single file at `tools/morris-helpers.R` provides the helper
functions used by every paper's Results-generation chunks:

```r
## morris_summary(replicates, by, beta_true)
## Takes a data.table of per-replicate output and a vector of
## grouping columns; returns a per-cell Morris-style summary
## table with all performance measures and their MCSEs.

morris_summary <- function(replicates, by, beta_true) { ... }

## morris_table(summary_df)
## Renders the summary as a kable table with the canonical
## column structure above; suitable for inline use in an Rmd
## via knitr::kable() with longtable / booktabs styling.

morris_table <- function(summary_df, caption = NULL) { ... }

## morris_power_plot(summary_df, x_axis, facet = NULL)
## Power-by-axis plot with binomial 95% CI error bars.

morris_power_plot <- function(summary_df, x_axis,
                              facet = NULL) { ... }

## morris_bias_plot(summary_df, x_axis, facet = NULL)
## Bias-by-axis plot with bias-MCSE error bars.

morris_bias_plot <- function(summary_df, x_axis,
                             facet = NULL) { ... }

## morris_heatmap(summary_df, x_axis, y_axis, fill = 'power')
## Heatmap of any performance measure across two axes; cells
## annotated with point estimates and MCSEs.

morris_heatmap <- function(summary_df, x_axis, y_axis,
                           fill = 'power') { ... }
```

These functions are sourced once at the top of each paper's
`setup` chunk:

```r
source(here::here('tools', 'morris-helpers.R'))
```

and used directly in inline R chunks throughout the Results
section. Building this shared infrastructure is the first
deliverable of the migration; once it exists, the per-paper
work becomes mostly mechanical.

## Per-paper migration assessment

### Paper 01-dgp-mean-moderation-vs-mvn

**Status**: substantive simulation discussion exists, but uses
LaTeX `\section{}` numbering rather than markdown headings; no
dedicated `# Results` section; existing power tables lack MCSE
columns.

**Existing data**: simulation outputs are preserved at
`analysis/data/figure4/` (per the `analysis/scripts/figure4/`
driver). Per-replicate granularity present.

**Migration action**:
1. Reorganise the body sections so that the simulation results
   move from §3 and §6 into a dedicated `# Results` section
   following the canonical six-subsection structure.
2. Convert `\section{}` numbering to markdown `#`/`##` headings.
3. Recompute MCSE columns from preserved per-replicate output;
   no rerun needed.
4. Replace ad-hoc tables and figures with `morris_table()` and
   `morris_*_plot()` calls.

**Rerun needed?** **No** — per-replicate data preserved.

**Estimated effort**: 1 day of writing; ~4 hours for figure
regeneration and verification.

### Paper 02-carryover-sensitivity

**Status**: substantive Results section already in place, with
power, bias, type I error, coverage, and MCSE language
throughout. The closest paper to canonical-format compliance.

**Existing data**: per-replicate output preserved at
`analysis/data/02-grid-summary.rds` and
`analysis/data/02-sensitivity-summary.rds`. Confirmed
per-replicate granularity.

**Migration action**: minor reformatting to align column
structure with the canonical `morris_table()` output.
Replace existing custom tables and figures with the shared-
infrastructure helpers.

**Rerun needed?** **No** — per-replicate data preserved.

**Estimated effort**: 0.5 days reformatting + 4 hours for
table/figure regeneration.

### Paper 03-latent-class-mixture-application

**Status**: conceptual/methods paper with no Results section
yet. The simulation programme's Phase 1 results are in the
Discussion's 'Progress to date' subsection. Pre-registered
five-study programme committed in
`analysis/scripts/latent-class-mixture-application/00-ademp-pre-registration.md`.

**Existing data**: Phase 1 pilot results at
`analysis/scripts/latent-class-mixture-application/output/03-study5-pilot.rds`
(per-replicate). The pre-registered Studies 1-4 production runs
have not been executed.

**Migration action**:
1. Run the pre-registered Studies 1-5 production simulations
   per the ADEMP pre-registration document. Total cell count
   ~2,300 (Study 1 alone is now 1,152 cells with the new
   pb_class_correlation axis).
2. Move the Phase 1 pilot results from the Discussion to a
   new `# Results` section under the canonical structure.
3. Populate the canonical Results subsections from the
   production-run output.

**Rerun needed?** **Yes — partial run already complete (Study 5
pilot, 100 reps); Studies 1-4 not yet run.**

**Estimated compute**: ~80 hours on 16 cores parallel for
Studies 1-3; Study 4 ~8 hours; total programme ~90-100
core-days. Feasible over 2 weeks of background compute on a
workstation.

### Paper 04-treatment-main-effect

**Status**: nominal Results section exists but is a stub.
ADEMP pre-registration committed in
`analysis/scripts/treatment-main-effect/00-ademp-pre-registration.md`
covering M1 (primary comparison) and S1-S12 (sensitivity sweeps).

**Existing data**: simulation infrastructure in
`analysis/scripts/treatment-main-effect/`. Whether the
production cells have been run requires checking
`analysis/data/derived_data/sim_workspace.RData`.

**Migration action**:
1. Verify production-run output at
   `analysis/data/derived_data/sim_workspace.RData` is
   per-replicate (not just summary) and aligned with the
   pre-registration cell counts.
2. If yes: populate the canonical Results section from
   existing data.
3. If no: extend the driver to write per-replicate output and
   rerun M1 plus S1-S12.

**Rerun needed?** **Likely partial** — some sweeps may have
summary-only output; rerun affected ones with the canonical
schema.

**Estimated compute**: M1 alone is ~1 hour per cell × 8 cells =
8 core-hours; S1-S12 sweeps ~30 core-hours total; full programme
~40 core-hours, feasible overnight on 8 cores.

### Paper 05-nof1-design-sensitivity

**Status**: substantive Results section in place, fully
compliant on power and MCSE for the twelve sensitivity sweeps,
but **coverage absent**. Sweeps S1-S12 reported with MCSE error
bars but no coverage analysis.

**Existing data**: per-replicate output preserved at
`analysis/data/05-sensitivity-replicates.rds` (verify).

**Migration action**:
1. Add a coverage subsection per sweep, computing nominal-95%
   CI coverage with binomial MCSE.
2. Reformat existing tables and figures to use the shared
   helpers.

**Rerun needed?** **Probably no** — if per-replicate output
preserves CI bounds, coverage is computable retrospectively.
If only point estimates and SEs were preserved, a rerun may be
needed to record CI bounds. Verification step: check the rds
schema.

**Estimated effort**: 1 day if no rerun needed; 2 weeks of
background compute if S1-S12 must be rerun.

### Paper 06-component-decomposition

**Status**: methods paper with worked examples; the §6
simulation-study placeholder describes the planned simulation
but no Results section exists. Pre-registration committed in
`analysis/scripts/component-decomposition/00-ademp-pre-registration.md`
covering three studies (A, B, C).

**Existing data**: none — simulation has not been executed.

**Migration action**:
1. Implement the three studies (Study A: bias under
   three-component DGP; Study B: identifiability of full
   decomposition; Study C: subadditivity sensitivity).
2. Run the production simulations per pre-registration.
3. Populate a new `# Results` section under the canonical
   structure.

**Rerun needed?** **Yes — full programme not yet run.**

**Estimated compute**: Study A is ~120 hours on 16 cores;
Study B ~40 hours; Study C ~60 hours; total ~220 core-days.
Largest single budget in the programme; needs dedicated compute.

### Paper 07-gompertz-evaluation

**Status**: stub Results section. Pre-registration committed in
`analysis/scripts/gompertz-evaluation/00-ademp-pre-registration.md`
covering Study 0 (calibration), Study 1 (power across families),
and Study 2 (sloppiness diagnostic).

**Existing data**: none.

**Migration action**:
1. Run Study 0 (calibration sub-study).
2. Run Study 1 (4 families × 144 cells × 1000 reps).
3. Run Study 2 (sloppiness eigenvalue diagnostics).
4. Populate canonical Results section.

**Rerun needed?** **Yes — full programme not yet run.**

**Estimated compute**: ~50 core-days (Study 1 dominates).

### Paper 08-test-procedure-design-sensitivity

**Status**: stub Results section. Pre-registration committed in
`analysis/scripts/test-procedure-design-sensitivity/00-ademp-pre-registration.md`
covering Study 1 (test comparison), Study 2 (cycle-by-period
grid), Study 3 (joint design-by-test heatmap).

**Existing data**: none.

**Migration action**:
1. Implement the three test-procedure variants (RM-ANOVA,
   linear-mixed, GEE) and the cycle-by-period DGP grid in
   `implementations/simple/simulation.R`.
2. Run Studies 1-3 per pre-registration.
3. Populate canonical Results section.

**Rerun needed?** **Yes — full programme not yet run.**

**Estimated compute**: ~80 core-days (Study 2 has 240 cells).

### Paper 09-informative-dropout-by-design

**Status**: stub Results section. Pre-registration committed in
`analysis/scripts/informative-dropout-by-design/00-ademp-pre-registration.md`
covering four studies (power × dropout, two-mode bias,
'happy accident' decomposition, dropout-pattern misspecification).

**Existing data**: none.

**Migration action**:
1. Extend pmsimstats `R/censordata.R` infrastructure with
   the explicit dropout-pattern grid.
2. Run Studies 1-4 per pre-registration.
3. Populate canonical Results section.

**Rerun needed?** **Yes — full programme not yet run.**

**Estimated compute**: ~30 core-days (~370 cells total
across the four studies).

## Summary of compute budget

| Paper | Per-replicate data preserved? | Rerun? | Compute estimate |
|---|---|---|---|
| 01 | yes | no | format only |
| 02 | yes | no | format only |
| 03 | partial (Study 5 pilot) | yes (Studies 1-4) | ~90 core-days |
| 04 | likely partial | likely partial | ~40 core-hours |
| 05 | needs verification | possibly no | format if no; 2 weeks if yes |
| 06 | no | yes (full) | ~220 core-days |
| 07 | no | yes (full) | ~50 core-days |
| 08 | no | yes (full) | ~80 core-days |
| 09 | no | yes (full) | ~30 core-days |

**Total compute**: roughly 500-600 core-days at the upper bound.
On a 16-core workstation running 24/7, this is approximately
4-5 weeks of dedicated background compute. On a smaller machine
running 8 cores in non-blocking background, 2-3 months. Feasible
at the project's typical research-programme timescale.

## Sequencing

The migration is staged in five phases. Each phase is a discrete
deliverable; a paper that has reached the end of its phase is
fully canonical-format compliant.

### Phase 1: Build shared infrastructure (1-2 weeks)

- Implement `tools/morris-helpers.R` with the six helper
  functions (`morris_summary`, `morris_table`,
  `morris_power_plot`, `morris_bias_plot`, `morris_heatmap`,
  plus a setup helper that loads `kableExtra`, `ggplot2`, etc.).
- Document the per-replicate output schema in
  `analysis/data/README-replicate-schema.md`.
- Commit unit tests of the helpers against synthetic data.
- Verify the helpers can read existing per-replicate outputs
  from papers 01 and 02 (to confirm the schema is realistic).

**Deliverable**: a single source-able R file that any paper
can use; tested helpers; documented schema.

### Phase 2: Reformat papers 01, 02, 05 (1-2 weeks)

The three papers with substantive existing simulation outputs
get their Results sections reformatted using the Phase 1 helpers.
No reruns required for papers 01 and 02. Paper 05 may need a
rerun for coverage; resolved by verification before committing
to a rerun.

**Deliverable**: papers 01, 02, 05 with canonical-format Results
sections, populated from existing data.

### Phase 3: Populate paper 04 (1 week + overnight compute)

Verify paper 04's existing simulation output schema. Populate
the canonical Results section from per-replicate data; rerun
any sweeps that lack per-replicate granularity.

**Deliverable**: paper 04 with canonical-format Results section.

### Phase 4: Run and populate stub papers 07, 08, 09 (3-4 weeks)

Build the simulation drivers for papers 07, 08, 09. Run the
production cells per pre-registration. Populate canonical
Results sections.

Order: 09 (smallest compute) → 07 → 08 (largest compute).
Each can run as background while the next is being prepared.

**Deliverable**: papers 07, 08, 09 with canonical-format
Results sections.

### Phase 5: Run papers 03 and 06 (4-6 weeks)

Largest compute commitment of the programme. Execute the
five-study latent-class programme (03) and the three-study
component-decomposition programme (06).

These run in parallel during the writing phase; the Results
sections are populated as compute completes per study.

**Deliverable**: papers 03 and 06 with canonical-format
Results sections.

## Cross-cutting requirements

- **Pre-registration discipline.** Each Results section is
  populated from cells specified in the corresponding ADEMP
  pre-registration document. Any deviation goes into
  `02-deviations-log.md` for that paper before the production
  runs commit.
- **Per-replicate output mandatory.** Every driver writes
  per-replicate data to `analysis/data/<paper-slug>/`. No
  driver writes summary-only output as its primary checkpoint
  (summaries can be derived from per-replicate; the converse
  is impossible).
- **Reproducibility.** Each paper's Results section is
  generated from its per-replicate rds via the shared helpers
  on render. The Results section text is mostly fixed prose
  with inline R chunks calling the helpers; rerendering after
  a fresh simulation run regenerates the tables and figures
  automatically.
- **Cross-paper synthesis.** A planned Phase-6 cross-paper
  discussion document at
  `analysis/report/cross-paper-synthesis.md` will use the
  shared helpers' canonical output to compare results across
  papers without per-paper formatting differences obstructing
  the comparison.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Some papers' existing per-replicate data are missing or in non-canonical form | Verification step at the start of each paper's phase; rerun if necessary |
| Compute budget exceeds available capacity | Run papers in priority order; defer Phase 5 until Phases 2-4 complete |
| Helper functions need iteration once real papers exercise them | Phase 1 includes unit tests; helper API may evolve in early Phase 2 papers |
| ADEMP pre-registration parameters need revision after Phase 2 results inform expectations | Deviation log mechanism already in place per paper |
| Coverage analysis on rerun-resistant papers blocked | Document the gap and run dedicated short coverage sweeps where the full rerun is infeasible |

## Deliverables (consolidated)

By the end of the migration:

- `tools/morris-helpers.R` and unit tests in
  `tools/tests/test-morris-helpers.R`.
- `analysis/data/README-replicate-schema.md`.
- Each paper's Results section in canonical-format markdown,
  populated from per-replicate output via the shared helpers.
- Each paper's per-replicate output preserved at
  `analysis/data/<paper-slug>/<study-slug>-replicates.rds`.
- Each paper's ADEMP pre-registration committed (already in
  place after the previous turn) with `02-deviations-log.md`
  recording any in-flight changes.
- A cross-paper synthesis document at
  `analysis/report/cross-paper-synthesis.md` (Phase 6, optional).

---

*Source: `~/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/report/morris-results-migration-plan.md`*
