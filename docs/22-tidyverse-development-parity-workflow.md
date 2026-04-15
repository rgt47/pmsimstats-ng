---
title: "Workflow for tidyverse-side development without losing parity"
author: "pmsimstats team"
date: "2026-04-15"
---

# Workflow for tidyverse-side development without losing parity

This note captures the recommended workflow for editing
`implementations/tidyverse/` while preserving numerical parity with
`implementations/original-extended/`. The cross-implementation parity
test at `implementations/test-parity-extended-tidyverse.R` is the
safety net; this document describes how to integrate it into
day-to-day development rather than treat it as an after-the-fact
audit.

## 1. Baseline before you start

Run the full parity sweep once on the clean branch and save the
report as the frozen "this is what parity looks like before I touched
anything" artefact:

```bash
git checkout -b tidyverse-refactor
Rscript implementations/test-parity-extended-tidyverse.R
cp implementations/parity_report.rds implementations/parity_baseline.rds
```

Commit `parity_baseline.rds` with the branch. At the end of the
development cycle, diff the current `parity_report.rds` against it to
confirm the only differences (if any) are the ones you intended.

## 2. Pick a tight inner loop

The default mode runs in approximately 10 seconds on current
hardware. That is fast enough to run after every substantive edit.
Use `--quick` (approximately 2 seconds, 16 cells) while iterating on
a specific function and escalate to the full sweep before committing.

A git pre-commit hook can enforce this, or the habit can be built
manually:

```bash
# pre-commit hook (optional)
Rscript implementations/test-parity-extended-tidyverse.R --quick || \
  { echo "parity broken; run --full for details"; exit 1; }
```

## 3. Edit semantics versus edit shape

Parity is preserved if a change is "same output, different code." It
is broken if the change is "different output, possibly better." Be
deliberate about which is intended:

- **Pure refactors** (rename, split function, swap `dplyr` for
  `purrr`, extract helper): parity *must* hold. If it breaks, a bug
  has been introduced.
- **Feature additions** that extend behaviour without changing the
  existing path (new optional argument with default value, new
  analysis option): parity *must still hold* for the default-argument
  call paths the test exercises.
- **Behaviour changes** (different DGP parameterisation, different
  MVN-draw ordering for speed, different Sigma caching scheme):
  parity *will* break by design. Before making such a change, decide
  whether this also means breaking parity with `original-extended`
  (in which case the parity test becomes the wrong oracle) or whether
  the change should be applied to both implementations in lockstep.

In practice, most development falls into the first two categories.
Treat the third as a separate project with explicit scope.

## 4. When parity breaks mid-refactor

The test reports which layer broke. Diagnose in this order:

1. **Sigma failure:** the change likely touched `build_sigma_matrix`
   or `build_correlation_matrix`. The BM-BR correlation branch is
   the most common culprit, as the earlier `dgp_architecture`
   scoping bug demonstrated.
2. **Data failure with Sigma passing:** the change touched the MVN
   sampling path, the post-draw mean-moderation shift, or the
   participant/path assembly logic.
3. **Analysis failure with data passing:** the change touched
   formula construction, `Dbc` computation, or the `nlme::lme`
   invocation in `lme_analysis`.

The test prints the first failing cell; that cell is usually enough
to reproduce the divergence interactively. A one-liner that compares
a single numeric column across implementations is faster than
re-running the full sweep during debugging.

## 5. Scope edits to one concern at a time

The `--quick` run covers 16 cells in approximately 2 seconds.
Re-running after every approximately 20 lines of changes localises
bugs to the edit that introduced them. The alternative, batching 200
lines of changes before debugging a dozen parity failures at once, is
substantially slower.

## 6. Guardrails around parallelisation and RNG

The parity test calls `set.seed()` at the per-cell level. If a
refactor consumes random draws in a different order (for example,
changing `purrr::map` to `furrr::future_map` or reordering per-path
MVN calls), parity will break even though the code is semantically
correct. Two defences apply:

- Do not change RNG-consuming loops unless `original-extended` is
  updated in the same commit.
- Keep `n_cores = 1` for any parity-relevant comparison. Introduce
  parallelism only after functional parity is locked in.

## 7. Expand the parity test as coverage gaps appear

When a feature is added that the current grid does not exercise (for
example, `carryover_form = "weibull"` or an analysis option toggling
random slopes), add a row to the `edge_cells` tibble at the bottom of
the script. The test grows with the surface area it protects. A
one-line edit per new feature is inexpensive, and the habit
compounds.

## 8. When divergence is intentional

If the intent is to change tidyverse behaviour in a way that
`original-extended` should not follow (for example, fixing a bug that
lives only in the data.table implementation, or using the tidyverse
implementation as a testbed for a new feature), document the
divergence explicitly:

- Note in a commit message or a `DIVERGENCES.md` file *which cells*
  of the parity test are allowed to fail and *why*.
- Mark those cells in the script with an `expected_divergent` flag
  so the summary reports "N expected divergences, 0 unexpected
  failures" rather than treating them as regressions.
- Re-baseline after the intentional divergence lands:
  `cp parity_report.rds parity_baseline.rds`.

This keeps "expected" and "broken" visibly separated. The trap to
avoid is letting an intentional divergence mask a subsequent
unintentional one.

## 9. Before merging

Run the full sweep (not `--quick`) one final time. If the diff
between `parity_report.rds` and `parity_baseline.rds` is empty (or
matches documented expected divergences), the branch is safe to
merge. Otherwise, the diff identifies exactly which cells regressed,
which is a tighter feedback loop than "tests pass locally but
something feels off."

A one-liner for the diff:

```r
a <- readRDS('implementations/parity_baseline.rds')
b <- readRDS('implementations/parity_report.rds')
all.equal(a$grid, b$grid)  # TRUE, or a character vector of diffs
```

## Checklist

- [ ] `cp parity_report.rds parity_baseline.rds` at start of branch.
- [ ] `--quick` after every focused edit (~2 s).
- [ ] Full sweep before every commit (~10 s).
- [ ] If parity breaks, read the first failing cell; fix before
      proceeding.
- [ ] Never change RNG consumption order in a pure refactor.
- [ ] Add edge cells to the script when new features are added.
- [ ] Document any intentional divergence explicitly in
      `DIVERGENCES.md` or the commit message.

## Optional infrastructure to land

If the recommendations above prove useful in practice, three pieces
of supporting infrastructure can be promoted from recommendation to
code in the repository:

1. A pre-commit hook that runs `--quick` and blocks the commit on
   failure.
2. An `expected_divergent` annotation mechanism in the parity script
   so intentional divergences do not read as regressions.
3. A baseline-diff wrapper script (`Rscript
   implementations/parity-diff.R`) that reports which cells differ
   between the current report and the baseline.

---

*Rendered on 2026-04-15 at 13:00 PDT.*
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/22-tidyverse-development-parity-workflow.md*
