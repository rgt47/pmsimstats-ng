# Rerun of the stress test and covariance study under the refined process
*2026-10-08 15:16 PDT*

Author: pmsimstats team

**Purpose.** The stress test (`11-strawman-stress-test.R`) and the
covariance study (`08-covariance-study.R`) were run with the published
recursive carryover on the drug-response (BR) mean and, in the stress
test, a fixed near-equal path allocation. Paper 01 (Table `tab:spec`,
Appendix B.3) and `analysis/report/NOTATION.md` ("Named
data-generating processes") define the refined process with anchored
carryover and randomized path allocation, and `docs/50` records that
the refined process "has not itself been simulated". This note lists
the code changes, the commands for the rerun (locally and on the
simulation server through `simrun`), the outputs, and the documents to
update afterwards. No job has been submitted and nothing has been
committed.

## 1. What changed and why

**The defect.** Hendrickson's code (commit 58b32a9) and the drivers
computed the off-drug BR mean recursively,
`br[p] <- br[p] + br[p - 1] * 2^(-tsd[p] / t_half)`. Because `tsd` is
the cumulative time since discontinuation and `br[p - 1]` has already
been decayed, the elapsed time is counted twice from the second
off-drug visit of a run onward. The anchored form, which
`R/generateData.R` (lines 277 to 290, `last_on`) already uses, is
`mu_BR,last * 2^(-tsd / t_half)`.

Verified (`check-anchored.R` in the session scratchpad; Hybrid path 1,
on-drug pattern 11110010, mean at discontinuation 10.187):

| $t_{1/2}$ | visit (week, tsd) | anchored | recursive | ratio |
|---|---|---|---|---|
| 0.5 | week 11, tsd 1 | 2.5468 | 2.5468 | 1 |
| 0.5 | week 12, tsd 2 | 0.6367 | 0.1592 | 0.25 |
| 1.0 | week 11, tsd 1 | 5.0935 | 5.0935 | 1 |
| 1.0 | week 12, tsd 2 | 2.5468 | 1.2734 | 0.50 |

On path 3 (11100010) the third off-drug visit at $t_{1/2} = 1$ is
1.2242 anchored against 0.1530 recursive. In every case the anchored
means equal `mu_last * 2^(-tsd / t_half)` exactly, in both the 04 and
the 11 generator.

**Changes (all switches default to the refined process).**

- `hendrickson-problems/04-power-simulation.R`: `build()` gains an
  argument `carryover = 'anchored'` (the logic stays inside `build()`
  because other scripts parse that function alone); new flags
  `--carryover anchored|recursive` and `--out-dir`; anchored runs add
  `-anchored` to the cell directory and output names; workers default
  to `SIM_CORES` (8 otherwise).
- `carryover-closed-form/08-covariance-study.R`: passes the carryover
  choice to 04's `build()`; flags `--carryover` and `--out-dir`;
  anchored runs carry `-anchored` at the end of the run tag (for
  example `summary-alt250-null500-dbc-anchored.csv`); `--cores`
  defaults to `SIM_CORES`. Allocation stays fixed (18/18/17/17,
  35/35); randomizing it was not in scope.
- `carryover-closed-form/11-strawman-stress-test.R`: anchored carryover
  (applied with the cell's exponential or Weibull weight `w`), and
  randomized allocation: in every trial each participant is assigned to
  one of the four Hybrid paths uniformly at random (`alloc()`), redrawn
  in the rare case that a path receives fewer than 2 participants,
  because the MAR dropout model standardizes Y within a path (redraw
  probability 0.0021 at $N = 35$, below 0.0001 at $N = 70$; simulated).
  The positive-definite pre-check per path is unchanged. Flags
  `--carryover anchored|recursive`, `--allocation random|fixed`,
  `--out-dir`. Tags: defaults `-anchored`; `-anchored-fixed`,
  `-recursive-random`; `--carryover recursive --allocation fixed` gives
  the untagged names of the stored run.
- `12-strawman-stress-decision.R`: `--tag T` reads
  `replicates-full{T}.rds` (and the `-extra{T}`, `-extra2{T}` runs) and
  appends `T` to every output; `--run` and `--dir` support smoke tests.
  Without flags it reads and writes exactly the stored files.
- `09-covariance-study-tables.R`: `--tag`, `--reps ALT,NULL`, `--dir`;
  defaults unchanged.
- `13-symmetry-check.R`, `14-posthoc-arm-summary.R`: take the same DGP
  switches or tag, write tagged outputs.
- `hendrickson-problems/12-common-grid-contrast.R`: now passes
  `'recursive'` to `build()` explicitly, because it checks the 04
  means against 58b32a9's `generateData()`.

**Backward compatibility (verified).** A two-replicate run of
`11 --pilot 2 --extra --carryover recursive --allocation fixed`
reproduces the stored `replicates-pilot2-extra.rds` exactly (234 of
234 rows matched, maximum absolute difference in estimate and p-value
0).

**Other consumers of 04's `build()` (not rerun, flagged).** These read
`build()` from 04 and now get anchored means by default if rerun, and
they write untagged files, so a rerun would replace stored outputs
with anchored ones: `01-e9-closed-form.R` (study 1 of paper 01: the
closed-form mean of the paired-difference statistic uses `b$mu`),
`06-single-path-lme.R`, `07-time-adjustment.R`,
`10-dbc-preexposure-check.R`, `hybrid-design-primer/01-figures.R`, and
`compendium-mathematics/check-closed-form-extensions.R`. Scripts that
use only the correlation matrix, or only $t_{1/2} = 0$, are unaffected
(`11-pd-repair-count.R`, `13-summarize-contrast.R`,
`check-prognostic*.R`, `02-baseline-propositions.R`). Pass
`carryover = 'recursive'` in those calls to reproduce stored values.

## 2. Smoke tests (local, 4 workers)

| run | replicates | wall time | result |
|---|---|---|---|
| 11 main, `--pilot 2` | 82 | 116 s | 902 fits, 0 failed |
| 11 `--extra`, `--pilot 2` | 82 | 16 s | 246 fits, 0 failed |
| 11 `--extra2`, `--pilot 2` | 82 | 26 s | 164 fits, 0 failed |
| 12 on the smoke output | | 0.3 s | four `decision-*-ds-nosym-all14-pilot2-anchored.csv` |
| 04 `--analysis car1 --arms C`, 2 reps | 120 | 21 s | 100 fits, 0 failed (10 cells infeasible, as before); `summary-reps2-car1-anchored-C.csv` |
| 08 main, 2 + 2 reps | 32 | 35 s | 640 fits, 0 failed |
| 08 `--arms dbc`, 2 + 2 reps | 32 | 427 s | 640 fits, 0 failed |
| 08 `--arms followup`, 2 + 2 reps | 32 | 454 s | 736 fits, 0 failed |
| 09 on the 08 smoke output | | 0.8 s | `summary-alt2-null2-combined-anchored.csv`, `size-adjusted-power-alt2-null2-anchored.csv` |

In 08 (fixed allocation, same seeds) a recursive and an anchored
two-replicate run give identical mean estimates in the $t_{1/2} = 0$
cells and differ by up to 0.0021 in the $t_{1/2} = 1$ cells, as
expected, since only the off-drug means at $t_{1/2} > 0$ change.

The S estimate is identical replicate for replicate between the main,
`--extra` and `--extra2` smoke runs (82 pairs, difference 0), so
randomized allocation keeps the post hoc arms paired with the main
run. With two replicates no analysis survives the screens; that is
expected and only shows that 12 runs end to end.

## 3. Estimated run times

Per-replicate cost is taken from the stored full runs (8 local
workers on a 4 performance plus 4 efficiency core machine) and checked
against the smoke runs, which used 4 performance cores and so ran about
twice as fast per worker. Anchoring changes means only, so the fitting
cost should not change; randomized allocation adds a negligible draw.

| run | replicates | worker-seconds per replicate | 8 local workers | 3 cores (simrun default) | 8 cores |
|---|---|---|---|---|---|
| 11 main | 28,000 (41 cells: 1,000 per null, 500 per alt) | about 11 | about 11 h | about 30 h | about 11 h |
| 11 `--extra` | 28,000 | about 1.3 | about 1.3 h | about 3.5 h | about 1.3 h |
| 11 `--extra2` | 28,000 | about 2.1 | about 2.1 h | about 6 h | about 2.1 h |
| 08 main | 6,000 (16 cells, 250 alt / 500 null) | about 14 | about 3 h | about 8 h | about 3 h |
| 08 `--arms followup` | 6,000 | about 120 | about 25 h | about 67 h | about 25 h |
| 08 `--arms dbc` | 6,000 | about 110 | about 23 h | about 61 h | about 23 h |

These assume the server's cores run at the speed of this machine's
average worker (assumed, not measured: no job on mercury has run these
drivers; the earlier jobs ran `01-run-factorial.R` only). The 08 driver
defaults are 500 alternative and 1,000 null replicates; the stored
results cited in `docs/45` use 250 and 500, so the commands below pass
those values to keep the comparison like for like.

## 4. Commands

All commands run from the repository root
(`/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng`). Defaults give
the refined process, so no DGP flag is needed.

### 4.1 Locally

```bash
S=analysis/scripts/quick-sim/carryover-closed-form
O=analysis/data/quick-sim/carryover-closed-form

Rscript $S/11-strawman-stress-test.R --check
Rscript $S/11-strawman-stress-test.R --cores 8 \
  > $O/strawman-stress/run-full-anchored.log 2>&1
Rscript $S/11-strawman-stress-test.R --extra --cores 8 \
  > $O/strawman-stress/run-full-extra-anchored.log 2>&1
Rscript $S/11-strawman-stress-test.R --extra2 --cores 8 \
  > $O/strawman-stress/run-full-extra2-anchored.log 2>&1

Rscript $S/12-strawman-stress-decision.R --tag -anchored
Rscript $S/12-strawman-stress-decision.R --strawman DS --no-symmetry \
  --tag -anchored
Rscript $S/12-strawman-stress-decision.R --strawman DS --no-symmetry \
  --with-extras --tag -anchored
Rscript $S/14-posthoc-arm-summary.R extra -anchored
Rscript $S/14-posthoc-arm-summary.R extra2 -anchored

Rscript $S/08-covariance-study.R --reps-alt 250 --reps-null 500 \
  --chunk 50 --cores 8 > $O/covariance-study/run-alt250-null500-anchored.log 2>&1
Rscript $S/08-covariance-study.R --reps-alt 250 --reps-null 500 \
  --chunk 50 --cores 8 --arms followup \
  > $O/covariance-study/run-alt250-null500-followup-anchored.log 2>&1
Rscript $S/08-covariance-study.R --reps-alt 250 --reps-null 500 \
  --chunk 50 --cores 8 --arms dbc \
  > $O/covariance-study/run-alt250-null500-dbc-anchored.log 2>&1
Rscript $S/09-covariance-study-tables.R --tag -anchored
```

Optional: `Rscript $S/13-symmetry-check.R` writes
`strawman-stress/symmetry/summary-reps1000-anchored.csv`.

### 4.2 On the simulation server (`simrun`)

How submission works (inspected: `~/dotfiles/bin/simrun`,
`~/dotfiles/bin/simjob`, `~/dotfiles/share/simrun/vignette.md`, and the
three records in `analysis/remote_runs/job-*/provenance.yml`):

- `simrun submit [--cores N] [--mem G] [--label TXT] -- CMD...` pushes
  a commit to the host (`mercury`; the SSH alias comes from
  `~/.config/simrun/config`), and the host runs `CMD` in the project's
  image (`localhost/simrun-res--36-pmsimstats-ng--pmsimstats-ng:<hash>`,
  the hash of `Dockerfile` plus `renv.lock`; the earlier jobs reused
  `10504a44fb8d`) from that commit's tree, with the project's data
  store overlaid. The job runs the committed tree, not the working
  tree: **`simrun submit` refuses while tracked files have uncommitted
  changes**, and `--commit HEAD` runs the last commit without them. The
  driver changes listed above must therefore be committed first; the
  working tree currently also holds many unrelated modified tracked
  files, which must be committed or stashed too before a plain
  `simrun submit` is accepted (untracked files do not block it). The
  stored `strawman-stress/decision-*-ds-nosym-all14.csv` files that
  `docs/47` cites are untracked at present (`git status`, 2026-10-08)
  and are worth committing before the rerun so that the old decision
  is on record.
- `CMD` is run as an argument vector, not through a shell, so it
  cannot contain redirections or `&&`; its output goes to the job log
  (`simrun log ID`). One job per command is the simplest form.
- The job sets `SIM_CORES` to its core share; the changed drivers use
  it as their default worker count, so `--cores` need not be passed to
  the R script.
- Every file the run creates or changes is collected as a result.
  `simrun fetch ID` copies results, log and provenance to
  `analysis/remote_runs/job-ID-SHA/` (gitignored), and
  `simrun apply ID` copies the results into the project, refusing to
  overwrite a differing file. The new runs write only new names, so
  `apply` will not touch the stored results.
- Jobs do not see each other's outputs. The decision and table scripts
  (12, 14, 09) therefore run locally after `fetch` and `apply`; they
  take seconds.

Suggested submissions (cores and memory are recommendations; the
earlier jobs used the defaults of 3 cores and 12 GB; each worker holds
one trial and its fits, a few hundred MB, so about 2 GB per core is
ample):

```bash
S=analysis/scripts/quick-sim/carryover-closed-form

simrun submit --cores 2 --mem 4 --label anchored-check -- \
  Rscript $S/11-strawman-stress-test.R --check
simrun submit --cores 8 --mem 16 --label stress-anchored -- \
  Rscript $S/11-strawman-stress-test.R
simrun submit --cores 8 --mem 16 --label stress-anchored-extra -- \
  Rscript $S/11-strawman-stress-test.R --extra
simrun submit --cores 8 --mem 16 --label stress-anchored-extra2 -- \
  Rscript $S/11-strawman-stress-test.R --extra2
simrun submit --cores 8 --mem 16 --label cov-anchored -- \
  Rscript $S/08-covariance-study.R --reps-alt 250 --reps-null 500 --chunk 50
simrun submit --cores 8 --mem 16 --label cov-anchored-followup -- \
  Rscript $S/08-covariance-study.R --reps-alt 250 --reps-null 500 --chunk 50 \
  --arms followup
simrun submit --cores 8 --mem 16 --label cov-anchored-dbc -- \
  Rscript $S/08-covariance-study.R --reps-alt 250 --reps-null 500 --chunk 50 \
  --arms dbc
```

Then, for each job: `simrun watch ID`, `simrun fetch ID`,
`simrun apply ID`, and run the 12, 14 and 09 commands of Section 4.1.
The number of cores the host can give at once is not known from the
repository (unverified); if fewer than 8 per job are available, the
times in Section 3 scale inversely with the cores. Checkpoint chunks are
written inside the job's own run directory, so a failed job is not
expected to resume from them (inferred from `simjob`, which unpacks a
fresh tree per job; not tested); `simrun resubmit ID` queues the same
commit and command again.

## 5. Outputs

Stress test, in `analysis/data/quick-sim/carryover-closed-form/strawman-stress/`:

- `cells-full-anchored/`, `cells-full-extra-anchored/`,
  `cells-full-extra2-anchored/` (chunk files and `progress.log`);
- `replicates-full-anchored.rds`, `replicates-full-extra-anchored.rds`,
  `replicates-full-extra2-anchored.rds`;
- `summary-full-anchored.csv`, `summary-full-extra-anchored.csv`,
  `summary-full-extra2-anchored.csv`;
- from 12: `decision-{screen,sensitivity,power,maximin}-anchored.csv`,
  `...-ds-nosym-anchored.csv`, `...-ds-nosym-all14-anchored.csv`;
- from 14: `posthoc-extra-anchored-{size,power}.csv`,
  `posthoc-extra2-anchored-{size,power}.csv`.

Covariance study, in `.../covariance-study/`:

- `cells-alt250-null500-anchored/`, `...-followup-anchored/`,
  `...-dbc-anchored/`;
- `replicates-` and `summary-alt250-null500-anchored`,
  `-followup-anchored`, `-dbc-anchored`;
- from 09: `summary-alt250-null500-combined-anchored.csv`,
  `size-adjusted-power-alt250-null500-anchored.csv`.

## 6. Documents to update after the rerun

- `docs/45` (covariance study): size, power and $\kappa$ tables, and
  the statement of the data-generating process (now anchored carryover;
  allocation still fixed in 08).
- `docs/47` (stress-test results): screens, survivors, maximin
  shortfall, McNemar comparisons, power curves and the AIC selection
  shares; state the two DGP changes as a protocol change against
  `docs/46`.
- Paper 01: the Results paragraphs that quote stress-test numbers
  ("Consequences for the analysis" and the prognostic paragraphs near
  Appendix A.4), and the ADEMP table `tab:ademp`, study 3 (the
  data-generating mechanism row now matches `tab:spec`; the replicate
  row is unchanged). If `01-e9-closed-form.R` is rerun, study 1 values
  at $t_{1/2} > 0$ change as well (Section 1).
- `analysis/report/compendium-mathematics.tex`: the sections that cite
  stress-test rates and the stress-test construct (search for
  `stress`, `docs/48, Section 2` and `11-strawman-stress-test.R`; the
  note near line 6560 that the stress test did not run the refined
  process in full).
- `docs/50`, "Open items": remove the statement that the refined
  process has not been simulated, and replace the cited numbers.
- `docs/48` and `docs/49` where they quote stress-test rates.

## 7. Comparing old and new results

The cells, analyses and seeds are the same, but the trials are not
paired across the two runs: the means differ, and randomized
allocation consumes random numbers, so the old and new rejection rates
are independent estimates. Compare them as independent proportions
(difference over $\sqrt{\text{MCSE}_\text{old}^2 +
\text{MCSE}_\text{new}^2}$). A sketch:

```r
library(data.table)
d <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress'
old <- fread(file.path(d, 'summary-full.csv'))
new <- fread(file.path(d, 'summary-full-anchored.csv'))
k <- c('cell', 'label', 'cbm', 'hyp', 'analysis')
cmp <- merge(old[, c(k, 'reject', 'mcse'), with = FALSE],
             new[, c(k, 'reject', 'mcse'), with = FALSE],
             by = k, suffixes = c('_old', '_new'))
cmp[, z := (reject_new - reject_old) / sqrt(mcse_old^2 + mcse_new^2)]
head(cmp[order(-abs(z))], 20)
```

The decision files compare directly (`decision-screen-ds-nosym-all14.csv`
against `...-all14-anchored.csv`, and so on). To separate the two
changes, the stress test can also be run with
`--allocation fixed` (tag `-anchored-fixed`, anchored carryover only),
which isolates the effect of the carryover correction. Expected
direction (inferred, not tested): anchoring raises the off-drug BR
mean at the second and later off-drug visits, so the misfit of a
binary or $D_{bc}$ drug term at those visits changes most in the
$t_{1/2} = 0.5$ and 1 week cells; cells with $t_{1/2} = 0$ are
unchanged in their means.
