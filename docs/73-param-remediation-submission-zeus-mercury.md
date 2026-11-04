---
title: 'Parameter remediation for papers 01 and 02: submission plan on zeus and mercury'
author: 'pmsimstats team'
---

# Parameter remediation for papers 01 and 02: submission plan on zeus and mercury
*2026-11-03 17:06 PST*

**Superseded by `docs/74`** (2026-11-03 17:51), which extends this plan
to every affected paper and adds athena and joe as hosts. Section 2's
inference that zeus runs at a quarter of julia's per-core speed was
wrong; the measured benchmark is in `docs/74` Section 2 and the
correction below.

## 1. Purpose

This plan submits every simulation that papers 01 and 02 need under
the published parameter set (`docs/70` Phases 1 and 2, paper 10
excluded) to the lab hosts zeus and mercury. It replaces the host
choice of `docs/71` (SDSC Expanse); the paper 02 job list, the code
changes and the downstream steps of `docs/71` still apply. No
allocation is spent.

## 2. What the timing record shows

Wall-clock times below are measured on zeus, not projected, wherever a
stored job ran the same work (`simrun status`, `provenance.yml`).

- **Driver 19 on zeus.** Block P with one analysis took 5h26m (job
  40, G9) and 5h33m (job 53, G10); block F took 32 minutes (jobs 41
  and 54); all on 8 cores. Block P is 72,000 (cell, replicate) jobs,
  so zeus spends about 2.2 core-seconds per replicate and analysis.
- **The local timing mode underestimates an 8-worker zeus job about
  fourfold.** `19-g9-planning-and-fig4b.R --time` gives about 0.53
  core-seconds per replicate and analysis. Corrected 2026-11-03: the
  cause is not per-core speed. One zeus worker matches julia natively
  on light cells and is about 35 percent slower on the heaviest, but
  zeus keeps only about a third of its single-worker speed when 8
  workers run together (measured, `docs/74` Section 2). Projections
  from `--time` still need multiplying by about 4 for an 8-worker zeus
  job.
- **Consequence for `docs/71`.** Withdrawn: the earlier version of this
  bullet concluded that `docs/71`'s Expanse budget was about four times
  too high, from the per-core reading above. Expanse's 8- or 16-worker
  throughput on this code has not been measured, so the budget stands
  as an unmeasured estimate.
- **Mercury** ran gap-fill block C in 13.0 hours on 6 cores (job 49);
  the same block is about 40 zeus core-hours, so mercury with 6 cores
  delivers about 0.38 of zeus's 8-core throughput.

## 3. The jobs

Hours are zeus wall-clock on 8 cores. Paper 02 hours are the `docs/71`
projections, which come from the stored zeus jobs plus about 20
percent for G10 where it was not fitted before. Paper 01 hours assume
that fitting G9 and G10 in one job costs the sum of the two stored
single-analysis jobs (the data draw is shared and small).

| Label | Paper | Driver and part | Zeus hours |
|---|---|---|---|
| p01-pf | 01 | driver 19, blocks P and F, G9 and G10 | 11.9 |
| p01-n80 | 01 | driver 19, block N, G9 and G10 | 4.7 |
| p02-factorial | 02 | driver 01, covariance moderation only | 2.3 |
| p02-tier1 | 02 | driver 27 | 4.9 |
| p02-e0iso | 02 | driver 38 | 3.3 |
| p02-scen | 02 | driver 06, the ten cells paper 02 reads | 2.4 |
| p02-gf-dropout | 02 | driver 32, dropout cells, both conventions | 7.2 |
| p02-gf-size | 02 | driver 32, size cells | 2.4 |
| p02-gf-B | 02 | driver 32, block B | 6.1 |
| p02-ref-N | 02 | driver 36, block N | 3.5 |
| p02-ref-A | 02 | driver 36, block A | 0.7 |
| p02-ref-S | 02 | driver 36, block S | 0.2 |
| **Total** | | | **49.6** |

The full commands are in Section 6.

Block N's 4.7 hours scale block P's measured time by the ratio of the
two blocks in the `--time` projection (9.2 against 21.4); it has not
run on zeus before. Every command also takes `--params published`.

**Already done, not resubmitted.** The paired rerun of paper 01
Appendix B.4 (`26-co-dbc-fix-paired.R --params published`) ran locally
on julia on 2026-11-03 in 26 minutes, in the host R session rather
than the project container. Nothing is paired with it. Optional: rerun
it on mercury in the container for provenance (about 4 hours, Section
5), writing to a separate folder.

## 4. Host assignment and pairing

Each host runs one job at a time from its own queue. Zeus takes the
paired jobs and paper 01's critical path; mercury takes jobs that have
no pairing requirement.

- **Job 27 and the E0 isolation on one host** (`39-check-e0-isolation.R`
  reproduces part of job 27, and the ladder reads both): both on zeus.
- **Gap-fill blocks A and B on one host and one commit**
  (`33-assemble-gap-fill.R` builds one table from them): all three on
  zeus.
- **Referee blocks, the scenario cells and paper 01 block N** are
  unpaired with anything else: mercury.
- **Paper 01 blocks P and F** fit G9 and G10 to the same trials within
  one job, so the G9 and G10 comparison holds on any host: zeus, first,
  because paper 01 waits on it.
- **One commit and one image** for every job. A commit that touches
  `simulation-core.R` or a driver after the first submission
  invalidates the jobs already queued.

| Host | Queue order | Hours on the host |
|---|---|---|
| zeus (8 cores) | p01-pf, p02-factorial, p02-tier1, p02-e0iso, p02-gf-dropout, p02-gf-size, p02-gf-B | about 38 |
| mercury (6 cores) | check job, p02-ref-A, p01-n80, p02-ref-S, p02-ref-N, p02-scen | about 31 (11.5 zeus-hours at 0.38) |

Both queues finish in about a day and a half of wall-clock time if no
job hangs. Zeus alone would take about 50 hours.

## 5. Before submitting

1. **Commit.** simrun runs commits only. The commit holds
   `simulation-core.R`; drivers 01, 06, 19, 26, 27, 32, 36 and 38;
   `42-check-param-sets.R`; docs/71, docs/72 and this plan; and the
   paper 01 revision (`report.Rmd`, `references.bib`, the figure
   scripts and the files the referee pass added). Other uncommitted
   work in the repository (the paper 02 text and ladder edits,
   `spec-labels.R`) stays out unless the authors decide otherwise;
   with it uncommitted, every submit command needs `--commit HEAD`.
2. **Checks at the commit.** The parameter-set checks
   (`42-check-param-sets.R`, all checks, about 4 minutes) and the
   local smoke of every driver
   under `--params published` (scratchpad script `smoke-p02.sh` and
   the driver 19 smoke; under 10 minutes together). The tinytest suite
   has not been run; the package code is unchanged.
3. **Image.** `Dockerfile` and `renv.lock` are unchanged since
   2026-07-06, so the image both hosts used for jobs 40 to 56
   (`...pmsimstats-ng:10504a44fb8d`) is reused; no build and no
   `--ship`.
4. **Mercury check job.** Mercury jobs 41 to 47 failed at once (exit
   127) before jobs 48 and 49 succeeded. A one-minute job confirms the
   route before the real queue:

```bash
simrun --host sim-mercury-ts submit --commit HEAD --cores 1 --mem 2 \
  --label p0x-mercury-check -- \
  Rscript analysis/scripts/carryover-sensitivity/42-check-param-sets.R \
    --reps 0
```

## 6. Submission commands

Run from the repository root. Zeus is the default host, so its
commands carry no `--host`. Driver 19 and the paper 02 drivers size
their workers from `SIM_CORES`, which simrun sets from `--cores`.

The three variables below hold single paths only. The simrun options
are written out in every command: zsh, the lab's shell, does not split
an unquoted variable into words, so a variable holding several options
would reach simrun as one argument.

```bash
D=analysis/scripts/carryover-sensitivity
D19=analysis/scripts/quick-sim/carryover-closed-form
D19=$D19/19-g9-planning-and-fig4b.R
P=analysis/data/carryover-sensitivity/published
```

**Zeus, in this order** (each queues behind the previous one):

```bash
simrun submit --commit HEAD --cores 8 --mem 24 --label p01-pf -- \
  Rscript $D19 --params published --specs E7null_cr2,E7cr2 --block all
simrun submit --commit HEAD --cores 8 --mem 24 --label p02-factorial -- \
  Rscript $D/01-run-factorial.R --params published --arch mvn \
  --out-dir $P/factorial
simrun submit --commit HEAD --cores 8 --mem 24 --label p02-tier1 -- \
  Rscript $D/27-run-tier1-g-vs-f-e0.R --params published \
  --reps 500 --n-boot 1 --out-dir $P/tier1-g-vs-f-e0
simrun submit --commit HEAD --cores 8 --mem 24 --label p02-e0iso -- \
  Rscript $D/38-run-e0-isolation.R --params published --reps 1000 \
  --out-dir $P/e0-isolation
simrun submit --commit HEAD --cores 8 --mem 24 --label p02-gf-dropout -- \
  Rscript $D/32-run-g9-gap-fill.R --params published --block A \
  --sub A1-dropout --dropout-per both --out-dir $P/g9-gap-fill/dropout
simrun submit --commit HEAD --cores 8 --mem 24 --label p02-gf-size -- \
  Rscript $D/32-run-g9-gap-fill.R --params published --block A \
  --sub A2-size --out-dir $P/g9-gap-fill/size
simrun submit --commit HEAD --cores 8 --mem 24 --label p02-gf-B -- \
  Rscript $D/32-run-g9-gap-fill.R --params published --block B \
  --out-dir $P/g9-gap-fill/B
```

**Mercury, in this order**, after the check job succeeds:

```bash
simrun --host sim-mercury-ts submit --commit HEAD --cores 6 --mem 12 \
  --label p02-ref-A -- \
  Rscript $D/36-run-referee-variants.R --params published --block A \
  --out-dir $P/referee-variants/block-A
simrun --host sim-mercury-ts submit --commit HEAD --cores 6 --mem 12 \
  --label p01-n80 -- \
  Rscript $D19 --params published --specs E7null_cr2,E7cr2 --block N
simrun --host sim-mercury-ts submit --commit HEAD --cores 6 --mem 12 \
  --label p02-ref-S -- \
  Rscript $D/36-run-referee-variants.R --params published --block S \
  --out-dir $P/referee-variants/block-S
simrun --host sim-mercury-ts submit --commit HEAD --cores 6 --mem 12 \
  --label p02-ref-N -- \
  Rscript $D/36-run-referee-variants.R --params published --block N \
  --out-dir $P/referee-variants/block-N
simrun --host sim-mercury-ts submit --commit HEAD --cores 6 --mem 12 \
  --label p02-scen -- \
  Rscript $D/06-g-vs-f-comparison.R --params published \
  --reps 500 --n-boot 1 --cells paper02 --out-dir $P/g-vs-f-comparison
```

Driver 19 needs no `--out-dir`: under the published set it writes to
`analysis/data/quick-sim/carryover-closed-form/published-planning/`
(blocks P and F) and to its `block-N/` subfolder (block N), which is
where paper 01 reads them. The paper 02 folders follow the stored
layout under `published/`, with the three gap-fill jobs in separate
subfolders so that `simrun apply` cannot collide.

## 7. Monitoring

```bash
simrun status
simrun --host sim-mercury-ts status
simrun log ID -f
```

The phone notification reports each job's end. Zeus and mercury
ignore `--time`, and zeus has stalled before (jobs 20, 39, 47 and 52
ran 125 to 263 hours on work that takes hours). Check each job against
its expected time in Section 3 and cancel any job running at more than
twice its expected time (`simrun cancel ID`, or
`simrun --host sim-mercury-ts cancel ID`), then resubmit it.

## 8. Fetch, apply and check

For each finished job:

```bash
simrun fetch ID          # zeus
simrun apply ID
simrun --host sim-mercury-ts fetch ID       # mercury
simrun --host sim-mercury-ts apply ID
```

Then, for each job:

- `meta.rds` records `params$param_set == 'published'`, and the
  output folder's `PARAM_SET` file says `published`.
- The log reports every (cell, replicate) job present and the
  pre-flight eigenvalues positive.
- Null-cell rejection rates fall inside the two-MCSE band;
  non-convergence is comparable to the stored run's.
- Paper 01 block N: count the cells flagged `extrapolated` in
  `n80-g10.csv`; if they lie below the lowest multiplier, rerun block
  N with a wider `--n80-mult`.

## 9. After the runs

**Paper 01.** Render
`analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`. Every
[pending] marker should disappear. Read the computed G10 sentences
with the real numbers: they were tested only against four-replicate
smoke output.

**Paper 02** (`docs/71` Section 8). Repoint the report's `readRDS()`
paths to `published/`; give `33-assemble-gap-fill.R` the three new
gap-fill job folders (zeus job numbers, no host suffix); run
`39-check-e0-isolation.R` against the new tier-1 output; remove the
G10 backfill block; rerun the ladder, referee reanalysis and figure
scripts; state the parameter table in Methods; rewrite the Results
from the new numbers; lint and render.

## 10. Evidence status

- **Measured:** the zeus and mercury durations of the stored jobs; the
  driver 19 replicate count per block; the image tag reused on both
  hosts; both hosts' queues empty on 2026-11-03.
- **Inferred:** (withdrawn: zeus at a quarter of julia's per-core
  speed; see Section 2) mercury at 0.38 of zeus's throughput (one job); the cost of
  two analyses in one job as the sum of two single-analysis jobs; the
  G10 overhead in paper 02's combined jobs; block N's duration.
- **Not checked:** mercury's core count (its gate refuses shell
  commands; jobs 40 and 49 ran with 6 cores); memory per worker at the
  published parameters (12 GB for 6 workers sufficed for block C).
