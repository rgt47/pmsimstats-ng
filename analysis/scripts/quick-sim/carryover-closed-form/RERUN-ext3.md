# Stress-test extension 3: data-chosen analysis half-lives
*2026-10-09 16:30 PDT*

Author: pmsimstats team

**Purpose.** `docs/57` sets out procedures for choosing the analysis
half-life $t_{1/2}^a$ of the biomarker-by-drug interaction test from the
data, and illustrates them on four stored cells (250 to 500 trials
each, stored recursive-carryover setting). Its Section 7.3 asks for a
stress-test extension at full scale under the refined process. This
note describes that extension as implemented in
`11-strawman-stress-test.R --extra3`: what it tests, the analyses and
cells, the commands (local and `simrun`), the cost, the outputs, how
`12-strawman-stress-decision.R` reads them, and the documents to
update afterwards. Nothing has been committed or submitted. The study is
exploratory. Several of its arms were defined after earlier results
were seen, and they are labeled post hoc below.

## 1. What the extension tests

- **Size of the data-chosen rules at scale.** In the illustration,
  F1A0 removed only about half of F1A's excess, AIC-weight averaging
  and the bootstrap max-$|t|$ held size, and the post hoc guard F1G
  behaved as F1. Those rates carry an MCSE of 0.007 to 0.014. The
  extension estimates them over all 15 null cells at 1,000 trials each,
  plus new null cells.
- **Prognostic biomarker combined with a wrong half-life.** No stored
  cell combined the two (`docs/57`, Sections 3.3 and 6.3). The new
  cells 42 to 53 do.
- **Power where the fixed value is wrong.** The illustration's only
  alternative cell had the true half-life equal to the fixed 0.5 week.
  The extension keeps cells 31 and 32 ($t_{1/2} = 0$, 1 at
  $c_{bm} = 0.30$) and adds $t_{1/2} = 2$ at 0.30 and $t_{1/2} = 0$, 1, 2
  at 0.45.
- **The fixed value itself.** F1 is fitted at 0.5 week, the value of
  the existing results, and also at the pharmacological 0.2 week.

Every trial is new. The base seed is 57261009 in place of 20261007, so
the extension does not reuse the trials of the main run, the anchored
rerun or the illustration. Rules defined after those results were seen
(F1G, MA0) are therefore judged here on trials their definition did
not see. Within the extension, every analysis is fitted to the same
trials, so comparisons between arms are paired.

## 2. Analyses

All are fitted per trial on the F1 mean model
`Sx ~ bmc + t + X + bmc:X + De + bmc:De` with a random intercept,
`corCAR1(~ t | ptID)` and CR2/Satterthwaite inference for `bmc:X`,
over the grid $\mathcal T_a = \{0, 0.25, 0.5, 1, 2\}$ weeks (ML for
the log-likelihoods, REML for the tests). The grid fits, F1, F1A,
F1A0, MA, the Bonferroni maximum and the null-model resample are the
definitions of `analysis/scripts/data-chosen-half-life/57-common.R`.
The driver evaluates those definitions from that file, so they are the
same code, not a copy. F1G and MA0 follow `57-posthoc-rules.R` line
for line.

| Analysis | Definition | Status |
|---|---|---|
| `F1` | fixed $t_{1/2}^a = 0.5$ week | as in `docs/46`, `docs/47` |
| `F1h02` | fixed $t_{1/2}^a = 0.2$ week (pharmacological value) | new |
| `F1A` | AIC over the grid with the interaction in the model; negative control | as in `docs/47` |
| `F1A0` | AIC over the grid without the interaction; REML refit with it | `docs/57`, Section 4.2 |
| `F1G` | keep 0.5 unless $2\{\max_h\ell_0(h) - \ell_0(0.5)\} > 3.84$, then the F1A0 choice | **post hoc** (`docs/57`, Section 4.2) |
| `MA` | AIC-weight averaging over the grid, $\ell_1$ weights | `docs/57`, Section 4.8 |
| `MA0` | the same with $\ell_0$ weights | **post hoc** (`docs/57`, Section 6.2) |
| `MaxTBonf` | $\max_h\lvert T(h)\rvert$, Bonferroni over the five | `docs/57`, Section 4.7 |
| `MaxTMVN` | $\max_h\lvert T(h)\rvert$ referred to the maximum of a multivariate normal whose correlation comes from the participants' model-based influence contributions to $\hat\beta_{bm:D}(h)$ (`docs/58`, Section 9.2; `infl_X` of `58-common.R` with target `bmc:X`) | new |
| `MaxTMVT` | the same with a multivariate $t$ on the smallest Satterthwaite df, rounded down (`pmvt` needs an integer) | new |
| `F1B`, `F1h02B`, `F1AB`, `F1A0B`, `F1GB`, `MAB`, `MA0B`, `MaxTB` | warp-speed bootstrap (Giacomini, Politis and White 2013) of the corresponding statistic | `F1AB`, `MaxTB` as in `docs/57`, Section 4.5; the rest are references (`F1GB`, `MA0B` post hoc) |

Table: Analyses of the `--extra3` arm (column `posthoc` is TRUE for F1G, MA0, F1GB and MA0B)

**Warp-speed bootstrap.** For each trial, the null model FM0 is refitted
by REML at the F1A0 choice $\hat h_0$, and one resample $Y^*$ is drawn
from it (`sim_null`). Every procedure, including its own selection, is
rerun on $Y^*$, and the resampled statistic is stored as `tstar` on
that procedure's row. After the run, `warp_speed()` pools the `tstar`
of each cell and analysis and gives each trial
$p^* = (1 + \#\{T^* \ge \lvert T_{\text{obs}}\rvert\})/(n + 1)$. These
are the rows with suffix `B`. The cell's 95% quantile is stored as
`crit_ws`. `docs/57` rejected when $\lvert T\rvert$ exceeded that
quantile; the $p^*$ rule is the same up to ties at the quantile. The
bootstrap draw and the multivariate-normal integrals use a separate
stream, seeded $571000000 + 10000 \times \text{cell} +
\text{replicate}$. The trial stream is restored afterwards, so the
trials do not depend on the analyses.

**Not included** (`docs/57`, Section 7.3 lists them): cross-fitting
(single split, average, Cauchy), and the profile-interval coverage of
the best-fitting grid value. Davies' bound for a continuous $h$ is not
computed either; the grid maximum uses the multivariate normal and $t$
above. F1G keeps the 0.5-week default as defined in `docs/57`; a
version guarded at 0.2 week would need 0.2 in the grid.

## 3. Cells

All 41 cells of `docs/46` (ids 1 to 41, unchanged), plus:

| id | label | $c_{bm}$ | true $t_{1/2}$ | prognostic pathway | reps |
|---|---|---|---|---|---|
| 42 | `prog-TV-const-thalf0` | 0 | 0 | $c_{bm,TV} = 0.3$ | 1,000 |
| 43 | `prog-TV-grow-thalf0` | 0 | 0 | $b_{TV} = 0.3$ | 1,000 |
| 44 | `prog-PB-thalf0` | 0 | 0 | $c_{bm,PB} = 0.3$ | 1,000 |
| 45 to 47 | the same three | 0.30 | 0 | as 42 to 44 | 500 each |
| 48 to 50 | `prog-TV-const-thalf1`, `prog-TV-grow-thalf1`, `prog-PB-thalf1` | 0 | 1 | as 42 to 44 | 1,000 each |
| 51 to 53 | the same three | 0.30 | 1 | as 42 to 44 | 500 each |
| 54 | `thalf2` | 0 | 2 | none | 1,000 |
| 55 | `thalf2` | 0.30 | 2 | none | 500 |
| 56 | `thalf0` | 0.45 | 0 | none | 500 |
| 57 | `thalf1` | 0.45 | 1 | none | 500 |
| 58 | `thalf2` | 0.45 | 2 | none | 500 |

Table: Cells added by `--extra3` (all other settings as the base cell: AR(1) at $\rho = 0.7$, $N = 70$, exponential decay, no dropout)

Cells 42 to 53 are the 12 cells of `docs/57`, Section 7.3: the
constant and growing natural-course pathways and the placebo pathway,
each at true $t_{1/2} = 0$ and 1 week, null and $c_{bm} = 0.30$.
Cells 54 to 58 add wrong-half-life alternatives beyond cells 31 and 32,
and the null twin that `thalf2` needs for size-adjusted power. Cell
54 is a non-prognostic null cell and therefore enters the pooled size
screen. They were added at the request that commissioned this
extension and are not in `docs/57`'s count. `--ext3-cells docs57`
drops them, writes tag `-extra3-docs57`, and leaves cells 1 to 53 and
their seeds unchanged. Every new matrix passes the positive-definiteness
check; the smallest eigenvalue, 0.042, is in cell 56 (verified, `--check`).

Totals: 58 cells, 22 null (12 non-prognostic; 10 prognostic, counting
`prog-BL`) and 36 alternative; 40,000 trials (22,000 null at 1,000,
18,000 alternative at 500). With `--ext3-cells docs57`: 53 cells,
37,000 trials.

## 4. Verification (local, 2026-10-09)

All verified by running code. Scripts were in the session scratchpad
and are not part of the repository. Threads were pinned with
`VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

- **Smoke run.** `--extra3 --pilot 2 --cores 4`: 116 trials, 2,088
  rows, **0 failed fits** (no NA p-value in any analysis), 6.4 minutes
  of wall time. In 1 of the 116 bootstrap resamples (cell 25, replicate
  1), the REML refit at the half-life chosen on the resample failed, so
  F1A and F1A0 have no `tstar` there. Such draws are dropped from the
  cell's pooled distribution. The `B` rows are 0 everywhere at 2
  replicates, as they must be, since $p^* \ge 1/3$ with two draws.
- **Same values as `docs/57`'s code.** On 10 trials (cells 2, 20, 45,
  50 and 56, replicates 1 and 2), the trials were regenerated through
  `57-common.R`'s `load_driver()`, analyzed with its `analyze_trial()`
  and with the F1G/MA0 code of `57-posthoc-rules.R`, and compared with
  the arm's output. h, estimate and p-value agree exactly (maximum
  absolute difference 0) for F1, F1A, F1A0, MA, F1G, MA0 and the
  Bonferroni maximum. The bootstrap statistic of the same resample,
  scored by `57-common.R`'s `boot_stats()`, also agrees exactly for F1,
  F1A, F1A0 and the maximum.
- **Existing arms unchanged.** The pre-edit driver (a copy of the
  working-tree file before this change) and the edited one were run
  with `--pilot 2 --cores 4` for the main, `--extra` and `--extra2`
  arms (defaults: anchored, random). The replicate files are
  `identical()` (902, 246 and 164 rows) and the summary CSVs are
  byte-identical. The edited main run resumed one chunk after a worker
  died under load, so checkpoint resumption also reproduces. The
  decision script, old and new, run three ways (defaults;
  `--strawman DS --no-symmetry`; plus `--with-extras`) on copies of the
  stored replicates gives 12 byte-identical CSVs and identical logs.
- **Decision script on the extension.** `12 --ext3 --no-symmetry
  --run pilot2 --tag -anchored` runs end to end on the smoke output.
- **Not tested.** The `--ext3-cells docs57` path was run with
  `--check` only: it builds cells 1 to 53 and fits one trial; no
  smoke run was made with it. The `simrun` image, `Dockerfile` plus `renv.lock`, which
  includes `mvtnorm` and `clubSandwich`, was not exercised; all runs
  were in the host R 4.6.1.

## 5. Cost

**Load.** This machine (Apple M1, 4 performance plus 4 efficiency
cores) was heavily loaded by unrelated R processes throughout, with
load averages from 12 to 22 (`uptime`, 16:01 to 16:26). Measured times
are therefore inflated, and the unloaded estimate is calibrated
indirectly.

| Measurement | Value |
|---|---|
| smoke run, per trial per worker (4 workers, load 12 to 19) | 13.0 s mean, 12.4 s median |
| one worker, same 8 base-cell trials, `--extra2` analysis | 2.77 s mean |
| one worker, same trials, `--extra3` analysis | 14.0 s mean |
| ratio, extra3 to extra2 | 5.05 (trial ratios 3.6 to 7.8) |
| stored full `--extra2` run, unloaded (`RERUN-anchored.md`, Section 3) | about 2.1 s |
| **estimated unloaded cost of `--extra3`** | **about 10.6 s per trial** |

Table: Timing of the `--extra3` arm (verified measurements; the unloaded figure is derived)

The ratio is mechanically plausible. Each trial runs the full grid
twice, on the data and on one resample: 15 fits and 5 CR2 computations
each, plus F1 at 0.2 week, against about 7 fits and 2 CR2 computations
for `--extra2`. It is also close to `docs/57`'s 7.1 worker-seconds per
replicate for a similar procedure (inferred). Weighted by the cell mix
($N = 150$ costs about twice $N = 70$, $N = 35$ about 0.6 times; cell-mix
factor 1.02):

| Cell set | Trials | Worker-hours | 3 cores | 8 cores | 16 cores |
|---|---|---|---|---|---|
| all 58 cells (default) | 40,000 | about 120 | about 40 h | about 15 h | about 7.5 h |
| `--ext3-cells docs57` | 37,000 | about 113 | about 38 h | about 14 h | about 7 h |

Table: Full-run estimate (derived; assumes the server's cores run at this machine's unloaded speed, which has not been measured)

Treat these figures as accurate to within a factor of about 1.5. The
calibration rests on one paired measurement under load and on the
stored extra2 figure. This is above `docs/57`'s estimate of 9 hours on
8 workers, which assumed 7.1 seconds per trial.

**Recommended `simrun` resources:** `--cores 16 --mem 32` if the host
can give 16 cores to one job, otherwise `--cores 8 --mem 16`. Memory:
each worker holds one trial and its fits, a few hundred MB. The final
merge holds about 720,000 rows (40,000 trials times 18 analyses), on
the order of 200 MB (inferred), so 2 GB per core is ample. The host's
available cores per job are not known from the repository
(unverified).

## 6. Commands

From the repository root
(`/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng`). The refined
process (anchored carryover, randomized allocation) is the default, so
the files carry the tag `-anchored`.

### 6.1 Locally

```bash
S=analysis/scripts/quick-sim/carryover-closed-form
O=analysis/data/quick-sim/carryover-closed-form/strawman-stress
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

Rscript $S/11-strawman-stress-test.R --extra3 --check
Rscript $S/11-strawman-stress-test.R --extra3 --pilot 2 --cores 4
Rscript $S/11-strawman-stress-test.R --extra3 --cores 8 \
  > $O/run-full-extra3-anchored.log 2>&1

Rscript $S/12-strawman-stress-decision.R --ext3 --no-symmetry \
  --tag -anchored
```

The pilot writes `replicates-pilot2-extra3-anchored.rds` into the data
directory; delete it or pass `--out-dir` to keep it out.

### 6.2 On the simulation server (`simrun`)

As in `RERUN-anchored.md`, Section 4.2: `simrun submit` runs the
committed tree and refuses while tracked files have uncommitted
changes. Commit `11-strawman-stress-test.R`,
`12-strawman-stress-decision.R`, the scripts in
`analysis/scripts/data-chosen-half-life/` (the driver reads
`57-common.R` at run time, and that directory is untracked at present),
and this note. The job sets `SIM_CORES`, which the driver uses as its
worker count. The command is an argument vector, so the BLAS
thread variables cannot be set inline. In earlier jobs the container's
BLAS was not checked for threading (unverified); if it threads, prefix
the command with `env VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1`.

```bash
S=analysis/scripts/quick-sim/carryover-closed-form

simrun submit --cores 2 --mem 4 --label ext3-check -- \
  Rscript $S/11-strawman-stress-test.R --extra3 --check
simrun submit --cores 16 --mem 32 --label stress-ext3 -- \
  env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  Rscript $S/11-strawman-stress-test.R --extra3
```

Then `simrun watch ID`, `simrun fetch ID`, `simrun apply ID`, and run
the decision script of Section 6.1 locally. A failed job does not resume
from its chunks (`RERUN-anchored.md`, Section 4.2); `simrun resubmit
ID` starts it again.

## 7. Outputs

In `analysis/data/quick-sim/carryover-closed-form/strawman-stress/`
(with `--ext3-cells docs57`, `-extra3` becomes `-extra3-docs57`):

- `cells-full-extra3-anchored/`: chunk files and `progress.log`;
- `replicates-full-extra3-anchored.rds`: one row per trial and
  analysis, with columns `cell`, `rep`, `analysis`, `h` (half-life
  used; weighted mean for MA and MA0; argmax for the maxima), `est`,
  `se`, `t`, `df`, `p`, `switched` (F1G), `tstar` (resampled
  statistic, nominal rows), `crit_ws` (`B` rows), `posthoc`, `secs`
  (seconds per trial), and the cell settings;
- `summary-full-extra3-anchored.csv`: rejection rate, failures, mean
  estimate, empirical SD, mean SE, $\kappa$ and MCSE by cell and
  analysis.

From `12 --ext3 --no-symmetry --tag -anchored`, with suffix
`-f1-nosym-ext3-anchored`:

- `decision-screen`, `decision-screen41`, `decision-sensitivity`,
  `decision-power`, `decision-maximin`.

## 8. How the decision script reads them

`12-strawman-stress-decision.R --ext3` reads
`replicates-{run}-extra3{tag}.rds` on its own. It does not combine the
file with the main run, because the trials differ, and it refuses
`--with-extras` alongside `--ext3`. It judges the 18 analyses of
Section 2. The strawman defaults to F1 at 0.5 week. Every label that
begins with `prog-` is a prognostic cell (rule 2 of `docs/46`); this
includes the 6 new prognostic null cells and `prog-BL`. The nominal
screen is pooled over the 12 non-prognostic null cells (the 11 of
`docs/47` and the new `thalf2` cell). `decision-screen41` repeats it
over the 11 original ones for comparison with `docs/47`. Size-adjusted
power uses the null cell with the same label, which exists for every
new alternative. The selection table reports F1A, F1A0 and F1G. The
warp-speed rows are judged by their $p^*$. Their size is the
bootstrap test's size, and its MCSE is somewhat larger than the
binomial one the script prints (`docs/57`, Section 4.5). Without
`--ext3`, the script reads and writes exactly what it did before.

## 9. Documents to update after the run

- `docs/57`: Sections 1, 4.2 to 4.8 (status lines), 6 (state that the
  illustration is superseded at scale), 7 (decision table and
  provisional plan; whether F1G, MA0 or MaxTB is to be preferred),
  and 8 (limitations: refined process, power at a wrong half-life,
  prognostic cells with a wrong half-life, and the multivariate normal
  and $t$ calibration, now run).
- `docs/46`: a further amendment recording the extension as an
  exploratory post hoc arm on fresh seeds, its 17 cells, and the
  analyses judged.
- `docs/47`: a section on the extension (screens, prognostic cells,
  size-adjusted power in cells 31, 32 and 55 to 58, the selection
  shares), stating that its trials are independent of the main run.
- `docs/56`, Section 6: the half-life recommendation.
- `docs/58`, Section 9.2: the influence-based multivariate normal and
  $t$ calibration, as applied to the interaction rather than the main
  effect.
- Paper 01 and `analysis/report/compendium-mathematics.tex`, where they
  state how the analysis half-life is chosen or quote F1A's size.
- `docs/50`, "Open items".
