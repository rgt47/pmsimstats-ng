---
title: 'Rerunning every simulation affected by the parameter error: plan for four hosts'
author: 'pmsimstats team'
---

# Rerunning every simulation affected by the parameter error: plan for four hosts
*2026-11-03 17:51 PST*

## 1. Purpose and scope

This plan reruns every simulation whose reported results rest on a
parameter set other than Hendrickson et al.'s published values
(`docs/70`): the 'simcore2026' set of `simulation-core.R` and the 'raw
fixture' copied from the unit tests. It extends `docs/73` (papers 01
and 02 only) to every affected paper and adds two hosts, athena and
joe, which became simrun hosts on 2026-11-03. Where the two plans
differ, this one governs.

The inventory comes from four read-only audits on 2026-11-03, which read
each driver's parameter source, the files each manuscript reads, and
the stored runs' provenance. Their findings are summarized in Sections
3 and 4; each statement there is marked verified (read in code,
provenance or a stored file) or inferred.

**Affected and in this plan:** papers 01, 02, 03, 06, 08, 09 and 10;
the cr2-olbdc check (read by papers 02 and 09); and the prognostic
material leaving paper 02 for a companion paper (paper 15, whose
folder does not exist yet).

**Confirmed unaffected (verified):** papers 04 and 05 use their own,
stated parameter sets; papers 07, 11, 13 and 14 load the published
values from the package data; paper 12 inherits published values from
an archived paper 11 run; paper 01's closed form, study 2 and study 3
use the published values.

## 2. Hosts and measured capacity

Benchmarks on 2026-11-03 ran driver 19 at commit c866404 in the
project image (`...:10504a44fb8d`) on each host.

| Host | Hardware | One worker, heaviest cell (s per replicate) | 8-worker block F run | Throughput relative to zeus |
|---|---|---|---|---|
| zeus | Linux server, 12 threads | 2.05 | 6.5 min | 1.0 |
| athena | M1 Max, 10 cores, Docker 24 GB | 5.57 | 8.2 min | about 0.8 |
| joe | M2, 8 cores, Docker 20 GB | 4.90 | not run | about 0.5 (inferred) |
| julia (overnight only) | M1, 8 cores, Docker 12 GB | 5.47 | not run | about 0.4 (inferred) |
| mercury | Linux, 6 cores used | not run | not run | about 0.4 (job 40) |

- The Macs run the amd64 image under Rosetta, about 2.6 times slower
  per worker than zeus, but athena keeps about two thirds of its
  single-worker speed at 8 workers, against about a third on zeus. One
  maths-library thread per worker did not change zeus's result (7.4
  against 6.5 minutes), so the cause of zeus's poor scaling is unknown.
- An earlier inference, in `docs/73` and `docs/72`, that zeus runs at a
  quarter of julia's per-core speed was wrong: one zeus worker matches
  julia natively on light cells and is about 35 percent slower on the
  heaviest. The fourfold gap is a loss when 8 workers run together.
- **Pairing across hosts.** Results are not draw-for-draw identical
  between hosts, and the Macs' emulated runs are not identical to
  zeus's. Jobs that must be compared trial by trial run on the same
  host (Section 6).
- **Memory.** athena's Docker Desktop has 24 GB (jobs use 22) and
  joe's 20 GB (jobs use 18), both from 2026-11-03.

## 3. Inventory by paper

Hours are wall-clock on zeus at 8 cores ('zeus-hours') where a stored
zeus job did the same work, and local hours on julia (native, about 7
workers) where only a local run exists; local hours are not converted,
because the relation between julia's native multi-worker speed and the
container hosts was not measured.

### Ready now: the parameter switch reaches the driver

| Paper | Job | Driver | Cost | Status |
|---|---|---|---|---|
| 01 | p01-pf | 19, blocks P and F, G10 and G9 | 11.9 zeus-h | ready (`docs/73`) |
| 01 | p01-n80 | 19, block N | 4.7 zeus-h | ready |
| 02 | ten jobs | 01, 27, 38, 06, 32 (three), 36 (three) | 33.0 zeus-h | ready (`docs/73` Section 6) |
| 15 | p15-gapC | 32 `--block C` | about 13.5 h on mercury (job 40) plus G10 | ready; needs a new `--out-dir` (verified) |
| 15 | p15-refT | 36 `--block T` | 0.13 zeus-h (job 36) | ready |
| 15 | (c_pb = 0.3 cells) | 36 `--block N` | in paper 02's ref-N job | no extra job (verified) |

### Small code changes first (Tier B)

| Paper | Job | Driver | Cost | Change needed |
|---|---|---|---|---|
| 10 | p10 | `run-paper10-refined.sh` (nine drivers) | 1.8 zeus-h (job 55) | `p10_out_dir()` and the wrapper write to a folder per parameter set (`output/paper10-published/`), record `param_set()` in metadata; the Rmd reads that folder; export `PMSIM_PARAMS=published` in the wrapper. The kappa script's printed 'sigma_u = 10' and Rmd line 720 ('B_i drawn from N(0, 1)') change. |
| 15 | p15-stress | 28 | 5.5 zeus-h (job 28) | add `guard_out_dir()` and `param_meta()`; new `--out-dir` (it resumes from `trials/` checkpoints) |
| 15 | p15-re | 29 | 3.2 zeus-h (job 30) | same as 28; same host as 28 (shared seeds) |
| 03 | p03-null | `12-referee-rerun.R`, null block | 2.0 zeus-h (job 38) | add `guard_out_dir()` and `param_meta()`; new out-dir; standardize `bm` before `lme_analysis()` and `hlme(classmb = ~bm)` (they use the raw biomarker) |
| 03 | p03-alt | same, alternative block | unknown: job 39 was cancelled after 152 h | as p03-null; run `--time` first and size from it |
| 02, 09 | cr2-p02 | `cr2-olbdc-check/01-run-cr2-null.R`, p02 setting | about 0.4 local h | export `PMSIM_PARAMS=published` (its wrapper passes no `--params`); new `CHK_OUT` (the summarizer pools every file in the folder); the p02 baseline rows use the raw biomarker |

### Larger code changes first (Tier C)

These drivers hard-code the raw fixture (TV and BR SD 5, PB SD 2, the
TV mean on the BR curve, biomarker 0/1, baseline 70/10; verified) and
do not source `simulation-core.R`. Each needs a parameter switch, a
true value taken from the active set, a biomarker scale decision
(Section 5, D1), an output folder that cannot overwrite the stored
run, and worker sizing from `SIM_CORES`.

| Paper | Job | Driver | Stored cost | Notes |
|---|---|---|---|---|
| 06 | p06-s1 | `run-study-a-prod.R` (Study 1) | about 13 local h (7 workers) | true value hard-coded at -2.25 and typed in the text; the m_PB and m_TV axes weaken against published SDs (D6) |
| 06 | p06-s2 | `run-study-2-feasible.R` (corrected Study 2) | 4.3 zeus-h (job 49; ran 11 workers on an 8-core cgroup) | contamination bias grows about 3 times relative to the signal under published SDs (inferred) |
| 06 | p06-s4 | `run-study-b-balanced.R` (Study 4) | about 0.7 local h | `TRUE_BETA` hard-coded |
| 06 | p06-bounds | `05-feasibility-bounds.R` | seconds | feasibility unchanged; SD-dependent columns change |
| 08 | p08 | `01-study1-test-procedure.R` | about 0.2 local h | nearly every result table is typed: a rerun means re-narration; `02-mancl-t-reference.R` must follow the new prefix |
| 09 | p09 | `20-paired-dropout.R` | about 1.6 local h (500 null replicates; 2.2 at the current default of 1,000) | full published set, not only BR SD 8 (D3); -2.25 typed in 15 places; the informative-dropout mechanism may weaken (inferred) |
| 02, 09 | cr2-p09 | `cr2-olbdc-check`, p09 setting | about 0.4 local h | follows `20-paired-dropout.R`; its `iso` setting edits `rp` with data.table syntax, which breaks on a tibble; same host as p09's iso run |

### Not rerun (recommendation; Section 5, D2)

Runs the manuscripts already describe as superseded, withdrawn or for
the record: paper 03's pilot (`quick-sim-driver-parallel.R`) and old
production (`10-study5-production.R`); paper 06's June Study 2
(contamination) and Study 3 (recovery) and its prototype; paper 09's
prototype (`09-dropout-prototype-10min.R`, Section `sec-old`). They
stay as stored and are labeled in the text as raw-fixture runs.
Paper 10's pre-correction (2026-05-20) outputs and paper 15's
repair-artifact cells (prognostic strength 0.6) stay as labeled
historical simcore2026 exhibits (D4, D5).

## 4. Findings outside the parameter error

- **Paper 12 (verified).** It reads a July 2026 run of the paper 11
  factorial; the current stored summary no longer has the N = 35 rows,
  the E1 specification or `sd_beta`, so `report.Rmd` lines 322-327
  describe a file that no longer exists in that form, and the July run
  predates the October DGP changes (inferred). A rerun of the paper 11
  factorial cells paper 12 needs is a provenance repair, not a
  parameter repair; it is listed as optional.
- **Paper 13 (verified).** `report.Rmd` line 70 reads
  `analysis/data/13-two-component-grid.rds`, which does not exist; the
  driver writes `13-drop-pb-grid.rds`. Paper 13 renders with every
  simulated cell missing. A one-line path fix; no simulation.
- **Paper 06 (verified).** `run-study-2-feasible.R` ran 11 workers in
  an 8-core job, and the report's core-hours use 11.
- **`docs/70` (verified).** It names `nof1power` as paper 04's engine
  (it is paper 05's), lists `00-common.R` as behind paper 09's results
  (it is an unrun scaffold), describes a `SIGMA_BR` axis for paper 09
  (it is a scalar), and lists `09-calibration-substudy.R` as feeding
  paper 03 (its output is not read).

## 5. Decisions (recommended defaults in force unless changed)

- **D1. Biomarker scale.** Standardize the biomarker by the active
  set's mean and SD before every analysis, as `prepare_long_data()`
  now does. Tests are unchanged (the map is affine), the true slope
  stays $-c_{bm}\sigma_{BR}$ per SD (-3.6 at $c_{bm} = 0.45$), and the
  raw biomarker (mean 124) cannot make the interaction nearly collinear
  with the drug term. The alternative is to keep the raw scale and
  report -0.234 per mm Hg.
- **D2. Superseded and withdrawn runs.** Do not rerun; label them as
  raw-fixture runs (Section 3).
- **D3. Paper 09.** Move to the full published set, not only BR SD 8,
  so all papers share one parameter table.
- **D4. Paper 10's pre-correction comparison.** Keep the 2026-05-20
  runs as the historical arm, labeled simcore2026, and rerun only the
  refined arm; say that the comparison mixes the DGP correction with
  the parameter change, or rerun the historical arm under published as
  well (about 1.8 more zeus-h).
- **D5. Paper 15's repair-artifact cells.** Keep as a labeled
  historical exhibit; under today's stop-on-non-PD default they would
  only fail. Rerun the feasible prognostic cells under published.
- **D6. Paper 06's design axes.** Keep the m_PB and m_TV grids and the
  drivers' structure, adopt the full published table (TV rate 0.35,
  PB max 6.51) as the base, and report that the axes are now weaker
  relative to the noise.
- **D7. Paper 03's alternative block.** Time it with `--time` under
  published before submitting; if it projects above about 24 hours on
  any host, split it by cell.

## 6. Host assignment

Pairing rules: job 27 with the E0 isolation, and gap-fill blocks A and
B, each on one host (paper 02); drivers 28 and 29 on one host (shared
seeds); the nine paper 10 drivers in one job (the wrapper); the
cr2-olbdc `iso` run with paper 09's p09 setting on one host. Everything
else is unpaired.

**Wave 1 (ready now or after Tier B).** Expected hours on each host
(zeus-hours divided by the host's relative throughput):

| Host | Queue, in order | Hours |
|---|---|---|
| zeus (8 cores, 24 GB) | p01-pf, p02-tier1, p02-e0iso, p02-gf-dropout, p02-gf-size, p02-gf-B | about 36 |
| athena (8 cores, 22 GB) | p02-ref-A, p02-factorial, p02-scen, p02-ref-S, p02-ref-N, p15-refT, p10, p15-stress | about 22 |
| mercury (6 cores, 12 GB) | check job, p01-n80, p15-gapC | about 28 |
| joe (8 cores, 18 GB) | p03-null, then p03-alt after its timing run | about 4 plus p03-alt |

**Wave 2 (after Tier C code and its smoke tests).** joe takes paper 09
with both cr2-olbdc settings (one host) and paper 06 Study 4 and the
bounds; julia takes paper 08 overnight (6 cores, 10 GB; julia is the
working machine, so it runs only small, unpaired jobs submitted in the
evening); Study 1 and corrected Study 2 go to whichever of zeus and
athena frees first. Each Tier C driver is timed on its host before
submission, because only local timings exist.

With every queue full, wave 1 ends in about a day and a half; wave 2
adds about a day once its code is ready.

**Submitted 2026-11-03 18:30 PST at commit 9d7821e.** zeus jobs 60 to
65; athena jobs 4 to 11; mercury check job 50 (passed), then jobs 51
and 52; joe job 2 (timing run of paper 03's alternative block) and job
3 (null block). Driver 29 (p15-re) was dropped under D5: every
confounded cell in it is a strength-0.6 repair-artifact cell, so the
stored run stays as part of the historical exhibit.

## 7. Before any submission

1. **Code.** Tier B changes (Section 3), then the commit: the
   remediation files of `docs/73` Section 5 plus the Tier B drivers.
   Other uncommitted work stays out; every command uses `--commit HEAD`.
   Tier C changes are a second commit before wave 2; wave 1 jobs keep
   their commit.
2. **Checks at the commit.** `42-check-param-sets.R`; a smoke run of
   every wave 1 driver under the published set; for Tier B drivers, a
   check that a rerun into the stored folder is refused.
3. **Image.** Unchanged on all four hosts (`...:10504a44fb8d`); athena
   and joe hold the registry copy built from the same `Dockerfile` and
   `renv.lock`.
4. **Hosts.** athena and joe plugged in, logged in, lid open; Docker
   Desktop's Resource Saver off. Mercury's one-minute check job first.

## 8. Commands, wave 1

Paths are set once; every simrun option is written out (zsh does not
split an unquoted variable into words).

```bash
D=analysis/scripts/carryover-sensitivity
D19=analysis/scripts/quick-sim/carryover-closed-form
D19=$D19/19-g9-planning-and-fig4b.R
P=analysis/data/carryover-sensitivity/published
```

**zeus** (`docs/73` Section 6 gives the same commands):

```bash
simrun submit --commit HEAD --cores 8 --mem 24 --label p01-pf -- \
  Rscript $D19 --params published --specs E7null_cr2,E7cr2 --block all
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

**athena:**

```bash
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p02-ref-A -- \
  Rscript $D/36-run-referee-variants.R --params published --block A \
  --out-dir $P/referee-variants/block-A
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p02-factorial -- \
  Rscript $D/01-run-factorial.R --params published --arch mvn \
  --out-dir $P/factorial
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p02-scen -- \
  Rscript $D/06-g-vs-f-comparison.R --params published \
  --reps 500 --n-boot 1 --cells paper02 --out-dir $P/g-vs-f-comparison
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p02-ref-S -- \
  Rscript $D/36-run-referee-variants.R --params published --block S \
  --out-dir $P/referee-variants/block-S
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p02-ref-N -- \
  Rscript $D/36-run-referee-variants.R --params published --block N \
  --out-dir $P/referee-variants/block-N
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p15-refT -- \
  Rscript $D/36-run-referee-variants.R --params published --block T \
  --out-dir $P/referee-variants/block-T
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p10 -- \
  env PMSIM_PARAMS=published bash $D/run-paper10-refined.sh
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p15-stress -- \
  Rscript $D/28-run-g9-prognostic-stress.R --params published \
  --out-dir $P/g9-prognostic-stress
simrun --host athena submit --commit HEAD --cores 8 --mem 22 \
  --label p15-re -- \
  Rscript $D/29-run-g9-random-effects.R --params published \
  --out-dir $P/g9-random-effects
```

**mercury** (after the check job of `docs/73` Section 5):

```bash
simrun --host sim-mercury-ts submit --commit HEAD --cores 6 --mem 12 \
  --label p01-n80 -- \
  Rscript $D19 --params published --specs E7null_cr2,E7cr2 --block N
simrun --host sim-mercury-ts submit --commit HEAD --cores 6 --mem 12 \
  --label p15-gapC -- \
  Rscript $D/32-run-g9-gap-fill.R --params published --block C \
  --out-dir $P/g9-gap-fill/C
```

**joe:**

```bash
D3=analysis/scripts/latent-class-mixture-application
P3=analysis/data/latent-class-mixture-application/published
simrun --host joe submit --commit HEAD --cores 8 --mem 18 \
  --label p03-null -- \
  Rscript $D3/12-referee-rerun.R --params published --block null \
  --out-dir $P3/referee-rerun-null
```

The paper 03 command's block argument and the paper 10 and 15 commands
take their final form with the Tier B changes; check each driver's
usage line before submitting.

## 9. Monitoring and failure handling

Zeus and mercury ignore `--time`, and zeus has stalled before (jobs 20,
39, 47 and 52 ran 125 to 263 hours on work that takes hours). Each job
is checked against its expected hours (Section 6, or the timing run for
Tier C) and cancelled at twice that. A Mac job stops if the Mac sleeps
or Docker Desktop quits; it then shows as failed, and is resubmitted.

## 10. After the runs, per paper

- **01:** render; every [pending] marker disappears.
- **02:** `docs/71` Section 8 (repoint paths, assemble the gap-fill
  jobs, ladder, referee reanalysis, figures, Methods parameter table,
  new Results).
- **03:** copy the referee rerun into the folder the report reads (the
  null block from job 38 was never copied, verified), relabel the
  superseded runs, re-narrate.
- **06, 08, 09:** replace typed numbers with values read from the new
  outputs, then re-narrate; the -2.25 true value and the fixture SDs
  disappear from the text.
- **10:** point the Rmd at `output/paper10-published/`; label the
  historical arm.
- **15:** create the paper folder (split plan step 5) and read the
  published-set outputs.

## 11. Evidence status

- **Measured:** the host benchmarks of Section 2; the stored jobs'
  durations; the image tag on every host.
- **Verified (read in code, provenance or stored files):** each
  driver's parameter source, output path and the files each manuscript
  reads; the findings of Section 4.
- **Inferred:** joe's and mercury's 8-worker throughput; the G10
  overhead; costs given as local hours; the scientific consequences
  flagged in Section 3 (weaker contamination signal relative to noise,
  weaker informative dropout).
- **Not done:** the Tier B and Tier C code changes; the paper 03
  alternative block's timing.
