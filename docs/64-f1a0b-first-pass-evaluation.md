# F1A0B: A First-Pass Evaluation Against the docs/62 Criteria

*2026-10-22 17:15 PDT*

**Author.** pmsimstats team

**Purpose.** `docs/62` catalogs the `F`-family half-life-selection
methods and states that F1A0B -- the project's current primary
data-driven candidate -- "has not been tested in paper 02 or cited
by it," and lists nine criteria ("what best would mean") against
which no method has been fully checked. This note reports a first,
small-scale evaluation of F1A0B (and its F1/F1A/F1A0 relatives)
against three of those criteria, on paper 02's own driver and the
refined-process DGP, as a first step toward closing that gap. It is
exploratory, not pre-registered or confirmatory.

## What was built

Paper 02's driver
(`analysis/scripts/carryover-sensitivity/simulation-core.R`)
previously fit only the `G`-family's simple model
(`Sx ~ bm + t + X + bm:X`). This session added the richer model
`Sx ~ bm + t + Dbc + bm:Dbc + De + bm:De` (bm already centered at 0
by the DGP's baseline parameters, so it serves as `bmc` directly),
and the three half-life-selection rules that differ among the
`F`-family methods cataloged in `docs/62` Part 2:

- **F1**: $t_{1/2}^a$ fixed at the project's policy value, 0.2 week
  (`CLAUDE.md`, "Policy 2026-10-11").
- **F1A**: AIC over the grid $\{0, 0.2, 0.5, 1, 2\}$ weeks with the
  interaction already in the selection model (the liberal pattern).
- **F1A0**: AIC over the same grid on the null (no-interaction)
  model, refit with the interaction at the selected value.
- **F1A0B**: F1A0's test calibrated by a warp-speed parametric
  bootstrap -- one resample per Monte Carlo replicate from the null
  model refit at F1A0's own choice of $t_{1/2}^a$, pooled across the
  cell, following the convention in
  `analysis/scripts/data-chosen-half-life/57-progbl-check.R`.

This machinery is self-contained in `simulation-core.R` (functions
`grid_fits_f`, `grid_methods_f`, `sim_null_f`, `boot_stat_f1a0`,
`fit_f_family`, `simulate_cell_f`); it does not depend on the
separate stress-test driver
(`analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R`)
or its own F-family code
(`analysis/scripts/data-chosen-half-life/57-common.R`), though the
two implementations follow the same conventions and are not
cross-validated against each other.

**A design blocker, found and worked around.** The driver's `CO`
design preset has constant expectancy (`De`, 0.5 at every
timepoint), which makes `De` and `bm:De` exactly collinear with the
intercept and `bm`: the `F`-family model is not identifiable on `CO`
at all (*verified*: every fit returns `NA`). `OLBDC` (open-label
run-in followed by blinded discontinuation, expectancy 1.0 then
0.5) was used instead for the design-robustness arm. This is
documented in the code and should be kept in mind by anyone
extending this evaluation to a true crossover cell.

## Method

Three cells of the refined-process DGP (graded coupling, separable
AR(1) with $\rho = 0.7$, anchored carryover, RH's fixed path
allocation), $N = 70$, $t_{1/2} = 0.5$ week, exponential carryover
decay:

| Cell | Design | $c_{bm}$ | Other |
|---|---|---|---|
| `hybrid_null` / `hybrid_alt` | Hybrid | 0 / 0.45 | -- |
| `olbdc_null` / `olbdc_alt` | OL+BDC | 0 / 0.45 | -- |
| `progbl_null` | Hybrid | 0 | $c_{bl} = 0.3$ (prog-BL) |
| `progtvconst_null` | Hybrid | 0 | $c_{tv} = 0.3$ (prog-TV-const) |
| `progtvgrow_null` | Hybrid | 0 | $b_{tv} = 0.3$ (prog-TV-grow) |
| `progpb_null` | Hybrid | 0 | $c_{pb} = 0.3$ (prog-PB) |

500 replicates per cell, one warp-speed bootstrap resample per
replicate ($n_{boot} = 1$, the cheapest setting). Driver:
`analysis/scripts/carryover-sensitivity/05-f1a0b-evaluation.R`.
Output: `analysis/data/carryover-sensitivity/f1a0b-evaluation/`
(`all.rds` pooled per-replicate results, `summary.rds` the table
below, `trials/` per-replicate checkpoints). Run time: 26.3 minutes
for the first 2,500 (cell, replicate) jobs (five cells) on 6 cores,
plus a further 15.6 minutes for 1,500 more jobs (the three
additional prognostic pathways, added once `c.bm.pb` was mirrored
into `implementations/tidyverse` and plumbed through this driver --
previously package-only, a gap flagged earlier in this session).

## Results

Rejection rate at $\alpha = 0.05$: Type I error at the null cells,
power at the alternative cells. Monte Carlo standard error
$\approx \sqrt{0.05 \times 0.95 / 500} \approx 0.010$ near the null,
so about $\pm 0.02$ at two MCSE.

| Cell | F1 | F1A | F1A0 | F1A0B |
|---|---|---|---|---|
| Hybrid null | 0.042 | 0.040 | 0.040 | 0.030 |
| Hybrid alt ($c_{bm}=0.45$) | 0.950 | 0.957 | 0.954 | 0.944 |
| OL+BDC null | 0.046 | 0.046 | 0.040 | 0.034 |
| OL+BDC alt | 0.849 | 0.844 | 0.844 | 0.844 |
| prog-BL null | 0.042 | 0.040 | 0.038 | 0.028 |
| prog-TV-const null | 0.040 | 0.044 | 0.042 | 0.032 |
| prog-TV-grow null | 0.042 | 0.040 | 0.040 | 0.030 |
| prog-PB null | 0.042 | 0.042 | 0.042 | 0.030 |

*Verified* (computed by running the driver above in this session).

## Interpretation, against the docs/62 criteria

**1. Type I error (size).** No inflation anywhere. F1A0B is the
most conservative of the four at every null cell checked (0.028 to
0.034 against nominal 0.05), the expected direction for bootstrap
calibration and never liberal in this pass.

**2. Power.** No advantage for F1A0B over plain F1 at either
design. OL+BDC: 0.844 vs. 0.849, within noise. Hybrid: both near
ceiling (0.94 to 0.96) at this effect size, too saturated to
discriminate among the four methods; a smaller $c_{bm}$ would be
needed to compare power meaningfully there.

**5. Robustness to trial design.** Both designs checked hold
nominal size and show no collapse comparable to `G3`/`G4`'s CO
failure (`docs/62`, Part 1). This is not evidence about a true
crossover cell, since `CO` itself is unusable for this model (see
above); `OLBDC` is a different design, not a substitute measurement
of the same thing.

**7. Robustness to a prognostic biomarker.** This is the criterion
`docs/62` singles out as checked for no `G`-family specification at
all. All four prognostic pathways of `docs/57` Section 7.3 are now
run ($c_{bm} = 0$ throughout, so the biomarker drives none of the
true interaction, only the confound named): prog-BL (biomarker
correlated with baseline severity), prog-TV-const (correlated with
the natural-history trajectory), prog-TV-grow (an additive mean
shift of that trajectory) and prog-PB (correlated with the placebo
response). F1A0B holds nominal size at every one (0.028 to 0.032),
as do F1, F1A and F1A0, with the same qualitative pattern across all
four: F1A0B mildly conservative, the other three close to nominal.
No pathway shows inflation. This is the first check of this specific
robustness claim across all four named confounds, at one effect
size ($c_{tv}/b_{tv}/c_{bl}/c_{pb} = 0.3$) and one design (Hybrid).

**8. Computational cost.** 26.3 minutes for 2,500 (cell, replicate)
jobs on 6 cores is approximately 0.38 seconds per job amortized
across cores, or roughly 2.1 seconds serially per job (three nominal
fits plus one bootstrap resample), consistent with `docs/62`'s
reported cost for the data-driven F-family methods.

**9. Validation completeness.** This pass checks criteria 1, 2
(partially, confounded by ceiling effects at Hybrid), 5 (partially,
substituting OL+BDC for CO) and 7 (all four named prognostic
pathways, one effect size each, one design) for F1A0B specifically,
at 500 replicates per cell. Criteria 3 (half-life misspecification),
4 (decay-shape misspecification) and 6 (DGP architecture) remain
entirely untested for F1A0B, as does any direct, shared-grid
comparison against `G8`/`G9` (criterion 3 of `docs/62` Part 3's "how
the two families connect").

## Caveats, unsolicited

- Single effect size per design; the Hybrid alternative cell is
  saturated near the power ceiling and cannot discriminate among
  the four methods there.
- 500 replicates, one pass: the gaps between F1A0B and the other
  three methods (0.006 to 0.012 at the null cells) are within or
  close to Monte Carlo noise, not established differences.
- $n_{boot} = 1$ is the noisiest warp-speed setting; a larger
  per-replicate resample count (or more replicates) would tighten
  F1A0B's calibration estimate specifically.
- All four prognostic pathways were tested, but at one effect size
  each ($0.3$) and one design (Hybrid); misspecification of that
  effect size, or whether a prognostic confound and the true
  interaction co-occur, is untested.
- No comparison against the `G`-family (`G1`-`G9`) on the same
  cells; this note evaluates F1A0B in isolation, not head to head.
- This driver's F-family code and
  `analysis/scripts/data-chosen-half-life/57-common.R`'s are
  independent implementations of the same methods; they have not
  been cross-checked against each other on shared trials.

## References

- `docs/62-carryover-specification-and-half-life-method-catalog.md`
  -- the nine criteria and the F-family catalog this note tests
  against.
- `docs/57-data-chosen-carryover-half-life.md`, Section 7.3 -- the
  four prognostic pathways, of which only prog-BL was run here.
- `analysis/scripts/carryover-sensitivity/05-f1a0b-evaluation.R` --
  the driver for this evaluation.
- `analysis/scripts/carryover-sensitivity/simulation-core.R` -- the
  F-family machinery (`grid_fits_f`, `grid_methods_f`, `sim_null_f`,
  `boot_stat_f1a0`, `fit_f_family`, `simulate_cell_f`), the `c.bm.bl`
  / `c.bm.tv` / `b.bm.tv` prognostic-coupling additions it builds on,
  and the documented CO-design identifiability blocker.
- `analysis/scripts/data-chosen-half-life/57-common.R` -- the
  independent F-family implementation this note's code does not
  share or cross-validate against.
- `CLAUDE.md` (pmsimstats-ng project file), "Recommended analysis
  for the interaction (F1)" -- the 0.2-week policy value for F1.
