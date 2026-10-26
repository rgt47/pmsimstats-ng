# F1A0B: A First-Pass Evaluation Against the docs/62 Criteria

*2026-10-22 17:15 PDT*

**Author.** pmsimstats team

**Purpose.** `docs/62` catalogs the `F`-family half-life-selection
methods and states that F1A0B -- the project's current primary
data-driven candidate -- "has not been tested in paper 02 or cited
by it," and lists nine criteria ("what best would mean") against
which no method has been fully checked. This note reports a first,
small-scale evaluation of F1A0B (and its F1/F1A/F1A0 relatives)
against seven of those criteria, on paper 02's own driver and the
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
| `hybrid_alt_lo` | Hybrid | 0.25 | de-saturation arm, added 2026-10-22 |
| `olbdc_null` / `olbdc_alt` | OL+BDC | 0 / 0.45 | -- |
| `t1hifhi_null` / `_alt` | Hybrid | 0 / 0.25 | $t_{1/2}=2$ wk (crit. 3) |
| `shapemis_null` / `_alt` | Hybrid | 0/0.25 | Weibull DGP, exp fit (crit 4) |
| `meanmod_null` / `_alt` | Hybrid | 0 / 0.25 | mean-mod. DGP (crit. 6) |
| `densegrid_null` | Hybrid | 0 | AIC grid widened to 16 points |
| `vdensegrid_null` | Hybrid | 0 | AIC grid widened to 31 points |
| `smalln_null` | Hybrid | 0 | $N=30$ (halved) |
| `offgrid_null` | Hybrid | 0 | true $t_{1/2}=0.35$ wk, between grid points |
| `smalln_offgrid_dense_null` | Hybrid | 0 | $N=30$, $t_{1/2}=0.35$, 31 pts |
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
a further 15.6 minutes for 1,500 more jobs (the three additional
prognostic pathways, added once `c.bm.pb` was mirrored into
`implementations/tidyverse` and plumbed through this driver --
previously package-only, a gap flagged earlier in this session),
a further 5.6 minutes for 500 more jobs (`hybrid_alt_lo`, added to
de-saturate the Hybrid power comparison), a further 40.9 minutes
for 3,000 more jobs (the six criterion-3/4/6 cells added
2026-10-23), and a further 58.2 minutes for 1,000 more jobs (the
grid-density stress-test cells, denser grids cost more per
replicate than the standard 5-point grid, also added 2026-10-23),
and a further 25.8 minutes for 1,500 more jobs (`smalln_null`,
`offgrid_null`, `smalln_offgrid_dense_null`, added 2026-10-23 to push
the two conditions most plausibly able to break F1A0's size, which
the grid-density test alone did not exercise).

## Results

Rejection rate at $\alpha = 0.05$: Type I error at the null cells,
power at the alternative cells. Monte Carlo standard error
$\approx \sqrt{0.05 \times 0.95 / 500} \approx 0.010$ near the null,
so about $\pm 0.02$ at two MCSE.

| Cell | F1 | F1A | F1A0 | F1A0B |
|---|---|---|---|---|
| Hybrid null | 0.042 | 0.040 | 0.040 | 0.030 |
| Hybrid alt ($c_{bm}=0.45$) | 0.950 | 0.957 | 0.954 | 0.944 |
| Hybrid alt, lower effect ($c_{bm}=0.25$) | 0.513 | 0.521 | 0.511 | 0.477 |
| OL+BDC null | 0.046 | 0.046 | 0.040 | 0.034 |
| OL+BDC alt | 0.849 | 0.844 | 0.844 | 0.844 |
| prog-BL null | 0.042 | 0.040 | 0.038 | 0.028 |
| prog-TV-const null | 0.040 | 0.044 | 0.042 | 0.032 |
| prog-TV-grow null | 0.042 | 0.040 | 0.040 | 0.030 |
| prog-PB null | 0.042 | 0.042 | 0.042 | 0.030 |
| Half-life misspec. null ($t_{1/2}=2$ wk) | 0.042 | 0.058 | 0.054 | 0.060 |
| Half-life misspec. alt ($t_{1/2}=2$ wk) | 0.193 | 0.291 | 0.269 | 0.269 |
| Decay-shape misspec. null | 0.040 | 0.036 | 0.034 | 0.034 |
| Decay-shape misspec. alt ($c_{bm}=0.25$) | 0.629 | 0.627 | 0.623 | 0.615 |
| Mean-moderation null | 0.042 | 0.040 | 0.040 | 0.030 |
| Mean-moderation alt ($c_{bm}=0.25$) | 0.614 | 0.622 | 0.612 | 0.561 |
| Dense-grid null (16 points) | 0.042 | 0.042 | 0.040 | 0.038 |
| Very-dense-grid null (31 points) | 0.042 | 0.042 | 0.040 | 0.034 |
| Small-$N$ null ($N=30$) | 0.044 | 0.058 | 0.052 | 0.044 |
| Off-grid null ($t_{1/2}=0.35$ wk) | 0.042 | 0.038 | 0.038 | 0.036 |
| Small-$N$ + off-grid + dense-grid null | 0.044 | 0.056 | 0.052 | 0.044 |

*Verified* (computed by running the driver above in this session).

## Interpretation, against the docs/62 criteria

**1. Type I error (size).** No inflation anywhere. F1A0B is the
most conservative of the four at every null cell checked (0.028 to
0.034 against nominal 0.05), the expected direction for bootstrap
calibration and never liberal in this pass.

**2. Power.** No advantage for F1A0B over plain F1 at any design or
effect size checked. OL+BDC: 0.844 vs. 0.849, within noise. Hybrid
at $c_{bm}=0.45$: both near ceiling (0.94 to 0.96), too saturated to
discriminate among the four methods. A second Hybrid cell at
$c_{bm}=0.25$ (added 2026-10-22 to de-saturate this comparison) puts
power at 0.48 to 0.52, well clear of the ceiling: F1A0B is lowest
(0.477) and F1 highest (0.513), a 3.6-point gap. The Monte Carlo
standard error here is $\sqrt{0.5 \times 0.5/500} \approx 0.022$, so
the gap is just under 2 MCSE and not conclusive on its own, but it is
directionally consistent with every other cell in this table (F1A0B
is the lowest-rejecting method everywhere it appears, null or
alternative). The pattern is consistent with bootstrap calibration
costing a small amount of power as well as controlling size, as
expected, but this single de-saturated cell does not establish the
size of that cost.

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

**3. Robustness to half-life misspecification.** This is the
scenario the F1A0/F1A0B machinery exists for: when the true
data-generating half-life is far from F1's fixed policy value of
0.2 week, a correctly-specified adaptive method should recover
power that F1 leaves on the table. At $t_{1/2}=2$ weeks (a
"behavioral" scenario, `CLAUDE.md`), that is exactly what happens:
power is 0.193 for F1 against 0.291 for F1A, 0.269 for F1A0, and
0.269 for F1A0B -- a 7.6- to 9.8-point advantage over F1
(*verified*; the F1-vs-F1A gap is about 5 Monte Carlo standard
errors, $\mathrm{SE}\approx\sqrt{0.25\times0.75/500}\approx0.019$).
This is the first result in this evaluation where an adaptive
method clearly outperforms F1, and it is the one that most directly
supports using a data-chosen half-life at all.

The null cell at the same $t_{1/2}=2$ complicates that picture. F1's
size is nominal (0.042), but F1A, F1A0 and (notably) F1A0B are all
somewhat above nominal: 0.058, 0.054, and 0.060 respectively. The
F1A0B deviation (0.0604 vs. 0.05) is about 1.06 Monte Carlo standard
errors ($\mathrm{SE}\approx0.0098$) -- not statistically
distinguishable from nominal at this replicate count, but notable
because it is the only cell in this entire evaluation where F1A0B is
not the most conservative of the four methods; here it is
tied for the most liberal. Whether this is noise or a genuine small
inflation specific to a strongly misspecified half-life cannot be
resolved with 500 replicates and warrants a larger run before this
cell's power result is used to argue for F1A0B specifically (as
opposed to F1A/F1A0, which show the same power gain without this
borderline-liberal null behavior).

**3b. Grid-density stress test (does a denser AIC grid break
F1A0's size?).** F1A0 never showed Type I error inflation anywhere
in this evaluation's standard 5-point grid ($\{0, 0.2, 0.5, 1, 2\}$
weeks), which raised the question of whether that is because the
scenario the bootstrap correction is meant to guard against --
AIC overfitting the half-life under the null, when it has more
candidate values to choose from -- was never actually exercised.
Two cells test this directly: a 16-point grid and a 31-point grid
(0 to 3 weeks by 0.1), both at the null, with the true DGP half-life
(0.5 week) and everything else matching `hybrid_null`.

AIC genuinely explores the wider grid: F1A0 selects from 11
distinct half-life values at the 31-point grid (*verified*, from
the saved per-replicate data), clustering near the true DGP value
(0, 0.1, ..., 0.9, with a mode around 0.3-0.5, plus one outlier at
1.5), not collapsing to one or two points. Despite that, F1A0's size
stays at 0.040 on both the 16- and 31-point grids, matching its
5-point-grid value (0.038-0.042 elsewhere in this table) and showing
no trend toward inflation as the grid widens six-fold. F1A0B tracks
slightly below it (0.038 and 0.034), its usual conservative margin.
F1 and F1A, both unaffected by the grid widening (F1's $t_{1/2}^a$
is fixed; F1A's own size was already nominal at 0.042 throughout),
serve as the unchanged baseline.

This first attempt did not succeed at this $N=70$, this true
half-life, and these two grid densities. Two further cells push the
conditions that most plausibly would: `smalln_null` halves $N$ to 30
(less information per AIC candidate, more room to overfit), and
`offgrid_null` sets the true half-life to 0.35 week, between the
standard 5-point grid's 0.2 and 0.5, so AIC cannot land on the exact
truth even on the narrowest grid tried. A third cell combines both
with the 31-point grid (`smalln_offgrid_dense_null`), the most
adversarial condition in this evaluation.

Here the picture changes. At $N=30$, F1's nominal size stays at
0.044, but F1A rises to 0.058 and F1A0 to 0.052 -- both above
nominal, in the same direction as (though smaller than) the
half-life-misspecification null result above. F1A0B, the bootstrap-
calibrated method, pulls back to 0.044, exactly matching F1's size
and correcting the inflation seen in the uncalibrated F1A/F1A0
(*verified*, from the saved per-replicate data). At the off-grid
true half-life (0.35 week, $N=70$), no method shows inflation: F1
0.042, F1A 0.038, F1A0 0.038, F1A0B 0.036, all within noise of
nominal and of each other. The combined cell (small $N$, off-grid
truth, dense grid) reproduces the small-$N$ result almost exactly --
F1 0.044, F1A 0.056, F1A0 0.052, F1A0B 0.044 -- suggesting $N$, not
grid density or the true half-life's position relative to the grid,
is driving what signal there is.

None of these three deviations individually clears 2 Monte Carlo
standard errors ($\mathrm{SE}\approx0.0098$ at $n_{sim}=500$,
$p\approx0.05$): F1A's 0.058 and F1A0's 0.052 at $N=30$ are each
under 1 SE from nominal. But this is the first demonstration
anywhere in this evaluation of F1A0B's bootstrap calibration doing
exactly the job it is designed for -- correcting an uncalibrated
method's size back toward nominal -- rather than simply being the
most conservative of the four by a fixed margin, as it is everywhere
else in this table. Read together, the grid-density test (3b above)
and these two additional stress cells suggest that if F1A0's
uncalibrated size is vulnerable anywhere in this grid, it is to a
smaller sample size, not to a denser AIC grid or an off-grid true
half-life on their own -- and that vulnerability, where it shows up
at all, is small and not statistically established at this replicate
count.

**4. Robustness to decay-shape misspecification.** Fitting the
F-family model under an assumed exponential decay when the true DGP
decay is Weibull ($k=2$, accelerated washout) produces no detectable
Type I error inflation (0.034 to 0.040 across all four methods at
the null) and no power advantage for the adaptive methods: F1
(0.629), F1A (0.627), F1A0 (0.623) and F1A0B (0.615) are all close
together, with F1A0B again the lowest, consistent with its pattern
everywhere else in this evaluation except the half-life
misspecification alt cell above. At this one shape and one
half-life, misspecifying the decay *functional form* behaves
differently from misspecifying the *half-life value*: it does not
visibly cost F1 power, and it does not reward the adaptive methods.

**6. Robustness to DGP architecture.** Under the mean-moderation
architecture (`dgp_arch = 'mean_moderation'`) instead of the default
covariance-moderation (`mvn`), the null cell shows no inflation
(0.030 to 0.042, the same conservative-F1A0B pattern as under
`mvn`). The alternative cell shows F1 (0.614), F1A (0.622) and F1A0
(0.612) close together, with F1A0B noticeably lower at 0.561 -- a
6.1-point gap from F1A, about 2.8 Monte Carlo standard errors
($\mathrm{SE}\approx0.022$), the largest power deficit for F1A0B
relative to the other three methods anywhere in this evaluation.
Because $c_{bm}$ has a different scale under mean-moderation (a mean
shift rather than a correlation; `CLAUDE.md`, "Biomarker
parameterization"), this alt cell's absolute power is not
comparable to the `mvn` alt cells elsewhere in this table; only the
relative ordering among the four F-family methods, and the absence
of Type I error inflation, should be read from it.

**8. Computational cost.** 26.3 minutes for 2,500 (cell, replicate)
jobs on 6 cores is approximately 0.38 seconds per job amortized
across cores, or roughly 2.1 seconds serially per job (three nominal
fits plus one bootstrap resample), consistent with `docs/62`'s
reported cost for the data-driven F-family methods.

**9. Validation completeness.** This pass now checks all seven of
the docs/62 criteria that have an operational test on this driver:
1 (size), 2 (power, at two effect sizes per design, one
de-saturated), 3 (half-life misspecification, one scenario, plus a
grid-density stress test aimed specifically at F1A0's size that did
not find inflation), 4
(decay-shape misspecification, one scenario), 5 (trial design,
substituting OL+BDC for CO), 6 (DGP architecture, mean-moderation
only) and 7 (all four named prognostic pathways). Each is checked at
one effect size and 500 replicates, which is enough to see the
clearest effects (the criterion-3 power gain, the criterion-6 power
deficit) but not enough to resolve the borderline criterion-3 null
result. Criterion 8 (cost) is addressed below. No direct,
shared-grid comparison against `G8`/`G9` has been run (criterion 3
of `docs/62` Part 3's "how the two families connect"), and the dual-
channel architecture (paper 11's `combined` DGP) was not included in
the criterion-6 check.

## Caveats, unsolicited

- Single effect size per design at OL+BDC; the OL+BDC alternative
  cell is not saturated but was only checked at $c_{bm}=0.45$.
  Hybrid now has two effect sizes (0.45, saturated; 0.25,
  de-saturated), but the de-saturated gap (F1A0B 3.6 points below
  F1) is just under 2 MCSE and should be treated as suggestive, not
  established.
- 500 replicates, one pass: the gaps between F1A0B and the other
  three methods (0.006 to 0.012 at the null cells, 0.036 at the
  de-saturated Hybrid alternative) are within or close to Monte
  Carlo noise, not established differences.
- $n_{boot} = 1$ is the noisiest warp-speed setting; a larger
  per-replicate resample count (or more replicates) would tighten
  F1A0B's calibration estimate specifically.
- All four prognostic pathways were tested, but at one effect size
  each ($0.3$) and one design (Hybrid); misspecification of that
  effect size, or whether a prognostic confound and the true
  interaction co-occur, is untested.
- Update, 2026-10-26: the co-occurrence gap above was run (job 28,
  `28-run-g9-prognostic-stress.R`) at strengths 0.3 and 0.6. The
  strength-0.3 cells used in this note are positive definite on
  every Hybrid path (verified by
  `30-check-prognostic-confound-pd.R`), so this note's prognostic
  results are unaffected. The strength-0.6 cells on prog-PB and
  prog-TV-const are not: their $\Sigma$ is non-positive-definite and
  was silently repaired under `stop_on_nonpd = FALSE`, which
  attenuates the generated interaction contrast by 33-45%, so the
  power loss every method shows there (F1A0B included) is a DGP
  artifact, not a confound cost. Details in `docs/60`'s 2026-10-26
  update. Any further prognostic cell should be run with
  `stop_on_nonpd = TRUE`.
- No comparison against the `G`-family (`G1`-`G9`) on the same
  cells; this note evaluates F1A0B in isolation, not head to head.
- The half-life misspecification null result (criterion 3) needs a
  larger run before it can be read as either "fine" or "a real
  liberal bias": at 500 replicates, F1A0B's 0.060 at a 0.05 nominal
  level is not statistically distinguishable from nominal, but it is
  the only cell in this evaluation where F1A0B is not the most
  conservative method, and that reversal co-occurs with its only
  clear power advantage over F1. Both deserve a rerun at a larger
  replicate count before either is trusted.
- Criteria 3, 4 and 6 were each checked at exactly one scenario (one
  misspecified half-life, one Weibull shape, one architecture), not
  a grid; none of the three claims generalizes past the specific
  values tested.
- The grid-density stress test (3b) did not break F1A0's size at
  $N=70$ and the true DGP half-life used everywhere else in this
  table (0.5 week). A smaller $N$ (30) and an off-grid true half-life
  (0.35 week) were tried as follow-ups: the off-grid condition alone
  showed no signal, but $N=30$ showed a small, not individually
  significant rise in F1A/F1A0's size (0.052-0.058) that F1A0B's
  bootstrap correction pulled back to nominal (0.044) in both cells
  where it appeared. This is suggestive, not established, at 500
  replicates, and was not tried at any $N$ between 30 and 70, or at
  any smaller $N$ than 30.
- The dual-channel (`combined`) DGP architecture (paper 11) was not
  included in the criterion-6 check, only mean-moderation.
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
- `analysis/scripts/carryover-sensitivity/30-check-prognostic-confound-pd.R`
  -- positive-definiteness and population-contrast check for the
  prognostic-confound cells (2026-10-26).
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
