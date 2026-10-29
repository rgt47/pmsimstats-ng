# Protocol change log: 07-gompertz-evaluation
*2026-10-28 18:01 PDT*

Changes to the working ADEMP protocol (`00-ademp-protocol.md`)
for this exploratory study, recorded per the protocol's
change-log requirement. Entries are cumulative; the 2026-10-28
entries supersede the 2026-06-17 entries where they conflict.

## 2026-10-28: common random numbers (supersedes the 2026-06-17 run)

- **Seeding.** The 2026-06-17 driver
  (`02-faithful-trajectory-sweep.R`) gave each of its 16 cells its
  own seed stream, so a between-family contrast was a difference
  of two independent binomial proportions (SE about 0.019 at 1,000
  replicates). The referee report of 2026-10-28 showed that the
  0.039 mean-moderation spread it reported, and the ordering that
  placed Gompertz lowest, were Monte Carlo noise. The production
  driver is now `03-common-random-numbers-sweep.R`: the same seed
  for the same replicate index in every family and both
  architectures, set immediately before the multivariate-normal
  draw. Family contrasts are paired, and the paired MCSE of a
  power difference is sqrt(b + c)/n with b and c the discordant
  counts. Per-replicate decisions, p-values and estimates are
  saved.

- **Replicates.** 1,000 common seeds per family and architecture
  at c.bm = 0.45; 2,000 common seeds per family at c.bm = 0. The
  protocol's 5,000-replicate null pass is still not run.

- **Null cells.** At c.bm = 0 the mean-moderation and
  covariance-moderation DGPs coincide (no mean shift, no
  biomarker-response correlation), so the 2026-06-17 run's eight
  null cells were four DGPs run twice. The null is now run once per
  family, labeled architecture `both`, with a 20-seed check that
  the two architectures return identical fits.

- **Bias at the alternative.** The 2026-06-17 manuscript stated
  that no closed-form reference existed. The reference is
  -c_bm sigma_BR / sigma_bm = -0.2344 points per biomarker unit
  (the on-drug conditional slope of BR on the biomarker, and the
  visit-averaged target under trajectory scaling). Bias and its
  MCSE are now reported against it in every alternative cell.

- **Verification claim withdrawn.** The 2026-06-17 entry below
  claimed that under covariance moderation the families return
  "bit-for-bit identical rejection decisions" with matched seeds.
  No shipped code performed that check, and under 1,000 common
  seeds the families disagreed in 8 of them. The claim is
  withdrawn; the paired artifact now records the actual count.

- **Sample size.** N = 70 is the total across the design's two
  paths (35 per path), per NOTATION.md rule 4. The 2026-06-17
  entries below wrote N = 35; they are corrected in place.

- **Design name.** The executed design is the two-path OL+BDC
  design (visits at weeks 4, 8, 12, 16, 17, 18, 19, 20;
  discontinuation at week 18 or 19; no crossover). It is not the
  four-path Hybrid design, and it is not the 16-week design with a
  crossover block that the protocol's "Common infrastructure"
  section originally described. The protocol now records the
  executed design.

- **Comparators.** The von Bertalanffy curve, named as a
  comparator in earlier drafts and in the bibliography, was never
  run. It has been removed from the manuscript and the
  bibliography.

## 2026-06-17: scope of the executed study

- **Cell grid.** The protocol's Study 1 specified 144 cells
  (4 families x 3 sample sizes x 3 effect sizes x 4 carryover
  levels). The submission runs a focused factorial
  (4 families x 2 DGP architectures x 2 effect sizes) at N = 70
  total (35 per path), t_half = 0. The N, effect-size, and
  carryover axes are deferred to the 144-cell follow-up.

- **Architecture factor added.** The protocol's Study 1 grid
  did not cross DGP architecture; it varied the trajectory family
  with the other components fixed. The executed study adds
  architecture (mvn vs mean_moderation) as an explicit factor.
  This is an addition to, not a subset of, the protocol grid.

- **Replicate count for null cells.** The protocol specified
  n_reps = 5000 for type I error cells. The executed study uses
  fewer (see the 2026-10-28 entry). A 5000-replicate null pass
  remains future work.

- **Study 0 calibration.** The protocol specified an
  optim-based calibration matching the marginal on-drug
  Cor(B, BR) at week 8. The executed study calibrates each family
  to the Gompertz BR level at a week-5 anchor via uniroot or
  closed form, with the asymptotic maximum held constant. The
  week-8 correlation target is not used.

- **Study 2 (sloppiness).** Not reported in the present
  submission. The identifiability/eigenvalue diagnostic is
  deferred. The introduction's identifiability framing should be
  read as motivation for future work, not as a delivered result.

## Data-generating mechanism (interaction encoding), 2026-06-17

- **Trajectory-family routing.** An earlier prototype driver
  (`01-architecture-a-trajectory-sweep.R`) realized the
  alternative families by a post-hoc additive shift of the BR
  column on Gompertz-generated MVN data. Because the Gaussian
  covariance is mean-independent, that shift is equivalent to
  swapping the BR mean curve. The current drivers route the family
  through `buildSigma` via the `br_family` argument.

- **Trajectory-scaled mean moderation.** The prototype applied
  the mean-moderation shift as a flat per-timepoint shift of
  magnitude beta_bm * sigma_BR, independent of the response
  trajectory. The current drivers implement the multiplicative
  form of paper 01 (`01-dgp-mean-moderation-vs-mvn`),
  Y = b1 * D * (1 + b_bm * B), in which the biomarker-dependent
  component follows BR(t), normalized so the mean on-drug
  magnitude equals sigma_BR. This is the
  `moderation_scaling = "trajectory"` argument to `generateData`;
  the default `"constant"` reproduces the original behavior.

- **Verification (withdrawn 2026-10-28).** The original entry
  claimed bit-for-bit identical decisions across families under
  covariance moderation and a 0.039 spread under mean moderation.
  See the 2026-10-28 entry.

## Package changes, 2026-06-17

- `R/generateData.R`: `buildSigma` gains `br_family`, `br_p2`,
  `br_p3` (BR mean via `trajectoryShape`, default gompertz);
  `generateData` gains the same plus `moderation_scaling`
  (default `"constant"`). All defaults reproduce prior behavior.

---
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/scripts/gompertz-evaluation/02-deviations-log.md*
