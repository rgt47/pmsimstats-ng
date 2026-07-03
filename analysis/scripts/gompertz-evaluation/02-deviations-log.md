# Deviations log: 07-gompertz-evaluation
*2026-06-17 09:26 PDT*

Departures from the pre-registered ADEMP plan
(`00-ademp-pre-registration.md`), recorded per the plan's
deviations-log requirement.

## Scope of the executed study

- **Cell grid.** Pre-registration Study 1 specified 144 cells
  (4 families x 3 sample sizes x 3 effect sizes x 4 carryover
  levels). The submission runs a focused 16-cell factorial
  (4 families x 2 DGP architectures x 2 effect sizes) at N = 35,
  t_half = 0. The N, effect-size, and carryover axes are deferred
  to the 144-cell follow-up.

- **Architecture factor added.** The pre-registered Study 1 grid
  did not cross DGP architecture; it varied the trajectory family
  with the other components fixed. The executed study adds
  architecture (mvn vs mean_moderation) as an explicit factor,
  because the family's inferential consequence turns out to depend
  on it. This is an addition to, not a subset of, the
  pre-registered grid.

- **Replicate count for null cells.** Pre-registration specified
  n_reps = 5000 for type I error cells. The executed study uses
  n_reps = 1000 for all cells, including the null cells (cell MCSE
  on type I error approximately 0.005). A 5000-replicate null pass
  remains future work.

- **Study 0 calibration.** Pre-registration specified an
  optim-based calibration matching the marginal on-drug
  Cor(B, BR) at week 8. The executed study calibrates each family
  to the Gompertz BR level at a week-5 anchor via uniroot or
  closed form, with the asymptotic maximum held constant. The
  simpler single-anchor calibration is sufficient for the
  matched-marginal-effect comparison; the week-8 correlation
  target is not used.

- **Study 2 (sloppiness).** Not reported in the present
  submission. The identifiability/eigenvalue diagnostic is
  deferred. The introduction's identifiability framing should be
  read as motivation for future work, not as a delivered result.

## Data-generating mechanism (interaction encoding)

- **Trajectory-family routing.** An earlier prototype driver
  (`01-architecture-a-trajectory-sweep.R`) realised the
  alternative families by a post-hoc additive shift of the BR
  column on Gompertz-generated MVN data. Because the Gaussian
  covariance is mean-independent, that shift is mathematically
  equivalent to swapping the BR mean curve, so the family did
  drive the BR mean. The current driver
  (`02-faithful-trajectory-sweep.R`) instead routes the family
  through `buildSigma` via the new `br_family` argument, which is
  cleaner, carryover-safe, and verifiable.

- **Trajectory-scaled Architecture-A moderation (material).** The
  prototype applied the Architecture-A biomarker moderation as a
  flat per-timepoint shift of magnitude beta_bm * sigma_BR,
  independent of the response trajectory. Under that encoding the
  trajectory family could not affect the bm:Dbc interaction under
  either architecture, and the prototype's apparent
  family-invariance was an artefact of the encoding rather than a
  property of trajectory family. The current driver implements
  Architecture A in the multiplicative form of paper 01
  (`01-dgp-mean-moderation-vs-mvn`), Y = b1 * D * (1 + b_bm * B),
  in which the biomarker-dependent component rides the response
  trajectory BR(t), normalised so the mean on-drug magnitude
  equals sigma_BR (matched marginal effect size across families).
  This is implemented in the package via the new
  `moderation_scaling = "trajectory"` argument to `generateData`;
  the default `"constant"` reproduces the original behaviour
  exactly, so all other pmsimstats runs are unchanged.

- **Verification.** With matched random-number seeds the four
  families return bit-for-bit identical rejection decisions under
  Architecture B (the covariance channel is exactly
  mean-independent), confirming that the Architecture-B
  invariance is structural rather than a small-effect artefact.
  Under the faithful Architecture A the families separate over a
  0.039 power range at N = 35, c.bm = 0.45.

## Package changes

- `R/generateData.R`: `buildSigma` gains `br_family`, `br_p2`,
  `br_p3` (BR mean via `trajectoryShape`, default gompertz);
  `generateData` gains the same plus `moderation_scaling`
  (default `"constant"`). All defaults reproduce prior behaviour;
  backward-compatibility checked against the committed original
  cell powers.

---
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/scripts/gompertz-evaluation/02-deviations-log.md*
