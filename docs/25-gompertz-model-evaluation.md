# A critical evaluation of the Gompertz model in pmsimstats-ng

*2026-05-07 06:11 PDT*

**Author.** pmsimstats team

## 1. Motivation and scope

The pmsimstats-ng simulation framework uses a modified Gompertz
function as the parametric template for three latent response
components: the biological response (BR), the time-varying
non-drug component (TV), and the placebo-belief component (PB).
The function appears in `R/utilities.R` and in
`implementations/{original,original-extended,tidyverse}/R/`
collections; it is invoked by the simulation core at every cell
of every parameter grid. Because the assumed functional form of
the response trajectory is loaded into the simulated outcome
before any analysis-model question is asked, the choice of family
is not a peripheral implementation detail but a foundational
modelling decision. This document evaluates that decision against
the relevant biological and biostatistical literature, articulates
where the Gompertz form is well-justified for the present setting,
identifies where it imposes assumptions the framework does not
explicitly defend, and outlines the alternative families that
could be used in sensitivity analyses.

## 2. The Gompertz function

The Gompertz function in its original form is a double-exponential
sigmoid

$$
y(t) = a \cdot \exp\bigl(-b \cdot \exp(-c \cdot t)\bigr),
$$

where $a > 0$ is an asymptote, $b > 0$ controls displacement of
the inflection along the time axis, and $c > 0$ is a rate
parameter [@gompertz1825nature]. The function was introduced as a
demographic model of human mortality and was first re-purposed as
a growth curve by [@winsor1932gompertz], who observed that the
sigmoidal shape of biological size-by-time data could be captured
parsimoniously by the same three-parameter form. Subsequent work
established the Gompertz model as a standard description of tumor
growth [@laird1964dynamics; @nortonSimon1977tumor;
@norton1988gompertzian], bacterial population dynamics
[@zwietering1990modeling], and innovation diffusion. A recent
review by [@tjorve2017gompertz] catalogues parameterisations and
emphasises the non-orthogonality of $b$ and $c$, which
complicates identification when both are estimated from data.

In pharmacodynamic modelling, the dominant family is the Hill or
Emax model, $y = E_{\max} \cdot t^{n} / (EC_{50}^{n} + t^{n})$
[@holfordSheiner1981; @sheiner1989clinical; @holford2017pharmacodynamic],
which differs from the Gompertz in two respects relevant here:
the Hill model is symmetric about its inflection (Gompertz is
not), and the Hill model has a closed-form sigmoid centroid at
$EC_{50}$, which has a direct pharmacological interpretation
(half-maximal effect concentration). The Gompertz inflection
occurs at $t = \ln(b)/c$ and has no comparable interpretive
anchor in the drug-response literature.

## 3. The modified form used in pmsimstats-ng

`R/utilities.R` implements

```r
y <- maxr * exp(-disp * exp(-rate * t))
vert_offset <- maxr * exp(-disp * exp(-rate * 0))
y <- y - vert_offset
y <- y * maxr / (maxr - vert_offset)
```

This is the Gompertz function shifted vertically so that
$y(0) = 0$ and rescaled so that $y(\infty) \to$ `maxr`. The shift
plus rescale maps the family $\{$ Gompertz with arbitrary $y(0)$
and asymptote `maxr` $\}$ onto the family
$\{$ Origin-passing curves with asymptote `maxr` $\}$. The
adjustment is convenient for symptom-improvement modelling
because clinical outcome scales (CAPS, PCL, HAM-D, BDI) are
typically expressed as change-from-baseline so that 0 is the
correct value at study entry.

Two cautions about this implementation. First, the rescaling
factor $\text{maxr} / (\text{maxr} - \text{vert\_offset})$
becomes numerically unstable as `disp` shrinks because
`vert_offset` approaches `maxr`. The current code guards against
the degenerate `maxr == 0` case but not against a degenerate
`disp` close to zero. Second, the term 'modified Gompertz' in the
microbiology literature [@zwietering1990modeling] denotes a
different reparameterisation, in which a lag time $\lambda$ and a
maximum specific growth rate $\mu_m$ are explicit parameters; the
pmsimstats-ng usage is unrelated to that established
modification, and overloading the name in code comments and
manuscripts invites confusion.

## 4. How Hendrickson et al. (2020) use it

[@hendricksonOptimizing2020] use the modified Gompertz to encode
three components of the simulated symptom trajectory. The
biological response (BR) is a function of cumulative time on
drug; the time-varying non-drug component (TV) is a function of
calendar time since study entry; the placebo-belief component
(PB) is a function of cumulative time under positive expectancy.
Each component receives its own
$(\text{maxr}, \text{disp}, \text{rate})$ triple. The default
calibration is

| Component | maxr | disp | rate (per week) |
|---|---|---|---|
| BR | 10.99 | 5 | 0.42 |
| TV | 10.99 | 5 | 0.42 |
| PB |  6.51 | 5 | 0.35 |

producing curves that rise from zero at study entry, pass through
inflection near $t = \ln(5)/0.42 \approx 3.8$ weeks for BR and
TV, and approach their asymptotes by approximately week 12. The
calibration target is the prazosin-versus-placebo CAPS-difference
trajectory in the Raskind active-duty cohort
[@raskindTrialPrazosin2013; @raskindReductionNightmares2003;
@raskindParallelGroup2007], reproduced as the "Figure 4" reference
in the original code base. The BR-versus-TV separation in the
Hendrickson model is an architectural choice rather than a
mechanism testable from the published data: the Raskind trials
report only the total CAPS difference, not a decomposition into
biological-response and non-drug temporal components, so the BR
and TV maxima are calibrated jointly to a common observed target.
This identifiability constraint is acknowledged implicitly in the
pmsimstats-ng default parameter sheets but is not, to our
knowledge, evaluated empirically in the published literature.

## 5. Strengths in the symptom-trajectory setting

**Sigmoidal shape.** Drug-induced symptom improvement in chronic
neurobehavioral conditions is widely observed to be sigmoidal
rather than linear: there is a delay before measurable benefit,
followed by a phase of more rapid improvement, followed by a
plateau as the underlying pathophysiology stabilises [@holford2017pharmacodynamic].
A monotone three-parameter sigmoid is therefore a defensible
parametric template for the mean structure of such trajectories,
and the Gompertz is one of several families that can produce this
shape.

**Parsimony.** Three parameters per component is the minimum
needed to represent magnitude, timing, and rate, which is
attractive for a simulation framework that must enumerate many
parameter grids without combinatorial explosion. Richer families
(four-parameter logistic, Weibull cumulative, splines) introduce
identifiability problems when the simulation is calibrated to a
single observed trajectory.

**Computational tractability.** The Gompertz is closed-form,
infinitely differentiable, and free of the pole and near-zero
behaviours that plague the Hill function at $t = 0$. The
double-exponential form is numerically stable across the
parameter ranges used in the framework.

**Established precedent in adjacent fields.** The Gompertz has
been used for decades in tumor-growth modelling
[@laird1964dynamics; @nortonSimon1977tumor; @norton1988gompertzian]
and bacterial-growth modelling [@zwietering1990modeling], giving
the family a level of biological face validity that less commonly
used sigmoids do not enjoy.

## 6. Limitations and concerns

**Asymmetry that may not match the data.** The Gompertz is more
strongly skewed than its biological cousins (logistic, Hill,
Weibull cumulative). Its inflection occurs at $1/e$ of the
asymptote (approximately 37 percent), rather than at 50 percent
as in the logistic. Whether prazosin-induced CAPS improvement in
PTSD truly accelerates earlier than it decelerates is an
empirical question; the published Raskind trajectories
[@raskindParallelGroup2018] are not granular enough to settle the
question of skewness, and the Hendrickson framework does not
report sensitivity to the choice of sigmoid family.

**Monotonicity.** The Gompertz cannot accommodate response
followed by relapse, which is documented in PTSD pharmacotherapy
[@steckler2007ptsdpharma; @hoskins2015ptsd] and chronic
depression. Patients who improve early in treatment and then lose
benefit at week 8-10 cannot be simulated under any
$(\text{maxr}, \text{disp}, \text{rate})$ choice; the curve is
monotone non-decreasing by construction. For a simulation
framework whose downstream analyses include the OL+BDC design
(open-label titration followed by blinded discontinuation), the
inability to simulate within-treatment relapse is a meaningful
restriction. The pmsimstats-ng framework currently models
relapse only as off-drug carryover decay, which conflates two
distinct biological mechanisms (treatment escape during continued
dosing, versus pharmacological washout after discontinuation).

**Plateau is fixed.** The asymptote `maxr` is a population-level
parameter; there is no participant-level random effect on the
plateau. Heterogeneity of treatment effect at the asymptote, which
is the empirical phenomenon that motivates predictive-biomarker
work in the first place, is grafted onto the Gompertz only via
the biomarker-treatment interaction at the architecture layer
(Architecture B differential correlation, or Architecture A
mean-moderation), not as a property of the Gompertz family. The
implicit assumption is that all participants approach the same
ceiling at different rates rather than reaching different
ceilings; this is a strong assumption that is not, to our
knowledge, defended empirically.

**Non-orthogonal parameters.** The displacement parameter $b$
(`disp`) and the rate parameter $c$ (`rate`) are strongly
correlated in the Gompertz likelihood surface
[@tjorve2017gompertz]. In a calibration setting this is a
nuisance; in a simulation setting this manifests as
parameter-grid axes that are not really independent. A grid that
varies `disp` while holding `rate` fixed is also implicitly
varying the inflection time and the early-trajectory slope, which
may be undesirable for sensitivity analysis. The pmsimstats-ng
parameter-sensitivity sweeps in
`analysis/scripts/carryover-sensitivity/` and
`analysis/scripts/nof1-design-sensitivity/` do not currently
isolate these geometric features of the Gompertz curve.

**Vertical-offset rescaling can amplify noise.** The
$y - y(0)$ shift plus rescale produces a curve whose value at any
given $t$ is sensitive to small changes in `disp` because both
`vert_offset` and the rescale factor depend on it. The function
is well-behaved in the parameter regime used by
[@hendricksonOptimizing2020] but the framework does not currently
flag the regions where the rescale factor diverges.

**Departure from established pharmacodynamic practice.** The
dominant pharmacodynamic family is the Emax/Hill model
[@holfordSheiner1981; @sheiner1989clinical; @holford2017pharmacodynamic].
Reviewers familiar with the PK-PD tradition may ask why a Gompertz
was preferred to an Emax model whose half-maximal-time parameter
has a direct interpretation. The pmsimstats-ng default
calibration does not currently report a sensitivity analysis under
an Emax alternative, and the published Hendrickson paper does not
defend the Gompertz choice over a Hill or Emax form on biological
grounds.

**Misspecification relative to the analysis model.** The
analysis model in pmsimstats-ng is a linear mixed-effects model
(`nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc, ...)`) with a continuous
exposure-decayed predictor `Dbc`. Under the Gompertz DGP, the
true mean structure is a sum of three sigmoidal trajectories; the
analysis fits this with a linear-in-time term plus a
binary-or-decayed-exposure term. The mismatch is intentional, and
the simulation studies the consequences for power and Type I
error, but the mismatch does mean that the simulated power
estimates are conditioned on a specific functional form of the
DGP that is itself a modelling choice rather than an established
fact.

## 6b. Implications for the biomarker-treatment interaction test

The founding pmsimstats question is the test of the
biomarker-treatment interaction in the mixed-effects model
$\text{Sx} \sim \text{bm} + t + \text{Dbc} + \text{bm}:\text{Dbc}
+ (1|\text{ptID})$, with $H_0: \beta_{\text{bm}:\text{Dbc}} = 0$
[@hendricksonOptimizing2020]. The Gompertz-induced response
trajectory enters the simulated data through the BR component
mean. The seven limitations enumerated in §6 do not contribute
equally to the validity of this interaction test. Ranking them
by their impact on the bias, power, and coverage of
$\hat{\beta}_{\text{bm}:\text{Dbc}}$:

**Most critical for the interaction test.**

1. **Fixed plateau (no random effect on `maxr`).** Architecture
   B encodes the interaction signal in the second-moment BM-BR
   correlation; Architecture A encodes it in the first-moment
   mean shift. Both architectures presuppose that
   between-participant heterogeneity in the biomarker-driven
   treatment effect manifests in the *speed* of approach to a
   common ceiling (timing), not in the height of the ceiling
   itself (magnitude). This is a strong restriction that the
   biomarker-stratification literature implicitly contests:
   responder-versus-non-responder heterogeneity is typically
   conceptualised as different asymptotic benefit, not different
   approach rates [@kentPATH2020; @rothwell2005treating]. If the
   true biological mechanism is plateau-heterogeneous, the
   Gompertz DGP under-represents the interaction signal that the
   analysis model attempts to detect, and the simulated power
   estimates are conservative (biased downward). A reviewer who
   accepts the latent-subtype framing of Architecture B will
   reasonably ask why the response ceiling is treated as a
   population constant.

2. **Monotonicity (no within-treatment relapse).** The
   Gompertz cannot generate response-then-relapse trajectories.
   In real PTSD pharmacotherapy data
   [@hoskins2015ptsd; @steckler2007ptsdpharma] a non-trivial
   fraction of participants show early benefit followed by
   decay, often correlated with comorbid depression and chronic
   pain status. If the true symptom trajectory has this shape
   and the simulation does not, the linear-in-time analysis
   model is being asked to fit a smoother surface than it would
   encounter in practice. The interaction test is then
   evaluated under conditions kinder than the application; this
   biases the simulated power upward (optimistic). Combined
   with point (1), the directions of bias point opposite ways
   and partially offset, but neither offset is principled.

3. **Misspecification relative to the linear analysis model.**
   The analysis model fits `bm + t + Dbc + bm:Dbc`. The true BR
   mean is sigmoidal, not linear in $t$. The linear-in-$t$ term
   absorbs trajectory curvature into a single slope, and any
   residual curvature is pushed into the random intercept and
   the AR(1) residual structure. The interaction term
   `bm:Dbc` is the residual quantity that survives this
   absorption; it is therefore the term most exposed to the
   choice of sigmoid family. A Hill-shaped DGP (more symmetric
   inflection) and a Gompertz-shaped DGP (asymmetric inflection
   at 1/e of the asymptote) push different amounts of curvature
   into the linear approximation and may yield different
   apparent `bm:Dbc` standard errors at matched calibration
   targets. The framework currently does not report sensitivity
   under an alternative sigmoid family; this is the single most
   actionable Gompertz-related contribution to interaction-test
   validity.

**Moderate impact on the interaction test.**

4. **Asymmetric inflection.** The Gompertz inflection at
   $t = \ln(b)/c$ corresponds to $1/e \approx 0.37$ of the
   asymptote, in contrast to the logistic $0.5$. Under matched
   `maxr` and matched inflection time, the early-trajectory
   slope is steeper for Gompertz than logistic, while the late
   trajectory is flatter. Because the BM-BR correlation under
   Architecture B is anchored at on-drug timepoints with
   `tod > 0`, the interaction signal is concentrated in the
   period where the BR curve is most rapidly changing. A
   Gompertz BR therefore weights the interaction signal toward
   early-trajectory observations more heavily than a logistic BR
   would, and the simulated power has a different mix of
   contributions across timepoints. The qualitative conclusions
   are likely robust; the absolute power numbers are not.

5. **Non-orthogonal `disp` and `rate`.** The displacement and
   rate parameters are correlated in the Gompertz likelihood
   surface [@tjorve2017gompertz]. For interaction-test
   simulation, this matters when sensitivity sweeps vary one of
   these parameters. A sweep that walks `rate` while holding
   `disp` fixed implicitly co-varies the inflection time and
   the early slope; the resulting effect on the interaction
   power is a mixture of two geometric properties of the curve,
   not a clean single-axis sensitivity. The carryover-sensitivity
   and nof1-design-sensitivity papers do not currently
   decompose this. A cleaner reparameterisation (Tjørve's
   alternative form, or fixing the inflection time as the
   primary parameter and deriving `disp` as a function of it)
   would produce sensitivity sweeps with clearer interpretation.

**Low impact on the interaction test.**

6. **Vertical-offset rescaling instability.** A numerical
   concern at extreme parameter values. The default
   $(\text{disp}=5, \text{rate}=0.42)$ regime is far from the
   instability region; the interaction test is not affected at
   the published parameter settings. A guard in the code is
   prudent but the issue does not undermine current results.

7. **Naming overlap with Zwietering's modified Gompertz.** A
   documentation concern, not a statistical concern. The
   interaction test is unaffected.

8. **Departure from Emax/Hill pharmacodynamic tradition.** A
   reviewer-facing concern about defensibility of the choice. The
   interaction test is unaffected at fixed parameters; it would
   change under an Emax DGP only through the same mechanism as
   point (3) above.

**Summary.** The Gompertz-related issues that most warrant
attention before the next interaction-test manuscript revision
are the fixed plateau (point 1), the monotonicity assumption
(point 2), and the absence of a sensitivity sweep over sigmoid
family (point 3). Points 4 and 5 are second-order. Points 6
through 8 are housekeeping. A direct sensitivity sweep at fixed
$(\text{maxr}, \text{inflection time})$ across Gompertz, logistic
4PL, and Hill sigmoids would address point 3 and is the highest-
value next experiment for the founding interaction-test question.

## 7. Alternative parametric families to consider

| Family | Strength | Weakness | Reference |
|---|---|---|---|
| Logistic (4PL) | Symmetric inflection at 50% asymptote | Same monotonicity problem as Gompertz | [@richards1959flexible] |
| Hill / Emax | Closed-form pharmacological interpretation; standard in PK-PD | Pole at $t = 0$ for some forms; less natural for time courses than for dose-response | [@holfordSheiner1981] |
| Weibull cumulative | Tunable shape parameter, recovers exponential at $k=1$ | Three-parameter and identifiability is similar to Gompertz | [@weibull1951statistical] |
| Bi-exponential | Two-phase response (rapid plus slow) | Non-monotone rate of approach; estimation can be unstable | [@holford2017pharmacodynamic] |
| Spline / GAM mean structure | Non-parametric, captures relapse and treatment escape | Loses the parsimony argument; harder to calibrate to a single target trajectory | [@brumback1998] |
| Latent change score / GMM | Allows participant-level trajectory heterogeneity in shape | Heavy machinery; identifiability requires more timepoints than a typical N-of-1 trial provides | [@muthen2001growth] |
| Stochastic differential equation | Endogenous treatment-escape and relapse | Loses closed-form simulation; substantially heavier compute | [@oksendalSDE] |

## 8. Recommendations for pmsimstats-ng

1. **Document the Gompertz choice as a modelling assumption,
   not a mechanism.** The pmsimstats-ng manuscripts (notably the
   companion to this document at
   `analysis/report/01-dgp-mean-moderation-vs-mvn/` and the
   sensitivity papers at `analysis/report/02-carryover-sensitivity/`
   and `analysis/report/05-nof1-design-sensitivity/`) should state
   in the simulation-design section that the Gompertz is chosen
   for parsimony and tractability, that its sigmoidal shape is
   defensible but not the only defensible choice, and that the
   power estimates are conditioned on this functional form.

2. **Add a Gompertz-versus-alternative sensitivity sweep.** A
   direct test of robustness would replace the Gompertz BR with a
   logistic 4PL and a Hill-form sigmoid at matched
   asymptote-and-inflection coordinates and report the resulting
   power changes. The infrastructure for such a sweep already
   exists in the
   `analysis/scripts/{carryover,nof1-design}-sensitivity/`
   pipelines and would require an additional axis in the cell
   tibble plus a small change to `simulation-core.R` to dispatch
   on response family.

3. **Rename `modgompertz` or document the deviation from the
   Zwietering modified form.** Either the function name should be
   changed to reflect the actual modification (vertical-offset and
   rescale to origin) or the in-code roxygen documentation should
   make the deviation from Zwietering explicit. At present, a
   reader who recognises 'modified Gompertz' from the
   bacterial-growth literature will encounter a different
   parameterisation under the same name.

4. **Guard the rescale factor.** Add an explicit check on
   `(maxr - vert_offset)` for near-zero values and either issue a
   warning or fall back to the unmodified Gompertz at low `disp`.

5. **Consider participant-level variation in `maxr`.** A
   random-effects extension in which each simulated participant
   draws $\text{maxr}_i \sim N(\text{maxr}_0, \tau^2)$ would
   accommodate plateau heterogeneity without abandoning the
   Gompertz family. This would also bring the framework into
   better alignment with the latent-class and finite-mixture
   formulations developed in the companion paper at
   `analysis/report/03-latent-class-mixture-application/`.

6. **Acknowledge response-then-relapse as a scope exclusion.**
   Until the framework supports non-monotone response trajectories,
   manuscripts should state explicitly that the simulation does
   not capture mid-treatment loss of benefit, and that the
   reported power estimates may be optimistic for clinical
   populations in which such loss is documented.

## 9. Summary

The Gompertz model is a defensible parametric template for the
simulated symptom trajectories in pmsimstats-ng. It is
parsimonious, numerically tractable, and has a long pedigree in
adjacent biological domains. It is also a strong assumption: it
imposes monotone non-decreasing trajectories with a fixed
population asymptote and a specific asymmetric inflection that may
not match the empirical PTSD or chronic-pain symptom-trajectory
literature. The current pmsimstats-ng manuscripts use the
Gompertz without a sensitivity comparison against the
established Emax / Hill alternative or against logistic and
Weibull-cumulative families. We recommend that future revisions
either defend the Gompertz choice empirically against a candidate
alternative or report a sensitivity sweep that demonstrates the
power conclusions are robust to the choice of sigmoid family. The
Gompertz alone should not bear the weight of a methodological
recommendation about trial design.

## References

- [@gompertz1825nature] Gompertz B. On the nature of the function
  expressive of the law of human mortality. *Philosophical
  Transactions of the Royal Society of London*. 1825;115:513-583.
- [@winsor1932gompertz] Winsor CP. The Gompertz curve as a growth
  curve. *Proceedings of the National Academy of Sciences*.
  1932;18(1):1-8.
- [@laird1964dynamics] Laird AK. Dynamics of tumor growth.
  *British Journal of Cancer*. 1964;13:490-502.
- [@nortonSimon1977tumor] Norton L, Simon R. Tumor size,
  sensitivity to therapy, and design of treatment schedules.
  *Cancer Treatment Reports*. 1977;61(7):1307-1317.
- [@norton1988gompertzian] Norton L. A Gompertzian model of human
  breast cancer growth. *Cancer Research*. 1988;48(24
  Pt 1):7067-7071.
- [@zwietering1990modeling] Zwietering MH, Jongenburger I, Rombouts
  FM, van 't Riet K. Modeling of the bacterial growth curve.
  *Applied and Environmental Microbiology*.
  1990;56(6):1875-1881.
- [@tjorve2017gompertz] Tjørve KMC, Tjørve E. The use of Gompertz
  models in growth analyses, and new Gompertz-model approach.
  *PLOS ONE*. 2017;12(6):e0178691.
- [@holfordSheiner1981] Holford NHG, Sheiner LB. Understanding the
  dose-effect relationship: clinical application of
  pharmacokinetic-pharmacodynamic models. *Clinical
  Pharmacokinetics*. 1981;6(6):429-453.
- [@sheiner1989clinical] Sheiner LB. Clinical pharmacology and the
  choice of model functions for predicting drug effects.
  *Clinical Pharmacology and Therapeutics*. 1989;46(6):605-615.
- [@holford2017pharmacodynamic] Holford NHG. Pharmacodynamic
  principles and target concentration intervention.
  *Translational and Clinical Pharmacology*. 2017;26(4):150-154.
- [@richards1959flexible] Richards FJ. A flexible growth function
  for empirical use. *Journal of Experimental Botany*.
  1959;10(2):290-301.
- [@weibull1951statistical] Weibull W. A statistical distribution
  function of wide applicability. *Journal of Applied Mechanics*.
  1951;18:293-297.
- [@brumback1998] Brumback BA, Rice JA. Smoothing spline models
  for the analysis of nested and crossed samples of curves.
  *Journal of the American Statistical Association*.
  1998;93(443):961-976.
- [@muthen2001growth] Muthén B. Latent variable mixture modeling.
  In: Marcoulides GA, Schumacker RE, eds. *New Developments and
  Techniques in Structural Equation Modeling*. Lawrence Erlbaum;
  2001:1-33.
- [@hendricksonOptimizing2020] Hendrickson RC, Thomas RG, Schork
  NJ, Raskind MA. Optimizing aggregated N-of-1 trial designs for
  predictive biomarker validation: statistical methods and
  practical considerations. *Frontiers in Digital Health*.
  2020;2:13.
- [@raskindReductionNightmares2003] Raskind MA, Peskind ER, Kanter
  ED, et al. Reduction of nightmares and other PTSD symptoms in
  combat veterans by prazosin: a placebo-controlled study.
  *American Journal of Psychiatry*. 2003;160(2):371-373.
- [@raskindParallelGroup2007] Raskind MA, Peskind ER, Hoff DJ, et
  al. A parallel group placebo controlled study of prazosin for
  trauma nightmares and sleep disturbance in combat veterans with
  post-traumatic stress disorder. *Biological Psychiatry*.
  2007;61(8):928-934.
- [@raskindTrialPrazosin2013] Raskind MA, Peterson K, Williams T,
  et al. A trial of prazosin for combat trauma PTSD with
  nightmares in active-duty soldiers returned from Iraq and
  Afghanistan. *American Journal of Psychiatry*.
  2013;170(9):1003-1010.
- [@raskindParallelGroup2018] Raskind MA, Peskind ER, Chow B, et
  al. Trial of prazosin for post-traumatic stress disorder in
  military veterans. *New England Journal of Medicine*.
  2018;378(6):507-517.
- [@steckler2007ptsdpharma] Steckler T, Risbrough V. Pharmacological
  treatment of PTSD: established and new approaches.
  *Neuropharmacology*. 2012;62(2):617-627.
- [@hoskins2015ptsd] Hoskins M, Pearce J, Bethell A, et al.
  Pharmacotherapy for post-traumatic stress disorder: systematic
  review and meta-analysis. *British Journal of Psychiatry*.
  2015;206(2):93-100.
- [@oksendalSDE] Øksendal B. *Stochastic Differential Equations:
  An Introduction with Applications*. 6th ed. Springer; 2003.

---

*Rendered on 2026-05-07 at 06:23 PDT.*<br>
*Source: ~/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/25-gompertz-model-evaluation.md*
