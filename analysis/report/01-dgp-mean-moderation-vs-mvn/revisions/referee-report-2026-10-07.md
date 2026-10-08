# Referee Report: Simulating Predictive-Biomarker Interactions in Aggregated N-of-1 Trials
*2026-10-07 09:38 PDT*

**Manuscript.** 'Simulating predictive-biomarker interactions in
aggregated N-of-1 trials: a refined data-generating process following
Hendrickson et al.' (`report.Rmd`, rendered `report.pdf`).

**Journal and article type.** *Trials* (BMC), Methodology.

**Recommendation.** Major revision.

**How the claims were checked.** Every numerical claim was traced to
the evidence documents named by the authors (`docs/34` to `docs/48`),
to the stress-test output (`decision-screen.csv` and
`summary-full.csv`), to the closed-form outputs under
`analysis/data/quick-sim/carryover-closed-form/`, and to the published
code at commit 58b32a9
(`analysis/data/quick-sim/hendrickson-58b32a9/R/`). Labels used below:
*verified* (recomputed or matched to a source file), *inspected* (read
in code), *inferred* (follows from verified quantities), *unverified*
(no source found). Line references are to `report.Rmd`.

---

## 1. Summary

The manuscript revisits the simulation framework of Hendrickson et al.
(2020, 'RH'), which represents a predictive biomarker as a correlation
$c_{bm}$ between the biomarker and the drug-response component of a
joint multivariate normal (MVN). Using a paired-difference
('on-minus-off') summary statistic whose moments are available in
closed form, the authors show that (i) RH's coupling rule, which gives
any visit with a nonzero drug-response mean the full correlation,
makes the coupling constant across on- and off-drug visits once any
carryover is present, so the interaction slope collapses to zero (or
to $1/(n_{\text{off}}+1)$ of its value when baseline is counted off
drug); (ii) coupling that decays with the remaining drug effect
('graded coupling') yields a smooth decline in power; (iii) under
compound symmetry the joint distribution exists only for $c_{bm}$
below about 0.256 in the Hybrid design without carryover, so RH's
effect sizes of 0.3 and 0.6 were generated from silently repaired
matrices; and (iv) a biomarker that is also prognostic inflates the
size of analyses that share the biomarker main effect between
baseline and post-baseline off-drug visits. The paper proposes
corrections to RH's code (Appendix B) and a refined data-generating
process (DGP) with AR(1) and separable correlation structures,
anchored carryover, and optional prognostic pathways.

## 2. Overall assessment

The diagnostic core of the paper is correct and useful. I verified
Propositions 1 and 2 algebraically, the Schur-complement ceiling, the
derivative in A.3, every entry of Table 2 (closed-form power), Table 3
(published-code patch), Table 4 (ceilings), the coupling-gap values
and the Gompertz landmarks. The quoted code in B.1 matches lines 131
to 138 of the 58b32a9 source exactly, and the replacement in B.1 is
correct given the published `buildtrialdesign()` semantics, in which
`tsd` is zero before first exposure. Identifying the step coupling
rule as the cause of RH's power collapse is a genuine and practically
important finding for anyone reusing `pmsimstats`.

The paper is not yet ready for publication, for five reasons that
are developed below. First, it frames covariance-based and mean-based
binding as two different encodings and credits RH with a 'first'.
Under an MVN the two are conditionally the same model and differ only
in residual covariance. Stating this would simplify the paper and
change several of its conclusions. Second, the abstract and
conclusions overstate what was proved (for example 'provably
monotone', 'power falls to the test size', 'exact explanations').
Third, the ceiling section contains a statement that contradicts its
own table, confounds a change of $\rho$ with a change of structure,
and omits the fact that the recommended AR(1) structures halve power
at a given $c_{bm}$. Fourth, several stress-test numbers are
misreported or reported selectively relative to the
pre-registration. Fifth, the 'refined DGP' that the title and third
aim promise is never fully specified in one place, is only partly
implemented, and is not checked against data.

---

## 3. Major issues

### M1. Covariance and mean binding are conditionally equivalent under the MVN; the framing, the 'first' claim and the comparison need revision

**Location.** Abstract, Background (lines 37 to 41); Background, 'The
framework of Hendrickson et al.' (lines 120 to 139); Methods,
'Covariance-based and mean-based binding' (lines 332 to 350);
Discussion, 'Covariance-based and mean-based binding' (lines 547 to
611); Results, 'Room for realistic effect sizes' (lines 493 to 497).

**Problem.** Let $b$ be the standardized biomarker, and suppose that
$\mathrm{Cor}(B, BR_t) = c_{bm} g_t$, with $B$ uncorrelated with $TV$,
$PB$ and $\mathrm{BL}$, as in RH's code. Under joint normality the
conditional law of the response block given $B$ is normal, with

$$
E[BR_t \mid b] = \mu_{BR,t} + c_{bm}\,\sigma_{BR}\, g_t\, b, \qquad
\mathrm{Cov}(R \mid b) = \Sigma_{RR} - c_{bm}^2 \sigma_{BR}^2\, g g^\top
\ \text{(in the BR block)},
$$

and all other conditional moments are unchanged. The covariance
encoding with coupling profile $g$ is therefore exactly a
mean-moderation model with exposure profile $h = g$ and slope
$\beta_{bm} = c_{bm}$, in which the residual BR covariance has been
reduced by the rank-one term $c_{bm}^2\sigma_{BR}^2 g g^\top$. Three
consequences follow, and the manuscript does not state any of them.

1. The result that 'mean-based binding that fades with the drug effect
   loses power almost identically' (Abstract; lines 594 to 603) is
   close to an identity rather than an empirical finding. Graded
   coupling and decayed mean moderation have the same conditional
   mean. The residual power difference (0.721 against 0.681 at
   $t_{1/2} = 0$) comes entirely from the convention that holds the
   marginal BR variance fixed under covariance binding and inflates it
   by $\beta_{bm}^2 h_t^2 \sigma_{BR}^2$ under mean binding. Listing this
   difference as a 'strength' of the covariance encoding (lines 565 to
   567) is therefore an artifact of parameterization. Matching the
   residual variances would make the two identical.
2. The ceiling is the requirement that the residual covariance
   $\Sigma_{RR} - c_{bm}^2\sigma_{BR}^2 g g^\top$ remain positive
   definite. This gives a sharper explanation than the one offered at
   lines 481 to 491. Under compound symmetry with $\rho = 0.8$, the
   within-person contrast directions of the BR block have eigenvalue
   near $1 - \rho = 0.2$. An on/off coupling profile loads on exactly
   those directions, which is why the ceiling is low whenever the
   profile has contrast and high (0.893) when it is constant.
3. The statement that effects as large as 0.6 'can be represented only
   by binding the biomarker in the mean' (lines 495 to 497) holds only
   when $\sigma_{BR}$ is fixed. The same conditional slope is
   representable under covariance binding with a larger $\sigma_{BR}$.
   The real constraint is on the share of BR variance the biomarker
   explains, $c_{bm}^2 g_t^2$.

Given this equivalence, the claim that RH were 'to our knowledge, the
first in the N-of-1 literature to encode the biomarker-by-treatment
interaction as a correlation' (lines 38 to 41 and 132 to 134)
distinguishes a parameterization rather than a model. It is also
unverified. Simulating a covariate correlated with a random
participant-by-treatment effect is a standard device in the
random-coefficient and joint-model literature, and the manuscript
reports no systematic search.

**Fix.** State the conditional equivalence as a proposition in the
Methods (with a short proof in Appendix A), and use it to organize the
comparison. The covariance and mean encodings are then two members of
one family, indexed by the exposure profile $h$ and the residual
covariance. Recast the 0.721 against 0.681 comparison as a residual
variance convention. Either remove the priority claim or replace it
with a narrower and documented statement (for example, 'RH
parameterized the moderation through the joint covariance rather than
through a mean term'), supported by a described literature search.

### M2. The abstract and conclusions overstate what was proved and shown

**Location.** Abstract, Results (lines 52 to 58); Results, 'Graded
coupling restores a smooth decline' (lines 435 to 436); Conclusions
(lines 652 to 659); A.3 (lines 769 to 772).

**Problem.**

- 'The interaction vanishes, and power falls to the test size at every
  positive half-life' (line 53 to 55). This holds for the post-baseline
  statistic (E9) only. For the baseline-inclusive statistic (E9b), which
  the authors use as the proxy for RH's analysis, power is 0.081. For RH's
  mixed model it is 0.08 to 0.11 (reconstruction) and 0.092 to 0.126
  (published code, Table 3). The slope does not vanish; it is
  $1/(n_{\text{off}}+1)$ of its no-carryover value (line 385 to 387).
- 'Restores a smooth, provably monotone decline' (line 56). In
  `docs/37` (Section 5.2), Proposition 3 is proved for the
  *path-stratified* statistic. For the pooled statistic with a common
  intercept, monotonicity was checked numerically on a 120-point grid,
  not proved. For RH's mixed model it is neither proved nor checked.
  The sentence then attaches 'provably' to the published-code result
  (0.13 to 0.64), whose graded column is not monotone (0.668, 0.638,
  0.656, 0.462). That non-monotonicity is within Monte Carlo error, but
  the juxtaposition misleads. The A.3 statement of Proposition 3 also
  omits the path-stratification qualifier.
- 'Both problems have exact explanations' (line 653). The explanation
  is exact for a summary statistic, and only approximate for RH's
  analysis (the Limitations concede this at lines 634 to 637).
- 'The refined process ... is a sound reference for simulation studies
  of these designs' (lines 657 to 659). The refined process was not
  validated against trial data, parts of it are unimplemented (see M9),
  and its parameter values are not estimated (lines 638 to 640).

**Fix.** Revise the abstract to: 'power falls to near the test size
(0.05 for the post-baseline statistic, 0.08 for the baseline-inclusive
statistic and RH's analysis)'. Restrict 'provably' to the
path-stratified statistic and say that the pooled and mixed-model
results were checked numerically. Replace 'exact explanations' with
'closed-form explanations for a summary statistic that tracks the
published analysis'. Soften the final sentence of the Conclusions to a
proposal.

### M3. The ceiling section contradicts its own table, confounds $\rho$ with structure, and omits the power cost of AR(1)

**Location.** Abstract (lines 59 to 62); Results, 'Room for realistic
effect sizes' (lines 476 to 506); Table 1, step 4 (line 299);
Recommendations 2 (lines 619 to 620).

**Problem.**

1. **A statement contradicts the table.** Lines 485 to 487 state that
   'AR(1) correlation alone therefore lowers the ceiling'. Table 4 shows
   the opposite in every design shown. At $\rho = 0.8$ the AR(1)
   column (0.291) exceeds compound symmetry (0.256), and at $\rho = 0.7$
   it does so as well (0.452 against 0.346). According to `docs/36`
   (Section 6.5), AR(1) raises the ceiling where the coupling has an
   on/off contrast and lowers it only where the coupling is constant
   (open-label; step-rule carryover cells). The sentence must be
   corrected.
2. **A change of $\rho$ is presented as a change of structure.** The
   abstract contrasts 0.256 (compound symmetry at RH's $\rho = 0.8$)
   with 'about 0.48' (separable AR(1) at $\rho = 0.7$). At equal
   $\rho = 0.8$ the separable minimum is 0.440 (Hybrid 0.440). At
   $\rho = 0.7$ most of the gain comes from AR(1) itself (0.346 to
   0.452). Separability adds only 0.022 (to 0.474). Lowering $\rho$
   from 0.8 to 0.7 raises even the compound-symmetry ceiling from 0.256
   to 0.346. The paper gives no justification for changing $\rho$.
3. **A feasibility claim fails at RH's $\rho$.** 'A biomarker correlation
   of 0.45 is feasible in every design under the separable form'
   (line 493) is false at $\rho = 0.8$ (minimum 0.440). It is true at
   $\rho = 0.7$, but there the non-separable AR(1) form (0.452) also
   meets it.
4. **The power cost is omitted.** At a fixed $c_{bm} = 0.25$, Hybrid E9
   power without carryover is 0.721 under compound symmetry, 0.326
   under AR(1) and 0.335 under separable AR(1) (`docs/37`, Section 5.8,
   verified). The contrast variance is 2.6 times larger in Hybrid and
   up to 5.9 times larger in CO. The 'room for realistic effect sizes'
   is therefore partly illusory: $c_{bm}$ does not mean the same thing
   across structures, and a reader who adopts Recommendation 2 at the
   same $c_{bm}$ will see power fall by more than half without being
   told why.
5. **AR(1) in weeks has its own implausibility.** With $\rho = 0.7$ per
   week and no person-level component, the correlation of a component
   between weeks 4 and 20 is $0.7^{16} \approx 0.003$. Pure AR(1) thus
   understates long-lag correlation about as badly as compound symmetry
   overstates it (Table 1, step 4). Symptom scales in PTSD trials
   typically show substantial correlation over months. The DGP in the
   stress test (base cell) has no random intercept in the generating
   covariance.

**Fix.** Correct the AR(1) sentence. Report ceilings at equal $\rho$
for all three structures, and justify $\rho = 0.7$ or keep 0.8. Report
power at a common $c_{bm}$ and at a common implied standardized effect
($\rho_{\Delta B}$) across structures. Add a structure with both a
person-level and a serial component, for example
$\phi + (1 - \phi)\rho^{|t-s|}$, and give its ceiling and power.
Ideally, calibrate $\rho$ and $\phi$ against repeated symptom scores
from a prazosin trial (Raskind 2013 or 2018).

### M4. Stress-test numbers are misreported or reported selectively, and deviations from the pre-registration are not disclosed

**Location.** Methods, 'Simulation studies' item 3 (lines 366 to 372);
Results, 'Consequences for the analysis' (lines 508 to 543).

**Problem.** Checked against `decision-screen.csv` and
`summary-full.csv` (verified):

| Manuscript claim | Source value | Status |
|---|---|---|
| 'The published random-intercept model rejected a true null 0.107 to 0.159' (lines 524 to 525) | R: 0.107 (TV constant), 0.053 (TV growing), 0.159 (PB) | **Incorrect range**; under a growing natural-course association R was nominal (0.053). Correct range 0.053 to 0.159 |
| 'An unstructured-covariance mixed model rejected up to 0.42' (lines 525 to 526) | A: 0.416, 0.195, 0.222; DA (unstructured with AIC-chosen half-life): **0.509**, 0.245, 0.256 | Selective; the worst unstructured analysis reached 0.51 |
| 'Analyses that free the biomarker's baseline association held the nominal level (0.040 to 0.063)' | ES 0.049/0.063/0.049; EA 0.049/0.057/0.050; F1 0.047/0.046/0.040 | Verified; 0.063 (ES, TV growing) is 1.9 binomial SE above 0.05 |
| Published analysis 0.069 on average, up to 0.081 at $N = 35$ | R pooled 0.0687, max 0.081 (N35) | Verified |
| Unstructured KR 0.087 at $N = 35$ | A max 0.087 (N35) | Verified |
| CR2 random-intercept models 0.036 to 0.069 | S min 0.036, C1 max 0.069 (N35) | Verified |
| 39 cells, 26,500 trials, 1,000 null replicates per cell | 39 cells; 26,500 fits per analysis; null cells 1,000 | Verified |

Further problems:

- **Omitted failure.** The strawman S, a CR2 random-intercept model
  with `corCAR1` residuals, which the text implicitly commends for size
  (lines 537 to 539), rejected 0.099 to 0.100 under the TV-constant and
  PB prognostic scenarios. DS rejected 0.112 to 0.115. The text names
  only R and A as failing under a prognostic biomarker.
- **Inconsistent standards.** R's pooled 0.069 is presented as
  anticonservative, while C1's maximum of 0.069 is presented as
  'staying' within bounds. With 1,000 replicates the binomial SE is
  0.0069, so 0.069 is 2.75 SE above nominal in either case.
- **Undisclosed deviation from the pre-registration.** `docs/46` was
  written on 2026-10-06 at 19:34 and the runs were made on 2026-10-06
  and 07. It is an internal document with no public timestamp. The
  pre-registered symmetry screen failed 9 of 22 tests through Monte
  Carlo noise, and the as-registered selection was C2 (flagged) or F1
  (`docs/48`, Section 4). The manuscript calls the study
  'pre-registered' but does not report the registered decision rule,
  the as-registered outcome, or the deviation.
- **No recommended analysis.** The paper never states which analysis
  it recommends. `docs/48` recommends F1, and the Recommendations
  (lines 615 to 630) say nothing about the analysis model.
- **Design and half-life left unstated.** The base half-life of 0.5
  week lies outside the authors' own 'pharmacological range' of 0.1 to
  0.2 week (lines 314 to 316), and the analyses' assumed half-life
  equals the true one in the base cell. The design (Hybrid) is not
  stated in the Results.

**Fix.** Correct the R range and report DA. Present a table of null
rejection rates (analyses by cells, with MCSE). Apply one size
criterion uniformly, for example the Bradley liberal interval or
$0.05 \pm 2$ MCSE. State the registration status accurately ('a
protocol fixed internally before the runs'), and report the
as-registered decision and the deviation. State the recommended
analysis and its power cost (ES against S: 0.364 against 0.514 under
AR(1); `docs/48`).

### M5. The prognostic-biomarker mechanism (A.4) is heuristic and does not explain the observed pattern

**Location.** Results, lines 514 to 529; Appendix A.4 (lines 787 to
802); Recommendation 5.

**Problem.**

- **The pattern is not explained.** A.4 asserts
  $\hat\beta_{bm:D} \approx (1 - f)s_{\text{prog}}$ with $f$ 'set by the
  weight the baseline row receives'. That is a description, not a
  derivation. It does not explain why R and S are nominal (0.053,
  0.049) when the natural-course association grows over time, although
  $s_0 = 0 \ne s_{\text{off}}$ there as well. The explanation in
  `docs/48` (timing of on-drug visits relative to a growing slope)
  implies that the bias depends on how drug status is confounded with
  time. In that case the claim that $B_i\mathbf{1}\{t = 0\}$ 'removes the
  bias' is exact only when the post-baseline prognostic slope is
  constant across visits. ES's 0.063 under the growing scenario is
  consistent with residual bias.
- **The PB case has a second mechanism.** Under PB-prognostic coupling
  the PB SD is $10e_t$, with $e_t = 1$ at open-label (all on-drug)
  visits and 0.5 at blinded visits. So $s_{\text{on}} \ne s_{\text{off}}$
  even without a predictive effect, because of expectancy, not only
  because of baseline sharing. The assumption
  $s_{\text{off}} = s_{\text{on}} = s_{\text{prog}}$ in A.4 is false for
  this scenario.
- **A plausible scenario is missing.** No scenario lets the biomarker
  correlate with baseline severity ($s_0 \ne 0$). Blood pressure and
  PTSD severity are plausibly associated at baseline, and the
  direction of the resulting bias in the shared-slope analyses could
  differ.
- **Prior literature is not cited.** The underlying issue, that
  treating baseline as a response imposes a constancy constraint, is
  established in the baseline-as-covariate literature (Liu et al.
  2009, which is in `references.bib` but not cited).

**Fix.** Either derive $f$ for a simple case, for example GLS weights
under compound symmetry with one baseline row, or present A.4
explicitly as a heuristic and remove 'removes the bias'. Distinguish
the expectancy mechanism for PB. Add a scenario with
$\mathrm{Cor}(B, \mathrm{BL}) \ne 0$. Cite Liu et al. (2009) and the
constrained longitudinal data analysis literature.

### M6. The estimand is not defined, so 'power loss under carryover' mixes DGP and analysis misspecification

**Location.** Methods, 'The published data-generating process' (lines
215 to 224) and 'Simulation studies' (lines 352 to 372); Results,
'The correction in the published code' (lines 456 to 474).

**Problem.** Under graded coupling the moderation persists at a decaying
level after the drug is stopped. RH's analysis uses a binary drug
indicator, so part of the decline in power attributed to 'a fading drug
effect' (line 461 to 462; `docs/37`) comes from misspecifying the
exposure in the analysis, not from the DGP. The paper never defines
the target quantity. Is it the on-drug moderation slope
$c_{bm}\sigma_{BR}/\sigma_{bm}$, the on-minus-off contrast, or the
coefficient of an exposure-weighted indicator? The analyses in the
stress test use $D_{bc}$ (exposure-weighted). Those in Table 2 and
Table 3 do not. The simulation studies are also not reported in the
ADEMP structure (Morris et al. 2019, which the authors cite), and MCSE
is given only for Table 3.

**Fix.** Add an ADEMP-structured Methods subsection that defines the
estimand(s) and states, for each table, which analysis targets which
estimand. Report Table 2 and Table 3 also for an analysis with a
correctly specified $D_{bc}$, to separate DGP-driven from
analysis-driven loss. Give MCSE (or maximum $|z|$ against the closed
form) for every simulated quantity.

### M7. The clinical calibration contradicts the effect sizes simulated

**Location.** Discussion, lines 581 to 584; Results, lines 476 to 497;
Recommendation 4.

**Problem.** The manuscript cites an observed effect of '14 points of
symptom change per 10 mm Hg' (Raskind 2016), but never carries out the
comparison it calls for. Using the manuscript's own values
($\sigma_{BR} = 8$, $\sigma_{bm} = 15.36$), the implied slope at
$c_{bm} = 0.25$ is $0.25 \times 8 / 15.36 = 0.130$ points per mm Hg,
or 1.3 points per 10 mm Hg. Matching 1.4 points per mm Hg would need
$c_{bm} = 1.4 \times 15.36 / 8 \approx 2.7$, which is impossible under
covariance binding at this $\sigma_{BR}$ (inferred). Either the cited
effect is on a different scale (CAPS total against RH's simulated
score, or a subgroup contrast), or the DGP's $\sigma_{BR}$ is far too
small, or the simulated effects are an order of magnitude smaller
than the motivating evidence. In each case the description of 0.45 as
a 'realistic' effect size (Results heading and line 493) is
unsupported.

**Fix.** Carry out the calibration explicitly in the Methods. Verify
the 14 points per 10 mm Hg figure and its scale against the source.
Report $c_{bm}$, the implied slope per SD and per 10 mm Hg, and the
implied $\rho_{\Delta B}$ for every scenario. Justify the range of
effect sizes studied against that calibration.

### M8. Appendix B: B.4 numbers are unsupported, B.5 omits a silent repair in the maintained package, and the patch run does not do what the abstract says

**Location.** Abstract (line 55 to 56); Appendix B.1 to B.5 (lines 804
to 907).

**Problem.**

- **B.4 numbers.** 'Restoring it raised the power of the
  exposure-weighted analysis in the crossover design from 0.51 to
  0.81' (lines 887 to 889). No source gives these numbers. In
  `docs/45` (Section 8) and `dbc-preexposure/summary-reps250.csv`, the
  values are 0.484 to 0.744 (compound symmetry, CR2), 0.44 to 0.708
  (model-based) and 0.256 to 0.556 (separable AR(1)). Manuscript 02
  itself reports 0.488 against 0.830 for the binary coding. `docs/45`
  also states that manuscript 02's pipeline 'has not been rerun with
  the correction'. B.4 concerns a companion manuscript rather than RH
  and is out of scope here.
- **B.5 omission.** B.5 does not disclose that the maintained package's
  sampler falls back to `make.positive.definite()` when the Cholesky
  factorization fails, whatever the value of `makePositiveDefinite`
  (`R/generateData.R`, lines 391 to 393; inspected; documented in
  `docs/36`, Section 6.7). This contradicts Recommendation 1 ('stop on
  non-positive-definite matrices') for users of the package that the
  paper offers as the implementation.
- **Patch run.** The abstract says that coupling 'in proportion to the
  remaining drug effect' raised power in the published code. In the
  patch run only the coupling lines were replaced, so the drug-response
  mean still follows RH's compounding recursion (lines 85 to 89 of the
  source). There the remaining mean at the second off-drug visit
  decays as $2^{-3/t_{1/2}}$, while the coupling decays as
  $2^{-2/t_{1/2}}$. In that run the coupling is therefore proportional
  to the anchored decay weight, not to RH's remaining mean.

**Fix.** Remove B.4 or replace its numbers with the verified ones and
the 'not rerun' caveat. Disclose the sampler fallback in B.5 and
remove it from the package before publication. Reword the abstract to
'coupling that decays with the stated carryover half-life'.

### M9. The 'refined DGP' (aim 3, title) is not fully specified, is partly unimplemented, and the sensitivity model is placed in the Methods without being evaluated

**Location.** Background aims (lines 167 to 170); Methods, 'Refinements'
and 'Covariance-based and mean-based binding' (lines 281 to 350);
Discussion (lines 605 to 611); B.5.

**Problem.** No single table or set of equations defines the refined
DGP, with all parameters, defaults, ranges and switches. Table 1 lists
changes but not values (for example $\rho$, $c_1$, the prognostic
couplings, measurement-error SDs, dropout mechanism and decay form).
According to B.5, the separable form, anchored carryover, the stop
rule, and the natural-course and measurement-error options exist only
in stress-test drivers. The 'sensitivity model' is introduced in the
Methods, but it is 'not yet implemented' (line 610) and no result is
reported for it. It is also the standard random-coefficient model
(random sensitivity to exposure with a covariate-by-exposure
interaction), and the manuscript presents it without citation as
'proposed for this compendium'.

**Fix.** Add a specification table, or an Additional file, that a
reader could implement without the repository. Move the sensitivity
model to the Discussion as future work, with citations to the
random-coefficient literature (for example Senn 2016; Araujo et al.
2016). Either implement the refinements in the package before
submission, or describe the code honestly as driver-level.

---

## 4. Minor issues

1. **Internal labels (lines 237 to 240, Table 2, Figures 1 and 2).**
   'E9' and 'E9b' are internal estimator labels. Use descriptive names,
   for example $\hat\beta_{\text{post}}$ and $\hat\beta_{\text{base}}$.
2. **Time on drug (lines 194 and 321 to 322).** The paper describes BR
   as driven by 'cumulative time on drug', yet states that the clock
   restarts at zero on re-exposure. The 58b32a9 code accumulates `tod`
   within a contiguous on-drug run only (`buildtrialdesign.R`, lines 95
   to 97; inspected). Use 'time since the current on-drug period
   began'.
3. **Proposition 2 is per path (lines 247 to 258).** In the Hybrid
   design $n_{\text{off}}$ is 3 in paths 1 and 2 and 4 in paths 3 and 4.
   State that the slopes are path-specific and give the pooled slope as
   the allocation-weighted average (pooled E9b floor $-0.029$).
4. **Proportionality to the coupling gap (lines 259 to 263 and 448 to
   450).** 'That correlation is proportional to the coupling gap' holds
   for E9 only; for E9b the slope is proportional to
   $1 - \frac{n_{\text{off}}}{n_{\text{off}}+1}\bar g_{\text{off}}$. The
   noncentrality is proportional to
   $\rho_{\Delta B}/\sqrt{1 - \rho_{\Delta B}^2}$, not to the gap, so
   'because the noncentrality is proportional to the gap' should read
   'approximately proportional for small $\rho_{\Delta B}$'.
5. **A.1 (lines 714 to 741).** Because the weights sum to zero,
   $m_1 = 0$, so every term containing $m_1$ vanishes. Simplify the
   expression. State that for E9b the baseline row has $S_0 = 0$.
   Proposition 1 holds only for unrepaired matrices: at $c_{bm} = 0.6$
   the repair changes response-block correlations by up to 0.097, and
   the repair depends on $t_{1/2}$ through the coupling (`docs/36`,
   Section 4.3).
6. **'whether 0.01 week or 10' (line 393).** Below about 0.012 week the
   recursive off-drug mean underflows to zero in double precision, so
   the published rule switches coupling off at late off-drug visits
   (`docs/37`, Section 5.1); the figure starts at 0.03 week for this
   reason. Replace with 'for every $t_{1/2}$ above about 0.02 week' or
   add the footnote.
7. **Two published-analysis estimates (lines 395 to 397 against Table
   3).** The reconstruction gives power 0.63 falling to 0.08 to 0.11.
   The published-code run gives 0.668 and 0.126 at $t_{1/2} = 0.1$, a
   difference that `docs/37` reports as about 2.4 combined MCSE. State
   both, and the seed difference.
8. **Characterizing the repair (lines 146 to 154; Abstract lines 60 to
   62).** At $c_{bm} = 0.3$ the repair changed correlations by at most
   0.007 (`docs/36`, Section 4.3). RH's 0.3 results are therefore
   essentially valid, and only the 0.6 results describe a different DGP.
   Say so, for fairness to RH. Also state the counting unit for '20 of
   54' (one matrix per design path, $c_{bm}$ and $t_{1/2}$), since
   `docs/36` also reports 40 of 162 under a different unit.
9. **Missing qualifier on the ceiling (Abstract lines 59 to 60).**
   'Under compound symmetry the joint distribution exists only for
   correlations below 0.256' needs 'without carryover, in the Hybrid
   and crossover designs'. Under the step rule with carryover the Hybrid
   ceiling is 0.893.
10. **'Share of drug-response variation the biomarker explains' (lines
    551 to 553).** This is $c_{bm}^2$ (or $c_{bm}^2 g_t^2$), not
    $c_{bm}$.
11. **Latent-class paragraph (lines 554 to 564).** 'No two-class mixture
    reaches the covariance encoding exactly, because a mixture always
    moves the conditional mean as well' is incorrect as worded: the
    MVN encoding also moves the conditional mean, linearly. The
    mixture's conditional mean is nonlinear in $B$, and its higher
    moments differ. 'Identifying it needs samples far larger than an
    N-of-1 program has' is unsupported. The identity cited
    ($c_{bm}^2$ as a product of between-class variance fractions) is
    exact for the covariance, whereas paper 03 says 'to leading order'.
    It rests on an unpublished companion paper; give the one-line
    derivation instead.
12. **'Full variance at off-drug visits' (lines 578 to 580).** This
    weakness is shared by RH's mean-moderation variant, in which BR
    also has SD 8 and mean 0 off drug. It belongs to the component
    construct, not to covariance binding.
13. **Decayed mean moderation (lines 335 to 337 and 594 to 598).** It
    is defined only in the Discussion, and $h_t$ ('the exposure') is
    ambiguous. Define both mean-moderation variants in the Methods and
    add the decayed variant to Table 2.
14. **Prazosin pharmacokinetics (lines 311 to 316).** 'Prazosin's plasma
    half-life is a few hours' is uncited. `docs/38` gives about 2.3 h
    (Hobbs et al. 1978). Explain why 0.1 to 0.2 week (17 to 34 h, about
    7 to 15 plasma half-lives) is the 'pharmacological range' for a
    symptom response.
15. **Line numbers in B.2 and B.3.** B.2 cites 'lines 145 to 149'. The
    block is lines 145 (comment) to 150, and the quoted 'published'
    snippet omits the enclosing `if(makePositiveDefinite){` at line 146.
    That guard matters: the published pipeline passes `TRUE`
    (`generateSimulatedResults.R`, line 138; inspected). B.3 cites
    'lines 85 to 89', but the replacement code begins with the
    `brmeans <- modgompertz(...)` assignment at line 83 and replaces
    the block through line 92. Give the exact replaced range in each
    case, and note that the B.1 block sits inside the `for(c in cl)`
    loop (it runs three times, harmlessly).
16. **B.3 correctness and testing.** The anchored code is correct given
    the published semantics (`tsd` is zero at on-drug visits and
    accumulates from the last on-drug visit; it is zero before first
    exposure because of `everondrug`). Two caveats should be stated.
    First, B.3 was never run: `docs/37` says its effect was 'not tested
    separately'. Second, combined with the restart of the Gompertz
    clock, it produces a downward jump in the mean at re-exposure when
    the residual exceeds $G(t_{od})$ at the first re-exposure visit.
    The `if(nP>1)` guard is redundant.
17. **Path allocation.** RH's code assigns participants to paths in
    fixed numbers (18, 18, 17, 17) rather than at random (`docs/37`,
    Section 6). Mention this and recommend randomized allocation in the
    refined DGP.
18. **Separable structure (lines 293 and 488 to 491).** The Kronecker
    form applies to the correlation, not the covariance, because the PB
    SD scales with $e_t$. It also imposes one $\rho$ on all three
    components, whereas RH allowed separate `c.tv`, `c.pb` and `c.br`.
    State both points.
19. **Figures.** The PNGs carry embedded titles and subtitles and use
    'E9b', 'lmer' and 'c_bm = beta_bm' as plain text. BMC requires
    titles and legends in the manuscript, not in the graphic. Figure 2's
    caption should say that it shows E9 only.
20. **Table 2 caption.** 'Monte Carlo estimates agree within simulation
    error' should give the replicate count and the maximum $|z|$
    (`docs/37`: 17 of 566 cell-variants beyond 2 SE).
21. **Stale companion file.** `appendix-b.Rmd` says it 'reproduces
    Appendix B of the manuscript', but it describes correlation
    structures with references to a Section 2.2.4 and 'Architecture B'
    that no longer exist in `report.Rmd`. Remove or retitle it before
    submission.
22. **Uncited background claims.** 'Closed-form power calculations exist
    only for simplified versions of these designs' (line 94) needs
    citations, for example Wang and Schork (2019) and Senn (2002).
23. **Satterthwaite (line 224).** Cite `lmerTest` (Kuznetsova et al.
    2017). The 58b32a9 vignette loads it (inspected).
24. **Heading mismatch (line 511 to 512).** 'Two results from the
    stress test bear on the DGP directly'. The second ('Size more
    generally') concerns analyses, not the DGP.
25. **Front matter.** The 'Article type ... Target journal' line (line
    32) belongs in the cover letter, not the manuscript.

---

## 5. Numerical claims that could not be verified

| Claim (quoted) | Location | Where I looked | Finding |
|---|---|---|---|
| 'with $m = 10.99$, $d = 5$, $r = 0.42$ and SD 8'; PB '$m = 6.51$, $d = 5$, $r = 0.35$'; TV 'the same curve as PB and SD 10' | lines 194 to 200 | 58b32a9 vignette (`extracted_rp` is loaded from package data, values not in source); `docs/37` (confirms SDs 10, $10e_t$, 8, 15.36 only) | Gompertz parameters unverified; the landmarks (3.8, 4.7, 9.2 weeks) are verified analytically *given* those parameters |
| 'Restoring it raised the power ... in the crossover design from 0.51 to 0.81' | lines 887 to 889 | `docs/45` Section 8; `dbc-preexposure/summary-reps250.csv` | Not found; sources give 0.484 to 0.744, 0.44 to 0.708, 0.256 to 0.556 |
| 'The published random-intercept model rejected a true null 0.107 to 0.159' | lines 524 to 525 | `decision-screen.csv`, `summary-full.csv` | Contradicted (TV-growing cell: 0.053) |
| 'An unstructured-covariance mixed model rejected up to 0.42' | lines 525 to 526 | same | True for A (0.416); DA reached 0.509 |
| '14 points of symptom change per 10 mm Hg' | line 584 | `references.bib` only (Raskind 2016 title) | Unverified; also inconsistent with simulated scale (M7) |
| 'Prazosin's plasma half-life is a few hours' | line 311 | `docs/38` (cites Hobbs 1978, 2.3 h) | Not verified against primary source; uncited in manuscript |
| 'Identifying it needs samples far larger than an N-of-1 program has' | lines 561 to 562 | paper 03 `report.Rmd` | No quantitative support found |
| 'no two-class mixture reaches the covariance encoding exactly' | lines 562 to 564 | paper 03, 'MVN approximation' section | Not stated there in this form (paper 03 says 'to leading order') |
| 'Higher pretreatment standing systolic blood pressure was associated with greater symptom reduction on prazosin but not on placebo' | lines 88 to 90 | `references.bib` | Consistent with the cited title; full text not checked |
| '5,000 Monte Carlo replicates per cell' | line 357 | directory `mc-cells-reps5000/`, Figure 1 footnote | Inferred from file names and figure, not from a results table |
| 'To our knowledge this was the first covariance-based encoding' | lines 132 to 134 | none provided | No literature search documented (see M1) |

All other numerical claims checked were verified, among them: 0.256
(Hybrid compound-symmetry ceiling); 20 of 54 (37%); 0.54 to 0.57;
0.74 to 0.12; 0.1% residual effect; every entry of Tables 2, 3 and 4;
the gap values 0.991, 0.906 and 0.754; the decayed mean-moderation
sequence; E9b mean agreement within 0.003; the pooling term under
0.7%; $-0.126$ and $-0.026$ to $-0.032$; 0.63 and 0.08 to 0.11; the
11-cell size figures (0.069, 0.081, 0.087, 0.036 to 0.069); and the
counts of 39 cells and 26,500 trials.

---

## 6. Missing literature

- **N-of-1 simulation with carryover and serial correlation.** Percha B
  et al. (2019), 'Designing robust N-of-1 studies for precision
  medicine: simulation study and design recommendations', *J Med
  Internet Res* 21:e12641. Wang Y and Schork NJ (2019), *Healthcare*
  (in `references.bib`, not cited). Schork NJ (2022), *Harvard Data
  Science Review* (in `references.bib`, not cited).
- **Carryover models in crossover designs.** Senn (2002) and Jones and
  Kenward (2014), both in `references.bib` and neither cited. Sturdevant
  and Lumley (2021) (in `references.bib`).
- **Baseline as response or covariate.** Liu et al. (2009) (in
  `references.bib`, not cited), central to M5.
- **Nearest positive-definite repair.** Higham NJ (2002), 'Computing the
  nearest correlation matrix', *IMA J Numer Anal* 22:329-343, on what
  the silent repair does.
- **Separable (Kronecker) covariance for multivariate repeated
  measures.** For example Galecki AT (1994), *Commun Stat Theory
  Methods* 23:3105-3119.
- **Pharmacodynamic carryover.** Dayneka et al. (1993) on indirect
  response models and Sheiner et al. (1979) on effect compartments,
  both used in `docs/38`; Hobbs et al. (1978) for prazosin
  pharmacokinetics.
- **Random-coefficient heterogeneity and moderators.** Senn et al.
  (2011) (in `references.bib`); Kraemer HC et al. (2002), *Arch Gen
  Psychiatry* 59:877-883, on moderator definitions; Ballman KV (2015),
  *J Clin Oncol* 33:3968-3971, on predictive against prognostic
  biomarkers.
- **Satterthwaite inference.** Kuznetsova A et al. (2017), `lmerTest`,
  *J Stat Softw* 82(13).
- **Previously used entries.** Zucker et al. (2010) and Duan et al.
  (2013) are in `references.bib` and were evidently intended for use.

---

## 7. Suggestions for the *Trials* Methodology format

1. **Structured abstract.** The current abstract is about 315 words,
   within the 350-word limit. Keep the four headings (Background,
   Methods, Results, Conclusions). Revise the Results per M2 and M3,
   and add the stress-test design (39 cells, Hybrid design) and the
   recommended analysis. Avoid internal labels.
2. **Section order.** Background, Methods, Results, Discussion,
   Conclusions are present. Move the sensitivity model from Methods to
   Discussion (M9). Add an ADEMP-structured 'Simulation study design'
   subsection to the Methods (M6), and a Results table for the stress
   test (M4).
3. **Appendices.** BMC journals publish appendices as Additional files.
   Convert Appendix A (derivations) and Appendix B (code) to Additional
   files 1 and 2. Give each a title and description in an 'Additional
   files' or 'Supplementary information' section, and refer to them in
   the text. Keep the main text readable by trialists, with the
   equivalence of M1 and the coupling-gap figure as the central
   exposition.
4. **Title page and authorship.** BMC requires named authors with
   affiliations and a corresponding author. 'pmsimstats team' is not
   acceptable unless it is listed as a group author with its members
   named. Remove the 'Article type ... Target journal' line.
5. **Declarations.** All seven required headings are present, but four
   are placeholders. Under **Competing interests**, disclose author
   overlap with RH if it exists. The RH entry in `references.bib` lists
   'Thomas, Ronald G.' as a co-author. A manuscript that critiques the
   authors' own earlier work should say so; this affects how readers
   and editors weigh the critique, and it bears on the 'to our
   knowledge' claim. **Availability of data and materials** gives
   repository-relative paths but no URL or persistent identifier for
   pmsimstats-ng. Archive a tagged release with a DOI (for example
   Zenodo), and give the full URL of RH commit 58b32a9.
6. **References.** Switch from `statistics-in-medicine.csl` to the BMC
   (Vancouver, numbered) style. Remove or justify citations of
   unpublished companion manuscripts (`pmsimstats-paper03`); BMC
   discourages citing unpublished work in support of claims.
7. **Figures and tables.** Supply figures as separate files without
   embedded titles. Put figure titles (at most 15 words) and legends in
   the manuscript. Give each table a title above and footnotes below
   (MCSE, replicate counts).
8. **Reporting.** Add line and page numbers to the submission PDF. Use
   the Morris et al. (2019) ADEMP checklist for all three simulation
   studies, and report MCSE throughout.
9. **Length and audience.** *Trials* readers are trialists. Lead with
   the practical message (the coupling rule removes the interaction;
   check feasibility; prefer analyses that free the baseline biomarker
   slope when the biomarker may be prognostic), and put the algebra in
   Additional file 1.
