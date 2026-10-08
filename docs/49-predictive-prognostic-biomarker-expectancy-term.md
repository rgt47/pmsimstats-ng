# Predictive or Also Prognostic? Whether to Include the Expectancy Term in the Hybrid-Design Interaction Model {.unlisted .unnumbered}
*2026-10-07 10:51 PDT*

**Author.** pmsimstats team

**Purpose.** A candidate biomarker such as standing blood pressure may
be *predictive*: it changes the effect of prazosin. It may also be
*prognostic*: it is associated with the symptom score whatever the
treatment, through the natural course, the placebo response or the
baseline severity.

This paper asks one practical question for the analysis of a Hybrid
N-of-1 trial: should the interaction model include the expectancy
term and its interaction with the biomarker, `De + bmc:De`?

It sets out where in the symptom model a biomarker can act, derives
what each candidate analysis assumes about it, gives the simulation
evidence from the stress test (`docs/46`, `docs/48`), and ends with a
decision framework and a recommendation. Results still running are
marked as pending.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`, `docs/45` and
`docs/48`. Symbols marked *local* are defined here.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $Y_{it}$ | symptom score (code `Sx`); $t = 0$ is baseline | canonical |
| $\mathrm{BL}_i$, $BR_{it}$, $PB_{it}$, $TV_{it}$ | baseline level; drug, expectancy and natural-course responses | canonical |
| $B_i$, $b_i$ | biomarker (code `bm`, centered `bmc`); standardized biomarker | canonical |
| $D_{bc,it}$ | exposure-decayed drug indicator (code `Dbc`) | canonical |
| $e_t$ | expectancy: 0 at baseline, 1 in the open-label phase, 0.5 when blinded (code `De`) | canonical (`docs/40`) |
| $\beta_B$, $\beta_{bm:D}$ | biomarker main effect; biomarker-by-drug interaction | local (`docs/41`); canonical |
| $\beta_{bm:e}$ | coefficient of `bmc:De`, the biomarker-by-expectancy term | local |
| $s_0$, $s_t$ | true slope of $Y$ on $B$ at baseline, and at post-baseline visit $t$ | local (`docs/48`) |
| $c_{bm}$ | predictive coupling: biomarker-BR correlation on drug | canonical |
| $c_{bm,TV}$, $c_{bm,PB}$, $c_{bm,BL}$ | prognostic couplings: biomarker correlation with TV, PB, BL | local |
| $\beta_{bm}^{TV}$ | biomarker moderation of the natural-course mean | canonical (component slope) |
| $\pi$ | power | canonical |

Table: Notation index and glossary of symbols used in this document

**Analyses** (all with a random intercept, `corCAR1` residuals, CR2
standard errors, linear time, baseline as a response row, and
$D_{bc,it}$ at an assumed half-life of 0.5 week):

- **S**: `Sx ~ bmc + t + Dbc + bmc:Dbc`.
- **ES**: S plus `bmc:bl0`, a separate biomarker slope at baseline.
- **F1**: S plus `De + bmc:De`, the expectancy and its interaction with
  the biomarker.
- **F2**: S plus `bmc:t`, a biomarker slope that changes linearly in
  time.

**Glossary.**

- **Predictive biomarker.** Associated with the size of the treatment
  effect [1].
- **Prognostic biomarker.** Associated with the outcome irrespective of
  treatment.
- **Constancy assumption.** The biomarker's association with the outcome
  is the same at baseline as at post-baseline off-drug visits
  ($s_0 = s_{\text{off}}$).
- **Size.** False-positive rate of the interaction test, nominal 0.05.

## 1. Summary

1. **Where the biomarker acts decides whether the interaction test is
   valid.** In the component model a biomarker can act on four
   components:
   - the drug response (predictive);
   - the natural course or the placebo response (prognostic,
     post-baseline only);
   - the baseline level (prognostic, at every visit).
2. **The analysis without the expectancy term (S) assumes the
   biomarker is purely predictive,** or prognostic only through the
   baseline level. It shares one biomarker slope between baseline and
   the post-baseline off-drug visits. When the biomarker is prognostic
   for the natural course or the placebo response, that slope is wrong
   at one set of visits, and the error is absorbed by the interaction.
3. **Adding `De + bmc:De` (F1) removes the bias in every prognostic
   scenario simulated.** The expectancy is 0 at baseline, so the term
   frees the baseline slope. It also lets the slope scale with
   expectancy, which is exactly the form a placebo-linked biomarker
   takes. In both scenarios the interaction is then identified from the
   blinded visits, where on-drug and off-drug visits share the same
   expectancy.
4. **Evidence** (verified, stress test, 1,000 null replicates per cell):

   | Biomarker linked to | S | F1 | ES |
   |---|---|---|---|
   | natural course, constant | 0.099 | 0.047 | 0.049 |
   | natural course, growing | 0.049 | 0.046 | 0.063 |
   | placebo response | 0.100 | 0.040 | 0.049 |
   | nothing (11 standard cells, pooled) | 0.048 | 0.051 | 0.048 |

   Table: Null rejection rates of S, F1 and ES by biomarker linkage, stress test

5. **The cost is power, and it depends on the correlation structure.**
   Under AR(1) data F1 gives up about 30% of S's power (0.364 against
   0.514 at $c_{bm} = 0.30$, size-adjusted). Under compound symmetry the
   loss is small (0.522 against 0.554 at $c_{bm} = 0.25$).
6. **The direct evidence on blood pressure is weak.** In the placebo arm
   of the trial that motivated blood pressure as a biomarker, baseline
   pressure showed no association with symptom change [2].
   That arm had 35 participants, so the evidence cannot rule out a
   prognostic association of the size that biases S.
7. **Recommendation.** Include the expectancy term in the primary
   analysis (F1), with S and ES as pre-specified sensitivity analyses.
   A material disagreement between F1 and S is itself evidence that the
   biomarker is prognostic. Reverse the roles only if independent data
   show that the biomarker is not prognostic.

## 2. Where a biomarker can act

The symptom score is $Y_{it} = \mathrm{BL}_i - (TV_{it} + PB_{it} +
BR_{it})$. Each component can carry an association with the biomarker.
The table gives the implied slope of $Y$ on $b_i$ at baseline and at
post-baseline visits.

| Pathway | Kind | $s_0$ (baseline) | $s_t$, post-baseline | Shape over visits |
|---|---|---|---|---|
| Drug response, $c_{bm}$ | predictive | 0 | $-c_{bm}\sigma_{BR}\,g_t$ | full on drug, decaying off drug |
| Natural course, $c_{bm,TV}$ | prognostic | 0 | $-c_{bm,TV}\sigma_{TV}$ | constant after baseline |
| Natural course, $\beta_{bm}^{TV}$ | prognostic | 0 | $-\beta_{bm}^{TV}\sigma_{TV}\,G_{TV}(t)/\max G_{TV}$ | grows with the natural course |
| Placebo response, $c_{bm,PB}$ | prognostic | 0 | $-c_{bm,PB}\cdot 10\,e_t$ | scales with expectancy |
| Baseline severity, $c_{bm,BL}$ | prognostic | $c_{bm,BL}\sigma_{BL}$ | $c_{bm,BL}\sigma_{BL}$ | constant at every visit, baseline included |

Table: Biomarker pathways by kind with implied baseline and post-baseline slopes

*Slopes per SD of biomarker, derived from the construct: $\sigma_{TV} = 10$,
the PB SD is $10e_t$, and $g_t$ is the graded coupling.*

Only the first row is the estimand. In every other row the
biomarker's association does not depend on drug, so a correct
analysis must report no interaction. Two features of the Hybrid design
make some of these pathways dangerous.

- **Baseline precedes every response component.** A prognostic
  association through the natural course or the placebo response is
  zero at baseline and nonzero afterwards. A baseline-severity
  association is the same at every visit.
- **Expectancy and drug are confounded in the open-label phase.** Every
  open-label visit is on drug and has $e_t = 1$, while every off-drug
  visit is blinded with $e_t = 0.5$. A slope that scales with
  expectancy is therefore larger on drug than off drug on average, even
  with no drug interaction.

## 3. What each analysis assumes

Write the biomarker slope the analysis fits at each kind of visit. With
$D_{bc} = 1$ on drug and close to 0 off drug:

| Visit | S | ES | F1 |
|---|---|---|---|
| baseline | $\beta_B$ | $\beta_B + \beta_0$ | $\beta_B$ |
| open label, on drug | $\beta_B + \beta_{bm:D}$ | $\beta_B + \beta_{bm:D}$ | $\beta_B + \beta_{bm:e} + \beta_{bm:D}$ |
| blinded, on drug | $\beta_B + \beta_{bm:D}$ | $\beta_B + \beta_{bm:D}$ | $\beta_B + 0.5\beta_{bm:e} + \beta_{bm:D}$ |
| blinded, off drug | $\beta_B$ | $\beta_B$ | $\beta_B + 0.5\beta_{bm:e}$ |

Table: Biomarker slope fitted by S, ES and F1 at each kind of visit

Here $\beta_0$ is the coefficient of `bmc:bl0`.

**S** ties baseline to the off-drug visits. Under any pathway with
$s_0 \ne s_{\text{off}}$, the fitted $\beta_B$ falls between the two
true slopes. The on-drug rows then differ from it by part of the
prognostic slope, which is attributed to $\beta_{bm:D}$
(`docs/48`, Section 2):
$\hat\beta_{bm:D} \approx (1 - f)\,s_{\text{prog}}$, with $f$ set by the
weight of the baseline row.

**ES** frees baseline but assumes one post-baseline off-drug slope. It
is exact when the prognostic slope is constant after baseline.

**F1** frees baseline through $e_0 = 0$, and lets the slope change with
expectancy. Under an expectancy-scaled pathway ($s_t \propto e_t$), the
`bmc:De` term absorbs it exactly. Under a constant post-baseline
pathway, F1 misfits the open-label and blinded slopes. But the misfit
is the same at blinded on-drug and blinded off-drug visits, so it
cancels in the blinded contrast that then identifies $\beta_{bm:D}$
(derived from the table). In both cases the interaction rests on the
comparison the Hybrid design was built to make clean: blinded on drug
against blinded off drug, at the same expectancy.

**Baseline severity** gives $s_0 = s_t$ at every visit, so the constancy
assumption holds and S should be unbiased. This is a prediction; the
simulation is running (Section 5).

## 4. What F1 gives up

F1 estimates two biomarker slopes that S constrains. Its interaction
then draws on less of the data:

- the open-label visits mostly inform $\beta_{bm:e}$ rather than
  $\beta_{bm:D}$;
- the baseline row no longer informs the off-drug slope.

Under AR(1) data the baseline row is a valuable observation of the
off-drug slope (`docs/45`, Section 6). Losing it costs about 30% of the
power. Under compound symmetry the participant effect already supplies
that information, and the loss is small.

Size-adjusted power (verified, stress test, 500 replicates per cell):

| Cell | S | F1 | ES |
|---|---|---|---|
| AR(1), $c_{bm} = 0.30$ | 0.514 | 0.364 | 0.364 |
| AR(1), $c_{bm} = 0.45$ | 0.808 | 0.672 | 0.656 |
| Compound symmetry, $c_{bm} = 0.25$ | 0.554 | 0.522 | 0.662 |
| Natural course prognostic, constant, $c_{bm} = 0.30$ | (biased) | 0.482 | 0.370 |
| Natural course prognostic, growing, $c_{bm} = 0.30$ | (nominal) 0.516 | 0.378 | 0.192 |
| Placebo prognostic, $c_{bm} = 0.30$ | (biased) | 0.446 | 0.498 |

Table: Size-adjusted power of S, F1 and ES by stress-test cell

Two observations:

- **ES and F1 are close when the biomarker is purely predictive.** F1
  keeps its power under every prognostic form, whereas ES loses most of
  it when the natural-course association grows over time.
- **S's power in the prognostic cells is not a fair comparison.** Its
  test is biased there; "size-adjusted" corrects the level, not the
  bias.

## 5. Pending results

Three post hoc runs on the same simulated trials are in progress
(`11-strawman-stress-test.R`):

1. **Baseline severity** ($c_{bm,BL} = 0.3$, null and alternative)
   under every analysis. Prediction: all analyses nominal, because the
   constancy assumption holds.
2. **F1 with an AIC-selected half-life** (F1A). It tests whether the
   robust model keeps its size, and how much power data-driven
   selection costs it (about 0.03 for S).
3. **F1's terms with phase-specific residual variances** (G1). It tests
   whether part of the power F1 gives up can be recovered while keeping
   its robustness.

These results will be added here when they are in. They are post hoc:
they can motivate a change in recommendation but not settle it.

## 6. Decision framework

**Question 1: is there independent evidence that the biomarker is not
prognostic?** For blood pressure, the only direct evidence is a null
placebo-arm association in 35 participants [2]. That
cannot exclude a prognostic slope large enough to bias S. A prognostic
coupling of 0.3 roughly doubles S's false-positive rate. Without
stronger evidence, treat the biomarker as possibly prognostic.

**Question 2: which pathway is plausible?**

- **Natural course or placebo response.** A noradrenergic marker could
  plausibly track spontaneous improvement or placebo response.
  F1 is robust to all three forms simulated.
- **Baseline severity only.** S is expected to be valid, but this does
  not remove the risk from the other pathways.

**Question 3: what is the cost of protection?** It depends on the
serial correlation the trial is likely to show. If pilot data suggest
little serial correlation beyond a participant effect, F1 costs little.
Under strong serial correlation it costs about 30% of the power. The
sample size should then be planned for F1, not for S.

**Question 4: how should disagreement be read?** Pre-specify S as a
sensitivity analysis. If F1 and S agree, the conclusion does not
depend on the constancy assumption. If they disagree materially, the
disagreement is itself evidence of a prognostic association. Report
F1, and report the estimate of `bmc:De` as an exploratory measure of
it.

## 7. Recommendation

1. **Primary analysis: F1.** Random intercept with `corCAR1` residuals,
   CR2 standard errors, linear time, baseline as a response row,
   $D_{bc}$ at a pre-specified half-life of 0.5 week, and
   `De + bmc:De`.
2. **Pre-specified sensitivity analyses:**
   - S (constancy assumed);
   - ES (baseline freed, expectancy not modeled);
   - binary $D$ in place of $D_{bc}$.
3. **Plan the sample size for F1,** with the serial-correlation
   assumption stated.
4. **Revisit when the pending runs are in.** F1 with an AIC-selected
   half-life may replace the fixed half-life on grounds of face
   validity, at a power cost to be measured. The phase-variance version
   of F1 may replace F1 if it holds size everywhere with materially more
   power. Either change would need a confirmatory run.

## 8. Limitations

- **One strength per prognostic form.** All prognostic forms were
  simulated at 0.3; how bias grows with the strength is not mapped.
- **Single pathways only.** Combinations, for example natural course
  and placebo together, were not simulated.
- **Linear, Gaussian couplings.** A biomarker that acts on the rate
  rather than the size of response is outside the construct.
- **One design.** The expectancy argument is specific to the Hybrid
  design's confounding of open-label and on-drug visits. In a fully
  blinded crossover the expectancy term would be constant and
  unnecessary.
- **The mechanism is heuristic.** The account of Section 3 explains
  the observed pattern but is not a derivation of the bias for each
  analysis.

## 9. Reproducibility

```bash
Rscript analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R \
  --cores 8 --chunk 25
Rscript analysis/scripts/quick-sim/carryover-closed-form/12-strawman-stress-decision.R \
  --strawman DS --no-symmetry
```

Outputs are in `analysis/data/quick-sim/carryover-closed-form/strawman-stress/`.

## 10. References

1. Kent DM, Paulus JK, van Klaveren D, et al. The Predictive Approaches
   to Treatment effect Heterogeneity (PATH) statement. *Annals of
   Internal Medicine* 2020; 172(1):35-45.
2. Raskind MA, Millard SP, Petrie EC, et al. Higher pretreatment blood
   pressure is associated with greater posttraumatic stress disorder
   symptom reduction in soldiers treated with prazosin. *Biological
   Psychiatry* 2016; 80(10):736-742. doi:10.1016/j.biopsych.2016.03.2108.
3. Liu GF, Lu K, Mogg R, Mallick M, Mehrotra DV. Should baseline be a
   covariate or dependent variable in analyses of change from baseline
   in clinical trials? *Statistics in Medicine* 2009; 28(20):2509-2530.
4. pmsimstats team. `docs/40`, `docs/41`, `docs/45`, `docs/46`,
   `docs/48`.
