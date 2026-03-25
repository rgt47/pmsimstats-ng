# Two Architectures for Simulating Biomarker-Treatment Interactions:
# Implications for Statistical Power Under Carryover

**Ronald G. Thomas, Ph.D.**
Department of Family Medicine and Public Health, UC San Diego

**Rebecca C. Hendrickson, M.D., Ph.D.**
VA Puget Sound Health Care System, Department of Psychiatry and
Behavioral Sciences, University of Washington

---

## Abstract

Simulation studies of predictive biomarker-treatment interactions
in clinical trials require a data-generating process (DGP) that
specifies how the biomarker moderates treatment response. We
identify two fundamentally different DGP architectures used in
practice and show that the choice between them determines whether
carryover effects reduce statistical power to detect the
interaction. Under *direct mean moderation*, the biomarker scales
the treatment effect in the population mean structure; power is
largely preserved under carryover because the proportional
relationship between biomarker and drug response is maintained at
every timepoint. Under *differential correlation* (the MVN
approach), the interaction emerges from treatment-state-dependent
correlation between biomarker and response in the joint
distribution; carryover erodes this differential correlation signal
and substantially reduces power. We formalize the mathematical
distinction, illustrate the consequences with simulation results
from an N-of-1 trial design framework, discuss the biological
assumptions underlying each architecture, and provide guidance for
selecting the appropriate DGP for a given clinical context.

---

## 1. Introduction

A central question in precision medicine trial design is whether a
baseline biomarker predicts differential treatment response -- that
is, whether the treatment effect is moderated by the biomarker
value. Simulation studies are the primary tool for evaluating the
statistical power of proposed trial designs to detect such
interactions. These simulations require a data-generating process
(DGP) that encodes the biomarker-treatment interaction mechanism.

The statistical literature on treatment effect heterogeneity
(Rothwell 2005; Kent et al. 2010) distinguishes between
*qualitative* interactions (the direction of the treatment effect
reverses across biomarker strata) and *quantitative* interactions
(the magnitude varies but the direction is consistent). Simulation
frameworks must choose how to generate this heterogeneity.

The dominant approach in the clinical trial simulation literature
is direct mean moderation: the biomarker enters the outcome model
as an interaction term that explicitly scales the treatment effect.
This approach is used in essentially all simulation frameworks for
biomarker-stratified designs (Simon 2010), adaptive enrichment
trials (Freidlin & Korn 2014), basket and umbrella trials (Renfro
& Sargent 2017), and platform trials. N-of-1 simulation studies
(Zucker et al. 1997; Araujo et al. 2016; Duan et al. 2013)
similarly use hierarchical random-effects models in which the
biomarker predicts the individual treatment effect through the
mean structure. We are not aware of a clinical trial simulation
framework other than Hendrickson et al. (2020) that generates the
biomarker-treatment interaction through differential correlation in
a multivariate normal joint distribution rather than through direct
mean moderation.

We encountered both approaches in the development of the
pmsimstats simulation framework for N-of-1 clinical trials
(Hendrickson et al. 2020). The original publication code uses the
MVN differential correlation approach -- an unusual choice that, as
we show, has substantive consequences for power estimation under
carryover. During a comprehensive audit, we discovered that a
parallel tidyverse reimplementation of the same simulation had
inadvertently adopted the standard direct mean moderation approach.
Both implementations produced plausible power estimates, but they
yielded qualitatively different conclusions about how carryover
effects influence statistical power. Tracing this discrepancy to
its root cause revealed the architectural difference described in
this paper, which, to our knowledge, has not been explicitly
discussed in the simulation methodology literature.

This paper formalizes the distinction, demonstrates its
consequences, and provides guidance for investigators designing
simulation studies for predictive biomarker-moderated treatment
effects.

---

## 2. Two DGP Architectures

### 2.1 Architecture A: Direct Mean Moderation

In this approach, the biomarker enters the mean structure of the
outcome directly. The treatment effect for participant $i$ at
timepoint $t$ is:

$$Y_{it} = \mu_0 + \beta_1 \cdot D_{it} \cdot (1 + \beta_{bm}
\cdot B_i) + \epsilon_{it}$$

where $D_{it}$ is the treatment indicator (or a continuous measure
of drug exposure), $B_i$ is the participant's biomarker value
(typically centered), and $\beta_{bm}$ is the biomarker moderation
parameter. The interaction is explicit in the mean: participants
with higher biomarker values have proportionally larger treatment
effects.

Key properties:

- The interaction is **deterministic** given the biomarker value
  and treatment status.
- The biomarker moderation $\beta_{bm}$ is a **population-level
  parameter** that applies uniformly to all participants.
- The signal for the interaction is in the **first moment**
  (conditional mean) of the outcome distribution.
- Adding noise, random effects, or carryover to $D_{it}$ scales
  the entire treatment effect (including its biomarker-dependent
  component) proportionally. The **ratio** of drug effect between
  high-biomarker and low-biomarker participants is preserved.

This architecture is used in the pedagogical simulation
(`simple/simulation.R`) and in the exploratory carryover analyses
(`simulation_carryover_spectrum.R`) of the pmsimstats project.

### 2.2 Architecture B: Differential Correlation (MVN)

In this approach, the biomarker and the biological response (BR)
component are jointly drawn from a multivariate normal distribution
with treatment-state-dependent correlation:

$$\begin{pmatrix} B_i \\ BR_{it} \end{pmatrix} \sim
\text{MVN}\left(\begin{pmatrix} \mu_B \\ \mu_{BR}(t, D_{it})
\end{pmatrix}, \Sigma(D_{it})\right)$$

where the covariance matrix $\Sigma$ has:

$$\text{Cor}(B_i, BR_{it}) = \begin{cases}
c_{bm} & \text{if } D_{it} = 1 \text{ (on drug)} \\
c_{bm} \cdot e^{-\lambda \cdot t_{sd}} & \text{if } D_{it} = 0
\text{ (off drug, with decay)} \\
0 & \text{if } D_{it} = 0 \text{ (off drug, no decay)}
\end{cases}$$

The interaction is not in the mean structure (the population mean
of $BR_{it}$ is the same for all participants with the same
treatment history). Instead, it emerges from the **conditional
distribution**: given a participant's biomarker value $B_i = b$,
the conditional expectation of their biological response is:

$$E[BR_{it} | B_i = b, D_{it} = 1] = \mu_{BR}(t) + c_{bm}
\cdot \frac{\sigma_{BR}}{\sigma_B} \cdot (b - \mu_B)$$

This conditional expectation has an interaction-like structure:
the biomarker value $b$ modulates the expected BR, but only when
the correlation $c_{bm}$ is nonzero (i.e., when on drug).

Key properties:

- The interaction is **probabilistic**: it manifests as a
  tendency, not a deterministic scaling.
- The biomarker moderation $c_{bm}$ is a **correlation
  parameter** that governs the strength of the association.
- The signal for the interaction is in the **second moment**
  (covariance structure) of the joint distribution.
- Carryover effects in the DGP inflate the BR means during
  off-drug periods, making them resemble on-drug periods. This
  reduces the **differential** in the correlation structure
  between treatment states, weakening the signal that the
  analysis model detects.

This architecture is used in the Hendrickson et al. (2020)
publication code and the audited `R/` package in pmsimstats-ng.

### 2.3 Mathematical Comparison

Under Architecture A, the expected outcome difference between
two participants with biomarker values $b_1$ and $b_2$, both on
drug, is:

$$E[Y | b_1, D=1] - E[Y | b_2, D=1] = \beta_1 \cdot \beta_{bm}
\cdot (b_1 - b_2) \cdot D$$

This difference is **constant** regardless of carryover: if
$D$ is replaced by $D \cdot \text{carryover}$, the difference
scales by the carryover factor, but the ratio is preserved.

Under Architecture B, the analogous quantity is:

$$E[BR | b_1, D=1] - E[BR | b_2, D=1] = c_{bm} \cdot
\frac{\sigma_{BR}}{\sigma_B} \cdot (b_1 - b_2)$$

During off-drug periods with carryover:

$$E[BR | b, D=0, t_{sd}] = \mu_{BR,\text{off}}(t_{sd}) +
c_{bm} \cdot e^{-\lambda t_{sd}} \cdot
\frac{\sigma_{BR}}{\sigma_B} \cdot (b - \mu_B)$$

The biomarker-dependent component decays exponentially with
time off drug. As $t_{sd}$ increases, the off-drug conditional
expectation converges to $\mu_{BR,\text{off}}$ for all
biomarker values, eliminating the differential signal.

---

## 3. Consequences for Statistical Power Under Carryover

### 3.1 Empirical demonstration

We ran identical simulation configurations under both
architectures using the pmsimstats-ng framework. The DGP
parameters matched as closely as possible between the two
approaches. The analysis model was the same in both cases:
`nlme::lme` with `corCAR1` and the `bm:Dbc` interaction term.

**Setup:**

- Trial designs: OL, CO, N-of-1 (Hybrid), OL+BDC
- Sample size: $N = 70$
- Biomarker correlation / moderation: $c_{bm} = 0.45$ /
  $\beta_{bm} = 0.45$
- DGP carryover half-life: $t_{1/2} \in \{0, 0.5, 1.0\}$ weeks
- Analysis: matched Dbc with same half-life
- Replicates: 50 per cell

**Results under Architecture A (direct mean moderation):**

| Design | $t_{1/2} = 0$ | $t_{1/2} = 0.5$ | $t_{1/2} = 1.0$ |
|--------|:------------:|:---------------:|:---------------:|
| CO | 0.38 | 0.36 | 0.34 |
| N-of-1 | 0.74 | 0.72 | 0.68 |
| OL+BDC | 0.62 | 0.58 | 0.54 |

Power declines modestly (~10-15% relative) with increasing
carryover.

**Results under Architecture B (differential correlation, MVN):**

| Design | $t_{1/2} = 0$ | $t_{1/2} = 0.5$ | $t_{1/2} = 1.0$ |
|--------|:------------:|:---------------:|:---------------:|
| CO | 0.40 | 0.40 | 0.32 |
| N-of-1 | 0.82 | 0.64 | 0.50 |
| OL+BDC | 0.62 | 0.38 | 0.24 |

Power declines substantially (~40-60% relative) with increasing
carryover.

### 3.2 Why the architectures diverge

**Under Architecture A**, carryover reduces the magnitude of the
drug exposure variable $D_{it}$ during off-drug periods, but the
biomarker moderation coefficient $\beta_{bm}$ is a fixed
population parameter. The interaction `bm * D` has reduced
variance (because $D$ has less contrast), but the signal-to-noise
ratio degrades only modestly because the noise structure is
independent of the carryover.

**Under Architecture B**, carryover operates on two levels
simultaneously:

1. **Mean blurring**: The BR means during off-drug periods are
   inflated by residual carryover, making them resemble on-drug
   BR means. This reduces the treatment contrast that the analysis
   model uses to distinguish drug effect from no-drug.

2. **Correlation erosion**: The BM-BR correlation during off-drug
   periods decays exponentially with $t_{sd}$. This is not just a
   statistical modeling choice -- it is a structural consequence
   of the DGP. The correlation between biomarker and BR reflects
   the degree to which the drug effect is active; as the drug
   washes out, the correlation weakens because the drug-mediated
   component of BR variance shrinks relative to the
   drug-independent component.

The analysis model's `bm:Dbc` interaction coefficient depends on
the covariance between `bm` and `Dbc`-weighted outcomes. Under
Architecture B, this covariance is being attacked from both sides:
the treatment contrast in $Dbc$ is compressed, and the
biomarker-BR correlation that generates the covariance signal is
simultaneously decaying.

---

## 4. Biological Assumptions

The choice between architectures is not merely a statistical
convenience -- it reflects different assumptions about the
biological mechanism by which the biomarker moderates treatment
response.

### 4.1 When Architecture A is appropriate

Architecture A (direct mean moderation) is appropriate when the
biomarker governs the **magnitude of the drug's biological
effect** for each individual participant. The drug effect for
participant $i$ is deterministically scaled by their biomarker
value: $\text{effect}_i = f(B_i) \cdot \text{drug}$.

**Clinical examples:**

- **Pharmacokinetic moderation (drug metabolism).** Blood pressure
  may influence the effective plasma concentration of prazosin
  through hemodynamic effects on hepatic clearance. Participants
  with higher resting blood pressure may achieve higher effective
  drug levels, producing a proportionally larger PTSD symptom
  reduction. The biomarker-drug relationship is mechanistic and
  deterministic: given the biomarker value, the drug effect is
  fully determined (up to random noise).

- **Receptor density moderation.** A biomarker that measures
  target receptor density (e.g., PET imaging of alpha-1
  adrenergic receptor availability) directly determines the
  pharmacodynamic response to an antagonist. More receptors
  means more drug effect, proportionally.

- **Genetic metabolizer status.** CYP2D6 metabolizer phenotype
  determines the rate of drug activation or clearance. Poor
  metabolizers of a prodrug receive less active compound,
  reducing their treatment effect proportionally.

- **Dose-response with biomarker-dependent effective dose.**
  When the biomarker determines the effective dose (through
  body composition, renal clearance, protein binding, etc.),
  the treatment effect scales with the biomarker through the
  dose-response curve.

In all these cases, the biomarker's moderating effect is
preserved under carryover because the residual drug effect during
off-drug periods maintains the same proportional relationship
with the biomarker. If participant A has twice the drug effect of
participant B when on drug, they will also have approximately
twice the residual effect during the washout period.

### 4.2 When Architecture B is appropriate

Architecture B (differential correlation) is appropriate when
the biomarker predicts **which participants are likely to
respond**, but the magnitude of any individual's response is not
deterministically governed by the biomarker. The biomarker and
the drug response are associated at the population level, but
the association is mediated by latent factors (disease subtype,
neurobiological heterogeneity, comorbidities) that the biomarker
imperfectly indexes.

**Clinical examples:**

- **Blood pressure as a PTSD subtype marker.** Elevated resting
  blood pressure may index a noradrenergic-predominant PTSD
  phenotype (Hendrickson et al. 2020) that is more likely to
  respond to prazosin (an alpha-1 adrenergic antagonist). The
  blood pressure does not *cause* the drug to work better --
  it correlates with an underlying neurobiological state that
  determines drug responsiveness. Two participants with
  identical blood pressure may have very different responses
  because blood pressure is an imperfect proxy for the
  noradrenergic state.

- **Baseline disease severity as a treatment predictor.** Trials
  in Alzheimer's disease, depression, and chronic pain often
  find that patients with more severe baseline symptoms show
  larger treatment effects (the 'floor effect' or 'regression
  to the mean' phenomenon). The baseline severity score
  correlates with treatment response, but the relationship is
  statistical (population-level tendency), not deterministic.

- **Inflammatory biomarkers in psychiatry.** C-reactive protein
  (CRP) or interleukin-6 levels may predict response to
  anti-inflammatory augmentation of antidepressants (Raison
  et al. 2013). High-CRP patients are more likely to respond,
  but the relationship is probabilistic -- many high-CRP
  patients do not respond, and some low-CRP patients do.

- **Genetic risk scores.** Polygenic risk scores predict disease
  trajectory and treatment response, but explain only a fraction
  of the variance. The remainder reflects unmeasured genetic,
  epigenetic, and environmental factors.

In these cases, carryover has a different implication: as the drug
washes out, the correlation between biomarker and response weakens
because the drug-mediated component of response variance
diminishes. The biomarker's predictive power is inherently tied to
the presence of active drug effect. Without drug, the biomarker is
just a biomarker -- it does not predict the non-drug component of
symptom change.

### 4.3 The prazosin-PTSD case

The motivating application for pmsimstats is the prazosin trial
by Raskind et al. (2013), analyzed by Murray et al. (in
preparation). The biomarker is resting systolic blood pressure
(SBP) and the outcome is CAPS (Clinician-Administered PTSD
Scale) score change.

The biological hypothesis is that elevated SBP indexes heightened
central noradrenergic tone, which characterizes a subtype of PTSD
that is preferentially responsive to alpha-1 adrenergic blockade
by prazosin. This is an Architecture B scenario: SBP is a proxy
for an underlying neurobiological state, not a direct determinant
of drug pharmacokinetics. Two patients with the same SBP may
differ in their actual noradrenergic tone, alpha-1 receptor
distribution, and downstream signaling, leading to different
treatment responses despite identical biomarker values.

The Hendrickson et al. (2020) DGP correctly uses Architecture B
(MVN differential correlation) for this clinical context. The
finding that carryover substantially reduces power under this
architecture is therefore clinically relevant: it suggests that
trial designs with longer off-drug periods will have meaningfully
less power to detect the SBP-prazosin interaction than designs
with shorter off-drug periods, and this effect is larger than
what a direct mean moderation simulation would predict.

---

## 5. Implications for Simulation Study Design

### 5.1 Choosing the appropriate architecture

The choice of DGP architecture should be guided by the biological
mechanism hypothesized for the biomarker-treatment interaction:

| Criterion | Architecture A | Architecture B |
|-----------|:-----------:|:-----------:|
| Biomarker role | Causal mediator | Statistical predictor |
| Mechanism | Deterministic scaling | Probabilistic association |
| Signal location | Mean structure | Covariance structure |
| Carryover impact on power | Modest (~10-15%) | Substantial (~40-60%) |
| Appropriate when | Biomarker determines effective dose or PK | Biomarker indexes a latent subtype |

### 5.2 Reporting requirements

Simulation studies evaluating biomarker-moderated treatment effects
should explicitly state:

1. Whether the biomarker-treatment interaction is generated through
   mean moderation or differential correlation.
2. The biological rationale for the chosen mechanism.
3. Whether the carryover model (if present) operates on the
   interaction signal consistently with the chosen architecture.
4. Sensitivity analyses under the alternative architecture, if the
   biological mechanism is uncertain.

### 5.3 When the choice is uncertain

If the biological mechanism is ambiguous (as it often is in early
biomarker development), investigators should:

1. Run the simulation under both architectures and report the
   range of power estimates.
2. Note that Architecture B produces more conservative (lower)
   power estimates under carryover, making it the safer default
   for trial design decisions.
3. Design the trial to minimize carryover (adequate washout
   periods) to reduce sensitivity to the DGP architecture choice.

---

## 6. Relationship to the Broader Literature

### 6.1 Treatment effect heterogeneity

The statistical literature on treatment effect heterogeneity
(Rothwell 2005; Kent et al. 2010; Varadhan et al. 2013)
distinguishes between *predictive* biomarkers (those that modify
the treatment effect) and *prognostic* biomarkers (those that
predict outcome regardless of treatment). Both architectures
generate predictive biomarker effects, but through different
statistical mechanisms. Architecture A corresponds to a
*parametric interaction model* where the effect modifier enters
the regression directly. Architecture B corresponds to a
*latent class* or *mixture model* perspective where the biomarker
is a noisy indicator of class membership.

### 6.2 Prevalence of each architecture

A survey of the simulation methodology literature reveals that
Architecture A (direct mean moderation) is the near-universal
standard. All simulation frameworks for biomarker-stratified
designs (Simon 2010), adaptive enrichment trials (Freidlin & Korn
2014), basket and umbrella trials (Renfro & Sargent 2017), and
platform trials use direct mean moderation with explicit
interaction terms. N-of-1 simulation studies (Zucker et al. 1997;
Araujo et al. 2016; Duan et al. 2013) use hierarchical
random-effects models in which individual treatment effects are
drawn from a population distribution -- a variant of Architecture
A where the biomarker predicts the random treatment effect
through the mean structure.

Architecture B (differential correlation via MVN) appears
primarily in latent variable and structural equation modeling
contexts (Rizopoulos 2012), Bayesian joint modeling frameworks,
and copula-based simulation studies. In the clinical trial
simulation literature specifically, the Hendrickson et al. (2020)
framework appears to be unique in using differential correlation
as the mechanism for biomarker-treatment interaction in an N-of-1
trial context.

This predominance of Architecture A has an important implication:
the power estimates reported in the vast majority of
biomarker-moderated trial simulation studies may be optimistic
for biomarkers that operate through an Architecture B mechanism,
because these studies do not account for the additional power
loss that carryover produces when the interaction signal resides
in the covariance structure rather than the mean structure.

### 6.3 Crossover and N-of-1 trial methodology

The crossover trial literature (Senn 2002; Jones & Kenward 2014)
addresses carryover extensively but does not distinguish between
architectures for biomarker-treatment interaction modeling. The
standard recommendation -- always include carryover in the
analysis model or use adequate washout -- applies to both
architectures but has different power implications under each.
Our finding that Architecture B is more sensitive to carryover
suggests that crossover designs with short washout periods may
be less suitable for detecting correlation-based biomarker
interactions than previously recognized.

### 6.4 Precision medicine trial design

Enrichment and biomarker-stratified designs (Simon 2010; Freidlin
& Korn 2014) use Architecture A exclusively for power
calculations. The power estimates from these calculations may be
optimistic for biomarkers that operate through an Architecture B
mechanism, particularly in designs with crossover components or
treatment switching.

---

## 7. Conclusions

1. The choice between direct mean moderation (Architecture A) and
   differential correlation (Architecture B) for generating
   biomarker-treatment interactions in clinical trial simulations
   has substantive consequences for power estimates under
   carryover.

2. Under Architecture A, carryover reduces power modestly
   (~10-15%) because the proportional biomarker-drug relationship
   is preserved during washout. Under Architecture B, carryover
   reduces power substantially (~40-60%) because the differential
   correlation signal erodes as the drug effect wanes.

3. The choice should be guided by the hypothesized biological
   mechanism: Architecture A for biomarkers that directly
   determine drug pharmacokinetics or pharmacodynamics;
   Architecture B for biomarkers that statistically predict
   treatment response through latent subtype membership.

4. For the prazosin-PTSD application with blood pressure as the
   predictive biomarker, Architecture B (MVN differential
   correlation) is the more appropriate model, and the resulting
   power sensitivity to carryover should inform trial design
   decisions regarding off-drug period duration.

5. Simulation studies for biomarker-moderated treatment effects
   should explicitly report their DGP architecture and consider
   sensitivity analyses under the alternative when the biological
   mechanism is uncertain.

---

## References

Araujo A, Julious S, Senn S. Understanding variation in sets of
N-of-1 trials. *PLoS ONE*. 2016;11(12):e0167167.

Duan N, Kravitz RL, Schmid CH. Single-patient (N-of-1) trials: a
pragmatic clinical decision methodology. *J Clin Epidemiol*.
2013;66(8):S21-S28.

Freidlin B, Korn EL. Biomarker enrichment strategies: matching
trial design to biomarker credentials. *Nat Rev Clin Oncol*.
2014;11(2):81-90.

Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
aggregated N-of-1 trial designs for predictive biomarker
validation: statistical methods and practical considerations.
*Front Digit Health*. 2020;2:13.

Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*.
3rd ed. Boca Raton, FL: CRC Press; 2014.

Kent DM, Rothwell PM, Ioannidis JPA, Altman DG, Hayward RA.
Assessing and reporting heterogeneity in treatment effects in
clinical trials: a proposal. *Trials*. 2010;11:85.

Raison CL, Rutherford RE, Woolwine BJ, et al. A randomized
controlled trial of the tumor necrosis factor antagonist
infliximab for treatment-resistant depression: the role of
baseline inflammatory biomarkers. *JAMA Psychiatry*.
2013;70(1):31-41.

Raskind MA, Peterson K, Williams T, et al. A trial of prazosin
for combat trauma PTSD with nightmares in active-duty soldiers
returned from Iraq and Afghanistan. *Am J Psychiatry*.
2013;170(9):1003-1010.

Rothwell PM. Treating individuals 5. Subgroup analysis in
randomised controlled trials: importance, indications, and
interpretation. *Lancet*. 2005;365(9454):176-186.

Senn S. *Cross-over Trials in Clinical Research*. 2nd ed.
Chichester, UK: John Wiley & Sons; 2002.

Simon R. Clinical trial designs for evaluating the medical
utility of prognostic and predictive biomarkers in oncology.
*Per Med*. 2010;7(1):33-47.

Varadhan R, Segal JB, Boyd CM, Wu AW, Weiss CO. A framework for
the analysis of heterogeneity of treatment effect in
patient-centered outcomes research. *J Clin Epidemiol*.
2013;66(8):818-825.

Relling MV, Evans WE. Pharmacogenomics in the clinic. *Nature*.
2015;526(7573):343-350.

Renfro LA, Sargent DJ. Statistical controversies in clinical
research: basket trials, umbrella trials, and other master
protocols. *Ann Oncol*. 2017;28(1):34-43.

Rizopoulos D. *Joint Models for Longitudinal and Time-to-Event
Data: With Applications in R*. Boca Raton, FL: CRC Press; 2012.

Zucker DR, Schmid CH, McIntosh MW, D'Agostino RB, Selker HP,
Lau J. Combining single patient (N-of-1) trials to estimate
population treatment effects and to evaluate individual patient
responses to treatment. *J Clin Epidemiol*. 1997;50(4):401-410.

---
*Rendered on 2026-03-25 at 09:27 PDT.*
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/dgp-architecture-white-paper.md*
