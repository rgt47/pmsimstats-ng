# Referee Report

*2026-06-12 09:57 PDT*

**Manuscript.** "Three-component decomposition of treatment response in
aggregated N-of-1 trials: pharmacological, expectation-driven, and
natural-history components."

**Target tier.** Statistics in Medicine / JRSS-C / Biometrical Journal
(applied biostatistical methodology).

**Recommendation.** Major revision.

---

## 1. Summary of the manuscript

The manuscript argues that treatment response in aggregated N-of-1
clinical trials should be modelled as the additive sum of three
causally distinct components: a pharmacological biological response
(BR), a placebo-belief response (PB), and a time-variant
natural-history component (TV). Each component is specified as a
modified Gompertz trajectory with participant-specific random effects;
the observed symptom score is the baseline minus the sum of the three
components plus AR(1) residual noise. The paper maps the
identifiability of each component onto trial-design phase contrasts
(open-label, blinded discontinuation, crossover), argues against
fitting a DGP-matched nonlinear model at trial-relevant sample sizes
in favour of a phase-augmented linear-mixed-effects approximation,
illustrates the inferential cost of lumping through worked numerical
examples, and pre-registers (but has not executed) a Monte Carlo
simulation programme calibrated to prazosin-PTSD reference parameters.
The stated central application is predictive-biomarker validation: the
claim that a biomarker must be shown to predict the BR component
specifically, not the lumped total, to qualify as a pharmacological
predictor.

## 2. Overall assessment and recommendation

The manuscript is clearly written, the clinical motivation is sound,
and the epistemic discipline around the unexecuted simulation is
exemplary: the abstract, Section 7, and the Discussion all correctly
label the empirical programme as pre-registered-but-pending, and the
pre-registration file exists and predates the pilot data. The
phase-contrast identifiability argument is the paper's strongest
methodological content.

However, the manuscript is **not yet at top-tier standard**, for three
reasons that jointly warrant major revision. First, there are two
**verified internal sign errors** in the formal content (Major
comments 1 and 2) that must be corrected before the formalisation can
be trusted. Second, the **novelty claim is overstated**: the additive
"natural-history + placebo + drug-action" decomposition by nonlinear
mixed-effects is the standard structure of the pharmacometric
disease-progression literature (Holford; Gomeni; Chen et al. 2021),
which the manuscript neither cites nor distinguishes itself from (Major
comment 3). Third, the **methodological positioning is dated**: the
paper handles carryover with a single exponential-decay scalar and
never states its estimand in ICH E9(R1) terms, while the closest
competitor — a distributed-lag model with autocorrelated errors in the
target journal (Liao et al. 2023) — is cited only as an example of
the lumped framing it shares, not engaged as a methodological rival
(Major comments 4 and 5). None of these is fatal; the core
contribution (design-based identification of the three components at
the aggregated-N-of-1 level) appears genuinely new and defensible once
re-scoped, but the paper needs substantial revision to earn the claim.

## 3. Significance and novelty

The clinical problem is real and the design-based identification
strategy is the right way to attack it. The difficulty is that the
*structural* decomposition the paper foregrounds as its contribution
is well-established prior art:

- **Holford (2015, *Br J Clin Pharmacol* 79:18-27)** consolidates the
  "clinical pharmacology = disease progression + drug action"
  tradition (foundational model Holford & Peace 1992, *PNAS*), in which
  the observed longitudinal outcome is the additive sum of natural
  disease progression, placebo effect, and drug effect, fitted by
  nonlinear mixed-effects with subject-level random effects. This is
  the manuscript's three-component additive structure.
- **Gomeni & Merlo-Pich and successors (*J Pharmacokinet Pharmacodyn*
  2007, 2009)** model the depression symptom trajectory as a
  disease/natural-history term plus a *separately parameterised
  placebo growth curve* (Weibull, inverse-Bateman), functionally
  analogous to the manuscript's Gompertz PB term.
- **Chen et al. (2021, *CPT Pharmacometrics Syst Pharmacol*)** is a
  recent worked example of exactly the additive baseline + linear
  progression + exponential-onset placebo + exponential-onset drug
  decomposition with subject random effects, differing only in the
  onset functional form and in not being N-of-1.
- **Muthén et al. (2009, *Stat Med* 28:3363-85)** separate drug from
  placebo response in longitudinal trajectories via growth-mixture
  modelling — same goal, different (latent-class) mechanism.
- **Kaptchuk et al. (2008, *BMJ* 336:999-1003)** already own the
  phrase "components of the placebo effect" for the within-placebo
  decomposition; the manuscript cites this but should be explicit that
  its BR-PB-TV partition is a *different* partition.

What survives as novel, on the searches run (Section 8), is the
**transplant of the structural decomposition into the aggregated /
single-trajectory N-of-1 setting with the three components identified
by trial-design phase contrasts** (open-label vs blinded vs crossover),
together with the explicit treatment of PB as a belief trajectory
distinct from natural history and modulated by an expectancy multiplier
$\eta(\text{phase})$. I located no paper performing the additive
three-component growth-curve decomposition at the N-of-1 level via
design-phase contrasts. This is the claim the paper should defend; the
generic "joint decomposition into pharmacological + placebo +
natural-history components" claim should be dropped or heavily
qualified.

**Bar for a top-tier methodology journal:** as re-scoped, the
design-based identification contribution is publishable in Statistics
in Medicine or Biometrical Journal *if* (a) the novelty is honestly
positioned against the pharmacometric literature above, (b) the
identifiability of three collinear growth curves from finite N-of-1
data is treated formally rather than asserted, and (c) the pending
simulation is actually run. It is not currently an Annals / JRSS-B
theoretical contribution and does not claim to be.

## 4. Major comments (correctness first)

**M1. Verified sign error in the subadditivity formalisation
(Section 3.2.1).** The generalised model is written
$Y_{it} = \mathrm{BL}_i - [BR + PB + TV + \gamma\, BR\cdot PB] +
\varepsilon_{it}$, so total improvement is
$S = \mathrm{BL} - Y = BR + PB + TV + \gamma\, BR\cdot PB$.
Subadditivity (combined improvement *less* than the sum of parts)
therefore requires $\gamma\, BR\cdot PB < 0$, i.e. $\gamma < 0$ (since
$BR, PB > 0$). The manuscript states the opposite: "$\gamma > 0$
encodes subadditivity ... $\gamma < 0$ encodes superadditivity." It
also defines subadditivity as $\partial Y / \partial(BR\cdot PB) < 0$
in two places (the bullet at the "additive model may be wrong"
paragraph and the orig paragraph). But
$\partial Y / \partial(BR\cdot PB) = -\gamma$, so subadditivity
($\gamma < 0$) gives $\partial Y / \partial(BR\cdot PB) > 0$; the
printed inequality describes superadditivity. *Why it matters:* the
entire subsection is built on subadditivity being the
empirically-supported case (Kaptchuk, Benedetti), and the downstream
biomarker-bias argument inherits the sign. *Remedy (verified
consistent):* keep the displayed equation; change "$\gamma > 0$ encodes
subadditivity" to "$\gamma < 0$" (and "$\gamma < 0$ encodes
superadditivity" to "$\gamma > 0$"), and change both
$\partial Y / \partial(BR\cdot PB) < 0$ statements to $> 0$. The
verbal biomarker-cancellation argument is correct under $\gamma < 0$
and needs no change. *(Applied in this revision.)*

**M2. Sign error in the PB heteroscedasticity claim (Section 3.2).**
The model sets $\mathrm{SD}(PB\mid\text{phase}) = \eta(\text{phase})
\cdot \sigma_{PB}^{\text{baseline}}$, so between-patient PB variance is
high in open-label ($\eta = 1$) and low under blinding
($\eta \approx 0.5$). The text states that ignoring this "inflates
standard errors in the open-label phase and deflates them in the
blinded phase." This is reversed. A homoscedastic fit uses a single
pooled variance intermediate between the two; in the high-variance
(open-label) phase the pooled value *understates* the truth, so the
model-based SE is too *small* (deflated/anti-conservative), and in the
low-variance (blinded) phase it is too *large* (inflated/conservative).
*Why it matters:* the claim is offered as a reason to model the
heteroscedasticity, and the stated direction would mislead a reader
about which phase's inference is anti-conservative. *Remedy:* swap
"inflates"/"deflates" in both the bullet and the orig paragraph.
*(Applied; confidence high for the phase-specific mean/contrast; the
author should confirm the direction holds for the specific coefficient
SE in the full GLS fit, which is the estimand the sentence names.)*

**M3. Overstated novelty; missing pharmacometric prior art (Sections
1 and 9).** The introduction repeatedly frames "explicit joint
decomposition ... simultaneous estimation of all three components" as
the methodological advance and states that "the published N-of-1
literature has, with rare exceptions, not adopted this
joint-decomposition framing." The latter is defensible; the former is
not, because the structural decomposition itself is standard
pharmacometrics (M3 references in Section 3 above). *Why it matters:* a
referee who knows the Holford/Gomeni literature will reject the novelty
claim as written. *Remedy:* (i) add a paragraph in Section 1 or 5
citing Holford 2015, Gomeni 2009, Chen et al. 2021, and Muthén et al.
2009, and distinguish the present method (single-trajectory N-of-1
identification by design-phase contrasts) from each (group-level
population PK/PD; latent-class mixtures); (ii) re-scope every novelty
sentence to the design-based N-of-1 identification, not the generic
decomposition.

**M4. Carryover is handled with a single exponential scalar; the
direct competitor is not engaged as a rival (Sections 1, 5).** Liao et
al. (2023, *Statistics in Medicine* 42:986-1004) model N-of-1
carryover via a Bayesian distributed-lag structure with autocorrelated
errors. The manuscript cites this paper only as an instance of the
lumped-response framing it shares, but a distributed-lag formulation is
a strictly richer carryover model than the
$D_{bc} = \exp(-\lambda t_{sd})$ scalar used here. *Why it matters:* in
the target journal, the obvious referee question is why a one-parameter
exponential decay is preferable to (or at least benchmarked against) a
distributed lag. *Remedy:* add a paragraph justifying the
exponential-decay choice and, ideally, include the distributed-lag
analysis as a comparator arm in the pre-registered simulation.
Daza (2018, APTE estimand) and Selukar/Daza et al. (2023, BMC Med Res
Methodol, G-estimation vs Bayesian vs linear under carryover) are the
other two carryover-specific references the paper should position
against.

**M5. No estimand statement; the interaction coefficient is never
formally defined as a target (Sections 5, 8).** The paper's inferential
object is the $bm{:}D_{bc}$ interaction, but it is never defined as an
estimand in the ICH E9(R1) (2019) sense (population, treatment,
endpoint, summary, intercurrent-event strategy), and dropout/carryover
are not assigned an intercurrent-event strategy. Relatedly, the
"pharmacological component" BR is treated as a well-defined causal
quantity without a formal identification statement; the estimands
addendum (principal-stratum or hypothetical strategies) is the natural
vocabulary for defining "the effect that would obtain absent the
belief-driven channel." *Why it matters:* defining a pharmacological
*component* as a causal estimand is exactly where a methodological
referee will press hardest, and the paper currently asserts rather than
identifies it. *Remedy:* add an explicit estimand statement for the
biomarker-interaction target and a sentence locating the BR component
within an estimands/principal-stratum framing. Engage Senn (2004, *BMJ*
329:966-8; 2018, *Nature* 563:619-21) on why separating a true
participant-level interaction from within-person variation is
non-trivial — this strikes directly at the paper's central premise and
is currently cited only for the participant-recovery caveat in the
Discussion, not for the headline interaction claim.

**M6. Identifiability of three collinear growth curves is asserted, not
demonstrated.** Section 4 gives the phase-contrast logic clearly, but
the three Gompertz trajectories share onset shapes and overlap in
calendar time, and the paper does not show that $(m_{BR}, m_{PB},
m_{TV})$ and their rate/displacement parameters are jointly identified
at trial-relevant $N$ beyond the marginal-maxima argument. The pilot
already exposed rank-deficiency of the full phase-augmented formula at
$N = 35$ (the $D_{bc}{:}\text{phase}$ term had to be dropped), which is
itself evidence the identification is fragile. *Remedy:* the
pre-registered identifiability deliverable (convergence rates, VIFs,
$\eta$-sensitivity) should be reported as actual results, and the
manuscript should state the minimal design (number of phases, off-drug
window length) under which the three components are jointly identified.

## 5. Minor comments

**m1. Results placed in the Discussion (Sections 7.1, 9).** The
prototype's actual numbers (power 0.95/0.93/0.83 vs 0.65/0.65/0.61;
coefficient means $-2.25$ vs $-2.21$) appear in the Discussion
subsection "Prototype Monte Carlo confirmation," while Section 7.1
defers them forward. Move the prototype results into Section 7 and
reserve the Discussion for their implications.

**m2. "Pilot" and "prototype" name the same 6-cell/100-replicate run.**
Unify the terminology.

**m3. Off-drug contrast wording (Section 4).** "All phases off drug
yield $TV$ alone, with no $BR$ contribution and $PB$ at the off-drug
expectation level" is self-contradictory; off-drug open-label still
carries residual $PB(\eta_{\text{off}})$, so the signal is
$TV + PB(\eta_{\text{off}})$, not "$TV$ alone." Reword.

**m4. Worked example holds $TV$ constant from week 8 to week 12**
(tables in Section 6), implying natural history plateaued over the
blinded window, while Section 3.3 states $TV$ accumulates with calendar
time. Add a one-clause acknowledgement of the simplification.

**m5. Zunhammer et al. (2018) magnitude.** The text says placebo
produces "large reductions in reported pain"; the source reports a
*moderate* behavioural effect (Hedges $g \approx -0.66$) and a very
small NPS effect ($g \approx -0.08$). Change "large" to "moderate."
*(Applied.)*

**m6. Hróbjartsson trial/condition counts.** "more than two hundred
trials and sixty conditions" are the figures of the 2010 Cochrane
review (202 trials, 60 conditions), not the 2001 NEJM paper (130
trials). The grouped citation `[@hrobjartsson2001; @hrobjartsson2004;
@hrobjartsson2010]` should attach the counts specifically to
`@hrobjartsson2010`. Two occurrences (Sections 1.2 and 3.3).

**m7. Bibliography hygiene — `wang2022comt`.** The author field reads
`... and Kaptchuk, Ted J. and others and Hall, Kathryn T.`. In
BibTeX/biblatex, `and others` must be the final element to render as
"et al."; mid-list it is parsed as a literal surname and corrupts the
rendered author list. *(Applied: moved `and others` to the end.)*

**m8. README path is stale (not the manuscript).** The paper's
`README.md` names the driver directory
`analysis/scripts/06-component-decomposition/`, but the manuscript and
the actual pre-registration file use
`analysis/scripts/component-decomposition/` (no numeric prefix). The
manuscript is correct; update the README.

**m9. CENT 2015 reporting standard.** The paper cites the CENT
extension; confirm it cites both `@vohra2015cent` and the E&E
companion `@shamseer2015cent` (it does). No change required; noted for
completeness.

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| §1, §5 novelty framing | Structural decomposition presented as novel; pharmacometric prior art absent | Holford NHG (2015) *Br J Clin Pharmacol* 79:18-27; Holford & Peace (1992) *PNAS* 89:11466-70 |
| §3.2 PB-as-growth-curve | Placebo modelled as a distinct parametric trajectory is prior art | Gomeni R, Merlo-Pich E (2009) *J Pharmacokinet Pharmacodyn* 36:?? (verify pages); Chen et al. (2021) *CPT Pharmacometrics Syst Pharmacol* |
| §3.2.1 drug-vs-placebo separation in trajectories | Alternative (latent-class) approach uncited | Muthén B et al. (2009) *Stat Med* 28:3363-85 |
| §1, §5 carryover model | Distributed-lag rival not engaged as a methodological alternative | Liao Z et al. (2023) *Stat Med* 42:986-1004; Daza EJ (2018) *Methods Inf Med* 57:e10-e23; Selukar S et al. (2023) *BMC Med Res Methodol* 23:203 |
| §5, §8 interaction estimand | No estimand definition | ICH E9(R1) Addendum (2019), EMA/CHMP |
| §8 biomarker-interaction premise | Heterogeneity-identifiability critique cited only late | Senn S (2004) *BMJ* 329:966-8; Senn S (2018) *Nature* 563:619-21 |
| §1 motivation | Canonical aggregated-N-of-1 precision-medicine citation absent | Schork NJ (2015) *Nature* 520:609-11 |
| §3.2 $\eta = 0.5$ assumption | "not directly measured" is defensible; nearest empirical anchor uncited | Open-hidden paradigm (Colloca, Benedetti, Amanzio); open-hidden discontinuation protocol NCT05191277 (PMC10286435) |
| §1 design comparison via simulation | Direct simulation-comparison precedents uncited as such | Blackston JW et al. (2019) *Healthcare* 7:137; Carrozzo E et al. (2025) *Biom J* 67:e70045 |

*All suggested entries verified to exist and to support the attached
claim via literature search (Section 8); page/volume details for Gomeni
2009 should be confirmed against the primary PDF before entry.*

## 7. Suggestions for strengthening the paper

1. **Run the pre-registered simulation.** The empirical core is absent;
   the paper currently rests on two worked arithmetic examples and a
   100-replicate structural check that inverted its own prediction
   (correctly explained as rank-deficiency, but still). Production
   results at the planned $N \in \{35,70,100,150\}$, with the
   bias-vs-power trade-off read off the coefficient mean as the pilot
   recommended, would move the paper from "framework" to "method."
2. **Add a distributed-lag comparator** (Liao 2023) to the carryover
   arm of the simulation, so the exponential-decay scalar is benchmarked
   rather than merely asserted.
3. **State the biomarker-interaction estimand** under ICH E9(R1) and
   give a one-paragraph principal-stratum/hypothetical identification
   of the BR component.
4. **Demonstrate joint identifiability** of the three growth curves
   (convergence, VIFs, $\eta$-sensitivity) as actual results, and state
   the minimal identifying design.
5. **Reproducibility:** the seed-handling and per-cell replicate counts
   are well documented; when the production run lands, report Monte
   Carlo SEs on every headline quantity (the pilot already computes
   MCSE $\approx 0.046$ at 100 reps), and version the driver under the
   path the manuscript names.

## 8. Scope of this review

**Verified (re-derived or executed):** the two sign errors (M1, M2)
were derived algebraically from the manuscript's own equations; the
worked-example arithmetic in Section 6 (open-label totals 11.0/10.0,
blinded totals 4.0/1.5, lumped drug effect 7.75, recovered mean BR
6.0) checks out; the existence and pre-pilot timestamp of the ADEMP
pre-registration file were confirmed on disk; the `wager2013nps` bib
entry was confirmed to point to the NEJM NPS paper (not the placebo
review), so that citation is correct.

**Inspected (read and judged):** notation consistency across sections;
the section cross-references (all consistent); the citation-to-claim
fit for the clinical literature, fact-checked against source abstracts
and PubMed metadata (claims for Wager 2013, Bingel 2011, Hall 2012,
Kirsch 2008, Kaptchuk 2008 are accurate; Zunhammer magnitude and
Hróbjartsson count-attribution need the tightenings in m5/m6); the
`wang2022comt` author-field corruption.

**Not checked:** I did not execute the pmsimstats-ng simulation code or
re-fit any model; the convergence and rank-deficiency claims from the
pilot are taken as reported. Page/volume numbers for several suggested
references (notably Gomeni 2009) are from search metadata, not
page-verified. I did not assess the companion pedagogy document or the
07-gompertz-evaluation paper.

**Literature searches run:** treatment-response decomposition into
pharmacological/placebo/natural-history components; pharmacometric
disease-progression + placebo + drug-action structural models (Holford,
Gomeni, Bateman/Weibull); growth-mixture decomposition of drug vs
placebo (Muthén, Tarpey, Boessen); Kaptchuk/Benedetti within-placebo
decomposition; balanced placebo design; mediation of expectancy vs
drug; post-2015 N-of-1 methodology (Bayesian aggregation, distributed
lag, G-estimation, causal N-of-1); ICH E9(R1) estimands; Senn on
individual-treatment-effect identifiability; PATH statement and
successors; CENT reporting standard; open-hidden discontinuation
quantification of residual expectancy. Closest prior art located:
Holford 2015 and Gomeni 2009 (structural decomposition);
Chen et al. 2021 (recent additive worked example);
Muthén et al. 2009 (latent-class alternative). No paper performing the
additive three-component growth-curve decomposition at the N-of-1
level via design-phase contrasts was found — an absence-of-evidence
finding, not proof of novelty.
