# Referee Report: Three-component decomposition of treatment response
in aggregated N-of-1 trials

*2026-06-13 08:27 PDT*

**Manuscript:** "Three-component decomposition of treatment response
in aggregated N-of-1 trials: pharmacological, expectation-driven, and
natural-history components"

**Journal orientation:** *Statistics in Medicine*

**Recommendation:** Major revision

---

## 1. Summary of the manuscript

The manuscript proposes a formal three-component decomposition of the
within-participant outcome trajectory in aggregated N-of-1 clinical
trials into a biological response component (BR, pharmacological), a
placebo-belief response component (PB, expectancy-driven), and a
time-variant natural-history component (TV). Each component is
modelled as a participant-specific modified Gompertz function, and
their sum, plus AR(1) residual noise, constitutes the observed
symptom trajectory. The paper maps identifiability of each component
onto trial-design contrasts supplied by the hybrid
open-label-plus-blinded-discontinuation design; argues that a
phase-augmented linear-mixed-effects model is the tractable
small-sample approximation; and advocates the decomposition as a
precondition for valid predictive-biomarker validation. A 100-replicate
pilot confirms the analysis pipeline runs end-to-end; the production
simulation (1,000 replicates per cell, $N \in \{35, 70, 100, 150\}$)
is described in detail but not yet executed. The paper claims
consolidation and extension rather than novel invention.

---

## 2. Overall assessment and recommendation

**Major revision.** The conceptual framework is coherent and the
application of a pharmacometric disease-progression decomposition to
the aggregated N-of-1 setting is a genuine and underexplored
contribution. The writing is precise and the epistemic hedging is
unusually disciplined for a framework paper. These are strengths.

However, the manuscript has one mathematical error and three
structural gaps that prevent acceptance in the current form.

First, the Gompertz ceiling interpretation is incorrect as stated
(Major Comment 1). Second, the identifiability argument conflates
design-contrast identification with temporal identification, and the
time-ordering gap is consequential because it undermines the central
analytical claim (Major Comment 2). Third, the characterisation of the
one-component analysis as producing 'biased' estimates conflates
estimand mismatch with statistical bias throughout the manuscript, and
the distinction matters for the paper's regulatory claims (Major
Comment 3). Fourth, and most critically for a *Statistics in Medicine*
submission, the empirical core of the manuscript -- the simulation
study -- has not been executed. The pilot is correctly described as a
structural check, but the paper in its current state cannot be
evaluated as an empirical methods paper; the bias-versus-power
trade-offs that are the paper's central quantitative claims are
explicitly described as conjectures. The manuscript should either be
reframed as a framework-and-pre-registration paper or held for
submission until the production simulation results are available.

The four open-question bullets in the Discussion are admirably honest
and are themselves a reason to require revision rather than rejection:
the authors know exactly what the paper is missing. The revision task
is to either supply what is missing or reframe the paper accurately
around what is present.

---

## 3. Significance and novelty

The authors correctly disclaim novelty for the decomposition itself.
The pharmacometric tradition of separate disease-progression, drug, and
placebo models is well established (Holford 2015; Gomeni 2009; Chen
2021, all cited). The paper's specific contribution is:

1. The explicit application of this structure to the aggregated N-of-1
   setting, where the components must be recovered from a single
   participant's trajectory rather than a parallel-arm average.
2. The formal mapping of each component's identifiability onto
   design-contrast conditions, including the hybrid blinded-
   discontinuation phase.
3. The phase-augmented LME as a degree-of-freedom-conserving
   small-$N$ proxy.
4. The argument that predictive-biomarker validation requires
   component decomposition as a logical precondition.

Contribution 1 and 4 are the most significant; they are not present
elsewhere in the N-of-1 methodology literature as reviewed here.
Contributions 2 and 3 are analytically sound but require the
production simulation to be quantitatively evaluated.

One paper directly relevant to the novelty claim was not cited:
Petkova, Tarpey, and Govindarajulu (2009, *International Journal of
Biostatistics*, 'Predicting potential placebo effect in drug treated
subjects') develops a linear-mixed-effects framework for decomposing
drug-treated outcomes into specific (pharmacological) and non-specific
(placebo) components in longitudinal trials. The parallel to the
present paper's Section 5 is close enough that the authors should
engage with it and characterise the difference explicitly.

For a top-tier journal, the contribution is plausible but the
manuscript must arrive with its simulation results. A framework-only
paper of this length would fit a specialty methods journal or a
short-form empirical letter to *Statistics in Medicine* better than a
full paper. If the simulation is completed and the mathematical issues
are resolved, a full *Statistics in Medicine* submission is
appropriate.

---

## 4. Major comments

**Major Comment 1. Mathematical error: the Gompertz ceiling
interpretation (Sections 3.1, 3.2; bullet block preceding Eq. 2,
and the calibration discussion of $m_{BR} \approx 11$)**

The modified Gompertz used throughout is:

$$
C(t) = m_C \exp(-d_C \exp(-r_C t)) - m_C \exp(-d_C).
$$

The offset $m_C \exp(-d_C)$ ensures $C(0) = 0$. The actual asymptote
as $t \to \infty$ is:

$$
C(\infty) = m_C \exp(0) - m_C \exp(-d_C) = m_C \bigl(1 - \exp(-d_C)\bigr),
$$

not $m_C$. The paper states in the bullet block (Section 2, notation):
'rescaled so that $C(0) = 0$ and $C(t) \to m_C$ as $t \to \infty$' --
this is incorrect. The paper then repeatedly describes $m_{BR}$ as 'the
biological ceiling, the largest reduction the drug can produce at
infinite exposure' and calibrates $m_{BR} \approx 11$ as 'eleven
points of nightmare-score reduction at saturation.'

For $d_{BR} \approx 5$, the error is small ($\exp(-5) \approx 0.0067$,
so the true asymptote is $\approx 0.993 \cdot m_{BR}$), but the
interpretation is wrong in principle. For smaller $d$ values, which
arise in sensitivity analyses or for TV trajectories with gradual
onset, the discrepancy is substantial. For $d = 1$, $C(\infty) = m_C
(1 - e^{-1}) \approx 0.632 m_C$, less than two-thirds of $m_C$.

**Why it matters.** The parameter $m_C$ is the regression coefficient
being modelled as a random effect; if readers interpret it as the
saturation value they will misread likelihood surfaces, prior
distributions, and identifiability conditions. Any participant-level
prediction built on $m_{BR}$ is miscalibrated.

**Remedy.** Either (a) correct the notation so that $m_C$ is
identified as the Gompertz scale parameter (not the asymptote) and the
correct asymptote $m_C (1 - e^{-d_C})$ is stated explicitly wherever
'ceiling' or 'maximum reduction' is discussed; or (b) reparameterise
the Gompertz so that the ceiling parameter is directly interpretable
as the asymptote (replacing $d_C$ with an explicit S-shape
parameterisation that forces the asymptote to be the stated maximum).
Option (b) is preferable for clarity.

---

**Major Comment 2. Temporal confounding gap in the identifiability
argument (Section 4)**

The paper claims that the open-label-versus-blinded contrast 'yields
$PB(\eta = 1) - PB(\eta = 0.5) \approx 0.5 \cdot PB$, which isolates
the belief component.' This argument treats the two phases as if the
Gompertz component $PB$ is at the same point in its trajectory in both.
It is not: the open-label phase is by design earlier in the trial
timeline than the blinded-discontinuation phase. At the time of the
blinded phase, $PB$ has been accumulating since trial entry via the
Gompertz function, so:

$$
PB_{\text{open}}(t_1) \;\neq\; PB_{\text{blinded}}(t_2) / 0.5
$$

unless $PB$ is at steady state in both phases (full Gompertz
saturation). The contrast conflates two changes simultaneously: the
drop in $\eta$ and the continued accumulation of the Gompertz
component over time. At $N = 35$ and a 12--16 week trial, steady-state
saturation of a Gompertz with $r_{PB}$ calibrated to the same scale as
$r_{BR} \approx 0.42$/week is unlikely for all participants.

The same issue affects the identifiability of $BR$: the 'blinded
on-drug versus blinded-placebo' contrast isolates $BR$ only if $TV$
and $PB$ are identical across the two allocation arms within the
blinded window. If participants are not randomly assigned to blinded
arms at the same timepoint within the Gompertz trajectories, or if
allocation is sequential rather than concurrent, the contrast is
confounded by differential Gompertz accumulation.

**Why it matters.** The paper's claim that the three components are
separately identified from the hybrid design is the central
methodological claim. If the identification argument holds only at
steady state, the design-contrast table in Section 4 requires explicit
steady-state assumptions that are currently absent, and the pre-
registered simulation must test identifiability as a function of how
close the trial timeline brings each component to steady state.

**Remedy.** (a) State explicitly that the open-label-versus-blinded
identification of $PB$ assumes either (i) $PB$ is at Gompertz
saturation in both phases, or (ii) the time argument of $PB$ is reset
at phase entry (i.e., $PB$ accumulates from zero at the start of each
phase). State which assumption the framework adopts and add it to the
bullet preceding the design-contrast table. (b) Add a pre-registered
sensitivity analysis to the Section 7 production simulation that varies
the degree of Gompertz saturation at the time of each phase transition
and quantifies identification loss as saturation decreases.

---

**Major Comment 3. 'Bias' versus estimand mismatch throughout
(Sections 1, 6, 8, and Conclusions)**

The manuscript repeatedly describes the one-component analysis as
producing 'biased' estimates of the pharmacological effect. For
example, Section 6: 'the one-component model has produced a single
number conflated with the placebo response'; the comparison table in
Section 6 reads 'Population mean drug effect? 7.75 points (biased).'

This framing is imprecise. The one-component analysis is not biased in
the statistical sense for the estimand it is targeting, which is the
average lumped treatment contrast. Its estimate of that contrast (7.75
points in the worked example) is consistent and converges to the true
value of $\mathbb{E}[BR + 0.5 \cdot PB]$ under the design's blinded
contrast, which is not $\mathbb{E}[BR]$. What the paper correctly
identifies is estimand mismatch: the one-component estimand is a
different quantity from the pharmacological component $\mathbb{E}[BR]$
that the precision-medicine question requires. This is the
'wrong-estimand' problem, not bias in a misspecified model.

The distinction matters clinically and regulatorily. The ICH E9(R1)
estimands framework (cited in Section 5) is precisely the machinery
for making this explicit: the paper should state that the
one-component analysis targets a valid estimand (the population-level
average treatment contrast under the blinded design), but that the
precision-medicine question requires a different estimand ($\beta^{BR}$
from the component regression) that the one-component model cannot
supply. Calling the one-component result 'biased' without this
qualification misrepresents what the one-component analysis does and
overstates the paper's critique of standard practice.

**Why it matters.** The paper's advocacy case for the decomposition
rests on showing that the one-component analysis gives the wrong
answer. If reviewers and practitioners understand 'wrong' as 'biased
estimator', they will question whether the paper has established this,
because for the one-component estimand it has not. The paper will be
more persuasive, not less, if it correctly states that the
one-component analysis answers a different question -- and then argues
that the precision-medicine question requires answering the right one.

**Remedy.** Replace 'biased' with 'estimand-mismatched' or 'targeting
the wrong estimand' throughout, and explicitly invoke the ICH E9(R1)
framework in Section 6 as well as Section 5 to characterise the
difference between the lumped estimand and $\beta^{BR}$.

---

**Major Comment 4. The empirical core is absent (Section 7)**

The manuscript's central quantitative claims -- that the phase-
augmented analysis reduces bias in the biomarker-treatment interaction
estimate, that the full decomposition is identifiable at specific
sample sizes, and that type I error is controlled -- are explicitly
described as forthcoming (Section 7.2: 'the production programme
described in the following section is the empirical core of this
manuscript and has not yet been executed'). The Discussion (Section
8.3) confirms: 'quantitative bias-versus-power trade-offs are
conjectures until that work is run.'

The pilot (100 replicates, 6 cells) reported in Section 7.1 ran in the
opposite direction from the paper's prediction: the phase-augmented
analysis was less powerful than the one-component analysis across the
entire $pb\_strength$ grid. The authors' explanation -- rank-deficiency
at $N = 35$ and coefficient mean rather than rejection rate as the
right summary -- is plausible, but it has not been tested in the
production simulation.

A *Statistics in Medicine* submission in which the primary empirical
claims are explicitly unconfirmed and the one executed pilot ran
counter to the prediction is not ready for peer review as a completed
methods paper.

**Remedy (two options).** (a) Execute the production simulation and
revise the manuscript to report results. This is strongly preferred.
(b) Reframe the manuscript explicitly as a framework-and-pre-
registration paper, remove claims that 'the decomposition reduces bias'
from the abstract and conclusions, and submit the results paper as a
follow-on. Under option (b), the abstract's 'Results' paragraph ('a
100-replicate pilot of the analysis machinery is reported') is
accurate; the Conclusions paragraph ('The three-component framework is
the principled basis...') makes claims that go beyond the presented
evidence and must be hedged as conjectures pending simulation.

---

## 5. Minor comments

**Minor 1. Cross-component covariance notation (Section 2, orig block
following the AR(1) residual definition).** The covariances are named
$c_{cf1t}$ and $c_{cfct}$ but these subscripts are not defined and
appear inconsistent internally. Introduce a clear subscript convention
(e.g., $\text{Cov}(C_1(t), C_2(t))$ for within-time cross-component
and $\text{Cov}(C_1(t), C_2(t'))$ for cross-time cross-component).

**Minor 2. The expectancy multiplier $\eta = 0.5$ is stated as 'a
deliberately central placeholder' (Section 3.2).** This is an honest
acknowledgement, but the design-contrast table in Section 4 uses this
specific value to characterise the blinded-phase $PB$ contribution.
The table should either carry a note that the listed components are
$\eta$-dependent, or the value should be replaced by $\eta$ in the
table to make the parameter dependence explicit.

**Minor 3. Subadditivity sign convention (Section 3.2.1, orig block
before and after the generalised model equation).** The manuscript
states '$\partial Y / \partial(BR \cdot PB) > 0$' as the condition for
subadditivity. Given the sign convention $Y = \mathrm{BL}_i - [BR +
PB + TV + \gamma \cdot BR \cdot PB] + \varepsilon$, this requires
$-\gamma > 0$, i.e., $\gamma < 0$. The text correctly states '$\gamma
< 0$ encodes subadditivity', but the derivative expression should be
written as $-\gamma > 0$ or equivalently as $\partial(-\Delta Y) /
\partial(BR \cdot PB) < 0$ to match the more natural framing in which
subadditivity means the total symptom improvement is less than the sum
of the parts. The current derivative expression requires the reader to
track the sign of $Y$ (symptom score, higher = worse) while
interpreting the interaction coefficient.

**Minor 4. The worked example holds $TV$ constant across phases
(Section 6, orig block for the blinded-phase table).** The text
acknowledges this ('for simplicity the table holds $TV$ at its week-8
value') but does not quantify the approximation error. For the prazosin
calibration, $TV$ continues to accumulate during weeks 9--12. If the
TV Gompertz rate is comparable to $r_{BR}$, the off-drug $TV$
increment over four weeks is non-negligible. State explicitly what
fraction of the total $BR$-vs.-$PB$ contrast the omitted $TV$
increment would represent.

**Minor 5. Rank-deficiency explanation (Sections 7.1, Discussion
8.3).** The rank-deficiency of the full phase-augmented formula at
$N = 35$ is described but the specific structural reason is not. Is
the rank-deficiency due to collinearity of $D_{bc}$ with phase
(because $D_{bc}$ is constant within each phase for many participants)?
Or due to small cell counts in one phase? Identifying the specific
source is necessary for the production specification to guarantee the
formula is not rank-deficient at larger $N$.

**Minor 6. Format: the rgt blocks throughout contain Lorem ipsum
placeholder text.** All `\begin{rgt}...\end{rgt}` environments are
unfilled. This is expected in a draft but must be resolved before
submission.

**Minor 7. Data availability statement.** The statement says data are
'openly available in the pmsimstats-ng project repository' but no URL,
DOI, or repository identifier is provided. For *Statistics in Medicine*
the repository must be publicly accessible and cited with a persistent
identifier at submission.

---

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| Section 1.2, claim that 'no published N-of-1 trial has employed a structural decomposition' | Strong claim requiring a systematic review; the referenced Li 2016 audit covers 1985--2013. A recent audit covering 2014--2024 would strengthen the claim or reveal exceptions. | Extend the Li 2016 audit or cite a more recent scoping review. |
| Section 5, fourth reason (estimand robustness), PATH citation | PATH [@kent2020path] focuses on statistical subgroup and risk-prediction methods; the claim that PATH 'defers the mechanism of treatment-effect heterogeneity to domain knowledge' is an interpretation that PATH does not state explicitly and that Kent 2018 BMJ [@kent2018bmj] addresses more directly. The citation is defensible but should be more precise. | Reframe to cite @kent2018bmj for the 'domain knowledge' point. |
| Section 3.2, BR-PB interaction subsection, 'subadditivity' citations [@kaptchuk2008components; @benedetti2014] | The Kaptchuk 2008 paper decomposes the placebo response itself, not a BR-PB interaction. The Benedetti 2014 review discusses placebo mechanisms generally, not specifically subadditivity. Neither directly documents the subadditivity of drug and placebo effects. The specific claim should cite the Bingel 2011 paper [@bingel2011] (already in the bibliography) and potentially Amanzio & Benedetti (1999, PNAS) on the partial-reversal-by-naloxone paradigm, which is the clearest evidence for non-additive drug-placebo interaction. | Replace or supplement with Bingel 2011 for the subadditivity claim; cite Amanzio & Benedetti (1999) if the naloxone-reversal evidence is invoked. |
| Section 1.2, 'most recent methodological advances retain this character' | Liao 2023 [@liao2023] is cited; Daza 2018 [@daza2018] and Gartner 2023 [@gartner2023] are cited later in Section 5 but not in Section 1.2. | Extend the Section 1.2 review to include Gartner 2023 for completeness. |
| Section 5, phase-augmented LME, 'distributed-lag analysis arm' acknowledged but not included in the pre-registration | The paper cites Liao 2023 as the reference for distributed-lag models. The Bayesian distributed-lag framework is a well-developed alternative; its omission from the produced pilot and the forthcoming simulation is noted by the authors themselves. | No new reference needed, but the pre-registration document should explicitly include the distributed-lag arm as a pre-specified comparison before the production run. |
| Uncited but directly relevant: Petkova, Tarpey, and Govindarajulu (2009) | 'Predicting potential placebo effect in drug treated subjects', *Int J Biostatistics* 5(1). This paper develops an LME decomposition of drug-treated outcomes into specific and non-specific (placebo) components and is the closest precursor in the statistical literature to Section 5. | Cite and discuss in the Introduction (Section 1.3, gap discussion) and Section 5. |

---

## 7. Suggestions for strengthening the paper

**S1. Execute the production simulation before resubmission.** The
framework's quantitative claims rest on the production simulation.
Without it, Sections 7.2 and the Conclusions are pre-registration
documents, not research findings. A *Statistics in Medicine* paper
requires the results to be present. Even a reduced production run
(500 replicates, $N \in \{70, 100\}$, the three primary analysis
strategies) would establish the central bias-versus-power finding.

**S2. Add a sensitivity analysis for $\eta$.** The expectancy
multiplier $\eta = 0.5$ is the single most consequential unvalidated
parameter in the framework. A simulation axis varying $\eta \in \{0.25,
0.5, 0.75\}$ for the blinded phase would characterise how sensitive
the $BR$ estimate is to this assumption. This is described as
'mandatory' in the Discussion but is absent from the pre-registration
specification. Add it to the production plan.

**S3. Resolve the Gompertz collinearity problem with a convergence
map.** The paper acknowledges (Section 4, identifiability, orig block)
that 'the symptom of near-collinearity is the rank-deficiency the
pilot of section 7 encountered'. A convergence map -- the proportion
of cells achieving successful NLME convergence as a function of $N$,
$r_C$, and time-per-phase -- would delineate the feasible operating
region of the full decomposition, and this is described as a
pre-registered deliverable. Ensure it is executed as the second
deliverable before resubmission.

**S4. Consider adding a Bayesian hierarchical implementation.** The
paper mentions `brms` as a feasible large-$N$ option but does not
present any Bayesian results. For the full decomposition, a Bayesian
fit with informative priors on the Gompertz parameters derived from the
prazosin-PTSD calibration would almost certainly outperform the NLME
convergence rate at $N = 35$--$70$. Including a Bayesian arm in the
simulation comparison would expand the paper's scope and address the
convergence-fragility concern more directly.

**S5. Include a brief formal identification argument.** Section 4
uses a design-contrast table to argue identifiability. A brief formal
argument -- e.g., specifying the observation model as a function of
$(m_{BR}, m_{PB}, m_{TV}, \eta)$ and showing that the system of
phase-by-allocation moment equations has a unique solution under the
stated assumptions -- would make the identifiability claim
mathematically rigorous. This argument is available and standard in
the longitudinal identification literature (e.g., by analogy with
parallel-group instrumental variables); the present text states it
only heuristically.

---

## 8. Scope of this review

**Verified (re-derived or directly checked):**

- Gompertz zero constraint $C(0) = 0$: verified algebraically.
- Gompertz ceiling $C(\infty) = m_C(1 - e^{-d_C}) \neq m_C$:
  verified algebraically; this is the basis for Major Comment 1.
- Worked-example arithmetic (Sections 6, 7.75-point average,
  $BR_A = 4$, $BR_B = 8$): verified by hand from the stated component
  values.
- $\mathrm{SD}(PB|\text{phase}) = \eta \cdot \sigma^{\text{baseline}}_{PB}$:
  verified as a consequence of the multiplicative $\eta$ scaling.
- Subadditivity sign convention with $\gamma < 0$: verified consistent
  with the stated partial derivative condition, though wording is
  imprecise (Minor Comment 3).

**Inspected (read and judged, not executed):**

- The temporal-confounding gap in the open-label-versus-blinded
  identification argument (Major Comment 2): identified by reasoning
  about Gompertz time accumulation across phases; not confirmed by
  simulation.
- Estimand-mismatch framing (Major Comment 3): assessed by comparison
  with ICH E9(R1) language and internal consistency of the manuscript.
- Cross-component covariance notation inconsistency (Minor Comment 1).
- Pilot finding of inverted power (less power for phase-augmented at
  $N = 35$): reported by the authors and not independently confirmed.

**Not checked:**

- The simulation code was not inspected.
- The ADEMP pre-registration file was not read.
- The pilot output files (`06-decomp-replicates.rds`,
  `06-decomp-summary.txt`) were not opened.
- The companion pedagogy document was not read.

**Literature searches conducted:**

- PubMed: 'N-of-1 trials placebo drug response decomposition' (0 hits
  on exact query; expanded search yielded no direct matches).
- PubMed: 'pharmacometric disease progression placebo drug separation
  longitudinal' (relevant hits: Chen 2021 and Goteti 2023, confirming
  that the pharmacometric literature has the component separation
  tradition but not in the N-of-1 setting).
- PubMed: 'treatment response components placebo pharmacological
  natural history' (Petkova 2009 identified as uncited but directly
  relevant; see Missing References table).
- PubMed: 'three component treatment effect decomposition clinical
  trial' (no direct matches beyond Petkova 2009).

---

*Prepared by external referee, 2026-06-13.*
