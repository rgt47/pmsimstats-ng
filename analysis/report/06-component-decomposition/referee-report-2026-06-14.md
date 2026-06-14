# Referee Report: Three-Component (BR-PB-TV) Decomposition of Treatment Response
*2026-06-14 09:40 PDT*

Manuscript: `analysis/report/06-component-decomposition/report.Rmd`
(pmsimstats-ng). Reviewed against the standard of a top methodological
journal (Statistics in Medicine / Clinical Trials / Biometrics), with
attention to the author's claim that the contribution applies to all
placebo-controlled longitudinal trials, not only aggregated N-of-1
designs.

---

## 1. Summary of the manuscript

The paper argues that treatment response in aggregated N-of-1 trials is
the sum of three causally distinct components -- pharmacological action
(BR), placebo-belief response (PB), and time-variant natural history
(TV) -- and should be analysed jointly rather than as a single lumped
treatment effect. It formalises each component as a modified Gompertz
trajectory with participant-specific random effects (Sections 2-3),
maps the identifiability of each component onto the phase contrasts of a
hybrid open-label-plus-blinded-discontinuation design (Section 4),
proposes a phase-augmented linear-mixed-effects approximation tractable
at trial-relevant sample sizes (Section 5), derives an
omitted-variable-bias identity for the one-component estimator and
illustrates it numerically (Section 6), and reports a pre-registered
Monte Carlo study (Section 7, "Study A," 250 replicates/cell). The
central application is predictive-biomarker validation (Section 8): a
biomarker that predicts the lumped response need not predict the
pharmacological component. The Discussion (Section 9) claims the
decomposition generalises to any longitudinal design with belief- and
exposure-varying contrasts.

## 2. Overall assessment and recommendation

**Recommendation: Major revision** (and, candidly, a recommendation to
target Statistics in Medicine / Clinical Trials rather than a
theory-first journal).

The framework is coherent, the central identity is correct, and the
question -- how much the conventional one-indicator analysis distorts
pharmacological inference -- is genuinely important. The manuscript is
also unusually honest: it reports a null result and explicitly corrects
its own motivating framing. But three issues block acceptance at top-tier
standard in the current draft, and the first two are serious:

1. **The empirical core does not yet support the headline claim.** The
   only executed simulation (Study A) is, by construction, a null on the
   innocuous axis; the one experiment that would demonstrate the
   biomarker-validation failure -- a biomarker coupled to PB -- is not
   run. As written, the paper's own evidence undercuts its Section 1
   motivation rather than confirming it.
2. **The broad-applicability claim is unsupported.** The manuscript
   conflates the existence of the three components in the
   data-generating process (true of essentially any trial) with their
   identifiability from a given design (true only of designs supplying
   within-subject belief x exposure contrasts). Standard
   parallel-group placebo-controlled longitudinal trials do not supply
   the contrast the participant-level decomposition needs.
3. **The novelty is framework/synthesis-level, not theoretical.** The
   key identity is elementary, and the trajectory decomposition has a
   substantial pharmacometric precedent the paper under-engages.

None of these is fatal. The framework, a completed simulation including
the decisive cell, a rigorous identifiability taxonomy across design
classes, and honest positioning would make a strong applied
methodological paper.

## 3. Significance and novelty

**What is genuinely new.** The specific contribution is (a) an explicit,
separately *identified* placebo-belief channel, pinned to a
within-subject belief contrast (the blinded/silent-discontinuation
phase) rather than inferred from a between-arm placebo comparison; and
(b) the framing of predictive-biomarker validation as recovery of the
biomarker-by-BR slope, with the lumped slope shown to be displaced by
the biomarker's covariance with PB and TV. The synthesis -- pulling
disease-progression modelling, placebo-component separation, and the
PATH treatment-effect-heterogeneity agenda into one decomposition for
the N-of-1 setting -- is useful and, to my knowledge, not assembled
elsewhere in exactly this form.

**What is not new, and is under-cited.**

- *The trajectory decomposition itself.* Additive
  disease-progression + placebo-response + drug-action modelling of
  longitudinal trial data is an established pharmacometric tradition
  (the paper cites Holford, Gomeni, Chen but treats the precedent as
  background). Recent disease-course-plus-placebo models in depression
  (NLMMRM), ALS, and Huntington's disease decompose trajectories into
  exactly disease-progression and placebo-response terms. The paper
  should engage this literature head-on and state precisely what its
  decomposition adds (the answer -- participant-level, design-identified
  PB -- is defensible, but must be argued against this body of work, not
  around it).
- *The within-subject belief contrast.* The "open-hidden" /
  balanced-placebo design that crosses actual treatment with
  *told*/believed treatment is the classical apparatus for separating
  drug from expectancy, and there is a current clinical literature using
  open-hidden and blinded-discontinuation designs for exactly this
  purpose (e.g. antidepressant-discontinuation open-hidden trials and
  open-label-placebo N-of-1 series). This is the empirical realisation
  of the manuscript's identifying contrast and is essentially absent
  from the references. Its omission weakens both the novelty positioning
  and the credibility of the identifiability argument.
- *The biomarker conflation insight.* That a lumped/marginal
  treatment-effect-by-covariate slope conflates predictive with
  prognostic (and, here, with placebo) signal is the core concern of the
  predictive-biomarker literature the paper cites only via PATH; it has
  been formalised under potential outcomes and in the
  treatment-by-covariate-interaction-is-descriptive literature. The
  manuscript's biomarker-validation failure is a sharpened special case,
  not a new phenomenon.

**Tier.** The Section 6.1 identity is a one-line consequence of linearity
of projection/covariance under an additive model (the paper says as
much). There is no new estimator, no asymptotic theory, no new
inferential procedure beyond a re-parameterised LME. This is a valuable
*applied methodological framework*, appropriate for Statistics in
Medicine or Clinical Trials, but it does not clear the
theory-contribution bar of JASA/JRSS-B/Biometrika.

## 4. Major comments

**M1. (Substance) The executed simulation is a null on the wrong axis,
and the decisive experiment is missing. [Sections 1.2, 6.1, 7.3]**
Study A varies component *magnitude* (`m_PB`, `m_TV`) with the biomarker
coupled to BR alone, so by the Section 6.1 identity the lumped slope is
unbiased by construction; the manuscript correctly reports mean absolute
bias of 0.040-0.012 within one MCSE of zero and concedes this "corrects
the looser claim of section 1." The consequence is that the paper's
motivating failure mode -- the one that justifies the whole framework --
is asserted via algebra and *contradicted* by the only data presented,
because the experiment that would exhibit it (biomarker-PB covariance
`beta_bm^PB != 0` at fixed `m_PB`) is explicitly deferred. For a
top-tier applied paper this is the central deliverable, not future work.
*Remedy:* run the contaminated-biomarker cell (vary `beta_bm^PB` and
`beta_bm^TV` directly), report bias of the one-component vs
phase-augmented vs full-decomposition analyses, and show the regime in
which decomposition recovers `beta_bm^BR` while the lumped analysis does
not. Until then the abstract and Section 1 should not be framed as if
the failure mode is demonstrated.

**M2. (Substance) The broad-applicability claim conflates DGP-level
existence with design-level identifiability. [Section 9 generalisability
paragraph; Section 4]** The Discussion asserts the decomposition applies
to "any longitudinal design" -- parallel-group RCTs with repeated
measures, crossover, observational cohorts. This is true of the
*components* but false of their *identifiability*. The participant-level
decomposition -- which is what predictive-biomarker validation requires,
since `beta_bm^BR` is a participant-level slope -- needs a within-subject
contrast in which belief varies independently of pharmacology. A
standard parallel-group placebo-controlled trial supplies a
*between-subject* placebo arm, which identifies the *population-mean* PB
but not a participant-level PB, and contains no within-subject
blinded-discontinuation contrast at all. The manuscript therefore does
*not*, as drafted, support application to "all placebo-controlled
longitudinal trials." *Remedy:* replace the generalisation paragraph
with an identifiability taxonomy by design class -- parallel-group (+/-
run-in), crossover, randomized-withdrawal/discontinuation, hybrid N-of-1,
open-hidden -- stating for each which components are identified and
whether recovery is population-level or participant-level. Then state the
broad claim precisely: the framework applies to longitudinal designs that
supply within-subject belief x exposure contrasts, of which standard
parallel-group placebo control is *not* one. A simulation in at least one
non-N-of-1 design (e.g. randomized withdrawal) would substantiate the
generalisation.

**M3. (Correctness/identifiability) PB is identified only up to the
expectancy multiplier eta, and eta is assumed, not estimated. [Sections
3.2, 4]** The manuscript states (correctly) that from a two-phase design
only the product `eta * m_PB` is identified, and that separating eta
requires an open-label-placebo phase; it nonetheless proceeds with
`eta ~= 0.5` in the blinded phase as a working value. This is a genuine
identification limitation, not a nuisance: without a belief manipulation
check (rarely collected in practice), eta is unfalsifiable, and every
PB-versus-BR attribution inherits its assumed value. *Remedy:* foreground
this in Section 4 (not only Section 3.2), give the identified set as a
function of eta, and tie the recommendation to designs/measurements that
pin eta down (open-label placebo phase, or direct belief/expectancy
measurement as in open-hidden designs). The Study A `eta`-mismatch
sensitivity cell (promised in Section 7.2) should be executed and
reported.

**M4. (Literature) Engage the open-hidden / balanced-placebo and
pharmacometric decomposition literatures directly. [Sections 1.3, 3.2]**
See Section 3 above and the reference table below. These are not
peripheral: the open-hidden discontinuation design is the empirical
instantiation of the paper's identifying contrast, and the
disease-progression-plus-placebo pharmacometric models are the closest
methodological precedent. Their near-absence makes the novelty claim look
larger than it is and leaves the identifiability argument without its
natural empirical support.

**M5. (Positioning) State the elementary nature of the identity and the
framework-level contribution honestly. [Sections 1.1, 6.1]** The 6.1
identity is correct and useful pedagogically, but it is a direct
consequence of additivity and linearity of covariance; presenting it as
a derived "result" overstates the theoretical content. *Remedy:* present
it as the immediate corollary it is, and let the contribution rest where
it actually lies -- the design-contrast identification map and the
biomarker-validation framing.

**M6. (Robustness) Trajectory-form sensitivity is deferred to a companion
paper. [Sections 3, 7]** A framework claiming broad applicability cannot
outsource the question of whether its conclusions survive replacing the
modified Gompertz with another saturating form. *Remedy:* summarise the
companion paper's (`07-gompertz-evaluation`) key sensitivity result
in-paper, ideally as one row of the simulation: does the
biomarker-by-BR recovery degrade under a misspecified component
trajectory?

**M7. (Assumption) Additivity is assumed by the identity and the main
analysis, while Section 3.3 documents plausible subadditivity. [Section
3.3]** The whole decomposition, and the linearity that makes 6.1 hold,
rests on BR + PB + TV additivity; Section 3.3 itself argues additivity
may fail (Senn scale-dependence; pharmacological-expectancy
subadditivity). The broad claim should state when additivity is tenable
and what a non-negligible interaction does to `plim beta_D` and the
biomarker slope. *Remedy:* add a sentence to 6.1 giving the
first-order correction term under the gamma-interaction model of 3.3, or
explicitly bound the regime in which additivity is assumed.

## 5. Minor comments

1. **Abstract vs Section 1 framing.** The abstract now reports the null
   honestly and states conditionality; Section 1 still opens with the
   three failure modes in declarative form before the "two points of
   precision" caveat (1.2). Tighten so a reader is not told a story that
   7.3 walks back. [Abstract; 1.2]
2. **Which section is "the empirical core."** Section 7.1 calls the
   Section 7.2 production programme "the empirical core," but the
   realised run is the 250-replicate Study A in 7.3, and 7.2 now
   describes deliverables largely not executed. Reconcile the
   cross-references. [7.1, 7.2, 7.3]
3. **`b_TV` definition.** In 6.1, `b_TV = 0` for a clean on/off contrast
   and strictly positive for a baseline-anchored pre-post analysis; the
   text is correct but the reader would benefit from stating up front
   which contrast the "lumped analysis" of later sections actually uses,
   since the failure-mode claims depend on it. [6.1]
4. **Notation.** Confirm `eta` is introduced once (Section 2 / Table of
   terms) before its use in 3.2 and 4; confirm `D_it` (continuous-decay
   drug state) vs the `1_{on drug}` phase-table shorthand are reconciled
   wherever both appear. [2, 4]
5. **Reproducibility.** Study A reports 250 replicates/cell; the
   pre-registration (Section 7.2) specified 1000 (bias) and 5000 (type
   I). State why the realised run used 250 and whether the MCSEs
   (0.013-0.05) are adequate for the bias claim, or raise the count.
6. **Type I / multiplicity.** The phase-augmented model adds `Dbc:phase`
   and `bm:Dbc:phase`; if both the lumped and phase-augmented analyses
   are run as a sequence, address the multiplicity of the attribution
   test. [5.2]

## 6. Missing and questionable references

| Claim / location | Issue | Suggested source |
|---|---|---|
| Identifying contrast = silent/blinded discontinuation (3.2, 4) | The empirical design realising this contrast is not cited | Open-hidden / hidden-vs-open antidepressant-discontinuation trials (e.g. open-hidden discontinuation protocols, 2023); open-label-placebo N-of-1 discontinuation series (2025) |
| eta / belief separated from pharmacology (3.2) | Classical apparatus uncited | Balanced placebo design (Ross et al.; Rohsenow & Marlatt) crossing administered x told treatment |
| "No published trial decomposes the response" (1.3) | Overstated vs pharmacometrics | Disease-progression + placebo-response models in depression (NLMMRM), ALS, and Huntington's disease-course models |
| Lumped biomarker slope conflates predictive/prognostic (1.2, 8) | Established under potential outcomes | Predictive-biomarker evaluation under potential outcomes / random effects (2015); treatment-by-covariate interaction as descriptive association |
| Additive disease-progression/placebo/drug tradition (1.1) | Cited as background only | Holford disease-progress framework + a modern placebo-response-model exemplar, engaged rather than listed |

(Full bibliographic completeness of `references.bib` was not audited
entry-by-entry; see scope.)

## 7. Suggestions for strengthening the paper

- **Run the decisive cell (highest priority).** A 3 x 3 grid over
  `beta_bm^PB` and `beta_bm^TV` at fixed component magnitudes, reporting
  bias of `hat-beta_bm^lumped` against `beta_bm^BR`, would convert the
  paper from "framework + null" to "framework + demonstrated failure mode
  + recovery." This single addition addresses M1 and most of the
  significance gap.
- **Add a non-N-of-1 design to the simulation** (randomized withdrawal
  or open-hidden parallel-group) to substantiate or appropriately bound
  the generalisation claim of M2.
- **Report the eta-sensitivity and trajectory-misspecification cells**
  already promised in 7.2; these are the robustness results a referee of
  a "broadly applicable" framework will require.
- **Replace the generalisation paragraph with the identifiability
  taxonomy** (M2); it is more useful and more defensible than the current
  blanket claim and directly answers the "all placebo-controlled trials"
  question the framework invites.
- **Consider a partial-identification framing for PB** (identified set as
  a function of eta), which is both more honest and more modern than a
  point assumption.

## 8. Scope of this review

- **Verified (re-derived):** the Section 6.1 probability-limit identity
  `plim hat-beta_D = m_BR + b_PB + b_TV` and the biomarker-slope
  decomposition `beta_bm^lumped = beta_bm^BR + w_PB beta_bm^PB + w_TV
  beta_bm^TV`. Both are correct and follow from additivity and linearity
  of covariance/projection; the cancellation of TV under a clean
  within-subject on/off contrast is correctly stated.
- **Inspected (read and judged, not executed):** the simulation
  specification and the Study A numerical results (bias 0.040-0.012;
  power 0.35 vs 0.59 at N=35, 0.62 vs 0.91 at N=70; CI coverage
  0.95-0.99; type I 0.020-0.060), the identifiability/phase-contrast
  argument of Section 4, the worked-example numbers of Section 6.2, and
  the component definitions of Section 3. I did not execute the
  simulation code or recompute these figures.
- **Not checked:** entry-by-entry bibliography completeness; the companion
  documents (`docs/24-...pedagogy`, `07-gompertz-evaluation`) beyond their
  role as cited support; numerical re-execution of any chunk.
- **Literature searches run:** (1) decomposition of treatment response
  into drug/placebo/natural-history in longitudinal mixed models; (2)
  predictive-vs-prognostic biomarker confounding and treatment-effect
  heterogeneity; (3) within-subject placebo identification via
  blinded/open-hidden discontinuation and expectancy. These surfaced
  overlapping pharmacometric decomposition work, the
  potential-outcomes predictive-biomarker literature, and current
  open-hidden / open-label-placebo discontinuation designs, all
  discussed above.

---

*Reviewer's note on the framing question.* On the specific question of
whether the manuscript supports application to *all* placebo-controlled
longitudinal trials: as written, it does not. It supports application to
longitudinal designs that supply a within-subject contrast in which
belief varies independently of drug exposure. Establishing the broader
claim requires the design-class identifiability taxonomy of M2 and at
least one demonstration outside the hybrid N-of-1 design. The components
are universal; their identification is not.
