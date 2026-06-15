# Referee Report (Second Round): Three-Component (BR-PB-TV) Decomposition
*2026-06-15 09:55 PDT*

Manuscript: `analysis/report/06-component-decomposition/report.Rmd`
(pmsimstats-ng). This is a second-round review of the revised draft;
the first-round report is `referee-report-2026-06-14.md`. Framing
retained: the authors position the work as broadly applicable to
placebo-controlled longitudinal trials.

---

## 1. Summary of the revision

The authors have responded substantively. Three changes dominate. (i)
The headline simulation (Study A) was re-run at the pre-registered
1,000 replicates per alternative cell, tightening the null on the
component-magnitude axis (mean absolute bias 0.016 at N=35 to 0.006 at
N=150). (ii) A new "Contaminated-biomarker pilot" subsection runs the
decisive experiment the first round demanded -- it varies the
biomarker-PB correlation c_bm,PB in {0, 0.15, 0.30, 0.45} and exhibits
the biomarker-validation failure (one-component bias rising to +0.51 at
c_bm,PB = 0.45, N=150, ~30 Monte Carlo SEs from zero). (iii) The
Discussion's broad-applicability paragraph was replaced by a
population-level vs participant-level identifiability taxonomy, and the
Section 6.1 identity, the eta-identification limitation, and the
additivity assumption were all sharpened. The result is a markedly
stronger and more honest paper.

## 2. Overall assessment and recommendation

**Recommendation: Major revision** -- but on a short and clear path to
acceptance at a journal such as Statistics in Medicine or Clinical
Trials. The two first-round blockers are essentially resolved, and what
remains is one genuine scientific gap (the framework's own remedy is not
yet demonstrated) plus completion and reproducibility work. This is no
longer a paper whose central claim is unsupported; it is a paper whose
diagnosis is now demonstrated and whose cure is still pending.

The single most important remaining point: the pilot shows that the
lumped analysis **is** biased under biomarker-PB coupling and that the
cheap remedy (the phase-augmented analysis) **provides no protection**;
but the full component decomposition that the paper advocates -- the one
that should recover the unbiased biomarker-by-BR slope -- is labelled
"Study B ... pending implementation." A framework paper must show that
its framework solves the problem it poses. As it stands the manuscript
demonstrates the disease and shows that the inexpensive treatment fails,
without yet showing that the expensive treatment it recommends works.

## 3. Status of the first-round major comments

- **M1 (central tension: failure mode unrun) -- substantially resolved.**
  The contaminated pilot demonstrates the failure mode unambiguously at
  N=150. Moreover the authors found, correctly, that the bias is driven
  by c_bm,PB x sigma_PB (the covariance of the biomarker with the
  *individual-level* PB amplitude) and not by the population mean m_PB --
  cells with m_PB = 0 and m_PB = 6 show identical bias. I verified this
  against the Section 6.1 identity: beta_bm^PB = Cov(m_PB,i, B)/Var(B)
  depends on the participant-level covariance, into which the population
  mean does not enter, so the mean shifts the intercept and the
  covariance shifts the slope. This is the right explanation, it
  retrospectively explains why Study A was a null on the wrong axis, and
  it is a genuinely useful refinement: the operative risk to a biomarker
  programme is correlation with PB *variance*, not the trial's average
  placebo response. Residual work is completion (below), not concept.

- **M2 (broad-applicability overclaim) -- resolved, and the broad claim
  is now correctly disavowed.** The taxonomy distinguishes population-
  from participant-level recovery and states that a parallel-group
  placebo-controlled trial identifies population-mean PB but not
  participant-level PB, so the biomarker-by-BR slope is not identified
  from that design. This is the correct position. See Section 4 below
  for the answer to the framing question.

- **M3-M7 and minors -- addressed.** The eta identified-set framing, the
  additivity caveat with its first-order interaction correction, the
  honest statement that the 6.1 identity is an elementary corollary, the
  imported Gompertz-robustness result from paper 07, the multiplicity
  note, and the "empirical core" cross-reference were all incorporated.

## 4. Is the broad claim now supportable?

**No, and the revised manuscript now says so itself.** This must be
stated plainly because the authors' framing ("applicable to all
placebo-controlled longitudinal trials") is in direct tension with their
own revised Discussion. The contaminated pilot sharpens exactly why: the
contamination enters through the *participant-level* covariance
Cov(m_PB,i, B), which requires participant-level PB amplitudes to be
identified, which requires a within-subject belief-by-exposure contrast.
A standard parallel-group placebo-controlled longitudinal trial does not
supply that contrast; it yields population-mean PB only. The genuinely
important and now-evidenced contribution is narrower and conditional:
*in designs that supply a within-subject belief contrast*, lumped
biomarker validation is biased precisely when the biomarker correlates
with individual-level PB or TV variance, and decomposition is the
mechanism that removes that bias. That is a valuable result. It is not
"all placebo-controlled trials," and the paper should keep resisting the
broader phrasing wherever it survives (e.g. any abstract or introduction
sentence that still reads as universal).

## 5. Significance and novelty (updated)

The increment over `hendrickson2020` (the program's foundational
"Optimizing aggregated N-of-1 designs for predictive biomarker
validation," already cited) is now clearer and better evidenced: the
explicit three-component decomposition, the omitted-variable-bias
identity for the lumped biomarker slope, and -- new in this round -- the
c_bm,PB x sigma_PB characterisation of when validation fails. The
contribution remains a methodological framework with a clean
identity-plus-simulation argument, not new estimation theory; the
appropriate venue is a strong applied-methods journal. The
sigma_PB-not-m_PB finding is the kind of crisp, counterintuitive,
practically actionable result that lifts the paper above a routine
"decompose your model" message.

## 6. Major comments (this round)

**R1. (Substance, highest priority) Demonstrate that the decomposition
recovers the unbiased slope. [Contaminated-biomarker pilot; "Study B"]**
The pilot establishes the problem and rules out the cheap fix. The paper
now needs Study B: under the same c_bm,PB contamination, fit the full
component decomposition and show that the recovered beta_bm^BR is
unbiased while the lumped and phase-augmented estimates are not. Without
this, the paper's central prescriptive claim -- decompose -- is
unsupported by its own evidence, and a sceptical reader can conclude only
that lumped validation is biased and that the paper's own tractable
alternative does not help. This is the one change I would insist on
before acceptance.

**R2. (Evidence strength) Promote the pilot to the production run.
[Contaminated-biomarker pilot, lines ~2905-2969]** The decisive cell is
100 replicates; the authors correctly flag that N=70 is uninformative
(MCSE ~0.06, comparable to the bias) and that type I error under
contamination is unrun. The N=150 result is solid at 30 MCSE, but the
headline of a top-tier paper should not rest on a single informative
sample size from a pilot. Run the pre-registered 1,000/5,000-replicate
grid, report bias across all N, and report type I error under
contamination.

**R3. (Reproducibility) Wire the reported numbers to the loaded data.
[setup chunk `load-sim-data`; pilot prose]** The chunk loads
`sim_contam` (study-contam-alt100-null100) but the pilot's numerical
claims (-2.31, -2.37, -2.08, -1.74, +0.51, the <=0.006 m_PB invariance)
are hardcoded in prose rather than computed inline from that object. The
same applies to the Study A figures. This is a live drift hazard: the
re-run already changed the Study A numbers once this revision cycle.
Compute the reported quantities with inline R from the summary objects so
the prose cannot fall out of sync with the data.

**R4. (Exposition) Make the abstract and Section 1 carry the
magnitude-vs-correlation distinction the pilot established.** The
abstract now mentions contamination, which is good, but the paper's
motivating language in places still risks conflating "large PB/TV" with
"biomarker coupled to PB/TV." The pilot's central lesson -- it is the
correlation with PB *variance*, not the placebo magnitude -- should be
stated crisply once in the abstract and once in Section 1, and every
residual "magnitude"-flavoured statement of the failure mode should be
corrected to the covariance axis.

## 7. Minor comments (this round)

1. **Bias sign interpretation.** The pilot explains the positive
   (toward-zero) bias by partial absorption into the bm main effect; this
   is plausible but is asserted, not shown. One sentence deriving the
   sign from the identity (the bm main-effect column competing with the
   bm:Dbc column for the same covariance) would close it. [lines ~2929-2934]
2. **"Study A ... at 1,000 replicates" vs pilot at 100.** With two
   simulation studies at different resolutions now in the paper, add a
   one-line table or sentence stating replicate counts, cells, and MCSE
   for each, so a reader is not left reconciling them across subsections.
3. **`sigma_PB = 2` provenance.** The individual-level PB SD that drives
   the contamination channel should be tied explicitly to the
   prazosin-PTSD calibration, since the whole effect size depends on it.
4. The `rgt` blocks remain one-sentence Lorem placeholders; the paper
   cannot be submitted until the author prose replaces them. (Noted as a
   production state, not a scientific defect.)

## 8. Missing and questionable references

No new missing references in this round. The closest prior work
(`hendrickson2020`) is cited and is correctly treated as the foundation
the paper extends. The first-round additions (`colloca2004`,
`rohsenow1981`) are present and appropriately placed. One optional
addition: when Study B is added, a sentence relating the
participant-level component-recovery identifiability to the individual-
treatment-effect cautions already cited (`senn2018`) would strengthen R1.

## 9. Suggestions for strengthening

- Run Study B (R1) and the production contaminated grid (R2) together;
  they are the remaining empirical core and convert a "major revision"
  into an acceptance.
- Consider reporting the bias as a function of c_bm,PB on a single
  figure (the paper currently has no figures); the monotone
  bias-vs-coupling curve is the paper's signature result and deserves a
  panel, alongside a recovery curve for the decomposition once Study B
  exists.
- State the practical rule the pilot licenses: a biomarker programme
  should screen candidate biomarkers for correlation with measured
  expectancy/placebo-responsiveness, because it is that correlation, not
  the trial's mean placebo response, that invalidates lumped validation.

## 10. Scope of this review

- **Verified (re-derived):** that the contamination bias scales with the
  participant-level covariance c_bm,PB x sigma_PB and not with the
  population mean m_PB, from the Section 6.1 identity; the consistency of
  the abstract's updated Study A figures with the body.
- **Inspected (read, not executed):** the contaminated-pilot numerical
  values (-2.31 to -1.74; +0.51 at 30 MCSE; m_PB invariance <= 0.006),
  which are reported in prose and which I did not recompute; the
  existence of the backing data files
  (`study-contam-alt100-null100/`, and an as-yet-unreferenced
  `study-contam-alt1000-null5000/` directory) on disk; the revised
  taxonomy, eta, and additivity passages.
- **Not checked:** numerical re-execution of either simulation; the
  contents of the `.rds` summaries; Study B (does not yet exist).
- **Literature searches run this round:** recent (2024-2025) N-of-1 /
  placebo-decomposition / biomarker work. The most on-topic hit was the
  program's own `hendrickson2020`; no new external work was found that
  threatens novelty or that the manuscript fails to engage.

---

*Reviewer's bottom line.* This revision did the hard thing: it ran the
experiment that could have falsified the paper's motivation, and the
experiment supported it -- while also correcting the paper's own account
of *why* (variance-coupling, not magnitude). The paper is now genuinely
promising. Two finite, well-defined runs (the decomposition-recovers-the-
slope study, and the production-resolution contaminated grid), plus
wiring the numbers to the data and holding the scope to what the taxonomy
licenses, would make it a clean accept at a strong applied-methods
journal. The "all placebo-controlled trials" framing should be retired in
favour of the precise, and now evidenced, conditional claim.
