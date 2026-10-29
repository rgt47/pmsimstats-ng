# Referee report: Three-component decomposition of treatment response
*2026-10-28 17:41 PDT*

**Manuscript.** `analysis/report/06-component-decomposition/report.Rmd`
(pmsimstats-ng), version dated 2026-10-08. Reviewed against the
standard of *Statistics in Medicine*. Prior reports in this directory
(2026-06-12, 06-13, 06-14, 06-15) were read; their status is tracked in
Section 4 and Section 7. `report-slim.Rmd` and `bullets.Rmd` exist in
the directory and were not reviewed.

**Evidence labels.** VERIFIED means I executed code against the stored
artifacts or re-derived the mathematics. INSPECTED means I read the
source or data and judged it without execution. SUSPECTED means the
evidence is indirect. UNVERIFIED means I did not check.

---

## 1. Summary

The manuscript proposes that treatment response in aggregated N-of-1
trials be written as the sum of a pharmacological component (BR), a
placebo-belief component (PB) and a natural-history component (TV),
each a modified Gompertz trajectory with participant-specific random
effects. It maps the identifiability of the components onto the phase
contrasts of the Hendrickson hybrid design, derives the omitted-
variable identity for the one-component treatment effect and for the
one-component biomarker slope, and reports three simulation studies.
Study A (1,000 replicates per alternative cell, 5,000 per null cell)
varies the magnitudes of PB and TV with the biomarker correlated with
BR alone and finds the one-component biomarker-treatment estimate
unbiased, with a phase-augmented analysis that is no less biased and
markedly less powerful. A contaminated-biomarker pilot (100 replicates
per cell) correlates the biomarker with PB and reports a one-component
bias that reaches $+0.51$ at $c_{bm,PB} = 0.45$ and $N = 300$, some 30
Monte Carlo standard errors from zero. Study B (1,000 replicates per
cell) shows that a belief-covariate decomposition recovers the
pharmacological slope on a balanced-placebo design but, the authors
state, not on the hybrid design. The paper concludes that
decomposition is useful only when the estimand is attribution or when
the biomarker may correlate with PB or TV, and that recovery is
conditional on a design that decouples drug exposure from belief.

## 2. Overall assessment and recommendation

**Recommendation: major revision.** The paper is substantially
stronger than the 2026-06-12 draft. The identity is correct, Study A is
clean, Study B on the balanced-placebo design is clean, and the
conditional framing of the Discussion is the right one. The
manuscript's own 2026-10-08 correction notes show a commendable
willingness to record what did not reproduce.

Three problems nonetheless prevent acceptance.

First, and most serious, the headline contamination result is an
artifact of silently repaired covariance matrices. I rebuilt the
covariance matrix for every cell of the contaminated-biomarker pilot
and of the hybrid-design recovery run. On the sixteen-visit two-path
design used for the $N = 300$ cells, the matrix is not positive
definite at every $c_{bm,PB} \ge 0.15$, and on the eight-visit design
it is not positive definite at $c_{bm,PB} = 0.45$. The June runs
passed `makePositiveDefinite = TRUE` through a code path that repaired
such matrices without warning, and the repair lowers the on-drug
biomarker-BR correlation from $0.45$ to between $0.44$ and $0.36$. The
true pharmacological slope in those cells is therefore not $-2.25$,
and the reported "bias" is in large part the gap between the nominal
and the repaired data-generating process. In the cells whose matrix is
positive definite the one-component bias is negative, not positive,
and its size agrees with what the Section 6.1 identity predicts. The
abstract, Section 7.4, the 2026-10-08 correction note on the sign, and
the hybrid-design half of Section 7.5 all need to be rebuilt on reruns
under the project's current refusal of non-PD matrices. (Major
comment 1 and 2.)

Second, several commitments the manuscript makes to its own reader are
not met in the stored code or data: the full nonlinear decomposition
is a stub that was never implemented, the convergence rate it quotes
for that fit is asserted without any run, the subadditivity study, the
$\eta$-mismatch grid and the dropout grid were never executed, and the
Reproducibility section describes seeds and paths that do not match
the drivers. (Major comments 3, 4 and 5.)

Third, the exposition has not caught up with the science. The
Discussion still says the simulation is forthcoming, the abstract
carries editorial correction notes in the body text, the cover letter
and README describe a different paper, the expectancy value assigned
to an "open-label placebo" phase differs in three places, and the
analysis that Study B recommends is the series' own F1 model without
being named as such. (Major comments 6 to 8.)

None of this is fatal. The identity-plus-design argument is sound, and
a corrected contamination study will very likely still show a
design-dependent bias whose sign and size the identity predicts. That
is a cleaner and more useful result than the present one.

## 3. Significance and novelty

The contribution is a methodological framework for the aggregated
N-of-1 setting, not new estimation theory. The prior reports have
already established, and the authors now acknowledge, that the
structural decomposition is standard pharmacometrics (Holford, Gomeni,
Chen) and that the Section 6.1 identity is an immediate corollary of
additivity. What survives as a contribution is (a) the design-contrast
identification map, (b) the biomarker-validation framing in terms of
the biomarker-by-BR slope, and (c) the result, first found in this
paper's revision cycle, that the one-component slope is displaced by
the biomarker's covariance with the participant-level PB amplitude and
not by the population-mean placebo response. Point (c) is the paper's
most useful finding and survives the artifact identified below, since
it follows from the identity and is confirmed in the positive-definite
cells. Point (d), that recovery requires a design that decouples
exposure from belief, is also likely to survive, but its hybrid-design
half rests on repaired matrices and must be re-established.

**Contemporariness.** The reference list runs to 2025 (Carrozzo) and
engages the carryover literature (Liao 2023, Daza 2018, Gaertner
2023). Two gaps remain from earlier rounds. Petkova, Tarpey and
Govindarajulu (2009, *Int J Biostat*, PMID 21556319), which
disentangles specific from non-specific response in drug-treated
subjects with a linear mixed model, and the companion Tarpey et al.
(2010, *JASA*, PMID 21494314), which applies it to an open-label
fluoxetine phase followed by blinded discontinuation, are the closest
statistical precedent for Sections 5 and 7.5 and are still not cited
(flagged 2026-06-13). The open-hidden antidepressant discontinuation
trial protocol of Meissner et al. (2023, *BMC Psychiatry*, PMID
37344789; NCT05191277) is the empirical realization of the Study B
design in a psychiatric population and is still not cited (flagged
2026-06-12 and 06-14). A balanced-placebo esketamine protocol (2023,
*Trials*, PMID 37568215) and a 2026 crossover balanced-placebo N-of-1
provocation protocol (PMID 42827697) show that the design the paper
recommends is being fielded, which strengthens rather than threatens
the paper. I found no published work that performs the three-component
decomposition at the N-of-1 level by design contrasts; the novelty
claim as narrowed in the current draft holds (absence of evidence,
searches listed in Section 8).

## 4. Major comments

Numbered in order of importance, correctness first.

### M1. The contamination cells were drawn from repaired non-PD matrices

The cells at $c_{bm,PB} \ge 0.15$ on the sixteen-visit design, and at
$0.45$ on both designs, were generated from covariance matrices that
were not positive definite and were repaired without warning.

**Location.** Abstract (Results paragraph); Section 7.4 throughout;
Section 7.5, the sentence "we verified that the belief-covariate
decomposition, a blinded-stratum drug contrast, and a Gompertz-basis
component decomposition all fail to recover $\beta_{bm}^{BR}$ on that
design"; Conclusions, third point.

**Problem.** I rebuilt the covariance matrix with `buildSigma()` at
the driver parameters ($c_{bm} = 0.45$, $\sigma_{BR} = 5$,
$\sigma_{PB} = 2$, $\rho = 0.7$, $c_1 = 0.1$, $c_\times = 0.05$,
$t_{1/2} = 1$) for each design and contamination level, with
`stopOnNonPD = FALSE` and `makePositiveDefinite = FALSE` so the raw
matrix is returned. Minimum eigenvalue of the positive-variance block
(VERIFIED, `scratchpad/pdcheck.R`):

| Design | $c_{bm,PB}$: 0 | 0.10 | 0.15 | 0.20 | 0.30 | 0.45 |
|---|---|---|---|---|---|---|
| eight-visit single path ($N = 70$) | 0.344 | 0.321 | 0.277 | 0.215 | 0.044 | $-0.307$ |
| sixteen-visit, path AB ($N = 300$) | 0.045 | 0.020 | $-0.023$ | $-0.084$ | $-0.260$ | $-0.635$ |
| sixteen-visit, path BA ($N = 300$) | 0.054 | 0.031 | $-0.011$ | $-0.072$ | $-0.247$ | $-0.624$ |
| balanced placebo (Study B) | 0.321 | 0.315 | 0.271 | 0.208 | 0.032 | $-0.348$ |

The drivers `run-study-contam.R` and `run-study-b-recovery.R` call
`generateData(..., makePositiveDefinite = TRUE)`. The package at the
time of the June runs (git history: `stopOnNonPD` was added 2026-10-21
and made the default 2026-10-26) repaired a non-PD matrix with
`corpcor::make.positive.definite(sigma, tol = 1e-3)` and emitted no
warning unless `verbose = TRUE`. Applying that repair to the raw
matrices above changes the on-drug biomarker-BR correlation, which is
the quantity that fixes the true slope (VERIFIED):

| Design | $c_{bm,PB}$ | repaired $\mathrm{Cor}(B, BR_{\text{on}})$ | implied true $\beta_{bm:D}$ | largest absolute change in any correlation |
|---|---|---|---|---|
| sixteen-visit AB | 0.15 | 0.445 | $-2.22$ | 0.006 |
| sixteen-visit AB | 0.20 | 0.433 | $-2.16$ | 0.019 |
| sixteen-visit AB | 0.30 | 0.404 | $-2.02$ | 0.062 |
| sixteen-visit AB | 0.45 | 0.359 | $-1.79$ | 0.169 |
| eight-visit | 0.45 | 0.399 | $-1.99$ | 0.112 |
| balanced | 0.45 | not run (excluded by the authors) | | |

The stored one-component means at $N = 300$ are $-2.37$, $-2.08$ and
$-1.74$ at $c_{bm,PB} = 0.15$, $0.30$ and $0.45$ (VERIFIED from
`study-contam-summary.rds`). Measured against the repaired truth
rather than $-2.25$, the biases are about $-0.15$, $-0.06$ and
$+0.05$. The "30 MCSE" result is therefore the distance between the
nominal and the repaired data-generating process, not a property of
the one-component estimator. The same applies to the hybrid-design
recovery run (`study-b-recovery-pilot/`, 1,000 replicates per cell),
whose cells at $c_{bm,PB} \ge 0.15$ use the same sixteen-visit design
and the same repair; its decomposition "failure" of $+0.82$ at $0.45$
is measured against a slope that the repaired DGP no longer has.

Study A and the balanced-placebo Study B are not affected: every
Study A cell is positive definite (VERIFIED over the full
$m_{PB} \times m_{TV}$ grid on all three paths; the smallest minimum
eigenvalue is $0.045$ on the sixteen-visit design), and Study B's grid
stops at $c_{bm,PB} = 0.30$ where the balanced matrix is still
positive definite, which is the feasibility bound the authors
correctly computed for that design. The balanced design does pass
through the old repair path because the covert visits have exactly
zero PB variance ($\eta = 0$), but the repair there only nudges a zero
variance to $10^{-3}$ and leaves the correlations unchanged (VERIFIED:
largest correlation change $0$ at $c_{bm,PB} \le 0.30$).

**Why it matters.** The abstract's quantitative headline, the sign
narrative of Section 7.4, the 2026-10-08 note that the contamination
bias "is toward zero", and the hybrid-design half of Section 7.5 all
rest on these cells. This is the same failure mode the project
documented for paper 02 (`docs/65`): a silently repaired matrix
reports a DGP attenuation as an analysis property. The manuscript's
own note ("the sign of the contamination bias ... depends on the
covariance construction, including any positive-definiteness repair")
identified the hazard without checking it.

**Remedy.** Rerun the contaminated-biomarker study and the hybrid
recovery study under the current default (`stopOnNonPD = TRUE`),
restricting each design to its feasible range: $c_{bm,PB} \le 0.30$ on
the eight-visit design and $c_{bm,PB} \le 0.10$ on the sixteen-visit
design, or lower $c_{bm}$ or the cross-component correlations until
the range widens. Report the feasibility bound for each design beside
its results, as Section 7.5 already does for the balanced design. Add
a sentence to Section 7 stating that no reported cell drew from a
repaired matrix. Until then, remove the $+0.51$ figure, the "30 MCSE"
sentence and the "scales monotonically" sentence from the abstract and
Section 7.4.

**Evidence.** VERIFIED (eigenvalues, repaired correlations, stored
means). The attribution of the residual bias to the repair, rather
than to a genuine estimator property, is INFERRED from the agreement
of the positive-definite cells with the identity (M2).

### M2. The identity predicts the sign and size of the contamination bias

The sign and size follow from the Section 6.1 identity and the
design's expectancy schedule; the positive-definite cells agree with
it, and they are not "uninformative".

**Location.** Section 7.4, paragraph beginning "At $N = 300$ the
one-component bias scales monotonically" and the 2026-10-08 note that
follows; paragraph beginning "The $N = 70$ results are uninformative
at 100 replicates".

**Problem.** Under the covariance-moderation DGP, $\mathrm{Cov}(B,
PB_{it}) = c_{bm,PB}\,\sigma_{PB}\,\eta_t$ at every visit, while
$\mathrm{Cov}(B, BR_{it}) = c_{bm}\,\sigma_{BR}$ on drug and decays off
drug. The one-component `bm:Dbc` coefficient is an on-minus-off
within-participant contrast of the biomarker slope, so the PB channel
displaces it by approximately
$$-\,c_{bm,PB}\,\sigma_{PB}\,(\bar\eta_{\text{on}} - \bar\eta_{\text{off}}),$$
where the bars are the design's mean expectancy over on-drug and
off-drug visits (sign negative because the components are reductions
in $Y$). This is the $w_{PB}\beta_{bm}^{PB}$ term of the identity with
the loading written out. On the eight-visit design
($\bar\eta_{\text{on}} \approx 0.92$, $\bar\eta_{\text{off}} = 0.5$)
the predicted displacement is about $-0.13$ at $c_{bm,PB} = 0.15$ and
$-0.25$ at $0.30$; the stored one-component biases in those
positive-definite cells are $-0.14$ and $-0.12$ at $0.15$ and $-0.37$
and $-0.43$ at $0.30$ (two $m_{PB}$ levels, MCSE about $0.06$;
VERIFIED). On the balanced design ($\bar\eta_{\text{on}} = 1/3$,
$\bar\eta_{\text{off}} = 1$) the predicted displacement is
$+0.13$, $+0.27$ and $+0.40$ at $c_{bm,PB} = 0.1$, $0.2$, $0.3$; the
stored one-component biases at $N = 150$ are $+0.16$, $+0.34$ and
$+0.48$ (VERIFIED). The sign is a design property: it is negative when
belief is higher on drug than off, positive when the reverse holds.
The manuscript's mechanism ("partially absorbed by the `bm` main
effect") is not needed and, as the 2026-10-08 note admits, was not
reproduced.

Two consequences. First, the $N = 70$ cells of the contamination
pilot are not uninformative: a bias of $-0.37$ at an MCSE of $0.057$
is six standard errors from zero and has the opposite sign to the
$N = 300$ cells. The manuscript labels the only cells that were
generated faithfully as uninformative and builds its narrative on the
cells that were not. Second, the sentence "At $c_{bm,PB} = 0.45$ ...
the one-component bias is approximately $0.51$ points ... representing
$23\%$ of the true pharmacological effect" should be withdrawn.

**Remedy.** Derive the displacement in Section 6.1 as a one-line
corollary, with the loading $w_{PB} = \bar\eta_{\text{on}} -
\bar\eta_{\text{off}}$ made explicit, and use it to predict each
design's bias before the simulation reports it. Replace the "toward
zero" narrative with the design-dependent sign. Report the $N = 70$
cells as the informative ones they are once the rerun of M1 is in
hand.

**Evidence.** VERIFIED (arithmetic and stored summaries). The
approximation ignores the carryover decay in $D_{bc}$ and the AR(1)
weighting of the mixed-model contrast, which is why the eight-visit
agreement is to within about one MCSE rather than exact.

### M3. The hybrid-design recovery claim is unsupported as reported

It is reported without numbers, cites an analysis that is not in the
stored data, and rests on the repaired cells of M1.

**Location.** Section 7.5, first paragraph; Conclusions.

**Problem.** The sentence "we verified that the belief-covariate
decomposition, a blinded-stratum drug contrast, and a Gompertz-basis
component decomposition all fail" is supported by
`study-b-recovery-pilot/` (eight cells, 1,000 replicates, four
analyses: `one_component`, `phase_augmented`, `decomposition`,
`decomposition_blinded`). No Gompertz-basis analysis exists in that
driver or its output (INSPECTED, `run-study-b-recovery.R` and the
summary). The run is named "pilot" but is at production resolution
and is the sole evidence for the paper's claim that the standard
design cannot recover the slope, which is half of the Conclusions'
"decisive practical lesson". Its numbers are not given, its directory
is not named in the manuscript, and its contaminated cells are the
repaired ones of M1. The uncontaminated and $c_{bm,PB} = 0.10$ cells,
had they been run, would be positive definite on this design.

**Remedy.** Rerun on the feasible range, report the run as a numbered
study with a table (bias and MCSE by analysis and $c_{bm,PB}$), name
its data directory, and either implement the Gompertz-basis analysis
or delete the claim. If the decomposition does recover the slope at
small $c_{bm,PB}$ on the hybrid design, as the near-zero stored bias
at $0.15$ ($+0.05$, which is within the repaired-truth gap) suggests
it may, the conclusion softens from "no decomposition recovers it" to
a quantitative statement about how fast recovery degrades with the
drug-belief correlation.

**Evidence.** VERIFIED (contents of the stored run); INSPECTED
(driver).

### M4. The full nonlinear decomposition was never implemented

The manuscript nonetheless quotes its convergence rate and promises
its diagnostics.

**Location.** Section 5.1, third reason ("At trial-relevant $N$, the
convergence rate of the component-matched fit is typically below
90\%"); Section 4, final paragraph ("The identifiability diagnostics
of section 7 (convergence rates of the component-matched fit,
variance-inflation factors ...)"); Section 7.2 (deliverables);
Section 3.3 ("A simulation study quantifying ... non-zero $\gamma$ is
forthcoming"); ADEMP protocol Study B and Study C.

**Problem.** `full_nonlinear_fit()` in `00-common.R` is a stub that
returns `NA` with `formula_dropped = 'nlme_stub'` (VERIFIED). No
driver fits a Gompertz component model. The "below 90\%" convergence
figure has no source. The protocol's Study B (identifiability of the
participant-level components) was replaced by a slope-recovery study
with the same letter, which the manuscript concedes in one clause
("the participant-level component recovery that motivated the original
identifiability aims of Study B is not addressed"). Study C
(`04-study-C.R`) has no output directory. The $\eta$-mismatch, TV-
magnitude and dropout grids of Section 7.2 were not run. The
production contaminated run (`study-contam-alt1000-null5000/`)
contains one cell.

**Why it matters.** A reader of Section 5 is told, as a fact, that
the component-matched fit converges poorly, and is told in Sections 4
and 7.2 that diagnostics mapping the recoverable region will be
reported. Neither is true of this manuscript. The Discussion's second
open question ("a sensitivity analysis varying $\eta$ ... is
mandatory") has been mandatory since the first draft.

**Remedy.** Delete the convergence-rate claim or replace it with a
citation or a run. Rewrite Sections 4 and 7.2 so that every promised
deliverable is either reported or explicitly listed as not done in
this paper. Rename the recovery study so that it does not share a
letter with the protocol's unexecuted identifiability study (NOTATION
Part 3 asks for numbered studies in any case). Either run the
$\eta$-sensitivity cell, which is cheap, or move the "mandatory"
language to a Limitations paragraph.

**Evidence.** VERIFIED (stub code, absent output directories).

### M5. The Reproducibility section does not describe the drivers

**Location.** Section "Reproducibility".

**Problem.** The section states that each driver sets
`set.seed(20260507)` and uses `parallel::mc.set.seed`, with
per-replicate seeds derived by `+ rep_idx`, and that per-replicate
outputs are at `analysis/data/quick-sim/06-decomp-replicates.rds`.
The drivers use study seeds 642 (Study A), 644 (contamination), 916
(Study B balanced) and 815 (hybrid recovery), seed each replicate with
`set.seed(study_seed + 1000L * cell_id + rep_idx)`, and do not call
`mc.set.seed` (VERIFIED). The quick-sim path holds only the six-cell
pilot; the three studies the paper reports live under
`analysis/data/derived_data/component-decomposition/`, which the
section does not name. The section also states the manuscript was
rendered inside the Docker container while the driver headers say
"host R, not Docker" (INSPECTED; which environment produced the stored
runs is UNVERIFIED). The stored `N` column was recoded in place on
2026-08-08 (`migrate-n-to-total.R`) without a rerun; that migration
should be disclosed where the data are described.

**Remedy.** Rewrite the section from the drivers: one row per study
giving driver file, seed, replicate count, output directory and the
date of the run. State the `N` migration. Add the one-line statement
requested in M1 that no reported cell drew from a repaired matrix.

**Evidence.** VERIFIED.

### M6. Stale and self-contradicting prose

**Location.** Discussion, second paragraph; Abstract; Section 1.1;
cover letter; README.

**Problem.**

- The Discussion says "the simulation study described in section 7 is
  forthcoming; the present manuscript reports the framework rather
  than its empirical evaluation, and the quantitative bias-versus-power
  trade-offs ... are conjectures until that work is run." Section 7
  reports three studies and the Conclusions cite them.
- The Abstract contains two parenthetical editorial notes ("corrected
  2026-10-08 ...") that belong in a change log, and ten of its
  sentences exceed the 30-word limit in `STYLE.md` (the longest is 68
  words; VERIFIED by word count).
- Section 1.1 opens "We begin with an opinionated methodological
  claim ... it should routinely be analyzed as such" and 1.3 closes
  that the decomposition is "a methodological priority rather than an
  academic refinement". The Discussion and Conclusions then argue the
  opposite: decomposition is counterproductive unless the biomarker
  is contaminated or the estimand is attribution. The 2026-06-14
  report asked for this tension to be removed; it has been noted
  ("Two points of precision ... should be carried through") but the
  opening claim still stands.
- The cover letter describes a 100-replicate pilot with the production
  program "forthcoming", $N \in \{35, 70, 100, 150\}$, and a Section 6
  simulation. The README says the driver directory "does not yet
  exist" (flagged 2026-06-12 as m8; still open).
- The 2026-10-08 correction notes inside Sections 3.1, 4, 5.2, 7.1,
  7.3, 7.4 and 7.5 should be absorbed into the text; a journal version
  cannot carry "Corrected 2026-10-08:" parentheticals.

**Remedy.** Rewrite the Discussion's open-questions paragraph; rewrite
the abstract to the STYLE rule; replace the Section 1.1 "routinely"
claim with the conditional claim the paper actually supports; update
the cover letter and README; fold every correction note into the prose
and move the audit trail to a change log.

**Evidence.** VERIFIED (text), VERIFIED (word counts).

### M7. Open-placebo expectancy defined three ways; deception undiscussed

The expectancy assigned to an open-placebo phase differs across the
table, the protocol and Study B, and the recommended design depends on
deception that the paper does not discuss.

**Location.** Table `tab:terms` ($\eta$ "between 0 and 0.5 in
open-label placebo"); ADEMP protocol ($\eta = 0$ for open-label
placebo); Section 7.5 and `run-study-b-balanced.R` (open-placebo
phase: off drug, $\eta = 1$, "believes treatment continues").

**Problem.** "Open-label placebo" in the Kaptchuk sense (patient told
it is placebo) and the "open placebo" phase of Study B (patient told
treatment continues while receiving placebo) are different
interventions with expectancies at opposite ends of the scale. The
paper uses the same phrase for both. The Study B design also contains
a covert on-drug phase ($\eta = 0$ while on drug). Both phases require
deception. The balanced-placebo literature the paper cites
(Rohsenow and Marlatt 1981) is explicit that two of its four arms
deceive participants, and the ethics literature (Miller and
Kaptchuk's "authorized deception", and the finding that safeguards
adequate for a balanced-placebo trial may reduce its ability to
answer the additivity question) is directly relevant to a paper whose
practical recommendation is that biomarker-validation trials "should
be fielded as a balanced-placebo or open-hidden design". The
open-hidden paradigm of Colloca et al. (2004) is feasible in inpatient
settings with concealed infusions; its translation to a multi-week
outpatient N-of-1 trial of an oral drug is not addressed.

**Remedy.** Fix the vocabulary: "open-label placebo" for disclosed
placebo, "covert discontinuation" or "concealed placebo substitution"
for the Study B phase, and give each its own $\eta$ in the table.
Add a paragraph on ethical feasibility citing authorized deception and
the Meissner et al. (2023) protocol as an existence proof that an
open-hidden discontinuation design can be fielded in a psychiatric
outpatient population, and state what expectancy that design actually
induces (its hidden arms tell participants there is a 50 percent
chance of discontinuation, which is the manuscript's $\eta = 0.5$, not
$\eta = 0$).

**Evidence.** INSPECTED (text, driver, protocol); literature from the
searches in Section 8.

### M8. The belief-covariate decomposition is the series' F1 analysis

Its hybrid-design verdict here conflicts with the series' verdict
elsewhere.

**Location.** Section 7.5; Section 5.2.

**Problem.** The Study B decomposition is
`Sx ~ bm + t + Dbc + De + bm:Dbc + bm:De`, which is the F1 model
(`Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De`) that the project
recommends for the interaction and reports as robust to a prognostic
biomarker on the Hybrid design (`docs/47` to `docs/49`, `docs/56`;
project CLAUDE.md). Paper 06 does not name it and reports that on the
hybrid design it "performs no better than, and sometimes worse than,
the one-component analysis". Part of the discrepancy is M1: the
project's own `docs/65` found that the strength-0.6 prognostic-
confound cells of paper 02 were a PD-repair artifact of exactly this
kind. Part is the design: paper 06's sixteen-visit hybrid variant has
$\mathrm{Cor}(D_{bc}, \eta) \approx 0.75$ against about $0.5$ for the
standard four-path Hybrid, as the manuscript's own 2026-10-08 note
computes, so its "hybrid" is a harder case than the series' Hybrid.
The deviations log of 2026-05-07 already proposed replacing `phase`
with the continuous expectancy `e` as option (b) for the aliasing
problem; that is the F1 model, and it would have removed the need for
the fallback formula that every phase-augmented result in Section 7
now uses.

**Remedy.** Name the decomposition as F1, cite the series documents
that evaluate it, and reconcile the verdicts explicitly: on which
design, at which drug-belief correlation and at which contamination
strength does F1 recover the slope. Consider replacing the
phase-augmented analysis (Section 5.2) with F1 throughout, since the
aliased `Dbc:phase` term has made the phase-augmented analysis a
different model from the one the paper describes.

**Evidence.** INSPECTED (formulas, series documents as summarized in
the project CLAUDE.md); the series-level robustness claims themselves
are UNVERIFIED here.

### M9. The sample-size axis confounds design with $N$

The sixteen-visit design also sits at the feasibility bound even when
uncontaminated.

**Location.** Section 7.2, final paragraph; Section 7.3.

**Problem.** The manuscript states that the step from $N = 70$ to
$N = 200$ changes the design (eight visits, one path, to sixteen
visits, two paths with crossover) and is "not attributable to sample
size alone". That is correctly disclosed, but the abstract and Section
7.3 still read the axis as a sample-size sweep ("$0.016$ at $N = 35$
and $0.006$ at $N = 300$"). Further, the sixteen-visit design's
smallest eigenvalue is $0.045$ with no contamination (M1 table), so it
is within $0.05$ of infeasibility at the reference parameters; any
sensitivity analysis on this design that raises a correlation will
cross the bound. The minimum eigenvalue of the eight-visit design is
$0.344$.

**Remedy.** Either run Study A at $N = 200$ and $300$ on the eight-
visit design as well, so that a clean sample-size column exists, or
present the two designs as separate arms. Report the feasibility
bound of each design with its parameters.

**Evidence.** VERIFIED (eigenvalues); INSPECTED (text).

### M10. The `Dbc:phase` aliasing was known on 2026-05-07

The deviations log documented it then; the manuscript carried a
contrary explanation until 2026-10-08.

**Location.** Section 4 ("the rank-deficiency the pilot ... encountered
when the full phase-augmented formula was fitted at $N = 35$");
Section 7.1; Section 7.3 ("the inverted ranking persists with the full
formula fitted and converged").

**Problem.** `02-deviations-log.md` states on 2026-05-07 that the rank
deficiency "is a property of the design matrix, not of the sample
size" and that at $N \ge 70$ the full formula fails on every
replicate. The manuscript nonetheless told the reader, through four
rounds of review, that the deficiency was an $N = 35$ symptom of
near-collinear trajectories and that the production run fitted the
full formula. The 2026-10-08 notes now correct this, but the
corrected text still frames the aliasing as a surprise. This is a
process finding rather than a scientific one, but it bears on how much
weight a referee can give to other statements about what was run.

**Remedy.** Absorb the corrections (M6), cite the deviations log, and
describe the phase-augmented analysis from the outset as the model
that was fitted.

**Evidence.** VERIFIED (deviations log, stored `formula_dropped`
counts: 99,995 of 100,000 Study A phase-augmented fits and all 1,800
contamination fits record `Dbc:phase` dropped).

## 5. Minor comments

Numbers are verified against the stored summaries unless stated.

1. **Numbers that reproduce.** Study A mean absolute bias ($0.0156$,
   $0.0157$, $0.0058$, $0.0058$), MCSEs ($0.027$, $0.019$, $0.0085$,
   $0.0069$), power ($0.585$ versus $0.350$ at $N = 35$; $0.904$ versus
   $0.623$ at $N = 70$), one-component coverage ($0.964$ to $0.992$),
   Type I error (one-component $0.0306$ to $0.0558$; phase-augmented
   $0.038$ to $0.049$), minimum convergence $0.998$, the Study B
   decomposition $z$-scores ($1.05$, $0.41$, $1.03$, $1.32$) and the
   quick-sim pilot figures all reproduce (VERIFIED).
2. **Numbers that do not quite.** Phase-augmented coverage runs to
   $0.987$, which rounds to $0.99$, not $0.98$. The $N = 200$ null
   cell is $1.8$ MCSE above nominal, not $1.9$. The $m_{PB}$
   invariance at $c_{bm,PB} = 0.30$ is $0.0067$, not "within
   $0.006$". The contamination sequence $-2.31, -2.37, -2.08, -1.74$ is
   the $m_{PB} = 6$ row; the $m_{PB} = 0$ row starts at $-2.24$. Say
   which row is quoted.
3. **Study B at $N = 70$, $c_{bm,PB} = 0.2$.** The decomposition bias
   is $-0.076$ at an MCSE of $0.029$, $2.6$ MCSE from zero. "The
   $N = 70$ cells show the same pattern with the decomposition
   unbiased" overstates; report the cell.
4. **Phase-augmented "indistinguishable" from one-component in the
   contamination pilot.** At $N = 300$, $c_{bm,PB} = 0.15$, $m_{PB} = 0$
   the two biases are $-0.029$ and $-0.111$ (MCSEs $0.030$ and
   $0.024$), about $2$ standard errors apart. "Indistinguishable at
   every cell" is not supported at pilot resolution.
5. **Notation rule 1.** Table `tab:terms` defines $D_{it}$ as the
   continuous-decay indicator. NOTATION.md reserves $D_{it}$ for the
   binary state and $D_{bc,it}$ for the decayed one, and Section 6.1
   then uses $D_{it}$ in the model display. `notation-lint.pl` is
   clean (VERIFIED) but does not check this rule.
6. **STYLE terminology.** The STYLE grep returns fifteen prose uses of
   "coupling" (lines 112, 113, 207, 260, 1506, 1547, 1590, 1606,
   1630, 1631, 1646, 1647) and "ceiling" in a non-feasibility sense.
   No Glossary section exists (STYLE Section 4). Study A and Study B
   are lettered; NOTATION Part 3 asks for numbers.
7. **US English.** "enrol" appears twice (Sections 1.2 and 3.3.1).
8. **Section 2.** The model display ends with a comma and the next
   sentence begins "Each component". The Gompertz display introduces
   $m_C$, $d_C$, $r_C$ without defining them until Section 3.1. The
   AR(1) residual is written $\rho^{|h|}$ on the visit index; the
   drivers space visits at $2.5$ or $1.25$ weeks and NOTATION uses
   calendar weeks. The cross-component covariances $c_{cf1t}$ and
   $c_{cfct}$ are still undefined (2026-06-13 Minor 1); NOTATION calls
   them $c_1$ and $c_\times$.
9. **$\sigma_{PB} = 2$.** The 2026-06-15 report asked that the
   individual-level PB standard deviation, which sets the scale of the
   whole contamination effect, be tied to the prazosin calibration. It
   is still introduced only inside Section 7.4.
10. **Replicate-count table.** Still absent (2026-06-15 Minor 2). With
    four runs at 100, 1,000 and 5,000 replicates the reader needs one.
11. **Figures.** The paper has none. The bias-versus-$c_{bm,PB}$ curve
    by design and analysis is the natural single figure.
12. **"Study A ... returns a result that does not support the
    framework's central bias prediction."** After the identity is
    stated, Study A is the predicted null; the sentence reads as a
    surprise the paper has already explained.
13. **Section 5.1, fourth reason**, says the one-component interaction
    "is the right marginal estimand for the validation question whether
    the moderation is biological or expectation-mediated"; Section 8
    says the opposite. Reconcile.
14. **Worked example.** The TV increment over weeks 9 to 12 is still
    unquantified (2026-06-13 Minor 4). The "After fitting, the model
    returns $\hat m^A_{BR} = 4.0$ ..." sentence describes an exact
    recovery by a model the paper never fits; mark it as the noiseless
    limit.
15. **Hróbjartsson.** The 202-trial count is now attached to the 2010
    review in Section 1.2 ("the last covering"), but Section 3.3.1 still
    groups the three citations before the count. Minor.
16. **Data availability.** No URL or persistent identifier
    (2026-06-13 Minor 7).
17. **`in_header: claudecode.tex`.** Inert per the project notes; remove
    for submission.
18. **Bibliography hygiene.** All 56 cited keys resolve and no entry is
    uncited (VERIFIED). `nurko2022` has `and others` before a named
    author, the same defect fixed in `wang2022comt` in June.
    `rohsenow1981` lacks a DOI. `senn2018` and `benedetti2014` lack
    DOIs. `iche9r1` has no URL.

## 6. Missing and questionable references

| Claim and location | Issue | Suggested source |
|---|---|---|
| Section 5, 7.5: linear mixed decomposition of drug-treated response into specific and non-specific parts | Closest statistical precedent uncited; flagged 2026-06-13 | Petkova, Tarpey, Govindarajulu (2009) *Int J Biostat* 5(1), PMID 21556319; Tarpey et al. (2010) *JASA* 105:968, PMID 21494314 (open-label phase then blinded discontinuation, the paper's own design) |
| Section 4, 7.5: open-hidden discontinuation as the identifying contrast | Empirical realization in a psychiatric population uncited; flagged 2026-06-12 and 06-14 | Meissner et al. (2023) *BMC Psychiatry* 23:457, PMID 37344789; NCT05191277 |
| Section 7.5: balanced-placebo design recommended for biomarker trials | Deception and feasibility undiscussed | Miller and Kaptchuk, authorized deception (*J R Soc Med* 2008, PMC2598073); Martin and Katz on authorized deception and placebo magnitude; the 2007 *BMC Health Serv Res* survey of antidepressant participants (doi in Section 8 sources) |
| Section 7.5: balanced-placebo design is being fielded | Contemporary instances strengthen the recommendation | Esketamine balanced-placebo fMRI protocol (2023) *Trials*, PMID 37568215; crossover balanced-placebo N-of-1 provocation protocol (2026) *Front Pharmacol*, PMID 42827697 |
| Section 1.3: "no published N-of-1 trial has employed a structural decomposition" | Audit basis is Li 2016 (1985 to 2013); recent aggregated N-of-1 trials should be acknowledged | Tamsulosin N-of-1 deprescribing RCT (2026) *JAMA Netw Open*, PMID 42406401, as an example of current practice with a single treatment indicator |
| Section 3.3.1: subadditivity evidence | Still cites Kaptchuk 2008 and Benedetti 2014, which do not document drug-placebo subadditivity (2026-06-13) | Bingel 2011 (already cited) for the interaction; Amanzio and Benedetti (1999) for naloxone reversal |
| Section 1 motivation | Schork (2015) *Nature* 520:609 suggested 2026-06-12; not added | Optional |

Bib entries present and correctly used, spot-checked against PubMed
metadata: `hendrickson2020`, `liao2023`, `colloca2004`, `rohsenow1981`,
`carrozzo2025`, `gaertnerCarryover2023`, `dazaCausalAnalysisSelftracked2018`,
`kaptchuk2008components`, `zunhammer2018nps`, `wager2013nps`
(INSPECTED; prior rounds verified the clinical citations' claim fit).

## 7. Suggestions for strengthening

1. Rerun the contamination and hybrid-recovery studies under
   `stopOnNonPD = TRUE` on each design's feasible range, and report
   the feasibility bound beside each design. Expect the sign of the
   one-component bias to follow the expectancy schedule as in M2.
   This converts the paper's weakest section into its clearest.
2. State the displacement $-c_{bm,PB}\,\sigma_{PB}\,
   (\bar\eta_{\text{on}} - \bar\eta_{\text{off}})$ in Section 6.1 and
   use it to predict, before each simulation, what the one-component
   estimator will do on each design. A paper that predicts its
   simulation is stronger than one that explains it afterwards.
3. Adopt F1 (biomarker by $D_{bc}$ and by $\eta$) as the decomposition
   throughout, name it, drop the aliased phase-augmented model, and
   reconcile the verdict with the series.
4. Run the $\eta$-sensitivity cell ($\eta_{\text{blinded}} \in
   \{0.25, 0.5, 0.75\}$) and a prognostic $c_{bm,TV}$ cell, since
   `c.bm.tv` now exists in the package; the TV half of the identity is
   currently untested.
5. Add one figure: bias of the one-component and F1 estimators versus
   $c_{bm,PB}$, one panel per design, with the identity's prediction
   drawn as a line.
6. Add a replicate-count table and a change log, and move every dated
   correction note into it.
7. Add a short ethics and feasibility paragraph for the recommended
   design (M7).

**Status of prior-round comments.** From 2026-06-12: M1 and M2 (sign
errors) applied; M3 (pharmacometric prior art) addressed; M4
(carryover rivals) addressed; M5 (estimand statement) addressed; M6
(joint identifiability) still asserted, not demonstrated (M4 above);
m8 (README) open; Schork 2015 not added. From 2026-06-13: Major 1
(Gompertz ceiling) fixed; Major 2 (steady-state assumption in the
phase contrast) not addressed; Major 3 (bias versus estimand mismatch)
not addressed, "biased" remains throughout; Major 4 (run the
simulation) done; Petkova not cited; Minors 1, 2, 4, 7 open. From
2026-06-14: M1 to M7 addressed as the 06-15 report records. From
2026-06-15: R1 (Study B) done on the balanced design; R2 (production
contamination run) not done, one cell exists; R3 (wire numbers to
data) done for Study B only, Study A and the contamination pilot
remain hard-coded; R4 done; Minor 1 (derive the sign) not done, and
M2 above supplies it; Minors 2 and 3 open; the suggested figure is
absent.

## 8. Scope of this review

**Verified (executed).** Reconstruction of every reported Study A,
contamination-pilot and Study B figure from the stored summaries
(`scratchpad/recompute.R`); the `formula_dropped` counts from
`study-a-all-reps.rds` and `study-contam-all-reps.rds`; the contents
of `study-b-recovery-pilot/`; positive definiteness of the covariance
matrix for every design and contamination level used, over the full
Study A grid, and the correlations the June repair would have produced
(`scratchpad/pdcheck.R`); the identity-based displacement arithmetic
of M2; the driver seeds and paths; the stub status of the nonlinear
fit; citation-key cross-check between the Rmd and the bib; the STYLE
grep; `perl tools/notation-lint.pl` (clean); abstract sentence word
counts; existence of the compendium sections and companion files the
manuscript cites.

**Inspected (read, not executed).** The four drivers and
`00-common.R`; the ADEMP protocol and deviations log; `buildSigma()`
and the repair lock in `R/generateData.R`; the git history of that
file; the prior referee reports; the cover letter, README and
`references.bib`; the project CLAUDE.md, STYLE.md and NOTATION.md.

**Not checked.** The compendium's derivations in the sections the
manuscript cites; the $\mathrm{Cor}(D, \eta)$ figures of the
2026-10-08 note (taken as reported); the series documents `docs/47`
to `docs/49`, `docs/56` and `docs/65` beyond their summary in the
project CLAUDE.md; the paper 07 Type I error figures imported in
Section 3.1; the clinical citations' claim fit beyond the spot checks
above; the rendered PDF; whether the stored runs were made on the host
or in the container; `report-slim.Rmd` and `bullets.Rmd`.

**Searches run (no identifier or email passed to any service).**
PubMed: Petkova Tarpey Govindarajulu placebo drug-treated; N-of-1 or
single-case with placebo response or expectancy and decomposition
(2019 onward); open-hidden, covert administration or balanced placebo
design with discontinuation or expectancy (2015 onward); placebo
response or disease progression model with drug effect, longitudinal,
pharmacometric (2019 onward); predictive biomarker or treatment-effect
heterogeneity with placebo response and confounding in crossover or
N-of-1 (2018 onward); placebo response with disease progression and
mixed models in depression (2020 onward). Web: open-hidden
antidepressant discontinuation protocol; authorized deception and the
balanced placebo design. Sources consulted:
[NCT05191277](https://clinicaltrials.gov/study/NCT05191277),
[Meissner et al. 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10286435/),
[Miller and Kaptchuk, deception in placebo research](https://pmc.ncbi.nlm.nih.gov/articles/PMC2598073),
[Martin and Katz, authorized deception](https://pmc.ncbi.nlm.nih.gov/articles/PMC1079669),
[BMC Health Serv Res 2007 survey](https://doaj.org/article/3cd6ab12fe384c82911eec39392fc4e8).
PubMed results cited by PMID above carry DOIs in the PubMed records;
the PubMed database is the source for all PMID-identified items.
