# Referee report: paper 13, two-component DGP without the placebo component
*2026-10-29 16:53 PDT*

**Manuscript.** 'Data generation machinery for N-of-1 clinical trial
simulations under a two-component response decomposition'
(`analysis/report/13-two-component-drop-pb/report.Rmd`, last committed
2026-10-09 in e077526a8, whose only changes are two 'Corrected
2026-10-08' notes; last substantive rewrite 2026-09-09 in ef8d58247).
The `report.pdf` and `report.tex` in the directory are dated
2026-09-09 08:52.

**Submitted as.** No cover letter or title page exists. The only venue
signal is `csl: statistics-in-medicine.csl`. Article type not stated.

**Reviewer.** pmsimstats team, internal referee pass.

**Evidence labels.** VERIFIED: I executed code and observed the result.
INSPECTED: I read the source and confirmed it. SUSPECTED: consistent
with what I read but not confirmed. UNVERIFIED: not checked.

**Scratch scripts** (not part of the repository):
`p13-check1.R` (every inline and hard-coded number against the stored
paper 13 grid), `p13-check2.R` (the paper 14 grid, MCSE of the
carryover costs, pooled Type I error, model-to-empirical SE ratios),
`p13-check3.R` (feasibility bounds of Appendix C.4 recomputed with the
current `buildSigma()`, under `c.cfct = 0.1` and `0.2`; positive
definiteness of every production cell with repair refused; the
coupling vectors; the placebo variance share), all in
`/private/tmp/claude-501/-Users-zenn-prj-res-36-pmsimstats-ng-pmsimstats-ng/12d860d5-0c62-4f77-a07f-554bceae2b3e/scratchpad/`.

## 1. Summary

The manuscript repeats a mean-moderation versus covariance-moderation
comparison under a reduced data-generating process (DGP) that keeps the
drug-response (BR) and natural-history (TV) components and removes the
placebo-response component (PB). It crosses two architectures, three
designs (CO, Hybrid, OL+BDC), $N \in \{35, 70\}$,
$c_{bm} \in \{0, 0.30, 0.45\}$ and $t_{1/2} \in \{0, 0.5, 1\}$ week
(108 cells, 500 replicates). It argues from the crossover and placebo
literature that dropping PB rather than TV is the conventional
reduction. It reports that the architecture contrast survives, that
the reference cell is compressed by a power ceiling, and that Type I
error is design-dependent. Appendix C derives the feasibility bound on
$c_{bm}$ as a Schur complement, with a closed-form spectrum under
compound symmetry, and shows that the component count barely moves the
bound.

The central defect is that the Results, Abstract, Conclusions and
rendered PDF do not describe this paper's simulation. Every reported
number, and several qualitative conclusions, come from the drop-TV run
that now belongs to paper 14. The Rmd also points at a data file that
does not exist, so it cannot be rendered. Against its own stored data
(`analysis/data/13-drop-pb-grid.rds`), four of the paper's stated
results are false. The stored run also predates three changes to the
project's DGP and inference. A rerun and a rewrite of Sections 3 to 5
are required.

## 2. Overall assessment and recommendation

**Recommendation: not submittable. Rerun under the current DGP, then
rewrite the Results and the claims built on them.**

The methods section, the literature argument of Section 1.1 and
Appendices B and C are largely sound. Appendix C's ceilings reproduce
exactly. The paper's own data, however, have never been written up.
Its framing also rests on a paper-01 comparison that no longer exists,
and on a reading of the two architectures that paper 11 withdrew on
2026-10-29. Its contribution as a stand-alone *Statistics in Medicine*
article is thin (Section 3). The most efficient path is probably a
merged paper 13/14 technical note or a supplement to paper 01. That is
a scoping decision for the authors.

## 3. Significance and novelty

- **The core result is predictable from first principles.** A latent
  component that is independent of the biomarker and enters the
  outcome additively changes only the residual covariance. It cannot
  move the expectation of $\hat\beta_{bm:D}$. The paper concedes this
  (Section 4.1, paper 06's omitted-variable identity). What remains is
  a variance effect on absolute power, which the paper itself says it
  cannot quantify without a matched three-component run (Section 3.4,
  not performed). INSPECTED.
- **Appendix C is standard algebra.** The Schur-complement bound and
  the spectrum of a block-exchangeable (Kronecker-sum) matrix are
  textbook results. The project already has the bound in closed form
  (`docs/34`, `docs/35`, `docs/65`, compendium). The genuinely new
  element is the observation that $\lambda_{\min}$ under compound
  symmetry is free of $n$ and $K$. INSPECTED.
- **The design-comparison use case is undercut by the paper's own
  Section 2.2.3.** With PB removed, the expectancy weight $\eta$ has no
  effect on the data. The Hybrid and OL+BDC designs exist to exploit
  the open-label versus blinded contrast. A DGP that cannot represent
  that contrast cannot rank them. Yet Section 4.2 recommends the
  two-component DGP for comparisons 'between trial designs'. This is
  an internal contradiction. INSPECTED.
- **Within the compendium.** Paper 14 (drop TV) has the same grid,
  driver and Section 4.5 text. Two papers presenting one sensitivity
  analysis under mirrored reductions will be read by a referee as
  salami slicing. INSPECTED.

## 4. Major comments

### 4.1 The reported results belong to paper 14's drop-TV run (correctness)

- **Location.** Abstract (Results, Conclusions); Sections 3.1 to 3.4;
  Section 5 items 1, 3 and 4; captions of the main tables; the
  rendered `report.pdf`/`report.tex`.
- **Problem.** The setup chunk reads
  `analysis/data/13-two-component-grid.rds`. That file was renamed to
  `14-drop-tv-grid.rds` in ef8d58247 (2026-09-09 10:42). Its metadata
  give `components = c('pb', 'br')`. INSPECTED (git) and VERIFIED
  (meta). Paper 13's own output is `13-drop-pb-grid.rds`, written
  2026-09-09 12:22 with `components = c('tv', 'br')`. VERIFIED. As the
  Rmd stands, `HAVE_SIM` is `FALSE`. The table chunks call `filter()`
  on `NULL`, so the document cannot knit. SUSPECTED (failure mode
  inferred from the code, not executed, since rendering was out of
  scope).
- **The hard-coded numbers match paper 14, not paper 13.** VERIFIED
  against both files (`p13-check1.R`, `p13-check2.R`):

  | Claim (location) | Text | Paper 14 grid | Paper 13 grid |
  |---|---|---|---|
  | Hybrid cost, covariance, $c_{bm}=0.30$, $N=70$ (Concl. 1) | 0.194 | 0.194 | 0.228 |
  | Hybrid cost, mean (Concl. 1) | 0.040 | 0.040 | 0.078 |
  | OL+BDC cost, covariance (Concl. 1) | 0.196 | 0.196 | 0.182 |
  | OL+BDC cost, mean (Concl. 1) | 0.078 | 0.078 | 0.086 |
  | Pooled Type I, CO (Concl. 4) | 0.061 | 0.061 | 0.065 |
  | Pooled Type I, OL+BDC (Concl. 4) | 0.021 | 0.021 | 0.033 |
  | Minimum power, $c_{bm}=0.45$, $N=70$ (Sec. 3.2, caption) | above 0.89 | 0.890 | 0.790 |

  The rendered `report.tex` (2026-09-09 08:52) states 'a reduced
  data-generating process carrying only BR and PB' and prints the
  paper-14 values. INSPECTED. No stamped render of the drop-PB text
  exists in `share/`. The last paper-13 render is
  `13-two-component-dgp-2026-09-09-0852-8fef8e4-wip.pdf`, which
  predates the drop-PB data. INSPECTED.
- **Four stated findings are false against paper 13's data.**
  VERIFIED.
  1. *'The reference cell is ceiling-bound ... every entry lies above
     0.89.'* Covariance moderation under OL+BDC falls to 0.790 at
     $t_{1/2} = 1$ (MCSE 0.018). The Hybrid cost at $c_{bm} = 0.45$ is
     0.110, about half the 0.228 at $c_{bm} = 0.30$. It is not 'roughly
     a sixth'. The 'factor of about six' compression (Sec. 3.2,
     Concl. 3) does not hold.
  2. *'Under CO neither architecture loses power ... point estimates
     move by less than one standard error.'* Mean moderation under CO
     loses 0.054 at $c_{bm} = 0.30$, $N = 70$ (unpaired SE of the
     difference 0.031, $z = 1.7$). Covariance moderation under CO loses
     0.050 at $c_{bm} = 0.45$ and 0.074 at $c_{bm} = 0.30$, $N = 35$.
  3. *'Mean moderation does not [lose power]' (Concl. 1); 'both
     mean-arm costs are within two or three [MCSE]' (Sec. 3.1).*
     The mean-arm costs are 0.078 (Hybrid, $z = 2.6$) and 0.086
     (OL+BDC, $z = 2.8$) on the unpaired SE of the difference. That is
     3.5 to 3.9 single-cell MCSEs.
  4. *'Comparison with the three-component setting: ... the ranking,
     the direction, the insensitivity of CO, and the design ordering
     are all as paper 01 reports them.'* Paper 01 no longer reports any
     of these (comment 4.3). No such check was possible.
- **What does survive.** The covariance-arm cost exceeds the mean-arm
  cost under Hybrid ($0.150 \pm 0.043$) and OL+BDC ($0.096 \pm 0.044$)
  at $c_{bm} = 0.30$, $N = 70$. The Type I pattern by design survives
  qualitatively (CO 0.065, $z = +4.8$; Hybrid 0.042, $z = -3.2$;
  OL+BDC 0.033, $z = -7.5$; 6,000 null replicates each). VERIFIED.
- **Remedy.** After the rerun of comment 4.2, point the setup chunk at
  the correct file. Remove every hard-coded result from the Abstract,
  the Conclusions and the captions. Rewrite Sections 3.1, 3.2 and 5
  from the new output. Re-render through `tools/render.sh`, and delete
  the stale `report.pdf`/`report.tex`.

### 4.2 The stored run is not the project's current DGP or inference

- **Location.** Section 2; Reproducibility;
  `analysis/scripts/two-component-dgp/01-run-grid.R`.
- **Cross-component correlation.** The driver sets `c.cf1t = .2`,
  `c.cfct = .1` (lines 125 and 142). This is the defect documented in
  `docs/61`. The extended DGP specifies one constant $c_1 = 0.2$ at
  every lag. INSPECTED. The change matters for Appendix C. Under
  `c.cfct = 0.2` the two-component bound rises from 0.469 to 0.487
  (Hybrid, $t_{1/2} = 0$) and from 0.464 to 0.481 (OL+BDC). VERIFIED
  (`p13-check3.R`). Its effect on power is UNVERIFIED. `docs/61`
  reports a larger minimum eigenvalue under the correction, so the
  residual covariance changes in every cell.
- **Positive-definite handling.** The driver passes
  `makePositiveDefinite = TRUE` and ran before the 2026-10-26 stop
  default. All 108 cells are positive definite under both `c.cfct`
  values with repair refused (0 failures over all paths). VERIFIED. No
  stored cell depended on a repair. Appendix C's motivation ('drivers
  in this program repair a non-positive-definite matrix silently') and
  C.5 ('rely on the simulation to complain, since it will not') are
  now false. Under `stopOnNonPD = TRUE` the simulation stops.
- **Elements that do match the extended DGP.** These are separable
  AR(1) at $\rho = 0.7$ in calendar weeks, the decaying correlation
  with zero before first exposure, and carryover decaying from the last
  on-drug mean. INSPECTED in `buildSigma()` at the data commit ef8d582,
  and unchanged since apart from the positive-definite branch (git
  diff). The coupling vectors were VERIFIED for the CO and Hybrid
  paths.
- **Inference.** The analysis is E2 at the true half-life, with a
  random intercept, `corCAR1` and model-based standard errors
  (`useDE = FALSE`). INSPECTED. The project now recommends CR2
  standard errors with Satterthwaite df (F1). The CR2 check of
  commit 8e810ca9d shows the model-based test conservative in OL+BDC
  (0.024 to 0.027) and CR2 calibrated (0.048 to 0.051). It reports
  ratios of mean model SE to empirical SD of 1.13 to 1.15
  (`analysis/data/cr2-olbdc-check/summary.csv`). INSPECTED. Paper 13's
  own summary already shows the same mechanism: that ratio is 0.91 to
  1.01 under CO and up to 1.20 (Hybrid) and 1.25 (OL+BDC). VERIFIED.
  Section 3.3's 'we did not diagnose the mechanism ... an untested
  conjecture' is therefore out of date on both counts.
- **Silent fallback.** `lme_analysis()` refits without `corCAR1` when
  the first fit errors, and counts the result as converged. INSPECTED.
  The reported non-convergence rate of zero says nothing about how
  often this happened. UNVERIFIED.
- **Random numbers.** Each cell uses `set.seed(SEED + i)`. Every
  across-half-life 'cost', which the paper calls the comparison that
  'carries the argument', is therefore unpaired. INSPECTED. Common
  random numbers across $t_{1/2}$ within a design and architecture
  would roughly halve the SE of the costs. That last point is
  SUSPECTED and depends on the cross-cell correlation of the
  rejection indicators.
- **Remedy.** Rerun with `c.cfct = 0.2`, the default `stopOnNonPD`,
  common seeds across half-lives, and both E2 and E2cr2. Record the
  CAR1 fallback count, and store bias, empirical SE, $\kappa$ and MCSE
  per cell. The local run took 6,871 s at 500 replicates (meta). At
  1,000 replicates on the simulation host this is a few hours.
  Recompute Appendix C.4 under `c.cfct = 0.2`.

### 4.3 The framing rests on a comparison that no longer exists

- **Location.** Abstract Background; Section 1, paragraph 2; Sections
  3.1, 3.4; README 'Differences from paper 01'.
- **Problem.** The paper presents itself as paper 01's question under
  fewer components ('paper 01's central finding, reproduced on a
  correlation matrix two-thirds the size'). Paper 01 was narrowed on
  2026-10-07 to the extended DGP specification under covariance
  moderation. Its current `report.Rmd` contains no occurrence of 'mean
  moderation'. INSPECTED. The project CLAUDE.md records that the
  matched-$N$ cross-architecture figure survives only in NOTATION.md
  and an archived draft. It also records that the comparison was
  flagged for non-comparable absolute-power baselines
  (`01-dgp-mean-moderation-vs-mvn/revisions/33-paper01-dgp-architecture-review.md`,
  Section 3.3). Paper 13 repeats that design. It compares absolute
  power costs across architectures at one nominal $c_{bm}$, from
  different baselines: at $t_{1/2} = 0$, Hybrid power is 0.722 versus
  0.664, and OL+BDC 0.640 versus 0.646. VERIFIED. A cost in absolute
  power depends on where on the power curve the baseline sits.
- **Section 2.3 is also stale.** It contrasts the carryover rule with
  'the one behind paper 01's published results'. Every current
  manuscript now uses carryover decaying from the last on-drug mean.
- **Remedy.** Restate the aim without paper 01 as its anchor. Report
  the carryover effect on a scale that does not depend on baseline
  power, such as the ratio of noncentrality parameters
  $(\beta/\mathrm{sd})_{t_{1/2}} / (\beta/\mathrm{sd})_0$ with MCSE.
  Alternatively, calibrate the two architectures to equal power at
  $t_{1/2} = 0$.

### 4.4 Appendix A misdescribes the two architectures

- **Location.** Appendix A.1, A.2; Section 4.5.
- **Problem.** A.2 states that under mean moderation 'two participants
  sharing a biomarker value share an expected treatment effect
  exactly, and any difference between their realized responses is
  measurement error'. In the code, BR is drawn with $\sigma_{BR} = 8$
  independently of the biomarker, and the shift
  $\beta_{bm} z_i \sigma_{BR}$ is added on top (`generateData()`).
  INSPECTED. Mean moderation therefore leaves the full BR variance as
  unexplained heterogeneity. Covariance moderation instead reduces
  $\mathrm{Var}(BR \mid b)$ by the factor $1 - c_{bm}^2$. The contrast
  in A.2 is backwards. Paper 11 (rewritten 2026-10-29, commit
  6ea37805e) now states that the two are 'two parameterizations of one
  Gaussian model, not two mechanisms'. Both set $E[BR \mid b]$, and
  only the correlation channel lowers $\mathrm{Var}(BR \mid b)$.
  INSPECTED. A.1's 'not two options among many, but the two components
  of a decomposition' conflicts with that reframing.
- **Unreported bias under mean moderation.** The mean-shift is applied
  at on-drug visits only, but E2's $D_{bc}$ assumes the interaction
  carries over. The mean-moderation estimate therefore drifts from its
  true value as $t_{1/2}$ grows. VERIFIED. At $c_{bm} = 0.30$,
  $N = 70$, the true value per raw biomarker unit is
  $-0.30 \times 8 / 15.36 = -0.156$. The mean estimate is $-0.176$
  (Hybrid) and $-0.186$ (OL+BDC) at $t_{1/2} = 1$, a 13 to 19 percent
  excess in magnitude. At $c_{bm} = 0.45$ it is 15 to 20 percent
  (MCSE of the mean estimate about 0.0036). Under covariance
  moderation the estimate stays at $-0.158$ to $-0.161$. The paper
  stores `beta` and `empirical_se` but reports neither. Its
  mean-moderation 'cost' therefore mixes a variance effect with a
  misspecification bias.
- **Remedy.** Rewrite A.1 and A.2 to agree with paper 11. Report bias,
  empirical SE and coverage alongside power for both architectures.
  State the true value of $\beta_{bm:D}$ in raw biomarker units, with
  the sign convention (NOTATION rule 5).

### 4.5 Drop-TV text left in a drop-PB paper

- **Location and problem.** INSPECTED in each case.
  - Abstract Conclusions: 'The TV component is inert ... and acts only
    as a variance term'. This paper removes PB.
  - Appendix A.1: 'Dropping TV removes rows from $\Sigma$'.
  - Appendix A.4: 'Removing TV therefore changes what the analysis is
    fitted to'.
  - The STATUS comment and README give the render path as
    `13-two-component-dgp/report.Rmd`. The directory is
    `13-two-component-drop-pb`.
  - The Reproducibility section names the data file as
    `13-two-component-grid.rds`.
- **Remedy.** Search the Rmd for 'TV' and 'two-component-dgp' and
  correct each instance against the drop-PB design.

### 4.6 The decaying correlation is misstated before first exposure

- **Location.** Section 2.2.2 ($\phi$); Appendix B.4 ($b_p$).
- **Problem.** Both definitions give $c_{bm} e^{-\lambda t_{sd}}$ at
  every off-drug visit. Before first exposure $t_{sd} = 0$, so the
  formula returns the full $c_{bm}$. The code sets 0 there. VERIFIED:
  CO path B gives $(0, 0, 0, 0, 0.45, \ldots)$. Appendix C.4 relies on
  the zero ('carry zero coupling at every half-life'). At
  $t_{1/2} = 0$ the code also sets the off-drug correlation to 0
  (`lambda_cor = 0` branch), which the formula leaves undefined.
- **Remedy.** Use NOTATION's $g_t$: 1 on drug, $\phi_t$ off drug
  after exposure, 0 before first exposure.

### 4.7 Missing bibliography entries

- **Problem.** Four cited keys are absent from the local
  `references.bib`: `schmid2022bayesian`, `hrobjartsson2001`,
  `hrobjartsson2004` and `hrobjartsson2010`, all in Section 1.1.
  VERIFIED (grep). The
  shared `analysis/report/references.bib` holds the Hróbjartsson
  entries and a Schmid and Yang entry under the key
  `schmidBayesianModelsNof12021`. INSPECTED. Section 1.1, the paper's
  main novelty argument, would therefore render with unresolved
  citations.
- **Remedy.** Copy the entries, reconcile the Schmid key and year, and
  check each against the source.

## 5. Minor comments

1. **MCSE.** The main tables carry no MCSE. Only the maximum (0.022)
   is quoted. VERIFIED. Each cost needs the SE of the difference
   (about 0.031 to 0.032 here). STYLE section 1 requires an MCSE or a
   replicate count beside every number.
2. **Type I error.** The check is present and correctly pooled
   (Section 3.3 table). Pooling over $N$ hides Hybrid at 0.037 for
   $N = 70$ against 0.046 for $N = 35$. VERIFIED. 500 null replicates
   per cell are too few to resolve a 0.01 departure. The 1,000-to-5,000
   replicate null cells used in the CR2 check are the precedent.
3. **'Monotone in off-drug information' (Section 3.3)** is asserted
   without a measure of off-drug information. Either define one or
   drop the sentence.
4. **The recursion check of Section 2.3** (three cells, 400
   replicates, at most 0.005 of power) has no script or stored output.
   UNVERIFIED. No matching script was found under `analysis/scripts/`.
   It is also moot, since every current driver uses the corrected
   rule.
5. **Notation (NOTATION.md).** `tools/notation-lint.pl` reports the
   file clean. VERIFIED. The following problems are not mechanized and
   need fixing by hand:
   - 'Architecture A/B' should be the named architectures (labels
     retired).
   - $c_{bm}$ is used as the mean-shift multiplier without the rule-2
     declaration.
   - $\mathrm{mg}(t; a, d, r)$ should be $G(t_{od}; m, d, r)$.
   - $N(0, \sigma_u^2)$ should be $\mathcal N$.
   - $\beta_t t_p$ should be $\beta_T w_t$.
   - In Appendix C, $A$, $B$ collide with the within-component block
     $A_c$, the compound-symmetry correlation should be $\rho_{cs}$,
     and $\lambda_{\min}$ collides with the decay rate $\lambda$.
   - No sign convention is stated (rule 5).
6. **STYLE.md.** The STYLE grep finds 51 hits of 'coupling',
   'ceiling', 'graded', 'anchored' or 'machinery' (the last in the
   title). VERIFIED. Further problems:
   - 'Covariance moderation' and 'mean moderation' should be the
     correlation and mean-shift encodings in prose.
   - There is no Glossary or Abbreviations section.
   - Abstract sentences exceed 30 words.
   - Semicolons join clauses throughout.
   - Evaluative phrasing recurs: 'striking', 'the price is', 'buys
     between 0.007 and 0.013', 'a narrower licence'.
7. **US English.** British spellings recur: 'modelled' and 'modelling'
   (10), 'favours' (3), 'licence', 'neighbours'. VERIFIED (grep).
8. **Embedded correction notes.** The two 'Corrected 2026-10-08' notes
   (Sec. 3.1; App. C.4) and the in-text repository paths belong in a
   change log, not a submitted paper. Their content was checked: the CO
   first off-drug value is $2^{-2.5} = 0.177$ and the second $0.031$.
   The bound for the drug-first path rises with $t_{1/2}$. VERIFIED.
9. **Appendix C.4 table.** It reproduces to three decimals under the
   stored `c.cfct = 0.1`, both the three- and the two-component
   columns. VERIFIED. The 'number of occasions' and 'index versus week
   gap' figures came from a reimplemented builder, as the provenance
   note says. UNVERIFIED.
10. **Placebo share.** 'Under 38%' is correct: 0.379 at open-label
    visits and 0.132 at blinded visits. VERIFIED. Section 2.2.3 should
    give both numbers.
11. **Designs differ by more than on/off pattern.** With $\eta$ inert,
    the designs still differ by calendar spacing, through the AR(1)
    week gaps and the TV mean. Sections 2.2.3 and 4.4 say 'only
    through their on-drug patterns'. INSPECTED.
12. **Exploratory framing.** Nothing describes the study as
    pre-registered or confirmatory. VERIFIED (grep). A one-line
    statement that the study is exploratory, with an ADEMP table, is
    expected (NOTATION 'Estimands and reporting').
13. **Front matter.** There is no cover letter, title page, data
    availability statement or article-type statement. The title
    ('Data generation machinery') does not describe the question.

## 6. Missing and questionable references

| Reference | Status | Relevance |
|---|---|---|
| Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. *Stat Med* 2019;38:2074-2102 | In bib, never cited | ADEMP structure and MCSE reporting |
| Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. *J Bus Econ Stat* 2018;36:672-683 | Missing | CR2 inference the rerun should use |
| Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS*. Springer, 2000 | Missing | `nlme`, `corCAR1` |
| Higham NJ. Computing the nearest correlation matrix. *IMA J Numer Anal* 2002;22:329-343 | Missing | Repair that Appendix C discusses |
| Rohsenow DJ, Marlatt GA. The balanced placebo design: methodological considerations. *Addict Behav* 1981;6:107-122 | Missing | Section 4.5 invokes balanced-placebo designs uncited |
| Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*, 3rd ed., 2014; Senn S. *Cross-over Trials in Clinical Research*, 2002 | In bib, never cited | Section 1.1's claim that period and drift terms are standard |
| Hróbjartsson and Gøtzsche 2001, 2004, 2010 | Cited, absent from local bib | Comment 4.7 |
| Schmid and Yang (Bayesian N-of-1 models) | Cited as `schmid2022bayesian`, absent; shared bib has a 2021 key | Comment 4.7; year to be checked |
| pmsimstats papers 01, 06, 11, 14 | Cited in prose without bib keys | 01 and 11 mis-described (comments 4.3 and 4.4) |

The reference details in this table are from memory and from the
repository bibliographies. No literature database or web search was
run in this pass, so page ranges and years are UNVERIFIED.

## 7. Suggestions for strengthening

1. Merge papers 13 and 14 into one short paper or supplement. Run one
   grid with three arms: full, drop PB and drop TV. Use common random
   numbers across arms. This also supplies the matched three-component
   run whose absence Section 3.4 and `02-deviations.md` record as a
   limitation.
2. State the variance result analytically. Removing an independent
   additive component changes the GLS information for $\beta_{bm:D}$ in
   a computable way. Report the predicted and simulated SE ratio by
   design. This would give Section 3.4 a number.
3. Report the carryover effect on a baseline-free scale (comment 4.3).
   Report bias and $\kappa$ by architecture (comment 4.4).
4. Restrict the recommendation in Section 4.2 to comparisons that do
   not involve blinding or expectancy. Say explicitly that design
   comparisons among open-label-first designs need the placebo
   component.

## 8. Scope of this review

- **Read in full.** `report.Rmd` (1,552 lines, including appendices),
  `README.md`, `01-run-grid.R`, `02-deviations.md`, STYLE.md, the
  relevant parts of NOTATION.md, `docs/61` (Sections 1 to 3),
  `buildSigma()`/`generateData()` at commit ef8d582, and the diff to
  HEAD.
- **Read in part.** `lme_analysis()` (formula and fallback), paper 01
  (title; search for 'mean moderation'), paper 11 (search for the
  one-Gaussian-law statements), paper 06 (search for the
  omitted-variable identity), the CR2 check summary, and `report.tex`
  (search for the numbers).
- **Executed.**
  - Every inline and hard-coded number, from both the paper 13 and
    paper 14 grids.
  - MCSE and $z$ of every carryover cost and architecture difference.
  - Pooled Type I error with $z$.
  - Model-to-empirical SE ratios.
  - Feasibility bounds for all designs, half-lives and both
    `c.cfct` values.
  - Positive definiteness of all production cells under
    `stopOnNonPD = TRUE`.
  - Coupling vectors and the placebo variance share.
  - The notation linter and the STYLE grep.
- **Not done.**
  - No simulation was run. The effect of `c.cfct = 0.2` and of CR2 on
    paper 13's power is unknown.
  - The document was not knitted. The knit failure in comment 4.1 is
    inferred from the code.
  - The CAR1-fallback frequency was not measured.
  - Paper 06's $+0.51$ bias figure and its 'PB or TV' conclusion
    (Section 4.3) were not checked against paper 06's current text.
  - Paper 14 was not reviewed beyond its data file.
  - No literature search was run.
- **No prior referee report exists for paper 13.** None was found in
  the manuscript directory, so there are no earlier comments to track.

## Readiness verdict

**Grade: D (needs new or rerun simulation).** The paper's own
simulation has never been written up. The text reports another paper's
run, and the stored run predates the `c.cfct` correction and the CR2
inference standard. The grade borders on E, because Sections 3 and 5
must be rewritten in full.

**Blockers, in order.**

1. Rerun the grid under the current DGP: `c.cfct = 0.2`,
   `stopOnNonPD = TRUE`, common seeds across half-lives, E2 and E2cr2,
   at least 1,000 replicates, with bias, empirical SE, $\kappa$, MCSE
   and the CAR1-fallback count stored.
2. Rewire the Rmd to the drop-PB output. Rewrite the Abstract, Section
   3, Section 5 and the captions from it. Remove all hard-coded
   results and the paper-14 conclusions (power ceiling, CO
   insensitivity, mean moderation 'does not' lose power).
3. Re-anchor the framing. Drop the 'paper 01's central finding'
   premise. Compare architectures on a baseline-free scale. Correct
   Appendix A to agree with paper 11's one-Gaussian-law statement, and
   report the mean-moderation bias under E2.
4. Decide on scope: merge with paper 14 or make it a supplement, and
   restrict the design-comparison recommendation given that $\eta$ is
   inert.
5. Copy-edit to submission standard. Fix the drop-TV remnants, the
   four missing bibliography keys, the pre-exposure case of the
   correlation, NOTATION, STYLE and US spelling, add the glossary,
   front matter and an exploratory statement, and update Appendix C
   for the positive-definite stop.

**Effort.** About 2 to 4 hours of compute on the simulation host for
blocker 1 (6,871 s locally at 500 replicates per cell before CR2). Then
about 2 to 3 working days for blockers 2, 3 and 5. Blocker 4 is a
decision. A merger with paper 14 would add about one day and remove a
second manuscript's worth of work.
