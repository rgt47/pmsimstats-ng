# Referee report: paper 14, two-component DGP without natural history
*2026-10-29 16:52 PDT*

**Manuscript.** 'Data generation machinery for N-of-1 clinical trial
simulations under a two-component response decomposition: dropping the
natural-history component'
(`analysis/report/14-two-component-drop-tv/report.Rmd`, drafted
2026-09-09, last substantive edit 2026-10-08 (two inserted correction
notes, committed in e077526a8 on 2026-10-09); rendered `report.pdf` of
2026-10-08).

**Submitted as.** No cover letter, title page or article type exists.
The bundled `statistics-in-medicine.csl` implies *Statistics in
Medicine*. The article type is not stated.

**Reviewer.** pmsimstats team, internal referee pass.

**Evidence labels.** VERIFIED: I executed code and observed the result.
INSPECTED: I read the source and confirmed it. SUSPECTED: consistent
with what I read but not confirmed. UNVERIFIED: not checked.

**Scratch scripts** (not part of the repository), all in
`/private/tmp/claude-501/-Users-zenn-prj-res-36-pmsimstats-ng-pmsimstats-ng/12d860d5-0c62-4f77-a07f-554bceae2b3e/scratchpad/`:
`p14-check1.R` (every inline number, carryover costs with MCSE, null
rejection rates), `p14-check2.R` (bias against the DGP's true value,
ratio of empirical to model standard error), `p14-check3.R`
(feasibility bounds under the run's and the current settings, positive
definiteness of every grid cell, Appendix C.2 eigenvalues).

## 1. Summary

The manuscript repeats a mean-shift versus correlation encoding
comparison under a reduced DGP that keeps the drug-response (BR) and
placebo-response (PB) components and drops the natural-history
component (TV). The grid is 2 encodings by 3 designs (CO, Hybrid,
OL+BDC) by $N \in \{35, 70\}$ by $c_{bm} \in \{0, 0.30, 0.45\}$ by
$t_{1/2} \in \{0, 0.5, 1\}$ week, 500 replicates per cell, with the
plain model `Sx ~ bm + t + Dbc + bm:Dbc`, random intercept, `corCAR1`
and model-based standard errors. The headline results are that the
correlation encoding loses 0.194 (Hybrid) and 0.196 (OL+BDC) of power
from $t_{1/2} = 0$ to 1 at $c_{bm} = 0.30$, $N = 70$, against 0.040
and 0.078 for the mean-shift encoding; that the reference cell
$c_{bm} = 0.45$ is ceiling-bound; and that null rejection is 0.061
(CO), 0.039 (Hybrid) and 0.021 (OL+BDC). The paper argues that
dropping TV is 'unconditionally safe' because no biomarker-TV parameter
exists, and Appendices B and C give the reduced correlation matrix and
a Schur-complement feasibility bound.

Every inline number reproduces from the stored file. The paper's
central argument, however, no longer holds. Since 2026-10-22 the
package has three biomarker-TV parameters, and the project has already
run biomarker-TV cells. The finding the paper claims to reproduce is
no longer in paper 01, and paper 11's current-DGP rerun now carries it.
Under the fixed analysis model the mean-shift estimator is biased away
from zero by up to 22 percent at $t_{1/2} > 0$. Part of the
mean-shift arm's apparent robustness to carryover is that bias. The
simulation uses `c.cfct = 0.1` (the `docs/61` defect), silent repair
enabled, and model-based inference. Several passages were copied from
paper 13 without being reversed, so the manuscript argues both sides
of its own choice.

## 2. Overall assessment and recommendation

**Recommendation: major revision with a rerun. Consider merging with
paper 13.**

The appendices are the most durable part. The Schur-complement bound
(C.1), the compound-symmetry spectrum (C.2) and the 2- versus
3-component bounds (C.4) are correct and reproduce exactly. The
main-text argument rests on three premises that have each lapsed since
the draft:

- the absence of a biomarker-TV parameter;
- a paper 01 result that paper 01 no longer reports;
- the claim that a fixed analysis model attributes every difference to
  the DGP.

The numbers were generated under a DGP that differs from the current
extended DGP in the cross-component correlation, and they were
analyzed with a test that is not calibrated by design. Papers 13 and
14 share one driver, one grid, one seed scheme and most of their text.
They report the same qualitative result (paper 13: correlation-encoding
costs 0.228 and 0.182, mean-shift costs 0.078 and 0.086, null rates
0.065, 0.042 and 0.033; VERIFIED from `13-drop-pb-grid.rds`). Two
journal papers would not survive a redundancy check. One paper, or a
section of paper 06, comparing both reductions against a matched
three-component run is the natural unit.

## 3. Significance and novelty

- **Within the field.** Whether a simulation DGP needs a separate
  natural-history component is a modeling-convention question. A
  statistical journal referee would ask what decision a reader makes
  differently after reading this paper. The answer offered (use two
  components for comparative work, three for design estimates) is
  plausible. However, it is supported by a single 108-cell grid with
  no matched three-component arm (Section 3.4 concedes this). The
  algebraic content (Appendix C) is standard: a Schur complement and
  the spectrum of a Kronecker sum. INSPECTED.
- **Within the compendium.** The covariance-feasibility material is
  in the compendium and `docs/34`, `docs/35` and `docs/65`. The
  omitted-variable identity belongs to paper 06. The encoding
  comparison under the current DGP now belongs to paper 11 (commit
  6ea37805e, common random numbers, current DGP). It finds the same
  direction for a different reason: 'only the correlation channel
  decays off drug.' The new contribution of paper 14 is therefore the
  statement that deleting TV changes no qualitative result. That
  statement is one table in a merged paper. INSPECTED.
- **Exploratory framing.** The manuscript makes no pre-registration
  or confirmatory claim. INSPECTED.

## 4. Major comments

### 4.1 The 'unconditional safety' premise is false for the current package

- **Location.** Abstract (Methods); Section 1.1; Section 2.2.3;
  Section 4.1; Section 4.3; Section 4.5 ('Safety'); Appendix B.1 and
  B.6; Conclusions 2; README.
- **Problem.** The paper's main argument for preferring this reduction
  is that 'there is no biomarker-TV counterpart anywhere in the
  package, the drivers, or any implementation collection.' That is no
  longer true. `R/generateData.R` now provides `c.bm.tv` (biomarker-TV
  correlation, commit ccb6a5060, 2026-10-22), `b.bm.tv` (biomarker-TV
  mean shift, 8d0d004fa, 2026-10-22) and `b.bm.tv.trend` (0635b0116,
  2026-10-27). INSPECTED. The 'untested half of paper 06's claim'
  (Section 4.3) has since been exercised:
  - paper 02's block S12 runs `prog-TV-const` and `prog-TV-grow` null
    cells;
  - `docs/64` reports their null rejection;
  - `docs/65` derives the feasibility bound for a biomarker-TV
    correlation and finds it identical to the biomarker-PB bound
    ('prog-PB and prog-TV-const have identical ceilings').

  All three are INSPECTED. The structural asymmetry between the two
  reductions has therefore disappeared. Both deleted blocks are now
  zero only by default.
- **Why it matters.** The paper's 'Safety' axis, half of its case for
  dropping TV rather than PB, and Conclusion 2 rest on this premise.
  Section 4.3 states that 'the safety of this reduction rests on a gap
  in the software.' That gap is now closed.
- **Remedy.** Restate both reductions as conditional on the deleted
  component's biomarker correlation being zero. Cite `docs/65` for the
  shared feasibility bound. Cite the biomarker-TV null cells for the
  empirical check. Drop the 'Safety' axis from Section 4.5 or recast
  it as a symmetric condition.

### 4.2 The mean-shift estimator is biased under the fixed analysis model

- **Location.** Section 3.1; Section 4.5; Conclusions 1; Appendix A.4
  ('every difference reported in Section 3 is attributable to the
  data-generating process alone').
- **Problem.** The mean-shift DGP adds $c_{bm} z_i \sigma_{BR}$ to BR
  at on-drug visits only (`R/generateData.R` at the run's commit
  955c604, the `mean_moderation` branch). INSPECTED. The analysis
  regresses on $D_{bc}$, which is positive off drug when $t_{1/2} > 0$.
  The analysis model is therefore correctly specified for the
  correlation encoding, whose off-drug correlation decays as $D_{bc}$
  does. It is misspecified for the mean-shift encoding. The true
  coefficient per raw biomarker unit is
  $-c_{bm}\sigma_{BR}/\sigma_{bm} = -0.30 \times 8/15.36 = -0.156$ at
  $c_{bm} = 0.30$. The stored mean estimates match it at $t_{1/2} = 0$
  in every cell ($|z| \le 1.8$). At $t_{1/2} = 1$ the mean-shift
  estimates are biased away from zero by 13.6 to 22.2 percent under
  Hybrid and OL+BDC ($z$ from $-5.9$ to $-16.6$ against the MCSE of
  the mean). The correlation-encoding estimates stay unbiased
  ($|z| \le 1.6$ in all 36 cells). VERIFIED (`p14-check2.R`).
- **Why it matters.** A biased-away estimate inflates power. Part of
  the mean-shift arm's small carryover cost is therefore the bias, not
  robustness. The comparison confounds the encodings' exposure profile
  with analysis misspecification for one arm. Paper 11's current text
  locates the divergence in the exposure profile ('only the
  correlation channel decays off drug'). CLAUDE.md
  ('Key architectural feature') gives the same reading. Appendix A.4's
  claim is false as stated.
- **Remedy.** Report bias and relative bias per cell with MCSE. Either
  analyze the mean-shift arm with the binary indicator (its correctly
  specified model), or define the estimand explicitly as on-drug
  moderation and report each arm's bias for it. Rewrite A.4 and A.2 to
  match paper 11's 'one Gaussian law' framing: both encodings set
  $E[BR \mid b]$, and only the correlation encoding lowers
  $\mathrm{Var}(BR \mid b)$.

### 4.3 The paper reproduces a finding that paper 01 no longer reports

- **Location.** Abstract (Background, Results); Introduction ('The
  companion paper 01 takes up the first choice ... shows the two are
  not interchangeable'); Section 3.1 ('paper 01's central finding,
  reproduced'); Section 3.4; Section 2.4 ('those of paper 01').
- **Problem.** Paper 01's current title is 'Simulating
  predictive-biomarker interactions in aggregated N-of-1 trials: an
  extended data-generating process building on Hendrickson et al.' Its
  text contains no mean-shift versus correlation comparison. A search
  for 'mean moderation' returns no hit, and the word 'moderation'
  occurs 7 times, all for the on-drug estimand. INSPECTED. CLAUDE.md
  records that the cross-architecture comparison was removed on
  2026-10-07 and that its figure 'is not currently narrated in any live
  manuscript'. Section 3.4 compares absolute power against 'paper 01's
  published figures'. Those figures now exist only in the archived
  pre-restatement draft.
- **Why it matters.** The paper's stated purpose is to test whether a
  paper 01 result survives the reduction. That result is no longer in
  paper 01.
- **Remedy.** Re-anchor on paper 11's rerun (current DGP, common
  random numbers), or on a matched three-component arm run under the
  same driver. The latter also resolves Section 3.4's admitted
  confound.

### 4.4 The DGP and analysis do not match the current extended DGP

- **Location.** Section 2; Appendix B.3; Appendix C; Reproducibility.
- **Problem.** The stored run (`14-drop-tv-grid.rds`, git 955c604,
  2026-09-08) and the driver
  (`analysis/scripts/two-component-dgp/01-run-grid.R`) use the
  following settings. INSPECTED except where marked.
  - `c.cf1t = 0.2`, `c.cfct = 0.1`. This is the defect of `docs/61`.
    The extended DGP sets one constant, $c_\times = c_1 = 0.2$.
  - `makePositiveDefinite = TRUE` under the pre-2026-10-26 code, which
    repairs silently. Impact check: every path matrix at
    $c_{bm} = 0.45$ is positive definite under the run's settings
    (smallest correlation-matrix eigenvalue 0.014, OL+BDC,
    $t_{1/2} = 0$) and under $c_\times = 0.2$ (0.047). No stored cell
    depended on a repair. VERIFIED (`p14-check3.R`). Appendix C's
    motivating sentence ('Simulation drivers in this program repair a
    non-positive-definite matrix silently') is now false.
  - Decaying correlation, a pre-exposure guard on both the
    correlation and $D_{bc}$, and carryover that decays from the last
    on-drug mean. All three match the extended DGP.
  - `useDE = FALSE`, raw (uncentered) biomarker and model-based
    standard errors. The project's recommended analysis is F1 (centered
    biomarker, `De` and `bmc:De` terms, CR2 with Satterthwaite df).
  - If `corCAR1` fails, `lme_analysis()` refits silently without the
    correlation structure, and the driver counts the fit as converged.
    The claim 'All 108 cells converged completely' therefore cannot be
    checked from the stored results, which carry no fallback flag.
    How often the fallback occurred is UNVERIFIED.
  - `set.seed(SEED + i)` per grid row. The $t_{1/2} = 0$ and
    $t_{1/2} = 1$ cells use different streams, so every 'cost' is an
    unpaired difference.
- **Impact on the bounds.** Under $c_\times = 0.2$ the two-component
  bounds at $t_{1/2} = 0$ rise from 0.469 to 0.487 (Hybrid) and from
  0.464 to 0.481 (OL+BDC). The reduction's gain is 0.007 to 0.009
  rather than 0.007 to 0.013. VERIFIED. Section C.3's statement that
  'the operative cross-component quantity is the difference
  $c_1 - c_\times$' describes the defect. Under the extended DGP that
  difference is zero, and the smallest compound-symmetry eigenvalue is
  $1 - \rho$.
- **Remedy.** Rerun under the current package defaults with
  $c_\times = 0.2$ and `stopOnNonPD = TRUE`. Add F1 with CR2 alongside
  the plain model. Use common seeds across half-lives, which gives
  paired costs. Record a fallback flag per replicate. The README's
  production estimate is about 75 minutes on a laptop; adding CR2
  roughly doubles it (SUSPECTED).

### 4.5 Design-dependent Type I error is reported but not addressed

- **Location.** Section 3.3; Abstract; Conclusions 4; Section 4.4.
- **Problem.** The pooled rejection rates reproduce: CO 0.061, Hybrid
  0.0388, OL+BDC 0.0208, 12 cells each. VERIFIED. By cell they range
  from 0.012 to 0.082. CO at $N = 70$ is liberal (0.067 pooled over 6
  cells). The ratio of empirical standard deviation to mean model
  standard error is 1.04 (CO), 0.94 (Hybrid) and 0.86 (OL+BDC) at the
  null. The model-based standard error therefore overstates
  uncertainty in the designs with open-label visits. VERIFIED. A
  plausible mechanism is that the PB standard deviation is 10 at
  open-label visits and 5 at blinded visits. The outcome is then
  heteroscedastic in a way that tracks drug status, and a homoscedastic
  working model misweights it. CO has uniform expectancy and is
  calibrated. SUSPECTED. The 'monotone in off-drug information'
  explanation in Section 3.3 is not supported by any analysis. The
  paper calls a CR2 correction 'a reasonable conjecture and an
  untested one'. The project has since tested it. Commit 8e810ca9d
  (2026-10-29) reports CR2 null rejection of 0.048 to 0.056 in OL+BDC
  under papers 02 and 09's drivers. The model-based test gives 0.027
  there. INSPECTED.
- **Why it matters.** Every cross-design statement in the paper is,
  in the authors' words, 'descriptive rather than equal-alpha'. A
  calibrated test makes those statements interpretable.
- **Remedy.** Report CR2 (Satterthwaite) rejection beside the
  model-based rate. Report $\kappa$ per design, with the STYLE
  definition (mean model standard error over empirical standard
  deviation). Cite 8e810ca9d's check. Put a Monte Carlo band beside
  every null rate.

### 4.6 Text copied from paper 13 contradicts this paper

- **Location and problem.** INSPECTED unless marked.
  - Introduction, line 174: 'reduced to two components, BR and TV,
    dropping the modelled placebo-belief term.' This paper keeps BR and
    PB.
  - Line 180: 'Section 1.1 sets out why dropping PB rather than TV is
    the reduction that aligns with established practice.' Section 1.1
    argues the opposite.
  - Section 3.4, line 602: 'PB carries a standard deviation of 10
    ... Removing it ...'. This paragraph explains the power gain by
    removing PB. It also claims the removed share is 'largest under
    designs with the most open-label time.' For TV removal the reverse
    holds. TV's share of component variance is 37.9 percent at
    open-label visits and 52.9 percent at blinded visits, because the
    PB standard deviation halves there. VERIFIED (arithmetic from
    `extracted_rp`).
  - Section 2.2.3 states that 'the share removed is the same under
    every design.' This is true of the absolute variance removed, not
    of the share. VERIFIED.
  - Section 4.2, line 681: the case for the reduction includes
    'conformity with the modelling conventions of the crossover and
    N-of-1 literature (Section 1.1).' This paper's Section 4.5 says the
    reduction departs from that convention.
  - Section 4.4 repeats two Limitations paragraphs verbatim: the
    calibration paragraph (lines 717 and 730) and the Type I paragraph
    (lines 721 and 741).
  - Driver header and `02-deviations.md` describe paper 13 only.
    README: 'Paper 13 compares', 'Paper 13's manuscript reads every
    quoted figure', 'Section 2.4' (the recursion is Section 2.3), and
    'the omitted-variable identity paper 13 relies on'.
- **Why it matters.** A referee reading the Introduction would take
  the paper to be the other paper.
- **Remedy.** Rewrite these passages from this paper's position. If
  the papers are merged (Section 2), most of them go away.

### 4.7 Overstated or unsupported quantitative statements

All items VERIFIED from the stored file unless marked.

- **'Mean moderation does not' lose power (Conclusions 1).** The
  mean-shift arm loses 0.078 under OL+BDC at $N = 70$ and at $N = 35$
  ($z = 2.7$ and 2.5 on the unpaired difference). That loss is net of
  the upward bias of comment 4.2.
- **'Both mean-arm costs are within two or three of [the MCSE]'.**
  This compares a difference with a single-cell MCSE. The MCSE of the
  unpaired difference is about 0.028, not about 0.022.
- **'Sixfold' understatement (Abstract, 3.2, Conclusions 3).** This
  holds for Hybrid (0.194 against 0.032). Under OL+BDC the cost falls
  from 0.196 to 0.098, a factor of 2.0.
- **CO: 'neither arm loses power to carryover'.** This holds at
  $N = 70$. At $N = 35$ the correlation arm loses 0.058
  ($c_{bm} = 0.30$, $z = 1.8$) and 0.056 ($c_{bm} = 0.45$, $z = 2.4$).
- **Section 2.3's recursion check.** This is 'at most 0.005 at
  $N = 70$, $c_{bm} = 0.45$, $t_{1/2} = 1$'. It was run in cells whose
  power is 0.89 to 0.998, which Section 3.2 calls ceiling-bound and
  uninformative. No artifact of that check is stored. UNVERIFIED.
- **'+0.51 at a biomarker-PB correlation of 0.45' (Section 4.3).**
  This attribution to paper 06 does not appear in paper 06's current
  text. Paper 06 now predicts a negative displacement on the Hybrid
  design. INSPECTED.
- **Appendix C.4 second and third observations.** These are the
  bound of 0.81 at $n = 4$ against 0.58 at $n = 12$, and 0.45 against
  0.35 for the index-gap exponent. The provenance paragraph says they
  come from 'a reimplemented matrix builder'. No script is cited or
  stored. UNVERIFIED.

### 4.8 The HAVE_SIM guard does not protect the document

- **Location.** Setup chunk and Section 3.
- **Problem.** If the RDS file is absent, `S` is `NULL`. The inline
  `max(S$mc_se_power[...])` in Section 3 and the three table chunks
  then fail, so the document does not 'knit either way'. INSPECTED.
- **Remedy.** Remove the guard, since the file exists, or guard every
  use. A production manuscript should fail loudly on missing data.

## 5. Minor comments

1. **Title.** 'Machinery' is on the STYLE avoid list. The title
   states a topic, not a finding.
2. **No glossary** (STYLE section 4). No Abbreviations section. CO,
   OL+BDC and Hybrid are not spelled out on first use.
3. **Retired labels.** 'Architecture A/B' appears in the section
   headings, the grid table and `arch_label()`. NOTATION Part 3
   retires the letters.
4. **STYLE terminology.** The STYLE grep (coupling, ceiling, graded,
   anchored, and so on) gives 45 hits, most in prose. 'Covariance
   moderation' and 'mean moderation' should be 'the correlation
   encoding' and 'the mean-shift encoding'. Use 'natural-history
   component' and 'placebo-response component' for TV and PB.
5. **Notation conflicts with NOTATION.md.** The linter reported the
   file clean but missed these:
   - $S_{x,ip}$ should be $Y_{it}$.
   - $p$ is used as the visit index. NOTATION reserves $t$ for visits
     and $p$ or $q$ for paths.
   - $\mathrm{mg}(t; a, d, r)$ should be $G(t_{od}; m, d, r)$.
   - $\phi$ is used for the correlation profile, including its on-drug
     value of 1. That is $g_t$. $\phi$ is the decay family.
   - $c_{bm}$ is used as a multiplier in Section 2.2.1. That breaks
     rule 2, which reserves $\beta_{bm}$ for that role.
   - $A$, $B$ in C.2 collide with $A_c$ in B.2. $\lambda$ is both the
     decay rate and an eigenvalue.
   - The sign convention (rule 5) is never stated. The true
     coefficients are negative, about $-0.156$ and $-0.234$ per raw
     biomarker unit.
6. **True value never stated.** The paper reports no bias, empirical
   standard error or coverage (ADEMP performance measures). The
   summary file holds `beta`, `empirical_se` and `mean_model_se`.
7. **MCSE.** Table 1 and Table 2 give power without MCSE. The 'Cost'
   column needs its MCSE (about 0.024 to 0.032).
8. **Prose against STYLE.** The draft has 41 semicolons, many
   sentences over 30 words, and evaluative phrasing ('worth stating',
   'striking feature', 'rather than pretending'). Two embedded
   'Corrected 2026-10-08' notes contradict nothing now but cannot
   appear in a submission.
9. **British spelling.** modelled (lines 174, 227, 726, 789),
   modelling (652, 681), favourable (694), favours (775, 785, 794),
   neighbours (1410).
10. **Section numbers typed into headings** ('\section{1.
    Introduction}') duplicate `number_sections: true`. The rendered
    PDF should be checked for double numbering. SUSPECTED.
11. **The Kronecker identity (C.2) is labeled 'Kronecker sum'.** It
    is $A \otimes I_n + B \otimes (J_n - I_n)$, which is not a
    Kronecker sum in the standard sense
    ($A \otimes I + I \otimes B$). The eigenvalues 7.7, 5.0, 0.5 and
    0.2 (multiplicity 14) reproduce. VERIFIED.
12. **Appendix B.5 worked example** uses $\rho_{br} = 0.6$,
    $\rho_{pb} = 0.4$ and $c_\times = 0.1$. None of these is the
    project value. Use $\rho = 0.7$, $c_1 = c_\times = 0.2$.
13. **Companion citations.** Papers 01, 02, 06 and 13 are named in
    prose but never cited with a key. The bib has no entry for 01 or
    13. The `pmsimstats-paper11` title ('A combined data-generating
    architecture ...') is stale. The current title begins 'A
    dual-channel data-generating process'. `references.bib` is
    byte-identical to paper 13's and has about 50 entries, 11 of them
    cited.
14. **Rendering.** `report.pdf` and `report.tex` sit in the directory.
    There is no staged copy in `share/` under `tools/render.sh`.
    SUSPECTED; only the directory was checked.

## 6. Missing and questionable references

| Reference | Status | Relevance |
|---|---|---|
| Jones and Kenward, *Design and Analysis of Cross-Over Trials*, 3rd ed., 2014 | In bib, not cited | Source for the 'time trend is standard' claim in Sections 1.1, 4.4 and 4.5 |
| Araujo, Julious and Senn, *PLoS One* 2016 | In bib, not cited | Mixed-model conventions for N-of-1 series |
| Zucker et al. 1997/2010 | In bib, not cited | N-of-1 combining models, with period terms |
| Morris, White and Crowther, *Stat Med* 2019 | In bib, not cited | ADEMP reporting of bias, MCSE and coverage |
| Pustejovsky and Tipton, *J Bus Econ Stat* 2018 (CR2) | Missing | Required once CR2 is added (comment 4.5) |
| Ernst and Resch, *BMJ* 1995, 'Concept of true and perceived placebo effects' | Missing; details UNVERIFIED | Separating placebo response from natural history |
| Hrobjartsson and Gotzsche, *NEJM* 2001 | Missing; details UNVERIFIED | No-treatment arms; separability claim in Section 4.5 |
| Colagiuri, *Clinical Trials* 2010 | Missing; details UNVERIFIED | Expectancy in blinded trials and balanced-placebo designs |
| `docs/64`, `docs/65`; paper 02 block S12 | Missing (internal) | Biomarker-TV cells and shared feasibility bound (comment 4.1) |
| pmsimstats paper 11 (current) | Cited with a stale title | Current-DGP encoding comparison (comments 4.2, 4.3) |
| Paper 06 '+0.51' figure | Questionable | Not in paper 06's current text |

The placebo-literature entries are from the reviewer's knowledge. They
were not checked against the originals. The claim in Section 4.5 that
'the placebo literature regards the belief and natural-history channels
as separable only under three-arm or balanced-placebo designs' has no
citation and needs one.

## 7. Suggestions for strengthening

1. Merge papers 13 and 14 into one paper. Use one driver run with
   three arms (all three components, BR+PB, BR+TV), common random
   numbers across arms and half-lives, and the current extended DGP.
   This removes the redundancy and the confound of Section 3.4, and it
   lets the paper state numerically what each component contributes.
2. Report a full ADEMP table: the true value, bias, empirical
   standard error, model and CR2 standard errors, $\kappa$, coverage,
   power with MCSE, and paired cost with paired MCSE.
3. Analyze each encoding under its correctly specified exposure
   coding as well as the common one. The difference between the two
   is the misspecification cost. It is reportable in its own right and
   connects to paper 11.
4. Add one prognostic cell per reduction: a biomarker-TV correlation
   for the full and BR+TV arms, and a biomarker-PB correlation for the
   full and BR+PB arms. Each must sit inside the `docs/65` bound. This
   shows directly that a reduction is safe only while the deleted
   block's biomarker correlation is zero.
5. Keep Appendix C, rewritten at $c_\times = c_1$. Under the
   separable form the component block is $A(\rho) \otimes R_K$, with
   eigenvalues equal to products. That gives a cleaner closed form than
   the compound-symmetry case and covers AR(1) (compendium; `docs/35`).
   INSPECTED in outline only.

## 8. Scope of this review

- **Read in full.** `report.Rmd` (all 1,510 lines and chunks),
  `README.md`, the driver `01-run-grid.R`, `02-deviations.md`,
  STYLE.md, and the referee report for paper 12 used as the format
  model.
- **Read in part.**
  - NOTATION.md: Parts 1 and 3, symbol rows.
  - `R/generateData.R`: current version, parameters and stop path.
  - `R/generateData.R` and `R/lme_analysis.R` at the run's commit
    955c604: means, correlation build, repair, mean-shift branch,
    $D_{bc}$ coding and fit fallback.
  - `docs/61` (header).
  - `docs/64`, `docs/65`, `docs/60` (grep hits for biomarker-TV).
  - Paper 06 (abstract, identity, future work).
  - Paper 01 (title; grep for the encodings).
  - Paper 11 (title, commit message).
  - Commit messages of 8e810ca9d and 6ea37805e.
- **Executed.**
  - Reproduction of every inline number in the Abstract, Sections 3.1
    to 3.3 and Conclusions from `14-drop-tv-grid.rds`. All matched.
  - MCSE of unpaired costs.
  - Pooled null rates by design, encoding and $N$.
  - Bias against the DGP's true value in all 72 alternative cells.
  - Ratio of empirical to model standard error.
  - Feasibility bounds for 2 and 3 components at $c_\times = 0.1$
    (reproduces Table C.4 exactly) and at $c_\times = 0.2$.
  - Positive definiteness of every path matrix at $c_{bm} = 0.45$.
  - The C.2 spectrum.
  - Paper 13's headline costs, for the redundancy comparison.
- **Not done.**
  - No new simulation was run.
  - The fallback-fit frequency was not measured, because the stored
    results carry no flag.
  - The Section 2.3 recursion check and the C.4 occasion-count and
    index-gap figures were not reproduced, because no script exists.
  - The PDF was not re-rendered.
  - The placebo-literature references were not read in the original.
  - Paper 13's text was read only through grep, for the copied
    passages.
- **No prior referee report exists for paper 14.** None was found in
  the manuscript directory. `docs/31` Part V summarizes the paper and
  repeats its 'safety' argument. That document needs the same
  correction as comment 4.1.

## Readiness verdict

**Grade: D (needs new or rerun simulation).**

The draft is complete and its numbers reproduce from the stored file.
The stored run, however, predates the `c.cfct` correction, uses
uncalibrated model-based inference, and pairs one encoding with a
misspecified analysis. Its central argument has been overtaken by code
added on 2026-10-22.

**Ordered blockers.**

1. Rewrite the premise. Biomarker-TV parameters (`c.bm.tv`,
   `b.bm.tv`, `b.bm.tv.trend`) now exist, so the 'unconditional
   safety' asymmetry with paper 13 is gone (comment 4.1). Decide
   whether papers 13 and 14 merge before any further work.
2. Rerun under the current extended DGP: $c_\times = 0.2$,
   `stopOnNonPD = TRUE`, a fallback flag, and common seeds across
   half-lives. Add F1 with CR2 inference and report $\kappa$ and
   calibrated null rates (comments 4.4, 4.5).
3. Fix the encoding comparison. Report bias and either analyze the
   mean-shift arm under its correct exposure coding or define the
   estimand. Re-anchor on paper 11 or on a matched three-component
   arm, not on paper 01 (comments 4.2, 4.3).
4. Remove the copied paper 13 passages, the duplicated Limitations,
   the incorrect variance-share statements and the README errors
   (comment 4.6). Correct the overstated statements of comment 4.7.
5. Bring the manuscript to STYLE and NOTATION: add a glossary,
   retire the A/B labels, add MCSE on every cost, add companion
   citation keys and the missing references.

**Effort estimate.**

- Rerun: about 1.5 to 3 hours of compute for 108 cells with an added
  CR2 fit. Add about 1 hour for a matched three-component arm. Merging
  with paper 13 doubles this.
- Driver and analysis changes: 0.5 to 1 day.
- Manuscript rewrite: 2 to 3 days separately, or 3 to 4 days for a
  merged paper.
- Appendix C rewrite at $c_\times = c_1$: 0.5 day.
