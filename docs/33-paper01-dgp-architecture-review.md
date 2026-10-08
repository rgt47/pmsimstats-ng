# Referee Review of Paper 01: Mean and Covariance Moderation Under Carryover {.unlisted .unnumbered}
*2026-09-30 08:06 PDT*

**Author:** pmsimstats team

**Scope.** This whitepaper reviews the three R Markdown sources in
`analysis/report/01-dgp-mean-moderation-vs-mvn/`: the manuscript
`report.Rmd` (2,292 lines), the standalone `appendices.Rmd`, and the
supplementary key-points summary `bullets.Rmd`.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## 1. Evidence basis

Each finding below carries one of four evidence labels.

- **Verified:** the claim was recomputed or checked against data or
  code output in this review.
- **Inspected:** the relevant source code or document was read and
  the claim confirmed by reading, without execution.
- **Inferred:** the claim follows from verified facts but was not
  checked directly.
- **Unverified:** the claim is a hypothesis requiring a new run.

The following checks were performed.

- All relative losses, within-architecture $z$ statistics, and
  difference-of-differences $z$ statistics in Section 3.1 were
  recomputed by hand from the reported proportions. All match the
  manuscript. (Verified.)
- The Section 3.1 power tables were compared against
  `analysis/data/quick-sim/01-dgp-summary.txt`. All 18 alternative
  cells and all 18 null cells match. (Verified.)
- The replicate-level file `analysis/data/quick-sim/01-dgp-replicates.rds`
  (36,000 rows, with `beta`, `betaSE`, `p`) was reanalyzed for the
  ratio of mean model-based SE to empirical SD, and for size-adjusted
  power. (Verified; the script is reproduced in Appendix 1.)
- `R/generateData.R`, `R/lme_analysis.R`, and the driver
  `analysis/scripts/quick-sim/01-dgp-prototype.R` were read, and the
  git history of the package source and the result files was
  examined. (Inspected.)
- The rendered `report.pdf` table of contents was extracted with
  `pdftotext`. (Verified.)
- `tools/notation-lint.pl` was run over all three files. All three
  pass. The linter does not mechanize NOTATION rules 2 and 5, which
  are assessed manually in Section 5. (Verified.)

No new simulation was run. The size-adjusted power figures in
Section 2.4 and the spacing-based design explanation in Section 2.5
are reanalyses of existing output and should be confirmed after the
re-run recommended in Section 7.

## 2. Critical findings

These findings bear directly on the headline results, the stated
mechanism, or both.

### 2.1 Section 3 results predate the carryover correction

**Evidence: verified (file dates and git history).**

The result files `01-dgp-summary.txt` and `01-dgp-replicates.rds`
are dated 2026-08-08. The carryover recursion in `buildSigma()` was
replaced by the anchored form in commit `955c604` of 2026-09-08
(`R/generateData.R:277-290`). Every number in Section 3 of the
manuscript was therefore produced under the recursive form that
Appendix B.6 describes as having been "corrected in September 2026".

Two passages are consequently inconsistent with each other and with
the reported results.

- Section 2.2.4 (lines 577-579) states that the carryover line "in
  the current package is character-for-character the published one".
  That was true of the code that produced Section 3, but is no longer
  true of the package.
- Appendix B.6 and B.7 describe `mean` and `covar` as using the
  anchored form, and do not disclose that the reported power figures
  were not produced under it.

The corrected form retains more residual effect at the second and
later off-drug occasions of each run (Appendix B.6's own worked
example gives $0.25\,\mu_{BR,1}$ against $0.125\,\mu_{BR,1}$). The
carryover-induced power losses are therefore likely to be larger
under the current code, particularly in the designs with consecutive
off-drug occasions (Hybrid and OL+BDC). (Inferred; the magnitude is
unverified.)

**Recommendation.** Re-run the Section 3.1 factorial on the current
package before any further revision of the Results, and state in the
Reproducibility section the commit that produced the reported
figures.

### 2.2 The conditional-variance-inflation mechanism is directionally inverted

**Evidence: inspected (code) and verified (algebra).**

Sections 2.2.3 and 3.2 attribute the architecture-specific loss to
the identity

$$
\mathrm{Var}(BR_{it} \mid B_i, D_{it}) =
\sigma_{BR}^2\bigl(1 - c_{bm}^2 D_{bc,it}^2\bigr),
$$

arguing that "as carryover erodes $D_{bc,it}$ toward zero off drug"
the conditional variance rises toward $\sigma_{BR}^2$ and adds noise
(lines 504-511, 841-847).

The direction of that argument is the reverse of what the
implementation does.

- At $t_{1/2} = 0$, off-drug `Dbc` is set to exactly zero
  (`R/lme_analysis.R:112-113`), and the off-drug biomarker-response
  correlation is zero because `lambda_cor` is set to 0 and the
  decay branch is skipped (`R/generateData.R:225-228, 350`).
- For $t_{1/2} > 0$, off-drug `Dbc` equals
  $(1/2)^{t_{sd}/t_{1/2}} > 0$.

Carryover therefore *raises* off-drug $D_{bc,it}$ from zero, and
under the stated identity the off-drug conditional variance *falls*
with carryover. Off-drug occasions are at the full unconditional
variance $\sigma_{BR}^2$ precisely in the no-carryover cell. The
qualitative statement in Section 3.2 that the mechanism "needs
sustained low-$D_{bc,it}$ time to accumulate" (line 907) has the
same problem: low-$D_{bc,it}$ time is maximal at $t_{1/2} = 0$.

The identity itself is correct. What fails is the claim that
carryover acts on power through it in the stated direction.

**Recommendation.** Withdraw the variance-inflation account as the
mechanism, or re-derive it as a statement about the whole design
(for example, about how the covariance-weighted information of the
`bm:Dbc` contrast changes as off-drug occasions acquire intermediate
exposure), and check the re-derivation against the replicate data
before it is restored.

### 2.3 The Section 3.3 empirical check is misread

**Evidence: verified (summary file and replicate reanalysis).**

Section 3.3 states that covariance moderation's signal-to-noise
erosion is "driven almost entirely by SE growth", and that "mean
moderation shows no comparable SE trend" (lines 961-965). The
signal-to-noise ratios quoted are arithmetically correct, being
`mean_beta / sd_beta` from the summary file, but the interpretation
does not hold.

| Design | Architecture | Empirical SD, $t_{1/2}=0 \to 1$ | Change | Mean $\hat\beta$, $t_{1/2}=0 \to 1$ |
|---|---|---|---|---|
| Hybrid | Mean moderation | 0.0782 to 0.0938 | +20.0% | -0.235 to -0.266 |
| Hybrid | Covariance moderation | 0.0708 to 0.0840 | +18.7% | -0.231 to -0.239 |
| OL+BDC | Mean moderation | 0.0807 to 0.0956 | +18.5% | -0.231 to -0.281 |
| OL+BDC | Covariance moderation | 0.0789 to 0.0904 | +14.6% | -0.232 to -0.234 |

Table: Empirical SD and mean $\hat\beta$ by design and architecture as carryover increases

Empirical SD grows as much or more under mean moderation. What
distinguishes the architectures is the numerator: under mean
moderation $|\hat\beta|$ grows with carryover by 13% (Hybrid) and
22% (OL+BDC), with MCSE near 0.003, so the change is far outside
Monte Carlo error. Under covariance moderation $\hat\beta$ is stable.

Three consequences follow.

- The statement in Section 3.2 that the interaction coefficient
  "stays essentially unbiased under both architectures" (line 850)
  is false for mean moderation.
- Mean moderation's apparent robustness to carryover is, at least in
  part, an estimand shift. As carryover grows, the fitted `bm:Dbc`
  coefficient targets a larger quantity, which offsets the SE growth
  and holds power up. This is the mechanical consequence of the
  DGP/analysis mismatch that Section 2.2.1 itself identifies (lines
  373-384): the mean-moderation DGP places no interaction signal
  off drug, while `Dbc` is nonzero off drug whenever carryover is
  present. (The direction of the shift is verified; the explanation
  is inferred.)
- The architecture contrast in Section 3.1 therefore compares one
  architecture whose target is fixed with another whose target moves
  with carryover. Power losses under the two are not losses of the
  same kind.

### 2.4 Test calibration accounts for a large share of the architecture gap

**Evidence: verified (replicate reanalysis); size-adjusted figures
are approximate.**

The ratio of mean model-based SE (`betaSE`) to empirical SD of
$\hat\beta$ is far from uniform across cells.

| Design | Mean moderation, $t_{1/2}=0 \to 1$ | Covariance moderation, $t_{1/2}=0 \to 1$ |
|---|---|---|
| CO | 0.95 to 0.95 | 0.97 to 1.00 |
| Hybrid | 1.02 to 1.04 | 1.12 to 1.16 |
| OL+BDC | 1.10 to 1.20 | 1.12 to 1.26 |

Table: Ratio of model-based SE to empirical SD by design and architecture

Values are for the $c_{bm} = 0.45$ cells. The null cells show the
same pattern (CO 0.92 to 0.98; OL+BDC 1.12 to 1.23).

The Type I error results are consistent with this.

- All six CO null cells lie above nominal (0.061 to 0.074; mean
  0.067). CO is anti-conservative.
- Hybrid and OL+BDC null cells lie between 0.010 and 0.053, and the
  OL+BDC covariance-moderation null rate falls from 0.032 to 0.010 as
  carryover increases.

Power was recomputed using, for each cell, the 95th percentile of
$|\hat\beta / \widehat{SE}|$ in the matching null cell as the
critical value. With 1,000 null replicates this quantile is itself
noisy, so the figures below indicate magnitude only.

| Design | Architecture | Nominal loss | Size-adjusted loss |
|---|---|---|---|
| CO | Mean moderation | 0.014 | 0.043 |
| CO | Covariance moderation | 0.022 | 0.054 |
| Hybrid | Mean moderation | 0.058 | 0.021 |
| Hybrid | Covariance moderation | 0.129 | 0.050 |
| OL+BDC | Mean moderation | 0.021 | 0.028 |
| OL+BDC | Covariance moderation | 0.238 | 0.125 |

Table: Nominal and size-adjusted power loss by design and architecture

The covariance-minus-mean gap falls from 0.071 to about 0.029 in
Hybrid and from 0.217 to about 0.097 in OL+BDC, and is about 0.011
in CO. Roughly half of the reported architecture gap in the two
non-CO designs is attributable to the model-based test becoming
progressively more conservative, not to a loss of information in
the data.

The following manuscript statements are contradicted.

- Type I error is "uniformly sub-nominal" (line 927). The maximum
  is 0.074, and every CO cell exceeds 0.05.
- The miscalibration is "common-mode with respect to the
  architecture contrast" (lines 931-933). It differs by design, and
  within Hybrid it differs by architecture (1.12 against 1.02 at
  $t_{1/2} = 0$).
- "The within-replicate model SE matches it at scale" (line 944).
  The ratio spans 0.92 to 1.26.
- The companion calibration study's "6 to 10 percent" overstatement
  (line 931) understates the OL+BDC overstatement observed here.

The architecture difference in SE calibration may itself be a real
and reportable finding. A homoscedastic residual model fitted to
covariance-moderation data, whose conditional variance differs
between on-drug and off-drug occasions, would be expected to
misstate its standard errors. That is an architecture-specific
effect, but it belongs to the analysis model, and it should be
reported as such rather than as irreducible information loss.
(Inferred.)

### 2.5 The design-based explanation contradicts the design specification

**Evidence: inspected (driver lines 65-103); the alternative
explanation is inferred.**

The designs used for Section 3 are as follows.

| Design | Occasions (weeks) | Spacing | Off-drug occasions per path |
|---|---|---|---|
| CO | 2.5, 5, ..., 20 (8 occasions) | 2.5 weeks | Path A: last 4 (10 weeks, after drug). Path B: first 4 (placebo-first, no carryover) |
| Hybrid | 4, 8, 9, 10, 11, 12, 16, 20 | 1 to 4 weeks | Paths A/B: weeks 11-12 and one of 16/20. Paths C/D: weeks 10-12 and one of 16/20 |
| OL+BDC | 4, 8, 12, 16, 17, 18, 19, 20 | 1 week in BDC | Path A: weeks 19-20. Path B: weeks 18-20 |

Table: Occasions, spacing and off-drug occasions per path for each design

The manuscript's account of design dependence (Abstract; Sections
3.1, 3.2 and 4.1; Conclusion 2) holds that OL+BDC "sustains its
off-drug window the longest of the three designs", that CO
"alternates quickly enough" for little variance inflation to act,
and that the OL+BDC blinded phase is "extended". None of these
holds.

- CO path A has the longest uninterrupted post-drug off period of
  any path in any design: four occasions over ten weeks.
- CO does not alternate quickly. Each participant switches once.
- The OL+BDC blinded phase is four weeks long, with two or three
  off-drug occasions.

A more parsimonious explanation is measurement spacing relative to
the carryover half-life. In CO the first off-drug measurement falls
2.5 weeks after discontinuation, when the residual fraction is
$2^{-2.5} \approx 0.18$ at $t_{1/2} = 1$ and $2^{-5} \approx 0.03$ at
$t_{1/2} = 0.5$. In Hybrid and OL+BDC the first off-drug measurement
falls one week after discontinuation, when the residual fraction is
0.5 at $t_{1/2} = 1$. CO's off-drug occasions are therefore nearly
uncontaminated at every half-life in the grid, which by itself
accounts for its insensitivity. (Inferred; the residual fractions
are verified by arithmetic.)

The competing explanation offered in the manuscript, that
"within-subject AR(1) correlation already carries most of the
signal" in CO, is not tested anywhere in the paper, and Section 3.2
concedes that it is "offered ... in place of that missing
derivation" (lines 901-902). The Abstract nonetheless presents it as
established.

**Recommendation.** Replace the uninterrupted-off-drug-time account
with an account based on the density of off-drug measurements
within one or two half-lives of discontinuation, and test it
directly, for example by varying CO spacing at fixed total duration.

### 2.6 The working-independence check is misreported

**Evidence: verified (comparison of the manuscript's own figures).**

Section 3.2 reports that summing $(D_{bc,it} - \bar D_{bc})^2$
predicts information losses of 36.8% (Hybrid), 4.6% (CO), and 51.2%
(OL+BDC), and states that this "correctly orders the three ... under
either architecture" (lines 877-882).

The observed mean-moderation power losses are 6.9% (Hybrid), 1.9%
(CO), and 2.8% (OL+BDC). The approximation places OL+BDC above
Hybrid, while the data place Hybrid well above OL+BDC. The ordering
therefore matches the covariance-moderation losses only. A predicted
51% information loss accompanied by a 2.8% power loss also
contradicts the statement that design-contrast erosion "alone
explains mean moderation's modest, design-dependent loss" (lines
836-837). Section 2.3 above supplies a likely reason: the
mean-moderation estimand grows with carryover and offsets the lost
contrast.

## 3. Major findings: internal consistency and estimand

### 3.1 The mean-moderation DGP is described two ways

**Evidence: inspected.**

Section 2.2.3 (lines 443-451) and Section 4.1 (lines 1007-1009)
describe a mean-moderation DGP in which the interaction scales with
exposure-decayed $D$, so that "the ratio ... is preserved at every
level of exposure" and "carryover keeps the same proportion off
drug". The implementation applies the shift at on-drug occasions
only, which Section 2.2.1 (lines 373-384) and Appendix B.6 state
correctly. Section 2.2.3's first paragraph and the Section 4.1
passage should be rewritten to match the implementation, and the
carryover-invariance claim for $\beta_1\beta_{bm}$ (line 485) should
be withdrawn in light of Section 2.3.

### 3.2 True value, units, and sign convention

**Evidence: inspected.**

- Section 2.3 defines the reference true value as
  $\theta_{\text{true}} = -c_{bm}\sigma_{BR}$ (line 691). Section 2.2.3
  gives the population slope as $c_{bm}\sigma_{BR}/\sigma_B$. The two
  differ by a factor of $\sigma_B$ and cannot both be the target of
  the unstandardized `bm:Dbc` coefficient.
- Bias relative to $\theta_{\text{true}}$ is listed as a performance
  measure but never reported.
- The sign convention required by NOTATION rule 5 (components are
  symptom reductions, so interaction coefficients are negative) is
  never stated. All DGP formulas in Section 2.2 are written with
  positive signs, while every reported estimate is negative.

### 3.3 Comparability of absolute power across architectures

**Evidence: inspected.**

Appendix B.8 shows that the two architectures have different
marginal variances at a common $c_{bm}$, and concludes that
"absolute power is not comparable between the two architectures"
and that Section 3 compares only within-architecture profiles
(lines 2280-2284). Section 2.2.1 states the opposite, that the
common parameter value "licenses running both" (line 368), and the
Abstract and Conclusions compare losses in absolute percentage
points computed from different baselines. The manuscript should
adopt one position. If cross-architecture comparison is retained,
calibrating the architectures on realized biomarker-response
correlation, or on power at $t_{1/2} = 0$, would make the
comparison defensible.

### 3.4 Section 2.2.4 and Appendix B disagree

**Evidence: inspected.**

- "Four differences" names different lists in the two places.
  Section 2.2.4 lists correlation structure, gate, carryover
  (described as unchanged), and analysis model. Appendix B.7 lists
  $A_c$, $C_{cc'}$, $b_p$, and carryover (described as changed).
  Appendix B.7 then states that the cross-factor form "is not listed
  among the four differences of Section 2.2.4" while listing it in
  its own table of four.
- Section 2.2.4 (lines 617-622) and Appendix B.2 state that compound
  symmetry restricts $c_{bm}$ more tightly than AR(1). Appendix B.7
  shows that for any $t_{1/2} > 0$ the ordering reverses, with the
  published implementation admitting the higher ceiling.

### 3.5 Unfinished material and reproducibility gaps

**Evidence: inspected.**

- Three placeholders remain in the manuscript: `[TO RE-VERIFY]` at
  line 623, `[TO RE-RUN]` at line 630, and `[TO RE-VERIFY]` at
  line 2178.
- The Reproducibility section points to the Hendrickson comparison
  arm, whose results have been withdrawn, and does not identify the
  primary driver (`analysis/scripts/quick-sim/01-dgp-prototype.R`),
  its seed (`20260509`, with per-cell seeds `100000 * i + idx`), or
  the result files.
- Line 285 states that the three designs "are defined in Section 3".
  They are not defined anywhere in the manuscript beyond their names
  and path counts; the Hybrid schedule appears only in Appendix B.2.
  A design table equivalent to the one in Section 2.5 of this review
  is needed, and the design explanation cannot be assessed by a
  reader without it.

### 3.6 Conflicting statements about recoverability

**Evidence: inspected.**

Section 4.3 describes the fading correlation as "lost information no
model can recover" (lines 1280-1281). Section 4.1 argues that a
class-aware analysis could in principle recover part of the loss
(lines 1155-1159), and Section 2.2.3 argues that the population
slope is unaffected. Section 4.3 also calls the shrinking exposure
variance "recoverable", although no analysis restores contrast that
the design does not provide. The passage should be reconciled with
Sections 2.2.3 and 4.1, and with Section 2.4 of this review, which
suggests that part of the loss is recoverable by better SE
calibration.

## 4. Minor findings

**Evidence: verified unless noted.**

- **Abstract, line 80.** "$z = 2.80$, $p = 0.005$" is attached to the
  within-architecture covariance-moderation loss. It is the
  difference-of-differences statistic. The within-architecture $z$
  for that loss is about 7.0, while the adjacent mean-moderation
  "$z = 3.34$" is within-architecture. The two parallel sentences
  report different statistics.
- **Line 703.** "The highest cell in the grid is 0.842" should read
  0.844 (mean moderation, Hybrid, $t_{1/2} = 0.5$).
- **Lines 780-782.** The step from 0.754 to 0.723 is not "within one
  MCSE". The SE of the difference is about 0.020, so the step is
  about 1.6 SE.
- **Appendix B.8, line 2248.** The conditional SD
  $8\sqrt{1 - 0.45^2} = 7.144$ should read 7.14, not 7.15. The
  statement that a simulation "reproduces each to three decimal
  places" should be checked against the corrected value.
- **Line 410.** The limit $\lambda = 0$ corresponds to an infinite
  half-life, not to a "no-carryover model variant".
- **Line 749 and elsewhere.** OL+BDC is called a "two-period design".
  It is a two-path design.
- **Lines 1324-1327.** The closing sentence of Section 4.4 is
  logically inverted. If the architecture choice matters least in
  parallel-group designs, assuming it harmless there is not wrong.
  The intended target is presumably within-subject designs.
- **Line 982.** "Essentially immune to carryover" conflicts with a
  statistically significant 6.9% loss in the focal Hybrid design
  ($z = 3.34$).
- **Line 189.** "As Section 4 shows" should refer to Section 3.2.
- **Uncited or unsupported claims.** The following need a citation,
  a supporting analysis, or removal.
  - $N = 70$ as "the smallest of the prazosin-PTSD-matched sample
    sizes used in the published methodological literature" (lines
    700-701).
  - "Methodological work on N-of-1 carryover sharing authors with
    Hendrickson et al." using mean moderation (lines 1301-1304).
  - "A companion paper on carryover-mitigation strategies" (lines
    1289-1290).
  - The "8-13%" and "40-60%" ranges from "earlier exploratory runs"
    (lines 751-753, 783-785, 1334-1336). These are unpublished and
    their removal is recommended.
  - Parallel-group enrichment designs giving "similar power under
    either architecture" (lines 1315-1316). This is untested and in
    tension with Appendix B.8.
- **Rendered PDF numbering.** The table of contents reads "1
  Abstract", "2 1. Introduction", "3.1 2.1 Notation and analysis
  model". The cause is `number_sections: true` combined with
  hand-numbered `\section{}` and `\subsection{}` titles. Either use
  `\section*{}` with manual numbers or remove the manual numbers.
- **Line 1272.** `\end{itemize} These requirements ...` begins a new
  paragraph on the closing line of the list; a blank line is needed.
- **Appendix B.5.2 (inspected).** The loop-order dependence of the
  cross-factor decay rate is a code defect at
  `R/generateData.R:334`, where `rho` is taken from the outer-loop
  factor and the later of the two writes survives. It should be
  fixed in the package rather than documented as an artifact.
- **Appendix B.7.** The smallest-eigenvalue expression
  $(1 - \rho_c) - (c_1 - c_{\times})$ assumes a common $\rho_c$
  across factors, which should be stated.
- **Appendix B.4 and B.5.3.** "At any positive half-life $b^{\text{orig}}$
  is constant" holds only when every off-drug occasion follows a
  drug exposure. It fails for placebo-first paths such as CO path B.

## 5. Notation

Assessed against `analysis/report/NOTATION.md`. The linter passes,
so the items below concern rules it does not mechanize.

- **$b$ is overloaded.** It denotes a realized biomarker value
  (Sections 2.2.2 and 2.2.3), the standardized biomarker $b_i$
  (Section 4.1, consistent with NOTATION), and the biomarker coupling
  vector (Appendix B). Sections 2.2.1 and B.6 use $z_i$ for the
  standardized biomarker where NOTATION specifies $b_i$. The coupling
  vector should be renamed (for example $\gamma$ or $r_{BM}$).
- **Biomarker moments.** $\mu_B$ and $\sigma_B$ (main text) and
  $\mu_{bm}$ and $\sigma_{bm}$ (Appendix B.6 notes, B.8) are used for
  the same quantities, as are $B_i$ and $\mathrm{BM}$.
- **Estimand label.** $\beta_{bm:D}$ (Section 2.3) and
  $\hat\beta_{bm:D_{bc}}$ (Section 3.3) denote the same coefficient.
  NOTATION specifies $\beta_{bm:D}$.
- **Shared $c_{bm}$ label.** The shared label is declared (line 360),
  which satisfies rule 2. However, Section 2.2.1's display uses
  $\beta_{bm}$ and Section 2.2.3 mixes $\beta_1\beta_{bm}$ with
  $c_{bm}$, so the declaration does not hold throughout.
- **Architecture labels.** "Architecture B" appears first at line
  547 and throughout Appendix B, but the main text never introduces
  the A/B/C labels given in the NOTATION table.
- **Time index.** $t$ is both the occasion subscript and the time
  covariate in `Sx ~ bm + t + Dbc + bm:Dbc`, and Appendix B switches
  to $p$ for occasions.

## 6. Companion files

### 6.1 `appendices.Rmd`

**Evidence: verified (diff).** The file is a verbatim copy of
`report.Rmd` lines 1449-2292; the diff differs by a single blank
line. It is untracked. Maintaining two copies guarantees drift,
which is already visible elsewhere in this directory. Either include
a single appendix source in both documents as a knitr child
document, or drop the standalone file and render the appendices
from `report.Rmd` alone. Every finding in Sections 3.4, 4 and 5 that
concerns Appendix A or B applies equally to this file.

### 6.2 `bullets.Rmd`

**Evidence: inspected.** The report describes this file as its
supplementary structured summary (lines 1438-1440), but it is out of
date in several respects.

- Its title and its description of the companion manuscript's title
  do not match the report's title, "Data Generation Machinery for
  N-of-1 Clinical Trial Simulations".
- Its section numbering does not follow the report. It places the
  architectures at 2.1 and 2.2 rather than 2.2.1 and 2.2.2, refers
  to `Dbc` as defined in "Section 3.1", and places the dual-channel
  architecture in "Section 2.2".
- Its Section 3.3 (lines 218-240) reports the $N = 140$ OL+BDC
  robustness check with specific figures. The report states that
  this check "is not reported here" (line 1430).
- Its Section 3.2 describes the superseded account in which the
  signal is "attacked from both sides" through correlation decay,
  not the report's current precision account (which is itself
  contested in Section 2.2 of this review).
- Its Section 3.1 heading describes the covariance-moderation loss
  as "absent at CO", against the reported 3.0%.
- It omits the `in_header: ../sim-preamble.tex` include used by the
  other two files.

The file should be regenerated from the report once the revisions
in Section 7 are complete, not edited in parallel.

## 7. Recommended order of work

1. Re-run the Section 3.1 factorial on the current package, recording
   the commit hash, and confirm convergence and null-cell rates.
2. Add to the reported performance measures the ratio of model-based
   SE to empirical SD, size-adjusted power, and bias against a true
   value whose units match the fitted coefficient (Sections 2.4 and
   3.2).
3. Rewrite the mechanism in Sections 2.2.3, 3.2 and 3.3 around what
   the data support: the estimand shift under mean moderation
   (Section 2.3), architecture- and design-dependent SE calibration
   (Section 2.4), and off-drug measurement density relative to
   $t_{1/2}$ (Section 2.5). Test the spacing explanation directly.
4. Add a design-specification table and complete the Reproducibility
   section (Section 3.5).
5. Resolve the internal inconsistencies of Sections 3.1, 3.3, 3.4 and
   3.6, the minor items of Section 4, and the notation items of
   Section 5.
6. Revise the Abstract, Discussion and Conclusions to match.
7. Consolidate `appendices.Rmd` into a single source and regenerate
   `bullets.Rmd`.

## 8. What this review did not do

- No new simulation was run. The effect of the carryover correction
  on Section 3 (Section 2.1) is inferred in direction and unknown in
  magnitude.
- The size-adjusted power figures use null-cell quantiles from 1,000
  replicates and were not bootstrapped; they indicate magnitude only.
- The spacing-based explanation of design dependence (Section 2.5)
  is a hypothesis consistent with the design specification and the
  residual fractions; it has not been tested.
- The Section 3.2 working-independence calculation was not rerun;
  only its reported figures were compared with the power table.
- The positive-definiteness ceiling of Appendix B.7, the
  two-million-draw check of Appendix B.8, and the Hendrickson
  comparison arm were not re-executed.
- Citations were not checked against `references.bib` or against the
  cited sources.
- Prose style was not reviewed beyond the specific wording problems
  listed in Section 4.

## Appendix 1. Reanalysis script

```r
library(data.table)
x <- readRDS(file.path(
  'analysis/data',
  'quick-sim/01-dgp-replicates.rds'
))
x[, z := beta / betaSE]
s <- x[, .(
  mean_beta = mean(beta),
  emp_sd = sd(beta),
  mean_model_se = mean(betaSE),
  se_ratio = mean(betaSE) / sd(beta),
  power = mean(p < 0.05)
), by = .(architecture, design, t1half, c.bm)]
setorder(s, architecture, design, c.bm, t1half)
print(s, digits = 3)

crit <- x[c.bm == 0, .(crit = quantile(abs(z), 0.95)),
  by = .(architecture, design, t1half)]
alt <- merge(x[c.bm == 0.45], crit,
  by = c('architecture', 'design', 't1half'))
adj <- alt[, .(
  nominal = mean(p < 0.05),
  size_adj = mean(abs(z) > crit[1])
), by = .(architecture, design, t1half)]
setorder(adj, architecture, design, t1half)
print(adj, digits = 3)
```

Run from the repository root with
`Rscript -e "source('path/to/script.R')"`.
