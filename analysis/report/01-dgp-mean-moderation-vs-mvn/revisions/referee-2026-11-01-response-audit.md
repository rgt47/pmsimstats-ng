# Audit of the response to the 2026-11-01 referee report on paper 01
*2026-11-05 10:38 PST*

Author: pmsimstats team

Scope: every item of
`analysis/report/01-dgp-mean-moderation-vs-mvn/referee-report-2026-11-01.md`
(Major comments M1 to M8, Minor comments 1 to 10, every row of
Section 6, the Section 7 suggestions), checked against the current
`report.Rmd` (3,071 lines; Appendices A, B and C are inline, lines
2494 to 3071) and the rendered `report.tex` and `report.pdf` of
commit 9c29ab3. The standalone `appendix-b.Rmd` was last changed on
2026-10-27, is not included by `report.Rmd`, and predates the
revision; it was not audited and should be deleted or regenerated so
that it cannot be mistaken for the current Appendix B. Line numbers
below refer to `report.Rmd`. The revision commit 9d7821eac was used
only to locate changes; every verdict rests on the current text.

Evidence labels used below:

- **Verified (data):** the rendered number was recomputed in R from
  the stored output file named.
- **Verified (run):** a script was executed for this audit.
- **Inspected:** the text or code was read; no recomputation.

## Summary

| Item | Verdict | Note |
|---|---|---|
| M1 | ADDRESSED | G10 rerun under published set; all G10 prose data-driven |
| M2 | PARTLY ADDRESSED | Anchored run primary; cause of changes not separated |
| M3 | ADDRESSED | Abstract qualified; pooled case not proved |
| M4 | ADDRESSED | Abstract states the no-carryover condition |
| M5 | ADDRESSED | G10 N80 searched in nine Hybrid cells |
| M6 | ADDRESSED | 4,000-replicate rerun (0.090) resolves the gap |
| M7 | PARTLY ADDRESSED | Strings and Europe PMC added; no WoS/Scopus |
| M8 | ADDRESSED | Non-replication and CSP 563 estimates reported |
| Minor 1 | ADDRESSED | Estimands table row split in three |
| Minor 2 | PARTLY ADDRESSED | Paper 01 consistent; companion's values unstated |
| Minor 3 | ADDRESSED | Underflow attributed to arithmetic |
| Minor 4 | ADDRESSED | rho = 0.7 called a convention |
| Minor 5 | ADDRESSED | A.5 derives q_ts = 0.976 |
| Minor 6 | ADDRESSED | Pooled band formula and interval stated |
| Minor 7 | ADDRESSED | Brookes condition stated in three places |
| Minor 8 | ADDRESSED | Figures relabeled |
| Minor 9 | NOT ADDRESSED | Declarations remain placeholders |
| Minor 10 | ADDRESSED | All listed tables read from data |
| S6: li2023 | ADDRESSED | Reworded as a contrast |
| S6: raskind 2003/13/18 | ADDRESSED | BP cites raskind2016bp only; 2018 null stated |
| S6: higham2002 | ADDRESSED | Now higham1988, description corrected |
| S6: sargent2005 | ADDRESSED | No longer cited for the factor of four |
| S6: brookes2004 extrapolation | ADDRESSED | Extrapolation stated |
| S6: zucker1997, senn2019 | ADDRESSED | Claim reworded to fit the sources |
| S6: liu2009baseline | ADDRESSED | "analogous to" |
| S6: pmsimstats-paper03 | ADDRESSED | Replaced by three mixture references |
| S6 missing: Grenet 2020 | ADDRESSED | Cited; false sentence corrected |
| S6 missing: Jones-Kenward, Senn 1992, Senn-Lambrou | NOT ADDRESSED | None cited |
| S6 missing: Senn 2004, Senn et al. 2011 | NOT ADDRESSED | Neither cited |
| S6 missing: Kraemer 2002 | NOT ADDRESSED | Not cited |
| S6 missing: Ding et al. 2019 | NOT ADDRESSED | Not cited |
| S6 missing: Haller-Ulm, Kent PATH | NOT ADDRESSED | In bib, uncited |
| S6 missing: post-2015 N-of-1 | NOT ADDRESSED | None of the seven cited |
| S6 missing: Van den Noortgate, Moeyaert | ADDRESSED | Cited, plus Moeyaert 2022 |
| S6 missing: Siepe, Pawel | NOT ADDRESSED | Not cited |
| S6 missing: Higham 2016, Ghosh-Henderson, Joe | NOT ADDRESSED | Not cited; block-repair alternative not discussed |
| S6 missing: Naik-Rao 2001 | NOT ADDRESSED | Not cited |
| S6 missing: Holford-Sheiner 1981 | NOT ADDRESSED | Not cited |
| S6 hygiene | ADDRESSED | All fixed except cosmetic tueller2022 key |
| S7: one parameter set, one generator | Partly adopted | Check script, not a unit test; three generators remain |
| S7: G10 over the half-life grid | Partly adopted | G10 at 0, 0.1, 0.2, 1 week only |
| S7: sensitivity model | Not adopted | Still a Discussion proposal |
| S7: calibration appendix | Partly adopted | CSP 563 values reported; nothing estimated |

Two defects not raised by the referee were found (Section "Further
findings"): an unresolved cross-reference that prints as
"reftab:planning" and "reftab:prognosticfirst" in two table captions
of the PDF, and an empty `git_sha` in the metadata of the G10 runs.

## Major comments

### M1. G10 simulated under a different DGP

Request: "rerun driver 19 with paper 01's stated parameters ... Then
rewrite the planning results, the G10 Figure 4B text, the Limitations
bullet and the Conclusions sentence. Add a check that the driver's
generator reproduces the closed-form moments of E9."

Location: Methods, "The original data-generating process" (l.404 to
412); Methods, "Simulation studies" (l.561 to 569); Stage 2, G10
Figure 4B paragraph and chunks (l.953 to 1106); planning section
(l.1791 to 1960); Table `tab:spec` row "Components" (l.1727);
Limitations bullet 1 (l.2219 to 2229); Conclusions (l.2272 to 2277);
Availability (l.2451 to 2460).

Verdict: ADDRESSED.

- The run exists and used the published set (verified, data):
  `published-planning/PARAM_SET` and `block-N/PARAM_SET` read
  `published`; `meta.rds` records `param_set = "published"`, TV and
  PB on the 6.51/0.35 curve, SDs 10, 10, 8, biomarker 124.3/15.4,
  baseline 83.1/18.5. The Rmd reads outputs only through `pub_ok()`,
  which refuses any folder not marked `published`, and the earlier
  `g10-planning/` and `g9-planning/` folders are not read.
- No `[pending]` placeholder remains in `report.tex` (verified, grep).
- Planning numbers (verified, data, from `planning-g10.csv`,
  `planning-g9.csv`, `planning-n80.csv`): null rate 0.046 to 0.053,
  MCSE at most 0.005; G10 above E9 by more than 2 MCSE in 23 of 27
  cells, within in 4, below in 0; Hybrid 0.25/0.2 week 0.357
  (MCSE 0.015) against 0.285; crossover 0.243 (0.014) against 0.202;
  at the E9 N, 0.844 to 0.936 in 20 uncapped cells and 0.601 to 0.915
  in 7 capped cells; no uncapped cell more than 2 MCSE below 0.80;
  nonconvergence at most 0.008; G9 versus G10 at most 0.040. All
  match the rendered text. The G10-over-E9 gap at the referee's cell
  is now 0.072, consistent with the referee's estimate of about 0.10
  on a common DGP.
- Figure 4B G10 numbers (verified, data, `fig4b-g10-cbm0.25.csv`):
  Hybrid 0.34, 0.38, 0.34; OL+BDC 0.34, 0.28, 0.24; crossover 0.20,
  0.19, 0.27; MCSE at most 0.031. The OL+BDC drop at 0.2 week is
  correctly classed as beyond 2 combined MCSE (difference 0.096
  against 0.081). The paragraph now says both figures use the
  original parameters and lists the remaining differences (structure,
  carryover rule, analysis, censoring rule).
- Generator check (verified, run): `42-check-param-sets.R --reps 1000`
  passed every check on 2026-11-05: E9 power 0.280 against 0.285
  (Hybrid), 0.199 against 0.210 (OL+BDC), 0.197 against 0.202
  (crossover), all within 2 MCSE; SD of change from baseline 19.02
  against 19.18 implied, 16.01 against 16.03. The run was in the host
  R session, not the container.
- Limitations and Conclusions now quote the computed ranges and no
  longer claim the confounded contrasts (inspected).

Remaining, outside the scope of paper 01's text:

- The referee's fourth bullet (paper 02 cites paper 01's "extended
  DGP" but ran under different response parameters) is not resolved
  from paper 01's side. Paper 02 (`02-carryover-sensitivity/report.Rmd`
  l.603) still says its DGP "is the extended DGP of the companion
  study" without stating its parameters (inspected). See Minor 2 for
  the proposed sentence in paper 01. The paper 02 rerun is planned in
  `docs/71`.
- `meta.rds` of both runs has an empty `git_sha` (git exit status
  128 inside the container), so the code version of the reported
  runs is not recorded. Proposed fix: record the commit hash in the
  job submission and write it into `meta.rds` from the host, or pass
  it to the container as an environment variable.

### M2. Stress-test table under the original carryover rule

Request: "report the anchored rerun as the primary table ... replace
the caption claim with the per-cell changes ... say which change is
responsible, or rerun with one change at a time."

Location: Methods, "Which form of the DGP each study used" (l.626 to
641); Stage 6 chunks `stressdata`, `prognostic` (l.1403 to 1517);
"Changes from the first run" (l.1540 to 1553); Table `tab:spec`
caption (l.1723); Appendix B.5 (l.3011 to 3022); Appendix C
(l.3054 to 3071).

Verdict: PARTLY ADDRESSED.

- Primary table is now the anchored run, read from
  `summary-full*-anchored.csv`; the first run moved to Appendix C
  (inspected; values verified, data). Original analysis, anchored:
  0.119, 0.077, 0.167, 0.069 (first run 0.107, 0.053, 0.159, 0.075);
  CAR(1) CR2 under the placebo pathway 0.114 (first run 0.100). These
  match the referee's figures.
- The caption claim "changed no analysis's verdict" is gone. The text
  now lists, computed from the data, every cell that crossed the bold
  threshold (the original analysis under the growing pathway, 0.053
  to 0.077) and every cell that crossed the single-cell band (three
  cells), and the largest change (0.031). Verified against the
  rendered text.
- Attribution: the paper states that both changes were made at once
  and that the stored outputs cannot separate them (l.636 to 639,
  l.1551 to 1553, l.3018 to 3019). This declines the referee's second
  option with a reason, but the reason is only that the run was not
  made.
- A residual inconsistency follows from it. Table `tab:spec` (l.1734)
  and B.5 define the extended DGP with fixed allocation, so the
  primary validity table is still not computed under the extended
  DGP as defined; it is computed under the extended DGP plus
  randomized allocation. The paper discloses this.

Proposed fix: run `11-strawman-stress-test.R --carryover anchored
--allocation fixed` for the four prognostic null cells and the
non-prognostic null cells (about 15,000 trials, the size of one
existing run). Report that run as the primary table, since it is the
extended DGP as defined, and attribute each of the two earlier
differences to the carryover rule or to the allocation.

### M3. Monotonicity stated too strongly

Request: "Either prove the pooled case ... or qualify the Abstract."

Location: Abstract, Results (l.69 to 74); Stage 2 (l.865 to 872);
Appendix A.3 (l.2583 to 2586); Conclusions (l.2266 to 2268).

Verdict: ADDRESSED (by qualification, not proof). The Abstract,
Stage 2, A.3 and Conclusions now say the decline is proved for the
path-stratified statistic, checked on a 120-point grid for the
pooled statistic and checked only at the simulated half-lives for
the mixed model (inspected). A.3 states why the pooled case is not
covered. No proof was attempted, which the referee allowed.

### M4. Unconditional existence claim in the Abstract

Request: "State the condition in the Abstract."

Location: Abstract (l.77 to 81); Stage 1 (l.672 to 711).

Verdict: ADDRESSED. The Abstract now reads "Without carryover, the
original compound-symmetry structure admits no joint distribution
for correlations above 0.256. Under the original rule with carryover
the bound is higher ... The original effect size of 0.6 therefore
relied on repaired matrices in part of the original grid" (inspected
in `report.tex` l.284 to 288). Table `tab:nonpd` (l.730 to 746, still
typed) carries the 12 of 40 count the referee cited.

### M5. Planning sample sizes for a statistic nobody analyzes with

Request: "search for the G10 N80 in at least the Hybrid cells ...
and present E9 as the analytical approximation it is, with its
measured bias."

Location: planning section, "The sample size for G10" and chunk
`n80g10` (l.1877 to 1894, 1953 to 1960, 2022 to 2043); Limitations
bullet 1; Conclusions.

Verdict: ADDRESSED.

- Driver 19 block N simulated G10 at six multiples (0.35 to 1.25) of
  the E9 N in all nine Hybrid cells, 500 replicates each, with a
  probit fit of rejection on log N (verified, data,
  `block-N/n80-g10.csv` and `meta.rds`). G10 N80 is 0.71 to 0.83
  times the E9 value and lies inside the simulated range in all nine
  cells (no extrapolation). The rendered table rounds up (193, 205,
  384, 104, 108, 193, ...), consistent with the stored N80_G.
- The text describes E9 as closed form for the paired-difference
  statistic and the ratio as "the measured difference between the
  paired-difference approximation and the analysis a trial would
  use" (l.1957 to 1960).

Residual, disclosed in Limitations: only the Hybrid design was
searched, and the N80 estimates carry Monte Carlo error that the
table does not show. Optional improvement: add a parametric-bootstrap
interval for N80_G from the probit fit, or the delta-method SE, as a
column of Table `tab:n80g10`.

### M6. Unexplained discrepancy between two implementations

Request: "explained or the comparison rerun with more replicates."

Location: setup chunk (l.31 to 40); Abstract (l.74 to 76); Stage 2
(l.751 to 755, 784 to 814); Table `tab:patch` caption (l.919);
Availability (l.2467 to 2468).

Verdict: ADDRESSED.

- Verified, data (`graded-patch/published-th0.1-reps4000.csv`): 16
  chunks, 4,000 replicates, no repairs, power 0.08975 (MCSE 0.0045),
  0.75 combined MCSE from the reconstruction's 0.080, as rendered
  ("within 0.8"). The Abstract now quotes the 4,000-replicate value
  (0.09) instead of 0.126.
- The 500-replicate stored value 0.126 lies 2.34 combined MCSE above
  the 4,000-replicate value of the same code. The paper calls this
  "consistent with a high draw of a 500-replicate run". Since both
  are the same code, chance is the only remaining explanation, but a
  2.3-SE draw (two-sided p of about 0.02) deserves its number in the
  text. Suggested edit (optional): state the 2.3 SE explicitly, and
  lead the first paragraph of Stage 2 (l.753 to 755) with the
  4,000-replicate value rather than 0.126.
- "Rerunning the stored cell with its original seed reproduces 0.126
  exactly" (l.811) was inspected only.

### M7. Novelty claim rests on one query

Request: "report the search strings and databases, extend to Web of
Science or Scopus and to the simulation-methods literature, and
soften the claim to 'in the N-of-1 literature we searched'."

Location: Background (l.107 to 108, 205 to 224, 296 to 298);
Discussion (l.2097 to 2098).

Verdict: PARTLY ADDRESSED.

- Done: full search strings and dates for PubMed (two searches) and
  Europe PMC are reported; the result counts are given; the paper
  notes that Grenet et al. were not retrieved, so the searches are
  not exhaustive; the claim is softened to "the N-of-1 literature we
  searched" in three places; single-case multilevel moderator models
  (Van den Noortgate and Onghena; Moeyaert 2014, 2022) are now
  discussed as mean-model precedents.
- Not done: Web of Science and Scopus were not searched (stated,
  l.219), and the simulation-methods literature (copula and NORTA
  joint generation) is acknowledged as not covered (l.221 to 224)
  but not engaged.

Proposed fix: run the 2026-11-03 string in Scopus (TITLE-ABS-KEY)
and Web of Science (TS=) and report counts; add one sentence placing
the correlation encoding against NORTA and copula generation of
covariates and outcomes from a target correlation matrix, citing
Ghosh and Henderson (2003), which also treats infeasible targets and
so connects to Stage 1. If the Scopus and WoS searches cannot be run,
say so as a limitation rather than leaving it only implied.

### M8. Motivating association did not replicate

Request: "report the failure to replicate in the Background and the
planning section; replace 'no parameter has been estimated' with
what CSP 563 does and does not identify ... and either run the CSP
563-calibrated sensitivity arm ... or state it as future work."

Location: Background (l.133 to 149); criterion 5 (l.328 to 334);
planning section (l.1760 to 1789); "Parameters by status"
(l.2045 to 2062); Discussion, "What remains uncalibrated" (l.2175
to 2189); Limitations bullet 4 (l.2238 to 2244).

Verdict: ADDRESSED.

- Background reports CSP 563's null primary result and the
  non-replication (+0.25 points per 10 mm Hg, 95% CI -2.72 to 3.21,
  279 participants) and recasts prazosin as a worked example.
- The planning section maps the CSP 563 interval to c_bm (-0.57 to
  0.48) and the paper's c_bm range to 1.4 to 2.5 points per 10 mm Hg.
  Both conversions were recomputed by hand for this audit
  (c_bm = -slope x 14.2 / 8) and agree.
- "What CSP 563 identifies" lists the baseline SD 16.0, biomarker SD
  14.2, r = 0.02 with baseline, placebo improvement 1.8 times, SD of
  change 17 to 24, and 0.975 per week with 10 percent noise, and the
  quantities not identified. These match `docs/69` (inspected; the
  reanalysis itself was not rerun).
- The sensitivity arm is stated as future work (l.1786 to 1789).

Residual (not raised by the referee, but a submission issue): the
non-replication is evidenced only by an internal repository document
(`docs/69`) on data that are not public. A journal will need either
a citable source or a short supplementary description of the
reanalysis (model, sample, covariates). Proposed fix: add a
supplementary note with the model and the 279-participant analysis
set, and cite it instead of `docs/69`.

## Minor comments

**Minor 1. Table `tab:estimands`.** Request: "Split the row."
Location: l.472 to 480. Verdict: ADDRESSED. The stress-test row is
now three rows: the original analysis (binary), fixed-half-life
analyses at 0.5 week, and AIC-chosen half-life analyses over 0, 0.25,
0.5, 1 and 2 weeks (inspected).

**Minor 2. TV curve and SD; does the companion paper use them?**
Request: "State whether the companion paper uses these values."
Location: Methods l.401 to 412; Availability l.2458 to 2460.
Verdict: PARTLY ADDRESSED. Paper 01 now states the TV curve and SD,
says every simulation in paper 01 uses the 58b32a9 values (verified
for G10, B.4 and the stress-test driver, which hard-codes SDs 10,
10 eta, 8), and notes that the earlier G10 runs "used other response
SDs and are not reported". It does not say whether paper 02 uses
these values, which was the request. Proposed fix: one sentence in
Methods after l.412, for example: "The companion paper's published
results were generated under a different response parameter set (TV
on the drug-response curve, SDs 6 and 2 eta_t); its rerun under the
values above is in progress." Remove the sentence once paper 02 is
rerun under the published set.

**Minor 3. Step-function threshold.** Location: Stage 2 l.778 to
781; Appendix B.1 l.2859 to 2866. Verdict: ADDRESSED. The underflow
is called "an artifact of double-precision arithmetic in the
original code, not a property of the DGP" and is explained in B.1
(inspected).

**Minor 4. Reason for rho = 0.7.** Location: Stage 4 l.1206 to 1212;
Table `tab:spec` l.1732. Verdict: ADDRESSED. The text calls 0.7 "a
convention only, not an estimate", states that 0.3 is already
feasible at 0.8 under the separable form (0.440), and adds that 0.7
has the largest separable bound of the four values. Verified (data,
`ceilings-by-rho-min.csv`): separable bounds 0.448, 0.462, 0.474,
0.440 at rho 0.5, 0.6, 0.7, 0.8.

**Minor 5. Derivation for Stage 5 property 4.** Location: Stage 5
l.1316 to 1328; Appendix A.5, "The summed response" (l.2767 to
2797). Verdict: ADDRESSED. A.5 derives Corr(S_t, S_s) = q_ts A_ts
with q_ts <= 1 by Cauchy-Schwarz, exact AR(1) within a phase, and
q_ts = 0.976 across the blinding boundary. The values 368, 257, 300
and 0.976 were recomputed by hand for this audit and agree.

**Minor 6. Interval used for the pooled band.** Location: l.1632 to
1635 and l.594 to 598. Verdict: ADDRESSED. The text gives
2 sqrt(0.05 x 0.95 / 11,000) = 0.004 and the band 0.046 to 0.054
(rendered; inspected).

**Minor 7. Brookes condition.** Location: Abstract l.44 to 47;
Background l.100 to 104; planning l.1905 to 1910. Verdict:
ADDRESSED. All three state the condition (interaction as large as
the main effect, two equal subgroups), and the planning section
flags the extension to a continuous slope as an extrapolation.

**Minor 8. Code labels in figures.** Location: Figures
`powerfig`, `gapfig`, `fig4b` (`16-covariance-rule-figures.R`,
l.25 to 26), `fig4bmc` (`16-fig4b-graded-heatmap.R`, l.22 to 24);
Glossary l.2307 to 2314. Verdict: ADDRESSED. The figure scripts map
`step` and `graded` to "Original rule (constant correlation after
first exposure)" and "Decaying correlation"; the regenerated
`37-fig3-power-heatmap-cov.png` was viewed and carries the new
labels (inspected). The Glossary keeps the code labels as
identifiers, which is appropriate.

**Minor 9. Declarations.** Location: l.2475 to 2485. Verdict: NOT
ADDRESSED. Competing interests, funding, contributions and
acknowledgements remain "To be completed by the authors". The
competing-interests text now says what must be disclosed (overlap
with the original authors, whether they were consulted) but does not
disclose it. This is an author action and cannot be completed from
the repository. Proposed fix: before submission, state the author
overlap with Hendrickson et al. (2020) explicitly and whether the
original authors reviewed the two extensions.

**Minor 10. Tables typed into the Rmd.** Location: chunks
`closedform` (l.821), `patch` (l.905), `rhostage` (l.1226),
`ceilings` (l.1353), `prognostic` (l.1507), `options` (l.1674),
`planning` (l.1962). Verdict: ADDRESSED for every table the referee
listed: each now reads its stored output and several assert row
counts with `stopifnot()` (inspected). Residual, not listed by the
referee: Table `tab:nonpd` (l.730 to 746) is still typed, as are
many prose numbers (Stage 2 bullets for Figure 4B, Stage 4 and 5
power values, the joint-ceiling margin of 0.006, the `docs/47`
size-adjusted powers). Spot checks for this audit (verified, data):
Stage 4 power 0.721 and 0.628 and Stage 5 power 0.335 and 0.290
(`power-by-rho.csv`); Figure 4B original-code bullets 0.66 to 0.08,
0.62, 0.68, 0.57, 0.12, 0.13, 0.60, 0.34, 0.33, 0.71, 0.78
(`fig4b-graded-cbm0.25-reps250.csv`); E9 planning N 272/278/480,
84/146 (`planning-n80.csv`); Abstract 0.64 and 0.67
(`graded-patch-reps500.csv`). All agree. Proposed fix: read
`tab:nonpd` from the output of `11-pd-repair-count.R`.

## Section 6. Citations

### Citations that did not support the attached claim

| Row | Verdict | Current text |
|---|---|---|
| li2023 | ADDRESSED | l.292 to 295: ICCs of modifier and outcome; modification in the mean model |
| raskind2003/2013/2018 | ADDRESSED | l.134 to 142: trials cited for "tested"; 2018 null stated; BP finding cites raskind2016bp alone |
| higham2002 | ADDRESSED | l.433 to 437: higham1988, nearest symmetric PSD; eigenvalues below 2e-3 raised |
| sargent2005 | ADDRESSED | l.100 to 101: cited only for the need for an interaction test |
| brookes2004 (planning) | ADDRESSED | l.1905 to 1910: binary-subgroup condition and extrapolation stated |
| zucker1997, senn2019samplesize | ADDRESSED | l.190 to 193: claim now "models and sample-size methods", which the sources support; simulators cite chen2014, wang2019 (l.231 to 232, 285) |
| liu2009baseline | ADDRESSED | l.2644 to 2646: "analogous to" |
| pmsimstats-paper03 | ADDRESSED | l.2115 to 2117: bauercurran2003, lubkemuthen2005, nylund2007 |

The higham1988 description was checked against the
`corpcor::make.positive.definite` algorithm (eigenvalues below
2 x tol raised to 2 x tol, no unit-diagonal constraint) from memory
of the package source, not re-read for this audit (inspected,
partial).

### Missing references

| Claim location | Verdict | Status in current text |
|---|---|---|
| Grenet et al. 2020 (l.131) | ADDRESSED | Cited l.162 to 164; "None of these four" corrects the false sentence |
| Jones and Kenward 2014; Senn 1992; Senn and Lambrou 1998 (carryover) | NOT ADDRESSED | None cited; `jones2014crossover` in bib, uncited |
| Senn 2004 BMJ; Senn, Rolfe, Julious 2011 (individual response) | NOT ADDRESSED | Neither cited; `senn2011individual` in bib, uncited |
| Kraemer et al. 2002 (moderator definition) | NOT ADDRESSED | Not cited, not in bib |
| Ding, Feller, Miratrix 2019 (c_bm^2 as explained variance) | NOT ADDRESSED | Discussion l.2104 to 2106 states the share without citation |
| Haller and Ulm 2018; Kent et al. 2020 (prognostic biomarker) | NOT ADDRESSED | Both in bib, uncited |
| Post-2015 N-of-1 methods (seven sources) | NOT ADDRESSED | None cited; `tang2020`, `tueller2022` in bib, uncited |
| Van den Noortgate and Onghena 2003; Moeyaert 2014 | ADDRESSED | Cited l.165 to 168, with moeyaert2022power |
| Siepe et al. 2024; Pawel et al. 2024 (simulation reporting, post hoc additions) | NOT ADDRESSED | Not cited; the exploratory-status paragraph (l.647 to 660) would carry them |
| Higham, Strabic, Sego 2016; Ghosh and Henderson 2003; Joe 2006 (repair versus stop) | NOT ADDRESSED | Not cited. Stage 1 (l.713 to 721) argues against an opaque repair, but not against a fixed-block repair that preserves c_bm, which was the referee's point |
| Naik and Rao 2001 (separable covariance) | NOT ADDRESSED | Only galecki1994 cited (l.287 to 289) |
| Holford and Sheiner 1981 (symptom versus plasma half-life) | NOT ADDRESSED | l.340 to 346 calls the range an assumption; no effect-compartment citation |

Proposed fix, in one pass: add the following, each to the sentence
the referee named.

- Jones and Kenward (2014) and Senn (1992) at the first mention of
  carryover in the Background and at Stage 3; Senn and Lambrou (1998)
  at Stage 3 on whether carryover is identifiable.
- Senn (2004) and Senn, Rolfe and Julious (2011) at l.190 to 193 and
  in the Discussion sensitivity-model paragraph (l.2158 to 2162).
- Kraemer et al. (2002) at the first definition of a predictive
  biomarker (l.198 to 200).
- Ding, Feller and Miratrix (2019) at l.2105 to 2106.
- Haller and Ulm (2018) and Kent et al. (2020) at the start of the
  Stage 6 prognostic paragraph (l.1391).
- At least Schmid and Yang (2022), Kravitz and Duan (2022) and
  Tueller et al. (2023) in the Background paragraph on N-of-1
  simulators (l.151 to 168).
- Siepe et al. (2024) beside Morris et al. at l.555 to 557, and
  Pawel, Kook and Reeve (2024) in the exploratory-status paragraph.
- Higham, Strabic and Sego (2016) at l.713 to 721, with one sentence
  saying why stopping is preferred to a fixed-block repair (for
  example: a block-preserving repair keeps c_bm but alters the
  response block, so it changes the noise the interaction is tested
  against, which Proposition 1 requires to be fixed). Ghosh and
  Henderson (2003) and Joe (2006) in the same paragraph or in M7's
  new sentence.
- Naik and Rao (2001) beside galecki1994 at l.289 and in A.5.
- Holford and Sheiner (1981) at l.343 to 346.

The DOIs the referee supplied should be re-checked when the entries
are added; this audit did not query Crossref.

### Bibliography hygiene

Verdict: ADDRESSED, one cosmetic residual (inspected,
`references.bib`).

- raskind2003 and raskind2013 now carry the DOIs the referee gave.
- sturdevant2021carryover: author "Sturdevant, S. Gwynn".
- Journal names: "PLoS ONE" throughout; "Healthcare" for both
  wang2019 and blackston2019.
- gartner2023 has `number = {1}`.
- `pmsimstats-paper11` is no longer in the bib.
- Residual: the key `tueller2022` still carries year 2023. The entry
  is uncited, so the PDF is unaffected; rename the key to
  `tueller2023` when it is cited (see the missing-references fix).

## Section 7. Suggestions

- **One parameter set, one generator.** Partly adopted. The
  parameter sets are now named and dispatched in
  `simulation-core.R`, output folders are guarded by `PARAM_SET`, and
  `42-check-param-sets.R` checks E9 power and the SD of change against
  the closed form (verified, run). But three generators remain (the
  closed-form `build()` in `04-power-simulation.R`, the stress-test
  driver's own generator with hard-coded values, and
  `simulation-core.R`), and the check is a standalone script, not a
  `tinytest` unit test run by `make test`.
- **Mixed-model power without the paired-difference proxy.** Partly
  adopted. G10 is simulated at 0, 0.2 and 1 week (planning) and 0,
  0.1 and 0.2 week (Figure 4B), not over the half-life grid, so the
  monotonicity claim for the mixed model remains "checked only at the
  simulated half-lives".
- **The sensitivity model.** Not adopted. It remains a proposal in
  the Discussion (l.2149 to 2173), "not yet implemented".
- **Calibration.** Partly adopted. CSP 563 estimates are reported in
  the planning section, but rho and sigma_BR are not estimated in an
  appendix, and the calibrated sensitivity arm is future work.

## Further findings not raised by the referee

1. **Broken cross-references in two table captions (verified, PDF
   text).** The captions of Table `tab:prognostic` (l.1513) and Table
   `tab:n80g10` (l.2042) render in `report.tex` as `\\ref{...}` and in
   the PDF as "reftab:prognosticfirst" and "reftab:planning", both in
   the table body and in the List of Tables. Proposed fix: change the
   cross-reference in those two R caption strings so that the LaTeX
   output contains a single backslash (for example, write the caption
   with `\\ref{tab:planning}` directly, since the output is LaTeX
   only), re-render through `tools/render.sh` and confirm with
   `pdftotext` that no "reftab" remains.
2. **Run provenance.** `meta.rds` for the published G10 runs has an
   empty `git_sha` (see M1).
3. **Abstract wording.** "With any carryover the interaction is then
   absent and power equals the test size" (l.67 to 68). Under the
   original rule E9b and the original analysis keep power of about
   0.08 to 0.09 (Stage 2; the 4,000-replicate rerun gives 0.090),
   because the baseline row keeps g_t = 0. Proposed edit: "power falls
   to near the test size".
4. **Rendered phrasing in the G10 Figure 4B paragraph.** The
   generated text reads "within Monte Carlo error of it at 0.1 week
   and lower, at 0.24 at 0.2 week in the OL+BDC design"
   (`report.tex` l.1284). Proposed fix in `f4_cmp()` (l.996 to 1000):
   generate "within Monte Carlo error of it at 0.1 week and lower
   (0.24) at 0.2 week".
5. **Stale standalone appendix.** `appendix-b.Rmd` (last changed
   2026-10-27) is not included by `report.Rmd` and predates the
   revision of Appendix B; delete it or mark it superseded.

## What this audit did not do

- The CSP 563 reanalysis, the 0.126 seed reproduction, the repair
  counts of `tab:nonpd`, and Appendices B.1 to B.3 code were inspected,
  not rerun.
- No literature database or Crossref query was run; the verdicts on
  missing references rest on whether the key appears in `report.Rmd`
  and `references.bib`.
- Paper 02 was read only at the sentence that cites paper 01's DGP.
- `report.pdf` was assumed current with `report.tex`; the two
  caption defects were confirmed in the PDF text.
