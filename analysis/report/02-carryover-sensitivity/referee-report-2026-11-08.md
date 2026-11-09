# Referee report: paper 02, carryover representation and inference for the biomarker-treatment interaction (predictive biomarker, Hybrid design)
*2026-11-08 11:58 PST*

Manuscript: `analysis/report/02-carryover-sensitivity/report.Rmd` as of
commit `d30613e63` (restricted to a predictive biomarker in
`eafe30762`; terminology aligned with paper 01 in `cf54870bc`).
Rendered output checked: `report.tex` and `report.pdf` of 2026-11-07
16:33, which postdate the last edit of the Rmd by one minute. Reviewed
to the standard of Statistics in Medicine; the authors target
*Trials*. Evidence labels: VERIFIED (re-derived or executed on the
tracked data), INSPECTED (read in source or rendered text and judged),
SUSPECTED (plausible, not confirmed). Line numbers refer to
`report.Rmd`.

## 1. Summary

The manuscript reports an exploratory simulation study of how to test
the biomarker-by-treatment interaction in an aggregated N-of-1 trial
when drug carryover is present. The data-generating process (DGP)
encodes the interaction as a correlation between a predictive-only
biomarker and the drug-response component that decays after
discontinuation (the "recommended DGP" of companion paper 01). The
paper builds an "analysis ladder": the starting analysis E0
(Hendrickson et al. 2020, `lmer`, no residual correlation), G1
(binary indicator with `corCAR1`), G6 (G1 with CR2 and Satterthwaite
degrees of freedom), G8 (exposure-weighted predictor at the true
half-life, an oracle) and G10 (half-life chosen by ML-AIC on the
model without the interaction, REML refit, CR2 test). G9, which
selects with the interaction in the model, is a side analysis. The
rungs are compared on paired trials in six "scenario families": a
30-cell base grid over three designs, misspecified half-life and
decay shape, serial correlation, patient-level dropout, a
patient-slope cell, and an independent-seed Hybrid run at up to
5,000 replicates. The authors recommend G10 for Hybrid trials with a
predictive-only biomarker at $N \le 70$.

## 2. Overall assessment and recommendation

Recommendation: MAJOR REVISION.

The revision is a substantial improvement on the 2026-10-27 version.
The selection step is now described correctly and a null-model
variant (G10) was implemented and run; the size cost of selecting on
the interaction model is reported with pooled and paired estimates;
the E0-to-G1 step is decomposed by an isolation run that cleanly
attributes E0's liberal size to its missing residual correlation; an
independent-seed run at 5,000 replicates per null cell now carries the
size claim; CR2 comparators (G6, G8) replace the model-based ones for
power; and the abstract is within the *Trials* limit. Most in-text
numbers are generated inline from tracked data, and the ones I
re-derived match.

The remaining problems are of two kinds. First, several statements
are not supported by the paper's own data: exposure weighting does
not "recover most of the loss" (it recovers 15% to 44% of it under
Hybrid and 19% to 25% under OL+BDC); the Limitations state that fixed
nonzero analysis half-lives were not evaluated, when they were run
and a fixed one-week half-life matches G10's power with no selection
step; and G10's coverage, called "not computed", is computable from
tracked files and shows under-coverage at $N = 70$. Second, the split
of the prognostic material left out-of-scope cells inside two
reported tables (the selection-frequency table is about two thirds
prognostic, repaired-matrix, mean-moderation and heavy-missingness
fits) and left several dangling counts. The DGP parameter values are
also never reported, and they are not those of the companion paper.
None of this is fatal. Most of it can be fixed from stored output
without new simulation, but the recommendation needs to be restated
against the fixed-half-life and G6 comparators before it can stand.

## 3. Significance and novelty

What is new, stated narrowly: under one correlation-encoded DGP, a
paired, stepwise evaluation of (i) residual correlation, (ii) CR2
inference, (iii) an exposure-weighted predictor and (iv) AIC selection
of its half-life with the tested term excluded, for the interaction
(not the main effect) in aggregated N-of-1 designs. The isolation of
E0's inflation to the missing residual correlation, and the
demonstration that excluding the tested term from the selection model
removes most of the post-selection size excess, are the two findings
most likely to transfer beyond this DGP.

Limits on significance. (a) Each component is standard; the
contribution is the assembly and its evaluation, which places the
paper closer to a *Trials* methodology article than to a Statistics
in Medicine methods paper. (b) After the split the recommendation
applies only to a biomarker known to be predictive only, which an
analyst cannot verify from trial data (Major comment M8). (c) The
value of rung 4 over a sensible fixed half-life is not shown and the
stored data suggest it is small (M3). (d) The analysis model is close
to the conditional mean of the DGP, so rung 3 partly shows that a
correctly specified regressor beats a misspecified one; the paper now
says this, which is to its credit.

Contemporariness. The bibliography now includes the 2019 to 2024
simulation literature the last report asked for (Gärtner 2023,
Percha 2019, Blackston 2019, Siepe 2024, Pawel 2024). It still omits
the Bayesian distributed-lag treatment of carryover in N-of-1 trials
(Liao et al., Statistics in Medicine 2023), the closest recent
competitor for modeling carryover as a decaying exposure, and the
small-sample CR2 literature for multilevel models (Section 6). My
searches found no 2023 to 2025 paper on the biomarker-interaction
estimand with carryover in aggregated N-of-1 designs, so the niche
appears open.

## 4. Major comments

Correctness first.

**M1. "Exposure weighting recovers most of the loss" is not what the
data show (VERIFIED, `tier1-g-vs-f-e0/summary.rds`).** Section
3.2.3 (lines 1858 to 1880), the Discussion (line 2753) and Conclusion
3 state that the exposure-weighted predictor recovers most of the
power that carryover removes under Hybrid and OL+BDC. In job 27 at
$c_{bm} = 0.45$, the CR2 binary indicator G6 loses, from
$t_{1/2} = 0$ to 1 week, 0.371 (Hybrid, $N = 35$), 0.165 (Hybrid,
$N = 70$), 0.477 (OL+BDC, $N = 35$) and 0.399 (OL+BDC, $N = 70$) of
power. The paired gain of G8 over G6 at $t_{1/2} = 1$ is 0.056,
0.073, 0.089 and 0.100. The oracle therefore recovers 15%, 44%, 19%
and 25% of the loss. The rendered sentence "The loss that rung 3
addresses is therefore up to half or more of the power" conflates the
size of the loss with the share recovered. Why it matters: the
correctly specified regressor at the true half-life is the best that
any exposure coding can do in this DGP, so most of the carryover loss
is information the data no longer contain, not something the analyst
can buy back. That is a central practical message for trial planners
(it argues for design changes such as longer off-drug spacing or
larger $N$, not for analysis) and the current wording inverts it.
Remedy: report the recovered share explicitly, rewrite the three
passages, and add one sentence on the design implication.

**M2. The DGP parameter values are not reported, and they are not the
companion paper's (INSPECTED: `simulation-core.R`,
`42-check-param-sets.R`, project memory note on the sim-core
parameters).** The Methods describe the structure of the DGP but give
no table of Gompertz parameters (maximum, rate, displacement),
component SDs, baseline mean and SD, biomarker distribution, or the
cross-component correlations other than $c_1$. The tracked results
were produced under the `simcore2026` set (response SDs 6, 2 and 8
for TV, PB and BR; baseline mean 70; biomarker mean 0), whereas the
companion paper 01 and Hendrickson's commit 58b32a9 use SDs 10, 10 and
8 and a biomarker mean of 124 (the `published` set; untracked
reruns under it now sit in `analysis/data/carryover-sensitivity/
published/`). The manuscript nonetheless says its DGP "is the
recommended DGP of the companion study" (line 602) and that its
numbers come from "a single parameter set calibrated to a prazosin
and PTSD trajectory" (line 2846). Why it matters: the D of ADEMP is
not reproducible from the paper; size and power under CR2 depend on
the variance ratios (a PB SD of 2 against 10 makes the placebo
pathway much weaker than in the starting DGP); and a reader
comparing with paper 01 will be
comparing different processes. The calibration claim is unsupported
for a set whose SDs were chosen for the simulation core. Remedy: add
an ADEMP table with every parameter value, state that the structure
is paper 01's recommended DGP but the parameter set differs and why,
and either rerun the headline Hybrid cells under the published set
(the reruns appear to exist) or drop "calibrated".

**M3. The fixed-half-life comparator was run and is not reported, and
the Limitations misstate this (VERIFIED,
`referee-variants/block-{N,A}/summary.rds`).** Line 2859 says
"Analyses at a fixed nonzero analysis half-life were not evaluated."
Blocks N and A store `E3fix_cr2`, the G8 estimator at fixed analysis
half-lives of 0.25 and 1 week, on the same trials as G8, G9 and G10
(the Reproducibility section, line 3199, acknowledges they exist and
are "not reported"). Under the null both are at nominal size (0.043
to 0.052). At the eight alternatives, the fixed one-week half-life
gives 0.273, 0.517, 0.723, 0.959 at true $t_{1/2} = 0.5$ and 0.117,
0.256, 0.370, 0.667 at true $t_{1/2} = 2$, against G10's 0.277,
0.526, 0.737, 0.960 and 0.108, 0.239, 0.362, 0.661. A fixed one-week
value is thus within about one point of G10 when the truth is half a
week and above G10 when it is two weeks, with no selection step and no
post-selection issue. A fixed 0.25-week value loses up to 15 points at
true $t_{1/2} = 2$. Why it matters: rung 4 is presented as the step
that "removes the need to know the half-life"; the stored data show
that a mid-range fixed value removes it almost as well, and the
project's own policy (CLAUDE.md, a fixed $t^a_{1/2} = 0.2$ week as the
validated fallback) is a fixed-value rule. The paper also reports that
AIC does not recover the half-life (Table `stage3-selection`), which
is consistent with the weight being insensitive to it. Remedy: report
the fixed-half-life rows beside G8, G9 and G10 in Tables
`stage3-variants-null` and `stage3-variants-alt`, correct line 2859,
and state what selection buys: protection against a badly chosen
fixed value (0.25 week when the truth is 2), at a cost of about one
point when the fixed value is reasonable. The recommendation may well
survive, but it must be argued against this comparator.

**M4. The prognostic split left out-of-scope cells in reported tables
and left dangling counts (VERIFIED).**

- Table `stage3-selection` (lines 2011 to 2029) reads the pooled rows
  of `referee-reanalysis/m9-selection-frequencies.csv`. At true
  $t_{1/2} = 0.5$ and $N = 35$ the pooled row (34,957 fits) contains
  about 8,000 prognostic-stress fits (16 cells), about 4,000
  random-effects fits from four prognostic cells at strength 0.6
  under repaired matrices, 1,498 mean-moderation fits (`B3-meanmod`)
  and 5,979 heavy-missingness fits (job 31); more than half of the
  row (about 19,500 fits) is out of scope. The in-text figures "56
  to 59 percent" and "29 to 33
  percent" (line 2063) come from these rows. Remedy: recompute from
  the predictive-only cells, or report the independent-seed run's
  selection file (`referee-variants/block-N/selection.rds`, which
  gives 0.59 to 0.64 at 0.25 week and 0.28 to 0.32 at 0.5 week).
- Table `stage3-coverage` (lines 2591 to 2618) filters out `C1-ladder`
  but keeps the 12 per-occasion heavy-missingness cells of job 31,
  which the paper elsewhere excludes; 21 cells per row instead of 15.
  The means barely move when they are dropped (G9: 0.938/0.929
  becomes 0.939/0.930), but the table should match the stated scope.
- Line 1332: "Across the eight Hybrid and OL+BDC null scenarios of
  Table `stage4-null`". The table now has four null rows
  (`hybrid_null`, `olbdc_null`, `t1hifhi_null`, `shapemis_null`); the
  eight included the removed prognostic cells.
- Lines 441 to 446 (Introduction) use the job-21 cell
  `progtvgrow_null` to show that a constant prognostic shift leaves
  G9's size unchanged (0.048 against 0.048). The cell's rows are
  bit-identical to `hybrid_null` across every method because the
  shift is absorbed exactly by the `bm` main effect and the seeds are
  shared; it demonstrates an algebraic identity, not a simulation
  result, and its stored name ("grow") contradicts the text. State
  the identity and drop the numbers, or label the cell honestly.
- Leftovers that are harmless but should go: the `cpb=` substitution
  in `tex_cell`, the redundant `c_pb == 0` filters on `rvN`, and the
  Reproducibility sentence on job 36 block T.

**M5. The recommended analysis's interval is not evaluated although
it can be, and its standard error is more optimistic than stated
(VERIFIED, `g10-backfill/all.rds`, `referee-variants/block-A`).**
Lines 2624, 2862 and 2826 say G10's coverage "was not computed" and
that its standard error is "a few percent too small" under the
alternative. `g10-backfill/all.rds` stores the G10 estimate, SE, $t$
and $p$ per replicate for every family-1 to family-4 trial, which is
all that the coverage reanalysis of `34-referee-reanalysis.R` needs.
Applying the same method (Satterthwaite df recovered from $t$ and
$p$ where $|t| \ge 0.5$, bounds otherwise; my reconstruction
reproduces the stored G9 bounds exactly) to the 30 predictive-only
alternative cells of the supplementary study gives mean G10 coverage
bounds of 0.946 to 0.960 at $N = 35$ and 0.934 to 0.945 at $N = 70$,
against G9's 0.939 to 0.953 and 0.930 to 0.940. So G10 under-covers
by about one to one and a half points at $N = 70$, roughly half of
G9's shortfall. In block A the G10 calibration ratio at the
alternative falls to 0.874 (true $t_{1/2} = 2$, $N = 35$,
$c_{bm} = 0.45$) and to 0.908 at $N = 70$; a 13% understatement is
not "a few percent". Remedy: compute and tabulate G10 coverage with
the existing script, report the range of $\kappa$ honestly, and add
one sentence on whether a selection-aware interval (for example the
parametric bootstrap of the project's F1A0B, or a fixed half-life per
M3) is needed for the interval as opposed to the test.

**M6. The size evidence is summarized with criteria that cannot
discriminate (VERIFIED numbers; INSPECTED method).**

(a) The scoreboard (Table `tab-ladder-score`) and the Results count
"families in the 0.04 to 0.06 band". G10 and G9 are both in band in
6 of 6 families; the oracle G8 is in 5 of 6 only because its family-4
rate is 0.0396 (MCSE about 0.0025), against G10's 0.0408. A count that
ranks the oracle below the recommended analysis on a 0.001 difference
is measuring noise. Family 5 is a single null cell of 1,000
replicates, yet it counts as one of the "six families".

(b) The family-level MCSEs treat pooled cells as independent (README
of `37-ladder-tables-figures.R`), while the cells of families 2 to 4
share random-number streams. Those MCSEs are therefore too small.
The paired SEs, by contrast, are the replicate-weighted mean of cell
SEs, an upper bound that is conservative (the caption is right on
this). The recommendation's key claim, that G10's size is "within one
paired standard error of the oracle in every family" (lines 2818 and
2960), uses that inflated SE as an equivalence criterion, which makes
equivalence easier to claim the noisier the comparison.

Remedy: replace band counts by the paired size difference to G8 with
a confidence interval, and state an equivalence margin in advance of
the reanalysis (for example 0.01) with a two one-sided test; use the
family S cells (independent seeds, 5,000 replicates) as the primary
size evidence, where the paired G10 minus G8 differences are -0.0012
to -0.0002 with paired SEs of 0.0014 to 0.0017, which supports the
claim well on its own; drop family 5 from size counts or call it a
single cell.

**M7. The headline numbers pool designs the recommendation excludes,
and G10's advantage over the feasible comparator is not uniform
(VERIFIED).** The abstract and Section 3.1 report base-grid figures
(E0 0.110, G1 0.038, exposure weighting +0.053) pooled over CO,
OL+BDC and Hybrid, while the recommendation is for Hybrid only and
G10 is the weakest G-family analysis under CO. In the scenario grids
G10 is below G6, the CR2 binary indicator, which also needs no
half-life, at 5 of 18 alternative cells (smallest margin $-0.010$);
the paper reports this in a subordinate clause. Against G6 the
practical question for a Hybrid analyst is "does exposure weighting
plus selection beat the simplest correctly sized test?", and the
answer is "on average by about three points, more when the half-life
is long or serial correlation strong, and not in about a quarter of
the scenario cells". Remedy: give Hybrid-specific numbers in the
abstract, report G10 minus G6 (paired) by design as the headline
power contrast, and keep G8 as the oracle reference.

**M8. Scope and practical value after the split (INSPECTED).** The
paper now recommends G10, which omits the design-expectancy terms,
"for a Hybrid trial with a predictive-only biomarker". Whether a
candidate biomarker is also prognostic is not knowable before the
trial and is hard to test within it. A reader who follows the paper
as written will run an analysis that the companion paper presumably
shows is not robust to a time-varying prognostic association, and
the paper gives no diagnostic or default for the common case of an
unknown role. Combined with the circularity of rung 3 (the analysis
model is close to the DGP's conditional mean), the single parameter
set (M2), one design and $N \le 70$, the recommendation is narrower
than the abstract's "we recommend G10" suggests. Remedy: state in the
abstract and the Recommendation section that the default for a
biomarker of unknown role is the companion paper's analysis, and
present G10 as the component of that analysis that this paper
justifies (the selection rule and the CR2 test); add one paragraph
in the Discussion on what would change under a mean-shift encoding
(already run in `B3-meanmod` and in the factorial's
`mean_moderation` slice, so a sentence on direction is cheap). The
study's exploratory status is stated correctly and should stay.

## 5. Minor comments

**m1.** Caption of Table `tab-spec-summary` (line 688): "G6 to G9 are
CR2 cluster-robust variants of G1, G2, G4 and G3" should read "of
G1, G2, G3 and G4" (G8 is the CR2 variant of G3, G9 of G4; setup
chunk `m_lab`). VERIFIED.

**m2.** Section 2.6 (lines 902 to 905): "The production grid holds 540
cells per encoding. These are the 162 configurations, a Weibull
$k = 1$ condition (54 cells), a linear decay form (216 cells), and a
fourth specification." The stored grid has 540 cells in total, 270
per encoding: 162 plus 54 Weibull $k = 1$ plus 54 linear; `E4` is a
specification, not a set of cells. VERIFIED (`02-grid-summary.rds`).

**m3.** The Introduction (line 456) says every result was produced
under the corrected code "with one labeled exception"; Section 2.8
item 2 says "No result from before this correction is reported." I
could not find a labeled exception in the current text. Reconcile.
INSPECTED.

**m4.** The OL+BDC check with 5,000 null datasets (lines 1627 to 1631)
reproduces: G6 rejects at 0.0486 ($N = 35$) and 0.0562 ($N = 70$) in
`analysis/data/cr2-olbdc-check/` (VERIFIED). Its data directory is not
listed in the Reproducibility section, and a parallel
`cr2-olbdc-check-published/` directory under the other parameter set
(0.050 and 0.054) sits beside it, so the source must be named. The
$N = 70$ rate is two MCSE above 0.05, which deserves a word. The
median CR2 Satterthwaite df is 12 at $N = 35$ and 24 at $N = 70$;
report the df, since they drive the power cost of CR2 at small $N$.

**m5.** Introduction, line 371: "a tenth to a sixth under the
crossover design"; the job-26 range is 7% to 16%. The rendered text
also shows "-0% to 5%" (line 1878) and "0.946 to 0.946" (line 2430).
VERIFIED.

**m6.** "The model-based tests are conservative at every rung"
(line 1909): the job-26 null range reaches 0.056, and in the
2,000-replicate cells G4's $\kappa$ is 0.99 to 1.00 at $N = 35$.
Write "on average". VERIFIED.

**m7.** Section 3.2.2 opens with "The model-based test rejected at
about 0.037 under the null at rung 3", a forward reference; and line
1975 hard-codes "Section 2.8" instead of `\@ref(sec-corrections)`.

**m8.** Discussion (lines 2769 to 2771): "on the base grid the G10
frequencies are within two points of those [of G9]". True when pooled
over the grid (largest difference 1.3 points); at the cell level the
two selections differ by up to 10 points (Hybrid, $N = 70$,
$c_{bm} = 0.45$, $t_{1/2} = 0.5$), as expected under the alternative.
Similarly, "the two selections choose the same half-life with nearly
the same frequencies" (line 2154) refers to marginal frequencies;
report the per-replicate agreement rate. VERIFIED.

**m9.** Blackston et al. (lines 2911 to 2914) are cited as showing
inflated Type I error when carryover is ignored, "consistent with the
E0 and G1 results here". It is not: G1 ignores carryover and is
conservative, and the isolation run attributes E0's inflation to the
missing residual correlation. The Introduction's own argument (line
348) is that carryover misspecification costs power, not validity,
for the interaction under the null. Blackston's estimand is the
treatment main effect, where carryover does bias. Rephrase as a
contrast. INSPECTED.

**m10.** Line 2896: "They [Schmid and Yang] note that carryover is
tractable only in aggregated analyses like ours." My search confirms
the paper discusses carryover extensions; I could not confirm this
specific claim from the abstract-level material. Check the wording
against the full text. SUSPECTED.

**m11.** Terminology. Commit `cf54870bc` aligned the paper with paper
01, but paper 01 calls E0 "the binary-indicator analysis" (19 uses)
and never "the starting analysis"; the glossary papers over this.
Choose one term across both papers. VERIFIED by search.

**m12.** Estimand (lines 651 to 656). The text calls $c_{bm}$ the
estimand, but every test is on the $bm{:}D_{bc}$ coefficient, whose
target is $\theta = -c_{bm}\sigma_{BR}$ with $\sigma_{BR} = 8$; the
value 8 appears only as "$\theta = -c_{bm} \cdot 8$" in a caption.
Under G9 and G10 the coefficient on $D_{bc}(\hat t^a_{1/2})$ is not
$\theta$ when the selected half-life differs from the truth, which
matters for the bias and coverage columns. "Power and Type I error
are on the native test-statistic scale" is unclear. State the
estimand per method.

**m13.** $\kappa$ still uses the mean SE (biased downward a few
percent at $N = 35$, as the text concedes). Use the root mean square
SE; it is a one-line change in the ladder script.

**m14.** Hard-coded numbers remain where inline R is available:
"0.062 to 0.064" (line 1337), "$\kappa$ 1.05 to 1.07, Type I 0.032 to
0.037" (line 1583), "0.91 to 0.93" and "0.98" (lines 1589 to 1590),
"56 to 59" and "29 to 33" (line 2063), "0.051 ... 0.058 to 0.059"
(line 2058), "77 ... 30" (caption, line 2003), "0.016 to 0.018",
"1.13 to 1.18", "1.29" (lines 2321 to 2322), "0.049 ... 0.056" (line
1629), "54,989 of 55,000" (line 3208). Those I checked are correct
today; the M4 example shows how hard-coded counts drift after a
restructuring.

**m15.** The isolation run shows that adding `corCAR1` moves the
binary-indicator estimate toward zero by 12% to 16% under carryover
(lines 1549 to 1559) and the text stops at "we did not examine this
weighting further". Since G1 is the reference for rung 3, one
sentence on why GLS down-weights the off-drug occasions that carry
the residual effect would help.

**m16.** The dropout comparison against "the no-dropout cell of the
independent Hybrid run" is across jobs, seeds and possibly hosts
(lines 2423 to 2430); the Limitations concede there is no rate-zero
cell under the same code. Add a rate-0 cell to the dropout block or
drop the before-and-after numbers.

**m17.** The slope-variance cell was not analyzed with a random slope
(line 2523). Since the cell exists and the variant code exists
(`G9_RANDOM_VARIANTS`), fitting G10 with `~1 + t | ptID` on the
stored trials would close the question cheaply.

**m18.** The REML refit after ML selection can fail separately from
the ML fits; the non-convergence column combines both. Report the
refit failure rate (persisting from the last report, m11).

**m19.** Title (34 words) is long for *Trials*; "under covariance
moderation" in the title uses the term the text replaces by
"correlation encoding". Suggest a title of under 20 words.

**m20.** Format for *Trials*: line numbers, table of contents and
lists of figures and tables should come out of the submission build;
the author line "pmsimstats team" must become named authors; the
Competing interests statement must name any author of Hendrickson et
al. (2020) who is an author here, since the paper evaluates that
work's analysis; the Availability statement points to a GitHub
repository whose public status I could not check and still lacks a
persistent identifier.

**m21.** `references.bib` holds three unused entries
(`jones2014crossover`, `senn2014four`, `vaida2005`); `senn2014four`
still carries a note attributing the paper to "Chen, Chen and Xia".
Remove it. Jones and Kenward (2014) could usefully be cited where
crossover carryover is introduced.

**m22.** The abstract (333 words in the body plus headings and
keywords) is at the limit; Hybrid-specific numbers (M7) can replace
the pooled ones without adding length.

**m23.** Wording that implies a planned or confirmatory design is
absent; the exploratory statement (lines 459 to 462) and "the ladder
was not planned in advance" (line 2719) are appropriate. "Specified
in a working protocol after those jobs" (line 930) is accurate. No
change needed; noted for completeness.

**m24.** Style is clean: no em dashes, US spelling, short declarative
sentences, and "coupling" appears only in the glossary as a code
label. Some sentences remain long (lines 2142 to 2168, 2324 to 2342).

## 6. Missing and questionable references

| Location | Issue | Suggested source |
|---|---|---|
| Introduction and Sec. 4.4 (carryover modeling) | The closest recent model of carryover as a decaying exposure in N-of-1 trials is absent | Liao Z, Qian M, Kronish IM, Cheung YK. Analysis of N-of-1 trials using Bayesian distributed lag model with autocorrelated errors. Stat Med 2023 (doi 10.1002/sim.9676; volume and pages to be checked) |
| Introduction (sample size, interaction power) | Brookes 2004 is for parallel trials; the N-of-1 sample-size literature with patient-by-treatment variance components is absent | Senn S. Sample size considerations for n-of-1 trials. Stat Methods Med Res 2019;28(2):372-383 (pages to be checked); Yang J, Steingrimsson JA, Schmid CH. Sample size calculations for n-of-1 trials (arXiv 2110.08970; check for a journal version) |
| Sec. 2.4 (CR2 with few clusters in multilevel models) | Only Bell and McCaffrey and Pustejovsky and Tipton are cited; recent simulation evidence for CR2 with Satterthwaite df in mixed models is relevant to $N = 35$ | Huang F, Li X. Using cluster-robust standard errors when analyzing group-randomized trials with few clusters (2021 or 2022; check year and journal); Huang F, Zhang B, Accounting for random slopes using cluster-robust standard errors in multilevel models (2025; check) |
| Sec. 2.4 (Satterthwaite df for CR2) | The df choice drives size and power at small $N$ (m4) | Imbens GW, Kolesár M. Robust standard errors in small samples: some practical advice. Rev Econ Stat 2016;98(4):701-712 (from memory; verify) |
| Sec. 4.4, Blackston et al. | Cited as consistent with E0 and G1; it concerns the main effect and contradicts the G1 result (m9) | Keep the citation; rephrase |
| Sec. 4.4, Schmid and Yang | Specific attribution not confirmed (m10) | Check full text: Harvard Data Science Review 2022, Special Issue 3 |
| Sec. 2.4.4 (AIC across LMM fits) | ML-AIC across candidates with differing fixed-effect columns is defensible; the conditional AIC literature is in the bib but uncited | Vaida F, Blanchard S. Conditional Akaike information for mixed-effects models. Biometrika 2005;92:351-370, or remove from bib |
| Introduction (crossover carryover) | Textbook treatment in the bib, uncited | Jones B, Kenward MG. Design and Analysis of Cross-Over Trials. 3rd ed. CRC; 2014 |
| Companion papers | Three unpublished companion manuscripts carry results the reader cannot check (calibration study for the working-covariance mismatch, line 758; paper 01 for the DGP; paper 15) | Provide preprint DOIs or a supplement |

References I checked against search results in this session: Gärtner
et al. 2023, Hendrickson et al. 2020, Liao et al. 2023, Schmid and
Yang 2022, Senn 2019, Percha et al. 2019. The remaining suggested
entries are from memory and need volume and page checks before use.

## 7. Suggestions

- Recompute Tables `stage3-selection` and `stage3-coverage` on the
  predictive-only cells, add G10 coverage, and add the fixed-half-life
  rows from blocks N and A. All of this uses stored files.
- Add an ADEMP table with every DGP parameter value and state which
  parameter set produced each result.
- Replace band counts with paired differences to the oracle and an
  equivalence test; lead the size evidence with family S.
- Report the recovered share of the carryover loss (M1) and a short
  design implication (where the off-drug occasions fall matters more
  than the regressor).
- Give a per-design decision table: analysis, size, power against G6
  and G8, interval coverage, and the evidence tier (independent-seed
  5,000 replicates, shared-stream 500 replicates, or none).
- If a rerun is affordable, put the Hybrid headline cells under the
  published parameter set and, for the default case of a biomarker of
  unknown role, show the cost of using the companion paper's analysis
  when the biomarker is in fact predictive only.

## 8. Scope of this review

Verified (read-only R, `Rscript --vanilla`, scripts in the session
scratchpad): the family-level ladder values and paired increments
(`ladder-predictive/*.csv`) behind the abstract and Section 3.1; the
G6 to G8 recovered share and the null table from
`tier1-g-vs-f-e0/summary.rds`; the four job-21 null rows and the
bit-identity of `progtvgrow_null` with `hybrid_null`; the
supplementary-study size block, scenario-grid null and $\kappa$
values (`g9-gap-fill/summary.rds`); the pooled G9 and G8 size and
discordant counts (`m2-pooled-size.csv`); the composition of the
selection-frequency rows (`m9`) and coverage rows (`m3`); the
independent-seed blocks N and A including the stored fixed-half-life
rows; the G10 backfill guard (54,989 of 55,000) and G9 versus G10
selection frequencies; a reconstruction of G10 coverage from
`g10-backfill/all.rds` (method validated by reproducing the stored G9
bounds); the factorial grid's cell counts and job-26 relative losses
(`02-grid-summary.rds`); the OL+BDC 5,000-replicate CR2 check; the
abstract word count; the citation keys against `references.bib`.

Inspected: the full Rmd including every chunk and inline expression;
the rendered `report.tex`; `simulation-core.R` (`select_half_life_aic`,
`cr2_extract`, the `E7null_cr2` and `E3fix_cr2` branches);
`37-ladder-README.md`; the header of `42-check-param-sets.R`;
`01-run-cr2-null.R`; relevant passages of paper 01's `report.Rmd`;
the 2026-10-27 report.

Not checked: the figures against the data; `38-run-e0-isolation.R`
and `41-run-g10-backfill.R` line by line; the E0 implementation
against Hendrickson's commit 58b32a9; `clubSandwich`'s CR2 theory for
`lme` with `corCAR1`; any rerun of a simulation; the companion papers'
results; the public status of the repository. No project file other
than this report was written; nothing was rendered or committed.

Searches run (web search; no personal identifiers sent): aggregated
N-of-1 carryover with biomarker interaction, 2023 to 2025; CR2 with
few clusters and misspecified covariance in mixed models; N-of-1
carryover and washout models 2021 to 2024; Schmid and Yang, Bayesian
models for N-of-1 trials; Liao, Qian, Kronish and Cheung, Bayesian
distributed lag model; Senn, sample size for N-of-1 trials.
Confidence: high for VERIFIED items; moderate for INSPECTED items;
low for the SUSPECTED attribution in m10.

## 9. Status of the 2026-10-27 report's comments

| Comment | Status | Note |
|---|---|---|
| M1 selection uses interaction model | Resolved | Text corrected; null-model variant G10 implemented and recommended |
| M2 selection costs size | Resolved | Pooled 0.0555 and paired +0.0059 reported |
| M3 G9 SE optimism under alternative | Partly resolved | Attributed to selection; G8 and G9 coverage reported; G10 coverage still missing (M5) |
| M4 shared streams, weak band evidence | Partly resolved | Independent-seed 5,000-replicate run added and shared streams disclosed; band counting persists in the scoreboard (M6) |
| M5 wrong counts | Mostly resolved | Earlier errors fixed; new ones from the split (M4, m2, m5) |
| M6 size and power mixed; E0 confounded | Resolved | G6 and G8 are comparators; isolation run separates E0's three differences; size-adjusted power for E0 still absent but E0 power is now explicitly not compared |
| M7 prognostic pathways | Moved | To paper 15; the vacuous `progtvgrow` cell is still used in the Introduction (M4) |
| M8 generality and circularity | Partly resolved | Circularity and Weibull analysis form now disclosed; single parameter set persists and is undisclosed in values (M2) |
| M9 selection noisy; fixed comparator | Partly resolved | Frequencies reported; fixed-half-life comparator run but suppressed and misdescribed (M3) |
| M10 heavy-missingness $\kappa$ | Resolved | Robust ratio reported |
| M11 dropout narrower than claim | Resolved in text | Conclusion qualified; informative dropout remains a stated limitation |
| M12 random effects | Moved or partly resolved | Slope-variance cell added; random-slope analysis not fitted (m17) |
| m1 notation $t^a_{1/2}$ | Resolved | |
| m2 $\kappa$ uses mean SE | Persists | Disclosed, not changed (m13) |
| m3 "overcorrects" wording | Resolved | |
| m4 paired SE | Resolved | |
| m5, m8, m9 | Moot | Material removed |
| m6 grid composition | Persists, now wrong | m2 |
| m7 data dictionary | Claimed resolved | Not checked |
| m10 CR2 nulls outside Hybrid | Persists | Disclosed |
| m11 REML refit failures | Persists | m18 |
| m12 cross-job comparisons | Resolved in text | Disclosed as unpaired |
| m13 linear form | Resolved | Disclosed |
| m14 *Trials* format | Mostly resolved | Abstract within limit, keywords and declarations added; DOI, authors and build options remain (m20) |
| m15 exploratory in abstract | Resolved | |
| m16 attributions | Resolved | Chen and Chen, Araujo et al. correct; stale bib entry remains (m21) |
