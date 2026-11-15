# Referee report, second round: paper 02, 'Testing a predictive biomarker under drug carryover in aggregated N-of-1 trials: a staged simulation study'
*2026-11-15 09:03 PST*

Manuscript: `analysis/report/02-carryover-sensitivity/report.Rmd`
(rendered `report.pdf`, 2026-11-14 17:10 PST, commit `7c27d7056`) and
`supplement.Rmd` (Tables S1 to S17, Figure S1). This is a second
reading after my report of 2026-11-14 (`referee-report-2026-11-14.md`)
and the revision in commits `7c218f474`, `3371920d0` and `7c27d7056`.
Target register: *Statistics in Medicine*.

## 1. Summary of the manuscript

The paper studies the test of a biomarker-by-treatment interaction in
aggregated N-of-1 trials when the drug has carryover. The interaction
is a correlation between the biomarker and the drug response, and
that correlation decays after the drug stops. Five analyses form a
ladder, A0 to A4. A0 is the published starting analysis. A1 adds a
`corCAR1` residual correlation. A2 adds a CR2 standard error. A3 uses
an exposure-weighted predictor at the true half-life. A4 chooses the
half-life by AIC on the model without the interaction. The revision
adds four things: a size check at the sample sizes of the sample-size
table (Table S17), coverage for families 1 and S (Table 7),
size-adjusted power for all five analyses (Table S16), and a
justification of the candidate grid. The recommendation is now for
the test only, not for estimation.

## 2. Overall assessment and recommendation

Recommendation: minor revision.

The revision answers every major comment of the first round with
data, not with argument. That is to the authors' credit. The new size
check is clean: A4 rejects at 0.040 to 0.058 in all 28 null cells.
The coverage computation confirms the problem I raised, and the
authors now say plainly that the A4 interval should not be reported
as a 95 percent interval. The size-adjusted comparison is honest. It
reports that A0 beats A4 under the crossover design, which the
authors could have left out.

In my view three problems remain, and two of them were introduced by
the revision itself. First, the Abstract now states a size range for
the Hybrid design that was checked only for OL+BDC. Second, the
Conclusions state the size-adjusted lead of A4 over A0 more firmly
than Section 3.5 does, and the crossover result is not carried into
the recommendation. Third, the calibration and coverage numbers were
added to Table 7 but not to Table 8, Figure 6 or the summary figures
of Section 3.4, so the paper now quotes different numbers for the
same quantity. None of these needs a new simulation. Each is a
matter of a few sentences or one table.

## 3. Significance and novelty

My assessment is unchanged. The paper is a careful simulation
comparison under one data-generating process. A4 combines known
pieces, so it is not a new method in the strict sense. The revision
now engages Gärtner et al. (2023), Liao et al. (2023) and Wang and
Schork (2019), and the novelty statement is accurate: none of them
tests a biomarker-by-treatment interaction. A search for 2024 to 2026
work on interaction tests under carryover in aggregated N-of-1 trials
found nothing that overlaps. The paper suits *Statistics in Medicine*
or *Pharmaceutical Statistics*.

One development deserves a sentence. The paper says a selection-aware
interval 'has not been evaluated'. There is now a literature on
inference after AIC selection, including in linear mixed models
(Charkhi and Claeskens 2018; Claeskens, Reluga and Sperlich 2021).
It would be helpful to cite it as the route to such an interval.

## 4. Status of the first-round major comments

| Comment | Status | Evidence |
|---|---|---|
| M1, sample sizes beyond the checked range | Resolved, but see new M1 | Table S17, 28 null cells, A4 at 0.040 to 0.058; OL+BDC checked at 480 and 1,012 |
| M2, A4 standard error and coverage | Resolved, but see new M3 | Table 7 now has families 1 and S; reference cell 0.891; recommendation limited to the test |
| M3, A0 excluded from power comparisons | Resolved, but see new M2 | Table S16; Section 3.5 |
| M4, candidate grid and estimand | Resolved | Section 2.4.6 explains the grid from the visit spacing and states that A4 is evaluated as a test |
| M5a, Limitations against Section 3.6 | Resolved | |
| M5b, conservatism of the two regressors | Partly resolved | Section 3.2.3 is fixed; the Introduction still says it 'applies equally' (minor 3) |
| M5c, 'at its edge' and 'half a point' | Partly resolved | Section 3.2.2 and Section 4.1 fixed; Conclusion 2 not (minor 5) |
| M5d, legacy names | Largely resolved | One trace remains (minor 6) |
| M5e, 'two incorrect ones' | Resolved | |
| M5f, unsupported paired-difference claim | Resolved | Removed |
| M5g, unlike comparisons | Resolved | Section 3.3.1 now compares A4 with A2 and A3; the Weibull cell is reported as showing no gain |
| M5h, Table 1 family 1 | Resolved | |
| M5i, selection independence | Resolved | Sections 2.4.6 and 4.2 |
| M5j, Reproducibility | Resolved in that section | Section 2.7 still says 'the two simulation hosts' (minor 7) |
| M5k, project jargon | Largely resolved | 'job', 'backfill' and one data path remain in the main text (minor 7) |

## 5. Major comments

**M1. The Abstract claims a size range for the Hybrid design that was
checked only for OL+BDC.** (Abstract, Conclusions; summary,
Recommendation.)

The Abstract says: 'For a Hybrid trial ... A4 is a correctly sized
test from $N = 35$ to about 1,000'. Table S17 does not support this
for Hybrid. Under this paper's parameter set Hybrid was checked at
152 and 300 only. Under the starting values it was checked at 300 and
600. The value 1,012 is a single OL+BDC run under the starting
values. I verified this from Table S17 and from the driver output
directories. The summary's recommendation has the same problem: 'it
holds its size up to $N = 1{,}012$ (Table S17)' sits in a paragraph
about Hybrid trials.

This matters because the Abstract is the sentence a reader will
quote. It is also unnecessary. Every Hybrid entry of Table 10 under
this paper's set is at most 292, inside the checked range.

I would suggest that the Abstract say 'from $N = 35$ to 300 under
this parameter set' or 'at every sample size checked'. The summary
should say 600 for Hybrid under the starting values, or drop the
number. Conclusion 7 is accurate as written ('at every $N$ checked,
up to 480 ... and 1,012'), because it does not name a design.

**M2. The size-adjusted comparison is stated more firmly in the
Conclusions than in the Results, and the crossover result is not
carried into the recommendation.** (Section 3.5, Section 4.2,
Conclusion 5, summary.)

Section 3.5 is careful. Under Hybrid A4 leads the size-adjusted A0 by
3.6 to 6.0 points, but no single difference exceeds two bootstrap
standard errors, and the text says so. Conclusion 5 and the summary
then say, without qualification, 'At equal size it also leads the
starting analysis A0 under Hybrid'. Four differences of the same sign,
each under two standard errors, are suggestive. They are not a
demonstrated lead. I checked the four ratios from Table S16: they are
1.9, 1.5, 0.8 and 1.7.

The crossover result points the other way, and more strongly. Under CO
A0 leads A4 by 1.4 to 11.6 points. A1 and A2 also lead A4 in every CO
cell of Table S16 (A2 by 0.6 to 3.5 points, A1 by up to 9.4). Section
4.2 still says that under CO 'the binary indicator is as good a
choice', and that the analyses 'are within a few points of one
another'. On the size-adjusted evidence the binary indicator with CR2
is the better choice under CO. That sentence should say so.

I would suggest two changes. First, Conclusion 5 and the summary
should carry the qualifier of Section 3.5, for example 'and at equal
size it is ahead of A0 in all four Hybrid cells, though no single
difference exceeds two standard errors'. A pooled difference over the
four Hybrid cells, with a bootstrap standard error, would give one
number that a reader can weigh. Second, Section 4.2 should recommend
A2 under the crossover design, citing Table S16.

**M3. The calibration and coverage results were updated in Table 7
but not elsewhere, so the paper quotes different numbers for the same
quantity.** (Section 3.2.4, Section 3.4, Table 8, Figure 6, Section
4.2, Limitations.)

a. **Table 8.** The caption says coverage is 'the mean of the nominal
   95 percent interval over families with coverage output'. Families 1
   and S now have coverage output (Table 7). The column still reads
   only families 2 to 4. I confirmed this in
   `37-ladder-tables-figures.R`, which reads `m3-coverage.csv` and
   `m5-coverage-g10.csv` only. With families 1 and S included the A4
   entry would be about 0.94, not 0.945. The difference is small, but
   it flatters the analysis the paper recommends.
b. **Two sets of numbers for families 2 to 4.** The first paragraph of
   Section 3.4 gives A4 coverage over families 2 to 4 as 0.940 to
   0.948 (one value per family, pooled over $N$). The next paragraph
   gives 0.946 to 0.955 at $N = 35$ and 0.934 to 0.941 at $N = 70$
   (Table 7, split by $N$). Both are correct for their aggregation.
   A reader will see two ranges for the same thing in consecutive
   paragraphs. The first can go, since Table 7 supersedes it.
c. **'A few percent'.** Section 3.2.4 says the selection leaves the
   standard error of A4 'a few percent too small under the
   alternative'. The same paragraph reports $\kappa = 0.87$ in family
   S, that is 13 percent. Section 3.4 and Section 4.2 give 18 percent
   at the reference cell.
d. **'Up to 18 percent in a single cell'.** This figure (Section 4.2,
   Limitations, summary) is the minimum of $\kappa$ over families 1
   and S. The sample-size runs have lower values. I computed $\kappa$
   for A4 over the 270 alternative cells of those runs. The minimum is
   0.79 (OL+BDC, this paper's set, $t_{1/2} = 1$, $c_{bm} = 0.4$,
   $N = 140$). A Hybrid cell reaches 0.80 (starting values,
   $t_{1/2} = 1$, $c_{bm} = 0.4$, $N = 128$). Two cells are below 0.82.
   The text should either say 'in the ladder families' or quote the
   overall minimum, about 21 percent.
e. **Figure 6.** The coverage panel shows families 2 to 4 only. The
   caption says so, which is honest. Adding families 1 and S would let
   the figure agree with Table 7.

The fix is mechanical. One source for coverage (the union of
`m3`, `m5` and `m6`) and one source for the minimum $\kappa$ would
keep the four places in step.

## 6. Minor comments

1. Limitations (Section 4.3): 'Above $N = 70$ size was checked only
   for A2 and A4, at two values of $N$ per parameter set'. OL+BDC now
   has three per set. This sentence was missed in the last revision.
2. Section 4.4, last paragraph: 'the inflation of A0 comes from its
   missing residual correlation (Section 2.4.2)'. The evidence is in
   Section 3.2.1.
3. Introduction, lines 205 to 208: the conservatism 'applies equally
   to every specification'. Sections 3.2.2 and 3.2.3 say the
   exposure-weighted predictor is slightly more conservative. Also,
   'therefore' appears in two consecutive sentences (my earlier minor
   comment 3).
4. Conclusion 5: 'Among the analyses fitted in every family that do
   not need the true half-life, A4 is the most powerful'. Only A3 and
   A4 were fitted in every family, so the set contains A4 alone. The
   sentence is true but empty. Possibly 'Among the analyses that do
   not need the true half-life' was meant.
5. Conclusion 2: CR2 'brings every analysis to within about half a
   point of nominal size'. That holds for pooled rates. Section 4.1
   now adds the cell range (0.026 to 0.074). The conclusion should too.
6. Section 2.4.5: 'Exposure-weighted is therefore a different model
   for each assumed half-life'. This reads as the old family name.
   'The exposure-weighted analysis' would read better.
7. Jargon and stale details in the main text. Section 2.7 refers to
   'the two simulation hosts'; the runs used at least five. 'The
   scenario job' (Sections 2.6 and 3.3.1), 'the backfill' (Section
   2.7) and the path `g9-gap-fill/job-31/summary.rds` (Section 3.3.2)
   belong in the Reproducibility section.
8. Sections 3.2.1 and 3.2.2 both say A0 is 'far outside the band
   under CO and OL+BDC and at its edge under Hybrid', in nearly the
   same words. Once is enough. Under Hybrid A0 is at 0.060, inside the
   500-replicate band.
9. Section 3.6: 'because the response components vary more' is
   plausible but was not decomposed (my earlier minor comment 9).
   'Serial correlation matters as much as the half-life' understates
   the paper's own numbers: at $\rho = 0.9$ the Hybrid $N$ falls to
   0.35 to 0.40 of its $\rho = 0.7$ value, while a 1-week half-life
   roughly doubles it (my earlier minor comment 10).
10. Section 2.6 and Section 3.6: the claim that the noncentrality is
    'close to proportional to $c_{bm}\sqrt{N}$' still has no argument
    or citation. Six of 36 curves fail the goodness-of-fit test, about
    three times the expected rate (my earlier minor comment 8).
11. Table 8: the A4 row still shows a 'worst cell' label next to a
    loss of 0.0 (my earlier minor comment 11).
12. Non-convergence in the sample-size runs is not reported. In the
    three OL+BDC cells at $N = 740$ under the starting values, 3.0 to
    3.6 percent of A4 fits and about 2 percent of A2 fits gave no
    $p$-value. The null cells at $N = 1{,}012$ had almost none. The
    power estimate drops the failures from the denominator. If they
    were counted as non-rejections, power in those cells would fall by
    at most about 3 points, and the effect on Table 10 would be small.
    The concentration at one $N$ and in both analyses suggests a cause
    other than the selection step. I have not diagnosed it. It would be
    helpful to report the rate in the caption of Table 10 and to look
    at those cells.
13. Table S17 prints 1012 without a thousands separator. Table 10
    prints 1010 the same way; the text uses '1,010' and '1,012'.
14. Bibliography: the journal name still appears as 'PLoS One',
    'PLOS ONE' and 'PLoS ONE' (references 3, 21 and 24). The three
    companion papers (references 2, 14 and 18) still have no stable
    identifier, and two of them carry load-bearing claims.

## 7. Missing and questionable references

| Location | Issue | Suggested source |
|---|---|---|
| Sections 3.4 and 4.2, 'a selection-aware interval has not been evaluated' | Post-AIC-selection inference is an established route and should be named. | Charkhi A, Claeskens G. Asymptotic post-selection inference for the Akaike information criterion. *Biometrika*. 2018;105(3):645-664. Claeskens G, Reluga K, Sperlich S. Post-selection inference for linear mixed model parameters using the conditional Akaike information criterion. arXiv:2109.10975 (2021). |
| Section 3.6 | N-of-1 sample size (first round; not added). | Senn S. Sample size considerations for n-of-1 trials. *Stat Methods Med Res*. 2019;28(2):372-383. |
| Section 2.4.4 | Satterthwaite and the Liang-Zeger estimator (first round; not added). | Satterthwaite FE. *Biometrics Bull*. 1946;2(6):110-114. Liang KY, Zeger SL. *Biometrika*. 1986;73(1):13-22. |
| Section 2.4.6, model averaging as an alternative to selection | Optional. | Hjort NL, Claeskens G. Frequentist model average estimators. *J Am Stat Assoc*. 2003;98(464):879-899. |

The two post-selection references were located by web search; I took
their volume and page details from memory and did not verify them on
Crossref. The other entries were suggested in the first round, and the
authors may reasonably decline the optional ones.

## 8. Suggestions for strengthening the paper

1. A pooled size-adjusted difference of A4 and A0 over the four
   Hybrid cells, with a bootstrap standard error (M2). This turns four
   suggestive results into one estimate.
2. A small evaluation of a selection-aware interval at the reference
   cell. The parametric bootstrap that repeats the selection is
   already named. The coverage gap at that cell (0.891 against 0.950)
   is large enough to show whether the bootstrap closes it.
3. One Hybrid null cell at $N = 600$ under this paper's parameter set
   would let the Abstract state a range on the paper's own set
   comparable with the starting values. This is optional; the current
   range covers every Hybrid entry of Table 10.
4. The suggestions of the first round on model averaging and on a
   departure from the analysis model's mean structure still apply.

## 9. Scope of this review

**Verified** (I ran a command or computed the value):

- The rendered text of `report.pdf` and `supplement.pdf` (2026-11-14
  17:09 and 17:10), read in full for the report and for Tables S12 to
  S17 of the supplement.
- Table S17 against the text of Section 3.6 and against the Abstract,
  including the designs and $N$ actually checked.
- The Section 3.5 claims against Table S16: the three ranges, the
  standard errors, the ratios of difference to standard error, and the
  A0 against A1 differences quoted in Section 3.2.1 (at most 8.0
  points).
- The Section 3.6 ratios against Tables 10 and S14 (A2 over A4: 1.14
  to 1.19 in Hybrid, 1.22 in OL+BDC; the starting values against this
  set: 1.87 to 2.13).
- The Table S15 power range 0.844 to 0.936.
- The minimum $\kappa$ of A4 and the non-convergence rates over the
  alternative cells of the sample-size runs, computed from their
  `summary.rds` files.
- Zero unresolved references in the report and the supplement.

**Inspected** (I read the source and judged it):

- The coverage and calibration code in `report.Rmd` (Section 3.4) and
  the coverage source of Table 8 in `37-ladder-tables-figures.R`.
- The sample-size and size-adjusted blocks of `shared-data.R`.
- `46-run-sample-size.R`, including the handling of failed fits.

**Not checked:**

- I did not rerun any simulation or recompute Table 7, Table S16 or
  the ladder CSV files from the per-replicate data.
- I did not check the figures beyond their captions and axis labels.
- I did not diagnose the non-convergence at $N = 740$ (minor 12).
- I did not reread Tables S1 to S11 in detail; they did not change in
  this revision.

**Literature.** I ran two web searches: (a) aggregated N-of-1 trials,
biomarker-by-treatment interaction and carryover, 2024 to 2025; (b)
post-selection inference after AIC selection in mixed models with
cluster-robust tests. The first found no new overlapping work. The
second found the Charkhi and Claeskens line of work, which is cited
above as a suggestion. I did not repeat the first round's DOI audit.
