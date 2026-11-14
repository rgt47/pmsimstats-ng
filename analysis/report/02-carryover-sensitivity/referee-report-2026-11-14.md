# Referee report: paper 02, 'Testing a predictive biomarker under drug carryover in aggregated N-of-1 trials: a staged simulation study'
*2026-11-14 11:10 PST*

Manuscript: `analysis/report/02-carryover-sensitivity/report.Rmd`
(rendered `report.pdf`, 2026-11-14 10:59 PST, commit `5dfc1eafa`) and
`supplement.Rmd` (Tables S1 to S15, Figure S1). Target register:
*Statistics in Medicine*.

## 1. Summary of the manuscript

The paper studies the test of a biomarker-by-treatment interaction in
aggregated N-of-1 trials when the drug has carryover. The data come
from a multivariate normal process. In it the interaction is a
correlation between the biomarker and the drug response, and that
correlation decays after the drug stops. The authors compare five
analyses arranged as a ladder. A0 is the published starting analysis.
A1 adds a `corCAR1` residual correlation. A2 adds a CR2 cluster-robust
standard error. A3 replaces the binary drug indicator by an
exposure-weighted predictor at the true half-life. A4 chooses the
half-life by AIC on the model without the interaction. The comparison
covers three designs and six scenario families. A final section gives
the sample size A4 needs for 80 percent power. The authors recommend
A4 for a Hybrid trial with a biomarker known to be predictive only.

## 2. Overall assessment and recommendation

Recommendation: major revision.

The study is careful and unusually candid. Every analysis is fitted to
the same trials, so the contrasts are paired. Monte Carlo error is
reported throughout. Two coding defects are disclosed, and their
effect on earlier conclusions is stated. The isolation of the A0 to
A1 step is a good piece of work. The reorganization around five
analyses has made the paper much easier to follow.

In my view three problems stand between the paper and acceptance.
First, the sample-size section reports $N$ well beyond the range in
which the size of A4 was checked. Second, the standard error of A4 is
too small under the alternative, by more than the text admits, and
coverage was not computed where the shortfall is worst. Third, A0 is
removed from every power comparison on the ground that it is liberal.
Yet in the Hybrid design, the design the paper recommends for, A0 is
only marginally liberal. A size-adjusted comparison is needed there.
None of these requires new theory. Two of them can be addressed from
data already stored. The third needs a small additional run.

## 3. Significance and novelty

The literature search (Section 8) found no prior study of a
biomarker-by-treatment interaction test with an exposure-decayed
indicator in aggregated N-of-1 trials. The novelty claim survives.
Three papers come close and should be engaged directly.

- Gärtner et al. (2023) use an exponential wash-in and wash-out term
  in a parametric model for aggregated N-of-1 trials. Their data
  include treatment effects that depend on covariates. The paper says
  only that their estimand is the treatment effect (lines 219 and
  2594). It would be helpful to state exactly how the two estimands
  differ.
- Liao, Qian, Kronish and Cheung (2023, *Stat Med* 42:2044-2060)
  estimate carryover from the data with a Bayesian distributed-lag
  model. That is the main alternative to choosing a parametric
  half-life, and it competes directly with A3 and A4.
- Wang and Schork (2019, *Healthcare* 7:84) give analytic power for
  crossover-based N-of-1 trials with serial correlation and carryover.
  It is the nearest precedent for the sample-size section.

The contribution is a careful simulation comparison under one
data-generating process with one parameter set. It is not a new method
in the strict sense: A4 combines known pieces. In my view it suits
*Statistics in Medicine* or *Pharmaceutical Statistics* once the major
comments are resolved. It does not reach the bar of a theory journal.
The paper is current on reporting practice (ADEMP, Siepe et al. 2024,
Pawel et al. 2024). It is less current on N-of-1 methodology. Senn
(2016, 2019), SPENT 2019 and the Bayesian aggregation literature
deserve at least a sentence each.

## 4. Major comments

**M1. Sample sizes are reported outside the range in which A4 was
validated.** (Section 3.6, Table 10, Abstract, Conclusion 7.)

The sample-size run has no null cells. I confirmed this in
`46-run-sample-size.R`: the grid uses $c_{bm} \in \{0.2, 0.3, 0.4\}$
only. Every size check in the paper is at $N = 35$ or 70. The
Limitations section says so plainly: 'behavior at larger samples was
not examined'. Yet Table 10 reports $N$ from 67 to 302 under the
paper's set and up to about 1,100 under the starting values. The
Abstract quotes 292.

This matters for two reasons. A sample size is only useful if the
test is correctly sized at that $N$. Also, the standard error of A4 is
too small under the alternative (M2). A too-small standard error
raises the rejection rate, so the $N$ for 80 percent power may be
optimistic. The size of the effect is unknown.

I would suggest two remedies. First, add null cells for A4 and A2 at
two or three values of $N$ in the sample-size grid (for example 150
and 300, Hybrid, $t_{1/2} \in \{0, 1\}$). Second, report the size
there beside Table 10. If the size holds, the Limitations sentence can
be revised. If it does not, the table should be restricted or
size-adjusted.

**M2. The standard error of A4 under the alternative is understated,
and coverage was not computed where it is worst.** (Sections 3.2.4,
3.4, 4.2.)

At the reference cell (Table 5, Hybrid, $N = 70$, $t_{1/2} = 1$) the
empirical SD of the A4 estimate is 1.28. Its mean reported standard
error is 1.05. That is $\kappa = 0.82$. The oracle A3 on the same
trials has 1.05 and 1.05. In other words, the selection step inflates
the sampling variance of the estimate by about half, and the CR2
standard error does not see it.

The text understates this in two places. Section 4.2 says the
standard error is too small 'by up to about 13 percent in single
cells'. That figure is computed from family S only. The reference cell
gives 18 percent. The Limitations section says the selection 'costs a
few percent of standard-error calibration'.

Coverage (Table 7) is computed only for families 2 to 4. Those are
the families with the mildest shortfall. Families 1 and S, which
contain the worst cells, have 'no coverage output'. The per-replicate
file of job 27 (`tier1-g-vs-f-e0/all.rds`) stores the estimate, the
standard error and the $p$-value. The Satterthwaite degrees of freedom
can therefore be recovered as in Table 7. The same holds for the
backfill refits.

I would suggest the following.

- Compute coverage for families 1 and S, and report it at the
  reference cell.
- Correct the '13 percent' figure, or compute it over all families.
- State in the recommendation that the A4 interval should not be used
  for estimation without a selection-aware correction. The paper
  already names one, a parametric bootstrap that repeats the
  selection. Possibly a small evaluation of it at the reference cell
  could be added.

**M3. A0 is excluded from the power comparisons in the very design
the paper recommends for.** (Sections 3.1, 3.2.1, Tables 5 and 8,
Figures 2 to 5.)

A0 is liberal under CO and OL+BDC. That is clear. Under Hybrid it is
not clear. In Table S2 (1,000 replicates per cell, pooled over $N$) A0
rejects at 0.061 (MCSE 0.005) under Hybrid. Its power there is 0.702,
against 0.631 for A1. The scoreboard (Table 8) puts A0 inside the 0.04
to 0.06 band in 3 of 4 families. Those families are Hybrid-only. The
paper's own word for A0 under Hybrid is 'marginal'.

A reader who uses the published analysis in a Hybrid trial will ask
whether A4 actually beats it. The paper cannot answer, because A0's
power is dashed everywhere. A size-adjusted comparison would answer
it. The paired null cells give each analysis an empirical critical
value. Power at that critical value can then be compared for all five
analyses on the same trials. I would suggest adding size-adjusted
power for the Hybrid cells of family 1 and the reference cell, at
least in the supplement.

**M4. The candidate grid and the estimand of A4 need justification.**
(Section 2.4.5, Table S4.)

The candidate set is $\{0.25, 0.5, 1, 2\}$ weeks. It excludes zero,
so A4 can never return the binary indicator when there is no
carryover. It also excludes values above 2. Family 2 places the truth
at 0.2 and at 4, outside the grid on both sides. The paper does not
say why the grid was chosen.

The selection also does not recover the half-life (Table S4). At the
reference cell A4 is biased by $-0.16$ at $t_{1/2} = 0$ and $+0.21$
at $t_{1/2} = 0.5$, while A3 is within 0.03. So A4 is a valid test,
but its coefficient does not target $\theta$ in general. Section 2.4
says this in one sentence. In my view the Discussion should say it
too, and the recommendation should be framed as a recommendation for
testing, not for estimating the interaction.

**M5. Internal inconsistencies that remain after the refactor.**

The relabeling to A0 to A4 is clean. I found no old G label in the
rendered text outside the mapping table and the glossary. Every
reference to Tables S1 to S15 points to the right supplement table,
and the report has no unresolved references. The following
inconsistencies remain.

a. **Limitations versus the sample-size section.** Section 4.3 says
   reruns under the starting values 'exist but are not reported
   here'. Section 3.6, Table 10, Table S14 and Figure S1 report them.
b. **Conservatism of the two regressors.** Section 3.2.2 says the
   working-covariance test is 'most conservative with the
   exposure-weighted predictor' (Type I 0.033 to 0.037 against 0.040
   to 0.042 in Table S3). Section 3.2.3 says 'the conservatism is the
   same for both regressors'. The data support the first statement.
c. **'At its edge'.** Section 3.2.2 calls A2 at OL+BDC, $N = 35$
   (0.074) 'at its edge' of the band $[0.031, 0.069]$. It is outside
   the band. The 5,000-replicate check explains it, so the wording can
   simply say so. Relatedly, Section 4.1 says CR2 brings 'every
   analysis to within about half a point of nominal size'. That holds
   for pooled rates only. Cell-level rates run from 0.026 to 0.074
   (Table 8).
d. **Legacy names.** Section 2.4.4 still uses 'Unadjusted' and 'the
   Exposure-weighted family' (rendered line 453). Section 2.4.5 opens
   with 'the Exposure-weighted formula'. These names are not defined
   anywhere in the current paper.
e. **'Two incorrect ones'.** Section 2.4.4 says the factorial grid
   'compares a correct model with two incorrect ones'. Section 2.6
   says the grid fits two analyses, one of them correct. This sentence
   appears to date from the version with a lag-adjusted analysis.
f. **Unsupported claim.** Section 4.4 says the paired-difference
   analysis 'trails the mixed-model analyses under every design' on
   the corrected grid. No table or number in the paper or supplement
   supports this. It should be shown or removed.
g. **Unlike comparisons in Section 3.3.1.** At the $t_{1/2} = 2$ cell
   A4 (CR2) is compared with A1 (working covariance, conservative).
   The like-for-like comparator is A2. The next sentence, on the
   Weibull $k = 2$ cell, gives the A4 power with no comparator, yet
   draws a conclusion about rung 3.
h. **Table 1, family 1.** The row lists $N$, $c_{bm}$ and $t_{1/2}$ as
   if fully crossed (54 cells). The grid has 30 cells: $c_{bm}$ varies
   at $t_{1/2} = 1$ and $t_{1/2}$ varies at $c_{bm} = 0.45$.
i. **Selection independence.** Section 2.4.5 says the selected
   half-life 'does not depend on the coefficient under test'. The
   selection does not use the coefficient. It does use the same data,
   so the two are not independent. Section 4.2 then gives this as the
   reason the size holds ('because the half-life is chosen without
   reference to the term under test'). That is an inference. The
   evidence is the paired size comparison.
j. **Reproducibility.** The section says the manuscript 'was rendered
   inside the project's zzcollab Docker container'. The current PDF
   was rendered with the host R installation
   (`/Library/Frameworks/R.framework/Versions/4.6`; I verified this
   from the render log). It also refers to 'the project's simulation
   host'. The runs used at least five hosts (zeus, mercury, athena,
   joe and julia). Both statements should be corrected.
k. **Project jargon.** 'Paper 01' (Glossary, starting analysis), 'job
   26', 'job 27', 'simrun job 42', 'the referee-variant runs' and
   'block S' appear in the main text. A journal reader cannot use
   them. They belong in the Reproducibility section only.

## 5. Minor comments

1. Abstract, Conclusions: 'A4 is correctly sized ... at $N \le 70$'
   sits next to a sample size of 292. See M1.
2. Line 155: the sentence on the companion paper runs on. It would
   read better as two sentences.
3. Lines 190-197: the argument that size does not depend on the
   exposure coding under the null is sound. 'Therefore' appears twice
   in two sentences.
4. Line 216: 'Earlier work under this encoding fixed the decay shape
   ... and evaluated only one analysis specification.' A citation
   would be helpful (presumably Hendrickson et al. 2020).
5. Lines 221-222: the factor of four from Brookes et al. (2004) is for
   a binary subgroup of equal size with an interaction as large as the
   main effect. The biomarker here is continuous. A qualifier would be
   helpful.
6. Line 250: `[-@bell2002bias]` suppresses the author, so the
   rendered sentence ends with a bare citation. Satterthwaite
   degrees of freedom with CR2 are better attributed to Pustejovsky
   and Tipton (2018) as well.
7. Section 2.4.4: the formula for $D_{bc}$ is the exponential form
   only. The factorial grid builds the predictor with the Weibull form
   in the Weibull cells. One sentence would cover it.
8. Sections 2.6 and 3.6: the noncentrality of the interaction test is
   said to be 'close to proportional to $c_{bm}\sqrt{N}$'. A short
   argument or a citation would be helpful. Six of the 36 curves fail
   the goodness-of-fit test. That is about three times the expected
   rate and deserves a sentence.
9. Section 3.6: the starting values need more $N$ 'because the
   response components vary more'. This is plausible but was not
   decomposed. 'Presumably because' would be more accurate.
10. Section 3.6: 'Serial correlation matters as much as the
    half-life.' At $\rho = 0.9$ the required $N$ falls to 0.35 to 0.40
    of its $\rho = 0.7$ value. That is a larger effect than the
    half-life. The sentence could say so directly.
11. Table 8: the A4 row shows a 'worst cell' label next to a loss of
    0.0. A dash would be clearer.
12. Section 3.5 and Table 8: the text quotes the power of the other
    analyses as 92.7 to 100.9 percent of A4. The range should exclude
    A0, consistently with the decision to dash its power.
13. Reproducibility, sample-size entry: 'the $N$ grid of Section 2.6'
    is not given in Section 2.6. It should be listed there or in the
    supplement.
14. The working-covariance test 'rejects at about 0.037' (Section
    3.2.2, Conclusion 2) is typed in. It matches the computed value
    today. Computing it inline would keep it in step with the data.
    The same applies to 'a fifth to two thirds' and '7 to 16 percent'
    in the Introduction, and '0.01 to 0.03' in the Limitations.

## 6. Missing and questionable references

| Location | Issue | Suggested source |
|---|---|---|
| Line 217 | Percha et al. (2019) is listed among simulations of *aggregated* N-of-1 trials. Its model is a single-patient time series. | Keep it at line 2603 only, or describe it accurately. |
| Lines 219, 2594 | Gärtner et al. (2023) simulate covariate-dependent treatment effects. The novelty statement should address this. | Gärtner T, et al. (2023), as cited. |
| Section 4.4 | The main data-driven alternative for carryover is missing. | Liao Z, Qian M, Kronish IM, Cheung YK. Analysis of N-of-1 trials using Bayesian distributed lag model with autocorrelated errors. *Stat Med*. 2023;42(13):2044-2060. doi:10.1002/sim.9676 |
| Section 3.6 | The analytic power precedent is missing. | Wang Y, Schork NJ. Power and design issues in crossover-based N-of-1 clinical trials with fixed data collection periods. *Healthcare*. 2019;7(3):84. doi:10.3390/healthcare7030084 |
| Section 3.6 | N-of-1 sample size. | Senn S. Sample size considerations for n-of-1 trials. *Stat Methods Med Res*. 2019;28(2):372-383. doi:10.1177/0962280217726801 |
| Line 680 | Liang-Zeger estimator and its small-cluster bias. | Liang KY, Zeger SL. *Biometrika*. 1986;73(1):13-22. Cameron AC, Miller DL. *J Hum Resour*. 2015;50(2):317-372. |
| Lines 249, 682 | Satterthwaite degrees of freedom. | Satterthwaite FE. *Biometrics Bull*. 1946;2(6):110-114. |
| Lines 541, 621, 734 | Software and REML comparability. | Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS*. Springer; 2000. Bates et al. (2015) for `lme4`; Kuznetsova et al. (2017) for `lmerTest`. |
| Line 2608 | Freeman (1989) concerns the two-stage procedure, and Jones and Kenward largely use first-order carryover models. Neither directly makes 'the case against a single lagged term'. | Senn S, Hildebrand H. *Stat Med*. 1991;10(9):1361-1374. Grizzle JE. *Biometrics*. 1965;21(2):467. |
| Line 2549 | MAR definition. | Rubin DB. *Biometrika*. 1976;63:581-592. |
| Companion papers (three keys) | They carry load-bearing claims (calibration mismatch at line 650, non-PD refusal at line 963). A referee cannot check them. | Deposit preprints and cite them with a stable identifier. |

Bibliography hygiene. These items come from the literature audit; I
did not recheck each one myself.

- `tang2020ttests` names the first author 'Jing'. Crossref gives Jillian Tang.
- `schmid2022bayesian` has a spurious `volume = {2022}`.
- `daza2018causal` lacks its DOI (10.3414/ME16-02-0044).
- `kenward1997` should give pages 983-997.
- `senn2014four` duplicates `chen2014` under a misleading key and is unused.
- `vaida2005` is unused.
- The journal name appears as 'PLoS ONE', 'PLoS One' and 'PLOS ONE'.
- The `note` fields recording verification may print under some styles.

## 7. Suggestions for strengthening the paper

1. Size-adjusted power for all five analyses in the Hybrid cells (see
   M3). This would settle the A0 question and also show how much of
   the CR2 gain is the correction of a conservative test.
2. A selection-aware interval for A4 at the reference cell. A
   parametric bootstrap that repeats the selection is already named
   in the text.
3. Model averaging over the candidate half-lives (Hjort and Claeskens
   2003) as an alternative to selecting one. It would address the
   inflated variance directly.
4. One scenario with a different mean structure. The Discussion
   concedes that the analysis model is close to the conditional mean
   of the data-generating process. A single departure would show
   whether the ranking survives.
5. Report the A4 sample sizes with the size check of M1 beside them,
   and flag cells in which the closed-form rule of the companion paper
   falls short (Table S15 already shows three).

## 8. Scope of this review

**Verified** (I ran a command or computed the value):

- The rendered numbers quoted above. I took them from the PDF text of
  `report.pdf` and `supplement.pdf` (rendered 2026-11-14 10:59).
- The supplement table order S1 to S15, and that every 'Table Sn'
  reference in the report names the right table.
- Zero unresolved references.
- The absence of null cells in `46-run-sample-size.R`.
- The cells behind the Abstract's Hybrid gain (all 40 are Hybrid).
- The render host for the Reproducibility claim.

**Inspected** (I read the source and judged it, without rerunning
anything):

- All 3,020 lines of `report.Rmd`, including every chunk.
- The relevant parts of `shared-data.R`.
- The rendered supplement tables S1 to S3.
- The equivalence, scoreboard, reference-cell and sample-size tables.

**Not checked:**

- I did not rerun any simulation or recompute any stored summary.
- I did not verify the ladder CSVs against the per-replicate data.
- I did not check the figures beyond their captions.
- I did not read the companion papers or the full texts of
  Gärtner et al., Schmid and Yang, or Siepe et al.

**Literature.** A delegated audit resolved all 21 DOIs on Crossref and
read the Europe PMC abstracts for Blackston, Percha and Gärtner. It ran
web searches on:

- power for crossover-based N-of-1 trials;
- biomarker-by-treatment moderation in N-of-1 simulations;
- carryover half-life estimation in N-of-1 trials;
- N-of-1 sample size;
- SPENT 2019;
- N-of-1 software;
- recent comparisons of crossover and N-of-1 meta-analysis.

Some lookups were rate-limited. Those references are marked as
recalled in the audit, and I have used them here only as suggestions.
