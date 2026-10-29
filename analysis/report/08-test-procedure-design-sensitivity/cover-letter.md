# Cover letter
*2026-10-28 18:30 PDT*

To the Editor,
Statistics in Medicine

Dear Editor:

We submit our manuscript 'Test-procedure choice for the
biomarker-treatment interaction in aggregated N-of-1 trials:
dichotomization, classical RM-ANOVA, mixed-effects, and generalized
estimating equations' for consideration as a research article.

Inference about the biomarker-by-treatment interaction in an
aggregated N-of-1 trial depends on the test procedure. The classical
repeated-measures ANOVA requires a dichotomized biomarker. Mixed
models and generalized estimating equations (GEE) do not. A direct
comparison of the two routes confounds the cost of dichotomization
with the cost of the test procedure. The manuscript separates them.

The manuscript reports three contributions. First, we derive the
closed-form strict RM-ANOVA $F$-statistic for the
biomarker-stratified treatment interaction in a balanced split-plot
design and show that it is a two-sample $t$-test on
within-participant phase differences, exact for any
within-participant covariance once each phase is averaged. Second,
we report a Monte Carlo study under the ADEMP framework (Morris,
White, and Crowther 2019) that compares five procedures: strict
RM-ANOVA on a median-split biomarker, a within-participant contrast
on the continuous biomarker, a linear mixed model with
continuous-time AR(1) residuals, GEE with the Liang-Zeger sandwich,
and GEE with the Mancl-DeRouen bias-corrected sandwich. The contrast
differs from RM-ANOVA only in keeping the biomarker continuous. The
study crosses the procedures with $\beta_{bm} \in \{0, 0.45\}$ and
$N \in \{25, 35, 70, 100\}$ at a two-path open-label plus blinded
discontinuation design, with 1000 replicates per power cell and 2000
per null cell, on the pmsimstats package data-generating process.
Third, we tabulate the small-sample calibration of the
continuous-biomarker procedures, including the calibration ratio
(model standard error over empirical standard error) for every arm
and cell.

The principal findings are as follows. Dichotomization accounts for
0.08 to 0.21 of absolute power and the test procedure for at most
0.04; the observed loss matches Cohen's $2/\pi$ efficiency factor.
GEE with the Liang-Zeger sandwich rejects at 7.4 to 9.7 percent under
the null at every sample size. The Mancl-DeRouen sandwich restores
nominal size only when paired with a $t$ reference on $K - p$
degrees of freedom (4.7 to 5.9 percent), at a power cost of 10 to 15
points relative to the uncorrected GEE at $N \le 35$. The linear
mixed model with model-based inference is conservative (size 2.1 to
3.4 percent) because its working covariance over-estimates the
standard error by 10 to 16 percent; a degrees-of-freedom correction
cannot repair this, and a cluster-robust standard error can. An LME
plus CR2 arm is specified in the driver and will be reported on
rerun.

The paper's scope is the test-procedure axis at a fixed design. A
cycle-by-period design grid is recorded as future work. The paper
complements companion manuscripts on the data-generating-process
architecture, on carryover and the recommended analysis, and on the
Type I calibration of the mixed-model interaction test.

This manuscript is original, has not been published elsewhere, and
is not under consideration by another journal. All authors have
approved the submission. The work received no specific funding and
the pmsimstats team declares no conflicts of interest. The data and
code that support the findings are available in the project
repository; the working protocol and the table of departures from
it are referenced in the manuscript. The study is exploratory, and
its findings should be confirmed in an independent study.

**Suggested reviewers** (to be confirmed by the corresponding
author):

<!-- Reviewer suggestions to be added at submission. -->

Thank you for considering our submission.

Sincerely,

Temperance Persons
On behalf of the pmsimstats team
tpersons@gmail.com
