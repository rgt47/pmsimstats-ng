# Cover letter
*2026-05-09 09:55 PDT*

To the Editor,
Statistics in Medicine

Dear Editor:

We are pleased to submit our manuscript titled 'Test-Procedure
Choices for the Biomarker-Treatment Interaction in Aggregated
N-of-1 Trials: Strict RM-ANOVA, Linear-Mixed, and GEE
Comparators' for consideration as a research article in
Statistics in Medicine.

Inference about the biomarker-by-treatment interaction in
aggregated N-of-1 trials depends on the analyst's choice of test
procedure. The classical literature offers strict
repeated-measures ANOVA, linear mixed-effects models, and
generalised estimating equations; each has well-characterised
operating properties in standard regimes but the comparative
behaviour at the trial-relevant sample sizes typical of N-of-1
biomarker-validation studies (N near 25 to 35) has not been
systematically reported. The relevant comparison must address
both type I error (the GEE small-sample inflation problem) and
power (the LME advantage over RM-ANOVA when within-subject
residuals are autocorrelated rather than compound-symmetric).

This manuscript reports two contributions. First, we derive the
closed-form strict-RM-ANOVA $F$-statistic for the
biomarker-stratified treatment-by-time interaction in a balanced
split-plot design with categorical predictors and the
sphericity assumption, and relate it analytically to the Wald
$t$-statistic on the corresponding linear-mixed coefficient.
Under sphericity and balance the strict-RM-ANOVA test is
identical to a two-sample $t$-test on the within-subject
treatment differences. Second, we conduct an ADEMP-compliant
Monte Carlo simulation study (Morris, White, and Crowther 2019)
crossing four test procedures (strict RM-ANOVA, linear-mixed
with corCAR1, GEE with naive sandwich, GEE with Mancl-DeRouen
small-sample correction) against $c_{bm} \in \{0, 0.45\}$ and
$N \in \{25, 35\}$, at 800 replicates per cell. The principal
findings are that (i) the linear-mixed analysis dominates
strict RM-ANOVA on power across both sample sizes by 0.20 to
0.28 absolute, (ii) GEE with the naive sandwich is
anti-conservative by approximately $2.4\times$ nominal at
$N = 25$, and (iii) the Mancl-DeRouen 2001 small-sample
sandwich correction restores nominal type I control at modest
power cost.

The paper aligns with Statistics in Medicine's editorial focus
on methodological developments for clinical trials. The scope
of the present submission is the test-procedure axis at the
prazosin/PTSD reference design; the joint cycle-by-period
design grid that the original abstract anticipated is deferred
to a follow-up. The paper complements companion submissions on
the data-generating-process architecture and on
carryover-mitigation analysis strategies; we anticipate
submitting those papers separately.

This manuscript is original, has not been published elsewhere,
and is not under consideration by another journal. All authors
have approved the submission.

We confirm that the work received no specific funding and that
the pmsimstats team declares no conflicts of interest. The data
and code that support the findings are openly available; the
ADEMP pre-registration is at
`analysis/scripts/test-procedure-design-sensitivity/00-ademp-pre-registration.md`.
The hand-rolled Mancl-DeRouen sandwich correction is documented
inline in the simulation driver.

**Suggested reviewers** (to be confirmed by the corresponding
author):

<!-- Reviewer suggestions to be added at submission. -->

Thank you for considering our submission.

Sincerely,

Temperance Persons
On behalf of the pmsimstats team
tpersons@gmail.com
