# Cover letter
*2026-05-09 09:55 PDT*

To the Editor,
Statistics in Medicine

Dear Editor:

We are pleased to submit our manuscript titled 'Comparative
Statistical Power of Aggregated N-of-1 and Parallel-Group Trial
Designs for Detecting Treatment Main Effects' for consideration
as a research article in Statistics in Medicine.

Aggregated N-of-1 trial designs and parallel-group randomised
controlled trials are the two principal methodological
alternatives for evaluating treatment effects under substantial
individual heterogeneity. The variance-component foundations of
each design have been characterised separately in the literature
(notably by Senn and colleagues), but a comprehensive
sample-size and power comparison across clinically realistic
effect-size and carryover regimes that admits Monte Carlo
standard-error reporting in the @morris2019using ADEMP
framework has not previously been reported. The absence of such
a comparison has limited the field's ability to make defensible
joint recommendations about design selection in
precision-medicine settings.

This manuscript reports an ADEMP-compliant Monte Carlo
simulation study comparing the aggregated N-of-1 hybrid design
to a parallel-group RCT at matched sample size $N = 35$, across
a grid of true-effect sizes ($0$, $-0.5$, $-1.0$, $-2.0$
nightmares per week), calibrated to the prazosin/PTSD reference
application. Each cell is run at 1{,}500 replicates;
performance measures (power, bias, empirical SE, 95% Wald
coverage) are reported with Monte Carlo standard errors. The
principal finding is that the N-of-1 hybrid dominates the
parallel-group RCT across the entire effect-size range
examined. At true effect $-2$, N-of-1 power was
$0.989 \pm 0.003$ versus parallel-group $0.189 \pm 0.010$. The
relative N-of-1-to-RCT power ratio peaked at approximately
$4.5\times$ at the moderate-effect cell, where parallel-group
RCTs are essentially powerless.

The paper aligns with Statistics in Medicine's editorial focus
on methodological developments for clinical trials and
contributes a quantitative basis for design selection in
personalised-medicine settings. It complements companion
submissions on the data-generating-process architecture for
biomarker-by-treatment interactions, on carryover-mitigation
analysis strategies, and on the hybrid N-of-1 design's
parameter sensitivity; we anticipate submitting those papers
separately.

This manuscript is original, has not been published elsewhere,
and is not under consideration by another journal. All authors
have approved the submission.

We confirm that the work received no specific funding and that
the pmsimstats team declares no conflicts of interest. The data
and code that support the findings are openly available; the
simulation workspace
(`analysis/data/derived_data/sim_workspace.RData`) and ADEMP
pre-registration (`analysis/scripts/treatment-main-effect/00-ademp-pre-registration.md`)
are committed alongside the manuscript.

**Suggested reviewers** (to be confirmed by the corresponding
author):

- [TODO: name 1, affiliation, email]
- [TODO: name 2, affiliation, email]
- [TODO: name 3, affiliation, email]

Thank you for considering our submission.

Sincerely,

Temperance Persons
On behalf of the pmsimstats team
tpersons@gmail.com
