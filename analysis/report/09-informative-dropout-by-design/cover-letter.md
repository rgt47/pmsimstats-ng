# Cover letter
*2026-05-09 09:55 PDT*

To the Editor,
Statistics in Medicine

Dear Editor:

We are pleased to submit our manuscript titled 'Informative
Dropout and Trial-Design Choice in Aggregated N-of-1
Biomarker-Validation Trials' for consideration as a research
article in Statistics in Medicine.

Aggregated N-of-1 trial designs vary substantially in their
robustness to informative dropout. The published methodological
literature on N-of-1 trials documents dropout as a generic
problem but has not engaged the coupling between the
trial-design family (open-label, open-label plus blinded
discontinuation, traditional crossover, hybrid) and the
inferential consequences of biased dropout for the
biomarker-treatment interaction estimand. The unaddressed
coupling is methodologically consequential: the bias and power
consequences of dropout differ by an order of magnitude across
design families even at a fixed overall dropout rate, and the
most powerful design under no dropout is not the most powerful
design under informative dropout.

This manuscript reports an ADEMP-compliant simulation study
(Morris, White, and Crowther 2019) on a 16-cell design-by-
dropout-pattern grid (4 designs × 4 dropout patterns), at 500
replicates per cell, calibrated to the prazosin/PTSD reference
application. Performance measures are reported with Monte Carlo
standard errors. The principal findings are that (i) the
hybrid Hendrickson design is the most dropout-robust of the four
candidates, losing less than half a percentage point of power
under 40% biased dropout, (ii) the OL+BDC design loses
approximately 5.6 percentage points and the traditional
crossover approximately 8.6 percentage points under the same
conditions, and (iii) the dominant determinant of the power loss
is the dropout *rate* rather than the dropout *mechanism*: the
biased-versus-MCAR contrast at fixed 25% dropout produces
within-design differences below the resolution of the present
simulation, while moving to 40% biased dropout produces clear
power losses in every design. The mean estimated coefficient
$\hat\beta_{bm:D_{bc}}$ is stable across dropout cells, so the
power degradation traces to standard-error inflation rather
than systematic estimand drift.

The paper aligns with Statistics in Medicine's editorial focus
on methodological developments for clinical trials. The
recommendation table in §6 maps three classes of dropout
mechanism to preferred N-of-1 design choices and is intended for
direct use by trial designers. The paper complements companion
submissions on the data-generating-process architecture, the
carryover-mitigation analysis strategies, and the hybrid
Hendrickson design's parameter sensitivity; we anticipate
submitting those papers separately.

This manuscript is original, has not been published elsewhere,
and is not under consideration by another journal. All authors
have approved the submission.

We confirm that the work received no specific funding and that
the pmsimstats team declares no conflicts of interest. The data
and code that support the findings are openly available; the
ADEMP pre-registration is at
`analysis/scripts/informative-dropout-by-design/00-ademp-pre-registration.md`.

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
