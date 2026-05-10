# Cover letter
*2026-05-09 09:50 PDT*

To the Editor,
Statistics in Medicine

Dear Editor:

We are pleased to submit our manuscript titled 'Power Sensitivity
of a Hybrid N-of-1 Design Under Carryover, Serial Correlation,
and Variance-Component Variation' for consideration as a research
article in Statistics in Medicine.

The hybrid N-of-1 design proposed by Hendrickson and colleagues
(2020), combining an open-label titration phase with a blinded
discontinuation phase and a brief crossover, is the most-studied
multi-cycle design for detecting the biomarker-by-treatment
interaction in aggregated N-of-1 trials and is increasingly
favoured in precision-medicine trial design. Its operating
characteristics under realistic carryover, serial correlation,
and between-patient and within-patient variance-component
variation have been examined piecewise across the literature; a
comprehensive sensitivity assessment within a single study has
not been previously reported.

This manuscript reports a 12-sensitivity-sweep simulation study
of the hybrid Hendrickson design's main-effect power, calibrated
to the prazosin-PTSD reference application that motivates the
pmsimstats research programme. The sweeps cross a representative
sample of designs (hybrid N-of-1, alternating A/B/A/B variant,
two-period crossover, parallel RCT) against effect size,
carryover half-life, decay form, between-patient and
within-patient variance components, cycle count, patient count,
AR(1) serial correlation, period length, biomarker-interaction
strength, and a misspecification factorial. Each sweep is
ADEMP-compliant (Morris, White, and Crowther 2019) and reports
performance measures with Monte Carlo standard errors. The
substantive findings reproduce the Hendrickson efficiency
advantage across the explored parameter ranges, quantify the
hybrid's robustness to within-patient and between-patient
variance variation, identify the autocorrelation regime in which
Type I error inflates for the parallel-RCT comparator, and
locate the carryover half-life at which the hybrid begins to
lose its saturation advantage.

The paper aligns with Statistics in Medicine's editorial focus on
methodological developments for clinical trials. It complements
companion submissions on the data-generating-process architecture
and on carryover-mitigation analysis strategies; we anticipate
submitting those papers separately.

This manuscript is original, has not been published elsewhere,
and is not under consideration by another journal. All authors
have approved the submission.

We confirm that the work received no specific funding and that
the pmsimstats team declares no conflicts of interest. The data
and code that support the findings are openly available; specific
paths are documented in the Reproducibility section of the
manuscript and an ADEMP pre-registration document is committed at
`analysis/scripts/nof1-design-sensitivity/00-ademp-pre-registration.md`.

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
