# Cover letter
*2026-05-09 09:50 PDT*

To the Editor,
Statistics in Medicine

Dear Editor:

We are pleased to submit our manuscript titled 'Robustness of
Carryover-Mitigation Analysis Strategies for Biomarker-Treatment
Interaction: A Factorial Simulation Study' for consideration as a
research article in Statistics in Medicine.

The biomarker-by-treatment interaction in aggregated N-of-1
clinical trials is the inferential target of an emerging
methodological literature on precision-medicine trial design. The
analyst has three principal options for accommodating carryover
effects in the within-subject treatment contrast: a binary
on-drug indicator (A1), an exposure-weighted continuous indicator
(A2), and a lagged-treatment specification (A3). Existing
methodological treatments assess these specifications in
isolation. The joint sensitivity of inference to
mis-specification at the data-generating and analyst layers,
including the functional form of carryover decay and the
analyst's assumed half-life, has not been comprehensively
characterised.

This manuscript reports an ADEMP-compliant factorial simulation
study (Morris, White, and Crowther 2019) crossing three
data-generating decay forms (linear, exponential, Weibull) with
the three analysis specifications across two architectures,
three trial designs, two sample sizes, and three
biomarker-effect levels (540 cells at 500 replicates per cell),
plus five Tier 2 sensitivity blocks. The principal finding is
that the exposure-weighted A2 specification consistently
dominates A1 and A3 by approximately 10 percentage points of
power at the prazosin-calibrated reference cell, with the
ranking robust to decay-form mis-specification and
analyst-versus-truth half-life mismatch. A high-precision rerun
at 600 replicates per cell on a 24-cell expansion confirms the
ranking with paired McNemar $p < 10^{-13}$ at the
highest-leverage cell. The work also addresses the comparison
against prior methodological work, including Hendrickson and
colleagues (2020), Jones and Kenward (2014), and Senn (2016).

The paper aligns with Statistics in Medicine's editorial focus on
methodological developments for clinical trials. It complements
companion submissions on the data-generating-process architecture
for biomarker-by-treatment interactions and on the N-of-1
design's hybrid configuration; we anticipate submitting those
papers separately.

This manuscript is original, has not been published elsewhere,
and is not under consideration by another journal. All authors
have approved the submission.

We confirm that the work received no specific funding and that
the pmsimstats team declares no conflicts of interest. The data
and code that support the findings are openly available; specific
paths are documented in the Reproducibility section of the
manuscript and an ADEMP pre-registration document is committed at
`analysis/scripts/carryover-sensitivity/00-ademp-pre-registration.md`.

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
