# Cover letter
*2026-08-20 11:05 PDT*

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
analyst has several principal options for accommodating carryover
effects in the within-subject treatment contrast: a binary
on-drug indicator (Unadjusted, G1), an exposure-weighted continuous
indicator committing to a parametric decay and half-life
(Exposure-weighted, G3), and a lagged just-off-drug term estimating
carryover non-parametrically (Lag-adjusted, G2), the classical
crossover device. Existing methodological treatments assess these
specifications in isolation. The joint sensitivity of inference to
mis-specification at the data-generating and analyst layers,
including the functional form of carryover decay and the
analyst's assumed half-life, has not been comprehensively
characterised.

This manuscript reports an ADEMP-compliant factorial simulation
study (Morris, White, and Crowther 2019) comparing nine
analysis-model specifications for the drug-exposure regressor, in
two tiers. The Tier 1 comparison, our primary evidence base, embeds
the three specifications above in a full factorial grid crossing two
data-generating decay forms (exponential, Weibull) with three trial
designs, two sample sizes, and three interaction strengths (216
cells, 500 replicates per cell). A Tier 2 extension adds six further
specifications (an AIC-selected half-life variant, a
paired-difference regression with no repeated-measures model, and
CR2 cluster-robust variants of four of the above) at targeted
reference cells, together with sensitivity blocks spanning
autocorrelation strength, analyst-versus-truth half-life mismatch,
dropout, biomarker-effect size, DGP decay shape, and cluster-robust
recalibration, all at 500 replicates per cell. The principal finding
is that the Exposure-weighted specification attains the highest
power only in the two designs with a blinded-discontinuation phase
(Hybrid and OL+BDC), leading by approximately 10 percentage points at
the prazosin-calibrated reference cell, but is markedly inferior to
Unadjusted under the classical crossover design (power 0.488 versus
0.830). Where Exposure-weighted leads, that lead is robust to
half-life mis-specification and to autocorrelation strength, and
widens under light-tailed decay-shape mis-specification, but narrows
to statistical insignificance under the two most heavy-tailed Weibull
decay shapes examined. The work also addresses the comparison
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

<!-- Reviewer suggestions to be added at submission. -->

Thank you for considering our submission.

Sincerely,

Temperance Persons
On behalf of the pmsimstats team
tpersons@gmail.com
