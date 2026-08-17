# Cover letter
*2026-08-16 16:33 PDT (revised to match the corrected manuscript
finding; see `docs/pub_review_remediation_2026-08-16.md`, series 07)*

To the Editor,
Statistics in Medicine

Dear Editor:

We are pleased to submit our manuscript titled 'Critical
Evaluation of the Modified Gompertz Response Function as a
Parametric Template for Treatment-Response Trajectories in
N-of-1 Simulation' for consideration as a Methods Note (or
Tutorial-track submission) in Statistics in Medicine. The
manuscript body is approximately 2,500 to 2,800 words, in
keeping with the journal's shorter-format venues.

The pmsimstats simulation framework that supports power
analyses for aggregated N-of-1 trials of biomarker-treatment
interactions adopts a modified Gompertz function as the
parametric template for each latent response trajectory. The
choice loads into the simulated data before any
analysis-model question is asked, yet the methodological
N-of-1 corpus is silent on family-misspecification
sensitivity. The present manuscript supplies that comparison.

This is a focused 16-cell sensitivity factorial executed at
1,000 replicates per cell, crossing four trajectory families
(modified Gompertz, symmetric logistic, hyperbolic-tangent,
piecewise-linear breakpoint) with two data-generating
architectures (mean moderation and MVN differential
correlation) and two biomarker effect-size levels at the
prazosin-PTSD reference cell. Performance measures follow the
Morris, White, and Crowther (2019) ADEMP standard with
binomial-proportion Monte Carlo standard errors. The headline
finding is architecture-conditional: under the
covariance-moderation (MVN) architecture, trajectory family is
exactly inert, as expected by construction, with matched
random-number seeds returning identical decisions across all
four families. Under the mean-moderation architecture, the four
families separate by a real, if modest, 0.039 spread in power
(0.743 for Gompertz to 0.782 for the symmetric logistic), with
the Gompertz default the lowest-power, mildly conservative,
family of the four. The Gompertz-logistic gap is approximately
two standard errors of the difference (SE approximately 0.019)
and is not attributable to simulation noise. An
earlier draft of this manuscript, produced before a
subject-invariant additive-shift implementation defect in the
family-manipulation code was corrected (see
`referee-report-2026-06-13.md`, item M1), had reported no
detectable separation under either architecture; that earlier,
pre-fix finding is superseded by the corrected result above.

We acknowledge that the present submission is a focused
factorial; the 144-cell extended grid (sample sizes, effect
sizes, DGP-side carryover crossed with family) anticipated in
the project pre-registration is deferred as future work. The
focused 16-cell factorial answers the central scientific
question with statistical clarity, and is the appropriate
scope for a Methods Note.

This manuscript is original, has not been published
elsewhere, and is not under consideration by another journal.
All authors have approved the submission. The work received
no specific funding; the pmsimstats team declares no
conflicts of interest. The data and code that support the
findings are openly available; specific paths are documented
in the Reproducibility section of the manuscript and an ADEMP
pre-registration document is committed at
`analysis/scripts/gompertz-evaluation/00-ademp-pre-registration.md`.

<!-- Reviewer suggestions to be added at submission. -->

Thank you for considering our submission.

Sincerely,

pmsimstats team
