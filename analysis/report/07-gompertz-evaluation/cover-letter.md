# Cover letter
*2026-10-28 18:01 PDT (revised to the common-random-numbers study;
see `referee-report-2026-10-28.md`)*

To the Editor,
Statistics in Medicine

Dear Editor:

We submit our manuscript 'Response-curve family and the
biomarker-treatment interaction test in N-of-1 trial simulation: a
common-random-numbers comparison of the modified Gompertz template
with three alternatives' for consideration as a Methods Note.

Simulation-based power analyses for aggregated N-of-1 trials need a
parametric template for the within-participant response trajectory.
The pmsimstats framework uses a modified Gompertz curve. The choice is
fixed before any analysis is run, and every reported power figure is
conditional on it. Existing N-of-1 simulators each fix one response
form; none compares forms. The present manuscript makes that
comparison.

We compared the modified Gompertz curve with the symmetric logistic,
the hyperbolic tangent and a piecewise-linear breakpoint model, each
matched at the asymptote and at a week-5 anchor, in the two-path
open-label plus blinded-discontinuation design with 70 participants
and no carryover, under both a mean-shift and a correlation encoding
of the biomarker-treatment interaction. The four families used common
random numbers, so every contrast is paired. We report paired power
differences with paired Monte Carlo standard errors, Type I error from
2,000 null replicates per family, and bias against the closed-form
interaction slope.

The finding is that families that coincide at the sampled visits give
identical decisions. The calibrated logistic equals the Gompertz curve
to within 0.12 points at every sampled visit and returned the same
decision in all 1,000 seeds under both encodings. The hyperbolic
tangent, which is the logistic at twice the rate, lowered power by
about 0.01 under both encodings. The correlation encoding is nearly
but not exactly inert to the family. The test ran below its nominal
size, at 0.034 against 0.05. The result is conditional on a visit schedule that
samples the rising phase of the response at two visits, which the
manuscript states as its main limitation.

An earlier version of this study, with independent seeds per cell,
reported a 0.039 power spread with Gompertz the lowest-power family.
That spread was Monte Carlo noise, and the manuscript says so. We
think the corrected result, and the paired design that exposed the
error, are useful to anyone comparing data-generating processes in a
simulation study.

The study is exploratory, and its findings should be confirmed in an
independent study. The 144-cell grid in the working protocol is
deferred.

This manuscript is original, has not been published elsewhere, and is
not under consideration by another journal. All authors have approved
the submission. The work received no specific funding; the authors
declare no conflicts of interest. The author list overlaps with that
of the framework evaluated. The data and code that support the
findings are in the project repository; specific paths are documented
in the Reproducibility section of the manuscript, and the ADEMP
working protocol is at
`analysis/scripts/gompertz-evaluation/00-ademp-protocol.md`.

<!-- Reviewer suggestions to be added at submission. -->

Thank you for considering our submission.

Sincerely,

pmsimstats team
