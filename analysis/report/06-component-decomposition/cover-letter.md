# Cover letter
*2026-05-09 09:55 PDT*

To the Editor,
Statistics in Medicine

Dear Editor:

We are pleased to submit our manuscript titled 'Three-Component
Decomposition of Treatment Response in Aggregated N-of-1 Trials:
Pharmacological, Expectation-Driven, and Natural-History
Components' for consideration as a research article in
Statistics in Medicine.

The treatment response in aggregated N-of-1 clinical trials is
not a single quantity but a sum of three causally distinct
contributions: pharmacological (BR), expectation-driven /
placebo-belief (PB), and time-variant natural-history (TV).
Routine analyses lump these components and report the total,
which is adequate for symptom-control endpoints but is not the
correct estimand when the inferential target is a predictive
biomarker. A biomarker that predicts the lumped total is not, in
general, a pharmacological predictor, and prescribing decisions
made on a pharmacological basis would mis-stratify patients if
the biomarker correlates strongly with placebo-belief or
natural-history components. The decomposition is therefore a
substantive methodological requirement for predictive-biomarker
validation in aggregated N-of-1 trials, not a parameterisation
convenience.

This manuscript formalises the three-component decomposition,
makes the trial-design contrasts that identify each component
explicit, supplies a phase-augmented linear-mixed analysis that
is identifiable at trial-relevant sample sizes, and clarifies
the target estimand for biomarker validation as the regression
of $BR$ on the biomarker rather than the lumped total. A
proof-of-concept simulation prototype (100 replicates per cell
on a 6-cell slice) confirms the analysis machinery runs
end-to-end on the three-component data-generating process. A
production simulation programme is specified in §6 of the
manuscript to quantify bias-versus-power trade-offs across
sample sizes $N \in \{35, 70, 100, 150\}$; the programme is
ADEMP pre-registered and is the subject of forthcoming
follow-up work.

The paper aligns with Statistics in Medicine's editorial focus
on methodological developments for clinical trials and provides
the formal framework through which predictive-biomarker claims
in aggregated N-of-1 trials can be evaluated rigorously. It
complements companion submissions on the data-generating-process
architecture, the carryover-mitigation analysis strategies, the
test-procedure choice, and the trajectory-template choice; we
anticipate submitting those papers separately.

This manuscript is original, has not been published elsewhere,
and is not under consideration by another journal. All authors
have approved the submission.

We confirm that the work received no specific funding and that
the pmsimstats team declares no conflicts of interest. The data
and code that support the findings are openly available; the
ADEMP pre-registration is committed at
`analysis/scripts/component-decomposition/00-ademp-pre-registration.md`.

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
