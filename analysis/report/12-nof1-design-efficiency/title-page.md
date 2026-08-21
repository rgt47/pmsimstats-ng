---
title: "Title page"
output:
  pdf_document:
    latex_engine: xelatex
fontsize: 11pt
geometry: "left=2.5cm,right=2.5cm,top=2.5cm,bottom=2.5cm"
---

# A Closed-Form Approximation to the Biomarker-Treatment Interaction Estimator in Aggregated N-of-1 Trials, and What It Implies for Design Efficiency Under a Mandatory Open-Label Lead-In

**Short communication.**

**Authors.** Temperance Persons, on behalf of the pmsimstats team.

\bigskip

Temperance Persons (corresponding author). The pmsimstats team
is a collaborative project; the full author roster, ORCIDs, and
affiliations will be recorded prior to submission.

\bigskip

**Affiliations.**

Affiliation details will accompany the resolved author roster.

\bigskip

**Corresponding author.**

Temperance Persons, tpersons@gmail.com. Postal address to be
recorded prior to submission.

\bigskip

**Keywords.** N-of-1 trials; biomarker-treatment interaction;
closed-form power approximation; carryover; trial design
efficiency; AR(1) correlation; aggregated N-of-1 design.

\bigskip

**Word count.** To be recorded prior to submission.

\bigskip

**Running title.** Closed-form interaction power in N-of-1 trials.

\bigskip

**Funding.** This work received no specific funding.

\bigskip

**Conflict of interest.** The pmsimstats team declares no conflicts of interest.

\bigskip

**Data availability statement.** The validation dataset is the
existing output of
`analysis/scripts/dgp-combined/01-run-combined-factorial.R`
(`analysis/data/dgp-combined/01-combined-summary.rds`), already
produced for the companion dual-channel-architecture manuscript.
The closed-form derivation is checked against it by
`analysis/scripts/nof1-design-efficiency/01-closed-form-validation.R`,
which reads that file and writes
`analysis/data/nof1-design-efficiency/closed-form-validation.{rds,csv}`.
No new Monte Carlo simulation was run for this manuscript. All
scripts are openly available in the pmsimstats-ng project
repository at https://github.com/[organisation]/pmsimstats-ng
(public-release organisation pending).

\bigskip

**Author contributions.** Author contributions will be recorded with the final author roster.

\bigskip

**Acknowledgements.** This manuscript synthesizes and extends
results from ten companion manuscripts in the same project
(cited throughout); its numerical contribution was checked
against Monte Carlo output already generated for the
dual-channel-architecture companion paper rather than from a
freshly commissioned simulation study, which is the appropriate
scope for a short communication.
