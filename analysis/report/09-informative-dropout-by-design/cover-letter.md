# Cover letter
*2026-05-09 09:55 PDT; revised 2026-10-29 11:35 PDT*

To the Editor,
Statistics in Medicine

Dear Editor:

We submit the manuscript 'Informative Dropout and Trial-Design Choice
in Aggregated N-of-1 Biomarker-Validation Trials: A Paired Simulation
at One Parameter Point' for consideration as a research article in
Statistics in Medicine.

The manuscript asks how participant dropout affects the test of a
biomarker-by-treatment interaction in three aggregated N-of-1 designs:
open-label titration with blinded discontinuation (OL+BDC), the
two-period crossover (CO) and the Hybrid design of Hendrickson et al.
(2020). Seven dropout conditions, from none to 40 percent, were
applied to the same 1,000 complete datasets per design. Losses and
coefficient shifts are therefore paired differences with small Monte
Carlo standard errors. The dropout mechanisms are defined exactly and
include worsening-driven, improvement-driven and symptom-return
hazards. Reporting follows the ADEMP framework of Morris, White and
Crowther (2019).

The principal findings, at one parameter point with N = 70, are as
follows. (i) The Hybrid and OL+BDC designs estimated the interaction
with the same precision. CO needed about 1.9 times the sample for the
same precision. (ii) Dropout of a given rate inflated the standard
error by nearly the same factor in all three designs (1.08 to 1.09 at
25 percent completely random dropout, 1.13 to 1.16 at 40 percent
magnitude-driven dropout). CO's larger absolute power loss followed
mainly from its lower starting power. (iii) Improvement-driven dropout
attenuated the estimate by up to 7 percent of its true value, most in
CO. Dropout did not change test size. (iv) An analysis half-life of
0.5 week biased the OL+BDC and Hybrid estimates by about 11 and 8
percent, because the generating interaction had no carryover. A
binary exposure indicator removed the bias.

The study is exploratory. No analysis plan was registered in advance,
and its findings should be confirmed in an independent study. It
covers one parameter point and one encoding of the interaction. It
does not include shared-parameter dropout or analysis-stage
sensitivity analyses. These limits are stated in the manuscript.

This manuscript is original, has not been published elsewhere, and is
not under consideration by another journal. All authors have approved
the submission.

The work received no specific funding, and the pmsimstats team
declares no conflicts of interest. The data and code that support the
findings are in the project repository. The driver is
`analysis/scripts/informative-dropout-by-design/20-paired-dropout.R`,
and the working ADEMP protocol is
`analysis/scripts/informative-dropout-by-design/00-ademp-protocol.md`.

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
