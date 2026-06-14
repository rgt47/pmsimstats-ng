# Referee Report
*2026-06-13 16:49 PDT*

**Manuscript:** Two Architectures for Simulating Biomarker-Treatment
Interactions: Implications for Statistical Power Under Carryover  
**Authors:** pmsimstats team  
**Target journal:** Statistics in Medicine (inferred from CSL and style)

---

## 1. Summary of the Manuscript

The manuscript formalises two data-generating process (DGP)
architectures for simulation studies of predictive biomarker-treatment
interactions in aggregated N-of-1 clinical trials. Architecture A
(direct mean moderation) places the biomarker's moderating effect in
the population mean structure; Architecture B (MVN differential
correlation, following Hendrickson et al. 2020) places it in the
treatment-state-dependent covariance between biomarker and response.
The authors demonstrate by Monte Carlo simulation that under carryover,
Architecture B loses substantially more power than Architecture A in
open-label and hybrid designs (relative loss 30.6% vs. 2.8% in the
reference OL+BDC cell), while the traditional crossover design shows
no detectable difference between architectures. A third architecture
(C, combined channels) is introduced and characterised through a
3x3 parameter grid. The paper provides guidance for choosing among
architectures and situates the work in the broader HTE, crossover
trial, and precision medicine design literatures.

---

## 2. Overall Assessment and Recommendation

**Recommendation: Major Revision**

The paper addresses a genuine gap. The distinction between a
biomarker interaction residing in the first moment versus the second
moment of the outcome distribution is consequential for power
planning, and the simulation literature has not made this explicit.
The main empirical finding (Architecture B is substantially more
carryover-sensitive than Architecture A in certain designs) is
internally consistent with the data I was able to verify, and the
framing against the PATH Statement and precision medicine design
literatures is appropriate.

However, three problems require attention before the paper can be
accepted. First, the Architecture C boundary validation claim is
demonstrably false: the boundary cells of the combined simulation
do not reproduce the Section 3.1 baselines, which undermines the
stated internal consistency argument. Second, the z-statistics that
appear throughout the Results section are never defined, leaving the
reader unable to verify or replicate them. Third, Section 6 presents
six alternative analysis strategies as conceptual proposals but
draws on no simulation evidence; the qualitative comparison table
is not appropriate for a methods paper of this type. Beyond these,
several minor issues need attention.

---

## 3. Significance and Novelty

The paper makes a genuine conceptual contribution by naming and
formalising the Architecture A vs. B distinction. To my knowledge,
no prior paper in the N-of-1 or aggregated N-of-1 simulation
literature has made this explicit; Hendrickson et al. (2020) uses
Architecture B without comparing it to Architecture A. The
observation that Wang and Schork (2019) and Schork (2022) share
senior authorship with Hendrickson yet use Architecture A in
their own simulation work makes the distinction practically
important: it is not a matter of laboratory convention but a
paper-specific choice that the field has not scrutinised.

The Architecture C contribution is thinner. The 3x3 parameter
panel provides useful descriptive results, but the super-additivity
finding reported in Section 3.4 is presented without any analytical
grounding. It is plausible that super-additivity follows from the
mechanics of combining two independent signal channels, but neither
a derivation nor a simulation that isolates the interaction from
the main effects is offered. As presented, the Section 3.4 results
read as descriptive rather than explanatory.

The claim that Architecture B is unique to Hendrickson et al. in
the N-of-1 trial context is plausible but rests on the authors'
survey of the literature rather than on a systematic search. The
manuscript should state explicitly what searches were conducted and
what was excluded.

---

## 4. Major Comments

**MC-1. Architecture C boundary validation failure (Section 3.4,
Table 1, data files).**

The manuscript states (Section 3.4, orig block following the
simulation setup paragraph): "The boundary cells (c_bm_a = 0.45,
c_bm_b = 0) and (c_bm_a = 0, c_bm_b = 0.45) reproduce Section 3.1's
Architecture A and Architecture B results respectively, providing the
internal validation gate." The Table 1 caption repeats this claim:
"The boundary column c_bm_b = 0 reproduces pure-Architecture A results;
the boundary row c_bm_a = 0 reproduces pure-Architecture B results."

Inspection of the data files reveals that these claims are false.
At (c_bm_a = 0.45, c_bm_b = 0), the Architecture C simulation
yields power = 0.584 for OL+BDC at t_1/2 = 1.0 week. The Section 3.1
Architecture A result at the matched cell (N = 35, c.bm = 0.45,
OL+BDC, t_1/2 = 1.0) is 0.730 -- a 14.6 percentage-point discrepancy
that far exceeds Monte Carlo sampling error (MCSE ~ 0.015). The
corresponding CO discrepancy is larger still: 0.280 (Architecture C
boundary) vs. 0.725 (Section 3.1 Architecture A). At (c_bm_a = 0,
c_bm_b = 0.45), the Architecture C OL+BDC boundary gives 0.309
vs. 0.539 for Section 3.1 Architecture B.

The most likely explanation is that the Architecture C simulation
driver sources simulation-core.R from the carryover-sensitivity
scripts, which may use different Gompertz or biomarker distribution
parameters from the Section 3.1 production driver. Whatever the
cause, the boundary validation claim in the paper is incorrect as
stated.

This has three consequences. (a) Table 1 cannot be interpreted
as spanning from "pure A" to "pure B" in the sense of Section 3.1;
it spans from one DGP to another DGP in the same software framework,
but at different effective signal strengths. (b) The Section 3.4
discussion of super-additivity is potentially artefactual: the
comparison should be within the Architecture C grid alone, not
between the Architecture C grid and the Section 3.1 baselines.
(c) The paper's claim that Architecture C recovers Architectures A
and B as boundary cases -- which is the key theoretical justification
for the Architecture C framework -- is unverified by the simulation.

**Remedy.** Either (a) re-run the Architecture C simulation using
parameter files identical to the Section 3.1 production driver and
confirm that the boundary cells match to within 3 pp; or (b)
remove the boundary-validation claim, describe the Architecture C
parameterization as internally consistent but not directly comparable
to Section 3.1, and reinterpret the Section 3.4 results accordingly.
Option (a) is strongly preferred.

---

**MC-2. Z-statistics undefined (throughout Results).**

Z-statistics appear in the abstract (z = 7.64, p < 1e-13; z = 3.96,
p < 1e-4), the Section 3.1 narrative, and the conclusions. No
definition is given anywhere.

I have inferred from the data that these are tests of the difference
in carryover-induced power loss between architectures. Specifically,
the test appears to be:

  z = (Delta_B - Delta_A) / SE_pooled

where Delta_B = power_B(t1half=0) - power_B(t1half=1.0) and
Delta_A is the corresponding Architecture A loss, both computed
from 1000-replicate proportions, and SE_pooled is derived from
the four independent binomial standard errors. Under this formula:

  Delta_B (OL+BDC) = 0.777 - 0.539 = 0.238
  Delta_A (OL+BDC) = 0.751 - 0.730 = 0.021
  SE_pooled = sqrt(0.777*0.223/1000 + 0.539*0.461/1000
                   + 0.751*0.249/1000 + 0.730*0.270/1000)
            = sqrt(0.001733 + 0.002485 + 0.001870 + 0.001971)/sqrt(10)
  -- computing properly --
  SE(Delta_B) = sqrt(p_B0*(1-p_B0)/n + p_B1*(1-p_B1)/n) = 0.02053
  SE(Delta_A) = sqrt(p_A0*(1-p_A0)/n + p_A1*(1-p_A1)/n) = 0.01960
  SE_pooled = sqrt(0.02053^2 + 0.01960^2) = 0.02839
  z = 0.217 / 0.02839 = 7.64 ✓

The same formula recovers z = 3.96 for the Hybrid design. This is
not obvious to a reader and should be stated explicitly in the
methods (Section 5.2 or a dedicated Statistical Methods subsection).

**Remedy.** Add a brief Statistical Methods paragraph defining the
test used for architecture comparisons. State whether this is a
two-sided or one-sided test. Note that the four proportions are
computed from independent Monte Carlo replications; if replicates
within a cell are not independent (shared seeds across cells), this
assumption should be verified.

---

**MC-3. Section 6 analysis strategies lack simulation support.**

Section 6 describes six alternative analysis strategies for
recovering power under Architecture B carryover. The comparative
table (Section 6.7) assigns qualitative expected benefits ("High",
"Moderate", "Unknown") with no simulation evidence. The closing
paragraph states "We regard this as a direction for future work."

For a methods paper in Statistics in Medicine, a section presenting
six analysis strategies as proposals requires at minimum preliminary
simulation evidence for the most promising candidates. The qualitative
table provides no basis for the reader to choose among the strategies.
Moreover, two of the strategies (on-drug-only restriction, weighted
analysis) are straightforward to implement and could be evaluated
using the existing simulation infrastructure in a single additional
factorial run.

**Remedy.** Either (a) run simulations for at least the two most
tractable strategies (on-drug only, weighted analysis) and present
quantitative power recovery estimates; or (b) re-frame Section 6
explicitly as a conceptual discussion, remove the comparison table
or replace it with a discussion of the relative plausibility of each
mechanism, and be clear in the abstract that the analysis-strategy
evaluation is outside the paper's scope. The current mixed framing
-- substantive strategies, substantive table, but no evidence -- is
misleading.

---

## 5. Minor Comments

**Mn-1. Super-additivity benchmark (Section 3.4).**

The orig block for the second bullets group in Section 3.4 states:
"power is 0.567, whereas the arithmetic mean of the corresponding
pure-architecture values (0.198 for pure A, 0.152 for pure B) is
only 0.175." The arithmetic mean of individual powers is not the
correct additivity baseline on the probability scale. The additive
benchmark (assuming independence of the two channels) would be
approximately p(c_bm_a, 0) + p(0, c_bm_b) - p(0, 0) = 0.198 +
0.152 - 0.025 = 0.325. The observed power 0.567 remains
super-additive relative to this corrected baseline (0.567 > 0.325),
so the qualitative conclusion stands. But the reported comparison to
0.175 overstates the degree of super-additivity. Replace with the
probability-scale additive benchmark.

Note: this comment assumes the Architecture C grid is internally
consistent (it appears to be), even if the boundary cells do not
match Section 3.1 (MC-1 above).

---

**Mn-2. Architecture B CO power non-monotone (Section 3.1 table).**

In the Section 3.1 production results, Architecture B (MVN) power
at the CO design, c.bm = 0.45 is:

  t_1/2 = 0.0 weeks: 0.745
  t_1/2 = 0.5 weeks: 0.754
  t_1/2 = 1.0 weeks: 0.723

The power at t_1/2 = 0.5 weeks is higher than at t_1/2 = 0, which
is inconsistent with the narrative that carryover erodes power under
Architecture B. The corresponding Architecture A values (0.739, 0.739,
0.725) are monotone. The non-monotonicity at CO is likely Monte Carlo
noise (MCSE ~ 0.014 at this cell), but it is unexplained and
noticeable given the paper's emphasis on Architecture B's carryover
sensitivity. The paper should acknowledge this as a sampling artifact
and note that the CO architecture-by-carryover contrast is
statistically indistinguishable from zero (which is the paper's
conclusion for the CO design), not a reversal of the expected pattern.

---

**Mn-3. Reproducibility section inaccuracies (Section 8 /
Reproducibility block).**

The Reproducibility section states "R 4.4.x inside the project's
zzcollab Docker container." The Dockerfile and CI configuration in
the repository use `rocker/tidyverse:4.6.0`, placing the actual
R version at 4.6.0. Update this claim.

The same section describes "8-core parallel::mclapply" for
parameter-grid sweeps. The Architecture C simulation driver uses
`furrr::future_map_dfr` with `plan(multicore)`, not `mclapply`
directly. The Section 3.1 driver may differ; the description should
be specific to each simulation.

The Architecture C data paths (`analysis/data/dgp-combined/`) are
not mentioned in the reproducibility section.

---

**Mn-4. Z-statistics for Section 3.1 within-architecture losses
(Section 3.3).**

Section 3.3 reports "Losses: 9.3 pp (Architecture B, p < 10^{-12})
vs. 1.7 pp (Architecture A, p ≈ 0.05); ordering matches N = 35."
The p-values here appear to be one-sample tests of whether the
N = 70 power loss (a difference of two proportions, each from 800
replicates) is significantly different from zero. These test
statistics are also undefined. The same definitional comment as
MC-2 applies.

---

**Mn-5. Varadhan et al. (2013) appears to be an orphaned reference.**

I did not find a citation to `@varadhan2013` in the manuscript text.
The entry appears in the .bib file and in the References section but
was not cited. Either cite it or remove it.

---

**Mn-6. Schork (2022) bibliographic entry incomplete.**

The Schork (2022) entry in references.bib contains only
`note = {Special Issue 3}` with no volume, number, or page fields.
The Harvard Data Science Review article at issue is:
"Accommodating serial correlation and sequential design elements in
personalized studies and aggregated personalized studies."
Harv Data Sci Rev. 2022;Special Issue 3. This is cited in the text
as a peer-reviewed article; the citation should conform to the
journal's reference style (issue/article identifier at minimum).

---

**Mn-7. Rizopoulos (2012) citation context.**

Section 7.2 cites Rizopoulos (2012) as an example of Architecture B
appearing in the literature. Joint models for longitudinal and
time-to-event data are a different application from biomarker-treatment
interaction in N-of-1 trials; the differential correlation structure
in those models relates longitudinal biomarkers to survival, not to
treatment response. The citation is not incorrect, but it is a stretch
that may mislead readers about how Architecture B relates to the joint
modeling tradition. A brief qualifier ("a context where
treatment-state-dependent covariance appears, albeit in a different
form from Architecture B here") would help.

---

**Mn-8. Mixed source-level markup.**

The manuscript source code uses raw LaTeX environments
(`\begin{bullets}`, `\begin{orig}`, `\begin{rgt}`) in sections
1-3 and pandoc fenced divs (`::: {.bullets data-latex=""}`) in
sections 4-7. Both render equivalently in the current xelatex build,
but the inconsistency in the source is confusing and could cause
problems if the build configuration changes. Standardise on one
markup style throughout.

---

## 6. Missing and Questionable References

| Location | Issue | Suggested action |
|---|---|---|
| Throughout | AR(1) residual correlation (corCAR1) is chosen without citation | Cite Pinheiro & Bates (2000) *Mixed-Effects Models in S and S-PLUS*, or the nlme package documentation with Jones & Kenward as context for temporal correlation in crossover models |
| §2, DGP description | Exponential carryover decay exp(-lambda * t_sd) is parameterised from pharmacokinetics but not cited | Cite a PK textbook or the original prazosin PK literature justifying this decay form |
| §4.2, finite mixture approximation | Architecture B is described as a "tractable second-moment approximation to a richer mixture model" without citing mixture model power analysis literature | McLachlan & Peel (2000) *Finite Mixture Models* or Muthén & Shedden (1999) on latent class analysis with treatment interactions |
| §4.3, prazosin-SBP mechanism | The neurobiological hypothesis (elevated SBP indexes noradrenergic tone) needs a mechanistic reference | A review of alpha-1 adrenergic signaling in PTSD, e.g., Raskind et al. (2018) supplement or a pharmacological reference |
| §7.2 | The claim that Architecture B appears in "Bayesian joint modeling frameworks" is asserted without a specific citation | Provide at least one specific reference |

---

## 7. Suggestions for Strengthening the Paper

1. **Resolve the boundary validation issue (MC-1) and add a
   validation table** to Section 3.4 that explicitly shows
   Architecture C boundary cell values alongside the Section 3.1
   baselines, whether or not they agree. If they differ, explain
   why. This transparency would itself be a contribution: it
   documents what "reproducing a pure architecture" means
   operationally when two simulation drivers are involved.

2. **Derive the super-additivity result analytically (or at least
   heuristically).** Under the assumption that the two channels
   contribute independently, the joint power should exceed the
   additive benchmark when the channels share a common signal path
   (the analysis model). A few lines of argument showing why the
   joint BM-BR correlation under combined DGP exceeds the sum of
   individual contributions would substantially strengthen Section
   3.4.

3. **Add a simulation for at least one Section 6 strategy.** The
   on-drug-only restriction (Section 6.1) is the simplest to
   implement within the existing framework. A single additional
   factorial run at N = 35 for the OL+BDC and Hybrid designs under
   Architecture B, comparing the A2 specification against
   on-drug-only, would provide the empirical anchor Section 6
   currently lacks.

4. **Report Monte Carlo standard errors explicitly in all tables.**
   Tables 1 and the Section 3.3 robustness table report point
   estimates only. MCSE is stated in the text (0.013-0.016) but not
   shown cell-by-cell. Adding MCSE columns (or parenthetical standard
   errors) would allow readers to assess statistical uncertainty
   directly from the tables, consistent with the ADEMP framework
   the paper endorses.

5. **Consider a power-versus-carryover-halflife curve figure** for
   the OL+BDC and Hybrid designs, showing Architecture A and B
   power as continuous functions of t_1/2 rather than at three
   discrete points. This would make the architectural divergence
   visually immediate and would be useful to applied trialists who
   want to read off their specific carryover duration.

6. **Document all literature searches conducted for Section 7.**
   State what databases were searched, with what queries, and over
   what date range. The claim that Hendrickson et al. is "unique"
   in using Architecture B for N-of-1 trials is important for the
   novelty argument and should rest on a systematic rather than
   informal survey.

---

## 8. Scope of this Review

**Verified** (re-derived or cross-checked against data files):

- All power values in the abstract and Section 3.1 against
  `analysis/data/quick-sim/01-dgp-summary.txt`.
- Z-statistics z = 7.64 and z = 3.96 via the formula inferred
  from the data (see MC-2); both reproduce to three significant
  figures.
- Relative-loss percentages (30.6%, 2.8%, 9.3 pp, 1.7 pp) against
  the data file.
- Section 3.3 N = 70 within-architecture losses (9.3 pp and 1.7 pp)
  are arithmetically consistent with the reported Table values.
- Architecture C boundary discrepancy: cross-checked
  `analysis/data/dgp-combined/02-power-grid-olbdc.csv` against
  `01-dgp-summary.txt`; discrepancies confirmed.
- Architecture C simulation driver (`01-run-combined-factorial.R`)
  and summarisation script (`02-summarise-combined.R`) read in full;
  the validation gate was confirmed to print but not compare values.

**Inspected** (read and judged plausible but not re-derived):

- Section 2 mathematical definitions of Architectures A, B, C.
- Section 3.2 mechanistic explanation (dual attack on Architecture B
  power) -- the argument is coherent and consistent with the
  simulation evidence.
- Section 3.3 N = 70 robustness table -- reported values appear
  consistent with ceiling saturation as argued.
- Section 4 biological discussion -- appropriate level of
  qualification; appropriately references Senn on variance artifacts.
- Section 5 decision table.
- Section 7 literature framing -- plausible but not systematically
  verified (see below).

**Not checked or suspected but unconfirmed**:

- The `simulation-core.R` function sourced by the Architecture C
  driver was not read; whether it uses different Gompertz parameters
  from the Section 3.1 driver is the most likely explanation for
  MC-1 but was not confirmed.
- The specific prazosin SBP mechanistic literature was not
  independently reviewed.
- Whether any post-2020 work on differential correlation in
  precision medicine simulation exists was not confirmed by
  literature search.

**Literature searches run**: None via external search engine. The
Section 7 literature framing was assessed against the references in
the manuscript's own bibliography plus general knowledge of the
field. Formal literature searches (PubMed, Google Scholar) for
competing work were not conducted; this is a scope limitation of
this review.

---

*Review completed 2026-06-13.*
