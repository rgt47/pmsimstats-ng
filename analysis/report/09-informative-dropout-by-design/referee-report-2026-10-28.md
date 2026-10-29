# Referee report: paper 09, informative dropout and design choice
*2026-10-28 22:29 PDT*

**Manuscript.** 'Informative dropout and trial-design choice in
aggregated N-of-1 biomarker-validation trials'
(`analysis/report/09-informative-dropout-by-design/report.Rmd`,
source dated 2026-10-08).

**Target journal considered.** Statistics in Medicine (the cover
letter's target); comments on suitability for Trials are given in
Section 3.

**Reviewer.** pmsimstats team (internal referee pass).

**Evidence labels.** VERIFIED: I executed code against the stored
artifacts or the current package and observed the result. INSPECTED:
I read the source or document and confirmed the point. SUSPECTED:
the evidence points this way but I did not establish it. UNVERIFIED:
not checked.

**Scratch scripts.** The verification script is
`/private/tmp/claude-501/-Users-zenn-prj-res-36-pmsimstats-ng-pmsimstats-ng-analysis-report-01-dgp-mean-moderation-vs-mvn/513155aa-ec85-4a1d-889a-5d57e8607d18/scratchpad/p09-check.R`
(run from the project root; about one minute; no long simulation).

## 1. Summary

The manuscript asks whether four aggregated N-of-1 designs (open
label, OL; open label plus blinded discontinuation, OL+BDC;
two-period crossover, CO; and the four-path Hybrid design of
Hendrickson et al. 2020) differ in robustness to informative dropout
for the biomarker-by-treatment interaction $\beta_{bm:D}$. It
proposes a missing-data formalization of the Hendrickson hazard
model, describes design-specific 'vulnerability windows', sets out a
four-study ADEMP protocol (372 cells), and reports a 16-cell
prototype: four designs by four dropout conditions (none, 25 percent
MCAR, 25 percent 'biased', 40 percent 'biased'), total $N = 70$,
mean-shift encoding with $\beta_{bm} = 0.45$, carryover half-life
0.5 week, 500 replicates per cell. The headline results are a power
ordering Hybrid (0.968) > OL+BDC (0.950) > CO (0.854) > OL (0.088)
that persists under dropout, a larger absolute power loss for CO
under 40 percent biased dropout (8.6 versus 5.8 and 5.6 points), a
claim that dropout rate dominates mechanism, and a claim that
dropout costs precision rather than accuracy.

The stored artifact reproduces every number in the Results power
and coefficient tables. The problems lie elsewhere. The dropout
mechanism that was simulated is not the one the paper formalizes
and names. The 'worsening-driven' condition is driven by the
squared magnitude of symptom change and in fact removes the most
improved participants slightly more often. The OL row tests a
different coefficient from the other three rows. The OL+BDC and
Hybrid estimates carry a 10 percent bias relative to the true
interaction in every cell, including the no-dropout cell. Several
sections of the Discussion, the Conclusions and the cover letter
still report the superseded unmatched run. Three of the four
contributions promised in the Introduction (the formal mechanism
study, the misspecification study, and the path decomposition) were
not carried out.

## 2. Overall assessment and recommendation

**Recommendation: reject in present form; encourage resubmission
as a substantially reworked paper.**

The question is relevant and the prototype is carefully matched on
total $N$, which corrects an earlier confound. The analysis code and
the stored output agree. These are real strengths. They do not
offset the following. (i) The paper's central conceptual object,
the dropout mechanism, is misdescribed. Section 2 formalizes
`censordata()`, which the driver never calls. The driver's own
mechanism is symmetric in the change from baseline and is not
worsening-driven. (ii) The design comparison mixes three effects:
information, a design-dependent estimation bias induced by a
mismatch between the generating and the analysis carryover, and a
design-dependent conservativeness of the model-based standard
error. (iii) The robustness ranking that the abstract and the
recommendation table rely on (CO 'clearly separated') is not
statistically supported: the difference in losses between CO and
either top design is about 3 points with a standard error of about
3 points. (iv) The manuscript is internally inconsistent, with
stale numbers from the unmatched run in the Discussion and the
Conclusions. (v) Literature attributions to Hendrickson et al. are
inaccurate in ways that change the paper's motivation. Each of these
is fixable, but together they require new simulations, a rewritten
Section 2, and a narrower set of claims. A single parameter point
also falls short of what Statistics in Medicine expects of a
design-comparison simulation study.

## 3. Significance and novelty

The topic is legitimate. Dropout in N-of-1 and multi-period
crossover designs is under-studied, and the estimand framework of
ICH E9(R1) has not been worked through for aggregated N-of-1
biomarker trials. A paper that did this well would interest
Statistics in Medicine, Trials or Clinical Trials.

The novelty claims are overstated, however.

- Hendrickson et al. (2020) already compared the same four designs,
  at the same visit schedules, under the same family of censoring
  patterns (their Figures 4, 6 and 7). The present prototype covers
  one point of their parameter space with a different dropout
  mechanism. 'Reproduce and extend Figure 4A' is not what was done
  (INSPECTED: driver and Hendrickson text).
- Matthews and Henderson (2013, Biostatistics 14:626-638) studied
  the robustness of MAR-based analysis of AB/BA crossover designs
  when dropout is non-ignorable and compared alternative sequence
  sets. This is the closest methodological antecedent for the
  'design-by-dropout coupling' and is not cited. The claim that
  'none of this machinery has been imported ... in a
  design-comparative way' does not hold for the crossover family.
- The paper claims an 'import' of pattern-mixture, selection-model
  and multiple-imputation machinery, but uses none of it. The only
  analysis is the likelihood-based mixed model under an implicit
  MAR assumption. No delta-adjusted, tipping-point, reference-based
  or shared-parameter analysis appears.

Contemporariness: the reference list stops at 2023 for N-of-1
material and at 2015 for estimands. Post-2015 sensitivity-analysis
practice (Cro et al. 2020; Permutt 2016) and recent N-of-1 design
work in Statistics in Medicine and Clinical Trials (sequential
monitoring of N-of-1 series, 2025) are absent.

Fit: as a single-scenario design comparison with practical
recommendations, the material is closer to a Trials methodology
article than to a Statistics in Medicine research article, but in
either venue the issues in Section 4 must be resolved first.

## 4. Major comments

Correctness items come first.

**M1. The dropout mechanism simulated is not the one Section 2
formalizes, and the manuscript says the opposite.**
*Location.* Section 2 (lines 329-389); Simulation design (lines
519-531: 'not by any dropout-specific code written for this
manuscript'); Reproducibility (line 926); title-page data statement.
*Problem.* Section 2 derives the MNAR property from
`R/censordata.R` (union of a flat and a biased Bernoulli draw,
$(\Delta_{Sx}+100)^{e_{b1}}$, patterns indexed by
$(\beta_0, \beta_1)$ such as 'balanced' and 'more-of-biased'). The
driver `analysis/scripts/quick-sim/09-dropout-prototype-10min.R`
never calls `censordata()`. It defines its own `apply_dropout()`
with per-visit hazard $\alpha\,(\beta_0 + \beta_1 z_{it}^2)$, where
$z_{it}$ is the change from baseline divided by its cross-sectional
standard deviation at that visit, $\beta_0 = 0.05$, $\beta_1 = 0.5$
fixed, and $\alpha$ solved by `uniroot` on each realized dataset so
that the expected ever-dropout fraction equals the target. None of
the five $(\beta_0, \beta_1)$ patterns of Section 2 and of the
protocol was run. The three dropout labels in the Results
(`MCAR_25`, `biased_25`, `biased_40`) do not appear in Section 2.
*Why it matters.* The paper's formal content (Section 2) and its
empirical content (Section 5) are about different mechanisms. The
statement that no dropout-specific code was written is false.
*Remedy.* Either run `censordata()` with the stated patterns, or
rewrite Section 2 to formalize `apply_dropout()` exactly (hazard,
standardization, calibration of $\alpha$ on the realized data,
monotone carry-forward), and delete the claim about existing code.
*Evidence.* INSPECTED (driver lines 132-206; `R/censordata.R`).

**M2. The 'biased' condition is not worsening-driven. It depends on
the squared magnitude of change and slightly favors removal of the
most improved participants.**
*Location.* Abstract ('worsening-driven dropout', three times);
Results; recommendation table rows 'Worsening-driven' and
'Improvement-driven'; Conclusions.
*Problem.* Because $z_{it}^2$ is symmetric and the change from
baseline is not centered, the hazard increases with $|\Delta|$. Mean
change is strongly negative (improvement; about -10 to -28 points
by visit in OL+BDC), so large $|\Delta|$ is mostly large
improvement. In a 4,000-participant OL+BDC path at the 25 percent
target, the ever-dropout fraction by quartile of mean improvement
was 0.224 (least improved), 0.236, 0.254, 0.283 (most improved).
*Why it matters.* Every statement that labels a result
'worsening-driven' is mislabeled. The recommendation table's
'Improvement-driven' row, and its rationale that CO 'loses signal
when responders drop early', rest on no simulated cell at all. The
weak, symmetric selection also explains why the coefficient is
stable (M6): selection on $z^2$ removes high and low responders
alike and largely cancels in the interaction.
*Remedy.* Relabel the condition as 'magnitude-of-change dropout'.
Add genuinely directional mechanisms: worsening-driven (hazard
increasing in $\Delta$), improvement-driven, and a mechanism tied to
the latent drug response or to off-drug symptom return in the
discontinuation phase, which is the clinical scenario the
Introduction motivates. Delete the improvement-driven row or support
it with cells.
*Evidence.* VERIFIED (script section 8; one large simulated path,
not a replicated estimate).

**M3. The 'biased' component of `censordata()` itself contradicts
its documentation and Section 2.**
*Location.* Section 2, lines 339-344 ('a large improvement lowers
the dropout probability').
*Problem.* In `generateData()` the columns `D_<visit>` equal
`rowSums` of the components, that is $BL - Sx$, the improvement.
`censordata()` uses $(D + 100)^{e_{b1}}$, which increases with
improvement. The same code is in `implementations/original`. The
flat component is proportional to `t_wk`, the block length, not to
elapsed time as Section 2 states.
*Why it matters.* If the authors adopt `censordata()` to fix M1,
its biased arm is improvement-driven, the reverse of the text and of
Hendrickson et al.'s verbal description. Hendrickson et al. describe
a hazard 'higher for participants experiencing worsening'. The code
does the opposite, which affects any claim to reproduce their
Figure 4A.
*Remedy.* State the sign explicitly, fix or document it, and report
the dropout-by-improvement gradient of whichever mechanism is used.
*Evidence.* INSPECTED (`R/generateData.R` lines 285-294;
`R/censordata.R` lines 51-62). The sign was not checked against
Hendrickson's commit 58b32a9 beyond the vendored `original` copy.

**M4. The OL row tests a different coefficient, and the stated
explanation (aliasing of $bm{:}D_{bc}$ with $bm$) is wrong.**
*Location.* Abstract lines 80-84; Results lines 594-615; Discussion
lines 793-801; recommendation table.
*Problem.* In `lme_analysis()`, a design in which no participant
varies in exposure after baseline sets `varInDb = FALSE` and fits
`Sx ~ bm + t + bm:t`. For OL the reported 'interaction' is the
biomarker-by-time slope (mean -0.09 points per week), not
$\beta_{bm:D}$. The baseline visit at $t = 0$ is in the data, so
$bm{:}D_{bc}$ would not be aliased in any case. The OL rejection
rate of 0.086 to 0.112 is power for a misspecified slope test under
the alternative, neither a test size nor an 'information floor'.
The Discussion also says the OL row 'sits at nominal alpha' and is
'an implicit null check', contradicting the Results.
*Why it matters.* The OL row is not comparable to the other three,
and both explanations offered for it are incorrect.
*Remedy.* Either drop OL from the comparison and say why (no
off-drug contrast, so $\beta_{bm:D}$ is not separable from the
natural-history and placebo components), or fit the same
$bm{:}D_{bc}$ model with the baseline row and report it. Remove all
'null check' language.
*Evidence.* VERIFIED (script section 6: fitted formula printed as
`Sx ~ bm + t + bm:t`).

**M5. OL+BDC and Hybrid estimates are biased by about 10 percent in
every cell, because the generating and analysis carryover models
differ. 'Informative dropout cost precision rather than accuracy'
is therefore stated against the wrong reference.**
*Location.* Abstract lines 91-98; Section 5.2; Conclusions point 3;
cover letter.
*Problem.* Under `dgp_architecture = 'mean_moderation'` with
`moderation_scaling = 'constant'`, the biomarker shift
$\beta_{bm} b_i \sigma_{BR}$ is applied only at on-drug visits and
has no carryover. The true $\beta_{bm:D}$ is therefore
$-0.45 \times 5 = -2.25$ points per biomarker unit. The analysis
builds $D_{bc}$ with $t_{1/2}^a = 0.5$ week, so off-drug visits a
week after discontinuation carry $D_{bc} = 0.25$ where the true
shift is zero. The stored means are -2.48 (OL+BDC) and -2.47
(Hybrid) with MCSE 0.028, a bias of about -0.23 (8 MCSE) already
in the no-dropout cell. CO is close to the truth (-2.24). In single
large-$N$ fits ($N = 4000$, no dropout), $t_{1/2}^a = 0$ gave
-2.25, -2.19 and -2.30 for OL+BDC, Hybrid and CO. $t_{1/2}^a = 0.5$
gave -2.72, -2.47 and -2.09.
*Why it matters.* The manuscript never states the true value and
never computes the planned 'bias relative to the known true value'
(protocol Study 2). 'No bias' means only 'no change across dropout
cells'. More importantly, the inflated $|\hat\beta|$ raises power
for OL+BDC and Hybrid relative to CO. Part of the headline ordering
is an artifact of analysis misspecification, not of design
information.
*Remedy.* Report the true value and bias against it. Align the
generating and analysis carryover for the interaction (for example,
apply carryover to the moderation shift, or use
`moderation_scaling = 'trajectory'`), or analyze at the project
policy value $t_{1/2}^a = 0.2$ and with an AIC-selected half-life
(paper 02's G9) and report all three.
*Evidence.* VERIFIED for the stored bias (script section 1) and for
the true value (code inspection plus arithmetic). The attribution to
the carryover mismatch is supported by single large-$N$ replicates
(standard errors 0.08 to 0.10), so it is likely but not established
by replicated simulation.

**M6. The 'model does not acknowledge the widening' claim is
contradicted by the stored output. The model-based standard error is
conservative, and more so for OL+BDC and Hybrid than for CO.**
*Location.* Conclusions lines 885-889; Section 5.2 lines 657-663.
*Problem.* The calibration ratio $\kappa$ (mean model-based SE over
empirical SD of $\hat\beta$) is 1.15 to 1.24 for OL+BDC, 1.09 to
1.13 for Hybrid, and 0.95 to 1.06 for CO. The model-based SE rises
by 13 to 14 percent from no dropout to 40 percent biased dropout,
which tracks or exceeds the 4 to 11 percent increase in the
empirical SD. The Conclusions claim the opposite. The Hybrid
empirical SDs quoted in Section 5.2 (0.411 to 0.454) are stale. The
stored values are 0.603 to 0.671. The 'last column' of the summary
file is the MCSE of the mean, not the empirical SE.
*Why it matters.* Design-dependent conservativeness changes power
comparisons across designs. It also reverses the inferential
warning in the Conclusions.
*Remedy.* Report $\kappa$ and 95 percent CI coverage (both planned
in the protocol) for every cell. Add CR2/Satterthwaite inference
(paper 02's G9) alongside the model-based test.
*Evidence.* VERIFIED (script sections 1 and 2).

**M7. The robustness ranking is not supported at 500 replicates.
The cells are unpaired, and absolute power loss is a
baseline-dependent scale.**
*Location.* Abstract lines 84-87 ('the crossover the only design
clearly separated on robustness'); Results lines 617-623;
recommendation table; Conclusions point 1.
*Problem.* Each cell has its own seeds
(`cell_seeds <- 100000L * i + seq_len(n_reps)`), so the no-dropout
and dropout cells of a design use different datasets. The losses
are 0.086 (SE 0.025) for CO, 0.058 (SE 0.015) for Hybrid and 0.056
(SE 0.017) for OL+BDC. The CO-minus-Hybrid difference in losses is
0.028 (SE 0.029, $z = 0.97$). The CO-minus-OL+BDC difference is
0.030 (SE 0.030, $z = 1.01$). Absolute power loss also depends on
where a design sits on its power curve. CO at 0.85 is on a steeper
part than Hybrid at 0.97, so equal information loss produces a
larger absolute loss for CO. On the SE scale the inflation is
nearly identical across designs (ratios 1.128 for CO, 1.138 for
OL+BDC, 1.139 for Hybrid).
*Why it matters.* The study's 'clearest discrimination' is within
Monte Carlo noise and is predicted by baseline power alone. The
sibling paper 07 had a headline that dissolved for the same reason.
*Remedy.* Apply each dropout mechanism to the same simulated
datasets (common random numbers within design), which makes the
loss a paired difference with far smaller MCSE. Compare designs on
a scale that does not depend on baseline power: the ratio of
standard errors, relative efficiency, or the $N$ needed to restore
80 or 90 percent power. Withdraw 'clearly separated' unless the
paired rerun supports it.
*Evidence.* VERIFIED (script section 3; driver seed line
inspected).

**M8. Internal inconsistencies: stale numbers from the unmatched
run remain in the Discussion, the Conclusions and the cover letter.**
*Location and problem.* VERIFIED against the stored summary:

| Location | Text says | Stored value |
|---|---|---|
| Results l. 703-704 | 25 to 40 percent biased costs Hybrid 0.6 points | 3.8 points |
| Results l. 663 | Hybrid empirical SE 0.411 to 0.454 | 0.603 to 0.671 |
| Discussion l. 805 | Hybrid loss about 0.004 | 0.058 |
| Discussion l. 808-809 | Hybrid biased_25 minus MCAR_25 = -0.038 | -0.018 |
| Discussion l. 811-812 | 25 to 40 percent costs 0.006 to 0.056 | 0.030 to 0.056 |
| Discussion l. 793-801 | OL 'sits at nominal alpha'; null check | 0.086 to 0.112; not a null (M4) |
| Conclusions l. 870-871 | 'clear power losses in every design except the hybrid' | Hybrid loss 3.8 (25 to 40) and 5.8 (none to 40) |
| Cover letter | Hybrid 'most dropout-robust', loses under half a point | 5.8 points; not separable from OL+BDC |
| Cover letter | 'most powerful design under no dropout is not the most powerful under dropout' | ordering unchanged in all four conditions |
| Cover letter | consequences 'differ by an order of magnitude' | losses 5.6 to 8.6 points |
| Abstract l. 92 vs Conclusions l. 880 | shift 'less than 0.08' vs 'below 0.04' | max 0.075 over all cells; 0.035 none-to-biased_40 |
| Results l. 549 | '8,000 fits in every cell' | 500 per cell, 8,000 in total |

The sensitivity CIs at lines 691-692 do not match $\pm 1.96$ MCSE
(computed [-2.2, 4.6] and [-3.4, 6.2]); the discrepancy is small.
*Remedy.* Regenerate every number from the summary file with inline
R, as the sibling papers now do, so the text cannot drift from the
data.
*Evidence.* VERIFIED (script sections 1 and 3).

**M9. The motivating claims attributed to Hendrickson et al. (2020)
are inaccurate. Two of them appear to be the protocol's own working
prediction.**
*Location.* Abstract lines 50-52; Introduction lines 128-134 and
239-243; Results line 565; Discussion lines 775-778; Conclusions
lines 860-863.
*Problem.* (a) The Introduction states that Hendrickson's Figure 4A
shows power losses 'exceeding 50 percent' and design rankings that
'reverse depending on the dropout pattern'. Hendrickson's text says
only that 'increased censoring lowered power across all trial
designs'. The two claims match, almost word for word, the 'working
prediction' in this paper's own ADEMP protocol
(`00-ademp-protocol.md`, Study 1). (b) Hendrickson et al. found the
N-of-1 (Hybrid) design had *lower* power than the traditional
crossover (abstract and Results). The present ordering (Hybrid >
OL+BDC > CO) is therefore not a reproduction of their ordering. The
repeated phrase 'reproduce the Hendrickson happy-accident
randomization-path ordering' is incorrect. (c) The 'happy accident'
in Hendrickson et al. refers to trial expense: participants with no
response, or with a response clearly lost on placebo, leave before
the later crossover blocks. It is not a power ordering, and it says
nothing about no-dropout power. (d) The OL+BDC negative bias with
$p < 0.0001$ in their Section 3.4 is a bias under the null
($c_{bm} = 0$) for most censoring patterns. It is not specific to
'symptom-worsening-driven dropout'. In their second, gold-standard
mode, only the OL design showed a censoring effect.
*Why it matters.* The paper's opening premise ('a first-order design
choice' with losses over 50 percent) is not supported by the cited
source. Once (a) is removed, the prototype's modest, design-uniform
losses no longer contradict any published finding.
*Remedy.* Quote Hendrickson et al. accurately, citing figure and
page. Do not present protocol hypotheses as published results.
Discuss why the present ordering differs from theirs (different
correlation parameters: $\rho = 0.7$ and cross-factor 0.1/0.05 here
against 0.8 and 0.2/0.1 there; carryover 0.5 week against 0 to 0.2;
a different analysis half-life; a different dropout mechanism).
*Evidence.* INSPECTED (Hendrickson PDF text,
`docs/references/hendrickson.pdf`; PubMed abstract, PMID 34713026;
figure 4A pixels not inspected, so the absence of a 50 percent loss
in the heatmap itself is UNVERIFIED).

**M10. Most of the stated contributions were not delivered, yet the
abstract and the Introduction describe them as done.**
*Location.* Abstract Methods ('reproduce and extend ... with wider
parameter ranges and sensitivity to dropout-pattern misspecification
at the analysis stage'); Section 1.4; Section 5.4 title.
*Problem.* One parameter point was run. No $N \in \{35, 70, 100\}$,
no $\beta_{bm} \in \{0, 0.3, 0.6\}$, no carryover sweep, no Type I
error cells for the powered designs, no gold-standard bias mode, no
path decomposition. Section 5.4, 'Sensitivity to dropout-pattern
misspecification', contains no misspecification. The same analysis
model is fitted in every cell, and the section compares two
generating mechanisms. Section 3's 'vulnerability windows' (for
example 'two discontinuation windows' for CO, 'largest set' for
Hybrid) are asserted, not tested. The paper's own results do not
bear them out (Hybrid and OL+BDC lose the same).
*Why it matters.* The abstract overstates what the evidence covers.
Statistics in Medicine requires the reported design to match the
executed design.
*Remedy.* Either run Studies 1 to 4 or retitle the work as a
prototype and rewrite the abstract and the contributions to match.
Misspecification at the analysis stage should mean fitting
different analyses (MAR mixed model, pattern-mixture or
delta-adjusted MI, a shared-parameter model) to the same data.
*Evidence.* INSPECTED.

**M11. The missing-data theory is incomplete, and the observed
stability of $\hat\beta$ is not explained by it.**
*Location.* Section 2, lines 351-367.
*Problem.* (a) The time-only component is MCAR, not merely
'ignorable'. (b) The MNAR argument is correct for the driver: the
hazard at visit $t$ depends on $Y_{it}$, and that value is then
deleted, which is a Diggle-Kenward selection on the current value.
But with serial correlation $0.7^{|w_t - w_s|}$ and weekly visits in
the discontinuation phase, $Y_{it}$ is well predicted by
$Y_{i,t-1}$. The mechanism is therefore close to MAR in its effect,
and the likelihood-based mixed model is nearly ignorable. This,
together with the symmetry in M2, is the likely reason for the small
coefficient shifts. It is a property of this mechanism, not of the
designs. (c) Likelihood ignorability under MAR also requires a
correctly specified full-data model. Here the fitted mean is linear
in $t$ with a random intercept and CAR(1), while the generating mean
is a sum of Gompertz curves with three AR(1) components. Under
misspecification, MAR dropout alone can bias $\hat\beta$. (d) The
statement that the violation 'scales with $\beta_1$' was not tested,
since $\beta_1$ was not varied.
*Why it matters.* The conclusion 'informative dropout cost
precision rather than accuracy' generalizes a result produced by a
weakly informative, near-MAR, symmetric mechanism. Shared-parameter
dropout, in which the hazard depends on the participant's latent
drug response (the quantity the biomarker predicts), is the case
the Introduction motivates. It is the case most likely to bias
$\beta_{bm:D}$, and it was not simulated.
*Remedy.* Add the conditions in (b) to (d) to Section 2, and add a
shared-parameter (random-effect-dependent) dropout condition to the
simulation (Wu and Carroll 1988; Follmann and Wu 1995).
*Evidence.* INSPECTED for (a), (b) and (c); (b)'s near-MAR claim
SUSPECTED (not quantified).

**M12. The analysis model and the dropout calibration differ from
project policy and from the text in ways that affect
interpretation.**
*Location.* Simulation design; Reproducibility.
*Problem.* (a) The fitted model is the historical
`Sx ~ bm + t + Dbc + bm:Dbc (+ De)` with model-based SE and an
analysis half-life equal to the generating half-life (0.5 week).
The project's current recommendation is F1 with CR2/Satterthwaite
inference, and paper 02 now recommends G9 (AIC-selected half-life,
CR2). (b) `lme_analysis()` silently refits without `corCAR1` if the
CAR(1) fit fails, and the driver does not store the warning field.
'100 percent convergence' therefore cannot distinguish CAR(1) fits
from fallback fits. (c) The dropout hazard is per visit, not per
unit time. Visit spacing differs by design: weekly in the Hybrid
discontinuation phase (weeks 8 to 12) and the OL+BDC
discontinuation phase (weeks 16 to 20), 2.5-weekly in OL and CO.
'Matched 25 percent dropout' therefore places dropout at different
calendar times and phases in different designs. (d) '25 percent
dropout' is the fraction of participants who ever drop. The
fraction of post-baseline observations lost is about 14 to 16
percent at the 25 percent target and 23 percent at 40 percent.
(e) The 'achieved dropout tracked targets closely' is true by
construction, because $\alpha$ is solved on each realized dataset.
It is not an outcome.
*Remedy.* Add G9 and the policy value $t_{1/2}^a = 0.2$ as analysis
arms. Record the fallback flag. Use a per-unit-time hazard or report
the timing of dropout by phase for each design. Report observation-
level missingness.
*Evidence.* (a), (b) and (e) INSPECTED; (d) VERIFIED on one large
path (script section 8); (c) INSPECTED.

**M13. The reproducibility statement is incorrect.**
*Location.* Reproducibility, lines 928-930.
*Problem.* The text says each driver sets `set.seed(20260507)`, uses
`parallel::mc.set.seed`, and derives seeds as master seed plus
`rep_idx`. The driver sets `set.seed(20260508)`, does not use
`mc.set.seed`, and seeds each replicate with
`100000 * cell_index + rep`. The protocol specifies a cell-descriptor
hash and a `01-base-seed.txt` file, neither of which exists. The
driver header and figure subtitle still say $N = 35$.
`analysis/scripts/informative-dropout-by-design/` contains only the
protocol, not the drivers the text locates there. The README says
the driver directory 'does not yet exist'.
*Why it matters.* A reader following the stated recipe cannot
reproduce the stored output. A sibling paper (22-zztable1) had a
fabricated seed claim flagged for the same reason.
*Remedy.* Describe the seeding that was actually used. Fix the
driver header and subtitle. Move or symlink the driver, or correct
the paths.
*Evidence.* INSPECTED.

**Checks that passed.** The stored summary file reproduces exactly
from the replicate file (all 16 cells, power, MCSE, mean and SD of
$\hat\beta$). VERIFIED. Every design is matched at total $N = 70$
with allocation 70; 35/35; 35/35; 18/18/17/17, consistent with
NOTATION.md rule 4. INSPECTED. The Hybrid design is the genuine
Hendrickson four-path schedule (visits at weeks 4, 8, 9, 10, 11,
12, 16, 20; paths A to D identical to `figure4/01-generate-power-
sweep.R`). OL, OL+BDC and CO also match Hendrickson's
specifications. INSPECTED. Every path covariance matrix used is
positive definite under the current package (minimum eigenvalue
0.168 to 1.0), and `stopOnNonPD = TRUE` accepts all nine. No stored
result depends on a repaired matrix. VERIFIED. The z statistics for
Hybrid versus OL+BDC (1.44, 3.59, 1.55, 0.85) are correct for
independent cells. VERIFIED.

## 5. Minor comments

1. STYLE.md: 'coupling' appears five times in prose (lines 45, 49,
   165, 219, 236), 'ceiling' twice (103, 632). 'Mean-moderation
   architecture' should be 'the mean-shift encoding'.
   'Architecture-B' (line 492) is a retired label. Many semicolons
   join clauses. Several abstract sentences exceed 30 words. Numbers
   are characterized without MCSE ('small', 'modestly'). There is
   no Glossary section. `perl tools/notation-lint.pl` is clean
   (VERIFIED).
2. 'OL+BDC' and 'CO' are used before definition in places. The
   Hybrid design is not defined on first use as STYLE.md requires.
3. Lines 140-143: 'the rate is relatively superficial' is
   contradicted by the paper's own conclusion that rate dominates.
4. Lines 168-182: the claim that 30 percent selective dropout 'can
   degrade power more than' 50 percent uniform dropout is
   unsupported by any cell.
5. Line 413 calls OL 'the robustness benchmark', which conflicts
   with M4.
6. Line 491: under the mean-shift encoding the analysis uses raw
   `bm` ($B_i$), not $b_i$. These coincide only because the
   biomarker is generated with mean 0 and SD 1. Say so.
7. Line 533 and the protocol use 'per design' for $N$, against
   NOTATION.md rule 4. The protocol also uses British spellings
   ('randomisation').
8. Lines 752-754: 'Wall clock was 421 s' is not recorded in any
   artifact. UNVERIFIED.
9. Section 5.3 should be retitled. It reports no decomposition.
10. Line 113: the 2,000-replicate figure is not derived. Give the
    calculation: the target difference, MCSE and power assumed.
11. 'Conflict of interest' and author placeholders ('Temperance
    Persons', '[TODO]' reviewers, 'github.com/[organisation]') are
    not submission-ready. Title-page keywords use 'design-by-dropout
    coupling' and 'Hendrickson hybrid'.
12. The figure `analysis/figures/quick-sim/09-dropout-power.pdf` is
    not included in the manuscript and its subtitle says $N = 35$.
13. `bullets.Rmd` exists in the directory. Per instructions it was
    not reviewed. The `claudecode.tex` header include is inert.
14. The CENT statement's bib year is 2015; the J Clin Epidemiol
    volume 76 issue is 2016.
15. US English: the title page has 'organisation'; the README and
    the protocol have 'programme', 'Formalise' and 'centre'.

## 6. Missing and questionable references

| Reference | Status | Comment |
|---|---|---|
| Hendrickson et al. 2020 | Cited; misattributed | See M9. Ordering, '50 percent' loss, ranking reversal and the 'happy accident' are misdescribed. |
| Stunnenberg et al. 2022; Natesan Batley et al. 2023 | Cited | Bibliographic data match PubMed (Stunnenberg PMID 34675101). That they rank dropout among 'the leading methodological weaknesses' is UNVERIFIED; give the specific finding. |
| Mallinckrodt et al. 2001, 2003; Prakash et al. 2008; Berger et al. 2003; Vohra et al. | Cited | Resolve on PubMed (PMIDs 11459446, 12729388, 18564199, 12689737, 26272792). VERIFIED metadata. |
| Morris, White and Crowther 2019 | Cited | DOI missing (10.1002/sim.8086). |
| ICH E9(R1) | Cited | Add the adoption date and URL. The paper names the estimand framework but never states an estimand with an intercurrent-event strategy for dropout. |
| Therneau and Grambsch 2000; Cox 1972; Kent et al. 2020 (two); Senn 1994; Senn 2002; Raskind et al. 2013 | In bib, not cited | Cite (Senn for crossover, Raskind for the prazosin parameters) or remove. |
| Rubin 1976, Biometrika 63:581-592 | Missing | Foundational MCAR/MAR/MNAR definitions. |
| Matthews and Henderson 2013, Biostatistics 14:626-638, doi:10.1093/biostatistics/kxt009 | Missing | Crossover designs under non-ignorable dropout. The closest antecedent to the paper's question. |
| Little 1993, JASA 88:125-134; Hedeker and Gibbons 1997, Psychol Methods 2:64-78 | Missing | Pattern-mixture models, which the paper claims to adapt. |
| Wu and Carroll 1988, Biometrics 44:175-188; Follmann and Wu 1995, Biometrics 51:151-168 | Missing | Shared-parameter informative dropout, the mechanism most relevant to a biomarker-predicted response (M11). |
| Carpenter, Roger and Kenward 2013, J Biopharm Stat 23:1352-1371; Cro et al. 2020, Stat Med 39:2815-2842, doi:10.1002/sim.8569 | Missing | Reference-based and delta-adjusted sensitivity analysis. |
| Ratitch, O'Kelly and Tosiello 2013, Pharm Stat 12:337-347; Permutt 2016, Stat Med 35:2876-2879, doi:10.1002/sim.6753 | Missing | Tipping-point and delta adjustment; regulatory view of sensitivity analysis. |
| National Research Council 2010; Little et al. 2012, NEJM 367:1355-1360 | Missing | Standard recommendations on missing data in trials. |
| Blackston et al. 2019, Healthcare 7:137, doi:10.3390/healthcare7040137 | Missing | Simulation comparison of aggregated N-of-1 against parallel and crossover designs, including selection bias. |
| Jiang, Arnold and Betensky 2025, Stat Med 44:e70177, doi:10.1002/sim.70177; Selukar et al. 2025, Clin Trials 22:257-266, doi:10.1177/17407745241304284 | Missing | Current N-of-1 design work in the target journal (early stopping is the planned analog of dropout). |
| Senn 2016 or Araujo, Julious and Senn 2016 | Missing | Variance components in series of N-of-1 trials. |

PubMed metadata for the added items were retrieved for Matthews and
Henderson, Cro et al., Permutt, Blackston et al., Jiang et al.,
Selukar et al. and Hendrickson et al. (INSPECTED). The other added
items are cited from the reviewer's knowledge and their pages should
be checked (UNVERIFIED).

## 7. Suggestions for strengthening

1. Define the estimand under ICH E9(R1). Dropout is an intercurrent
   event. State whether the target is the hypothetical
   $\beta_{bm:D}$ had participants remained, and choose analyses
   that target it.
2. Use common random numbers: simulate each complete dataset once
   and apply none, MCAR and each informative mechanism to it. Power
   loss and coefficient shift become paired quantities. The
   2,000-replicate requirement will probably fall substantially.
3. Cross three directional mechanisms (worsening, improvement, and
   symptom return after discontinuation) with a shared-parameter
   mechanism tied to the individual drug response. Index each by an
   interpretable selection strength, such as the odds ratio of
   dropout per SD of the driver.
4. Report bias against the true value, $\kappa$, coverage, power
   and the $N$ needed to restore 90 percent power. Compare designs
   on SE ratio or relative efficiency, not absolute power loss.
5. Add $\beta_{bm} = 0$ cells for every design to check size under
   each dropout mechanism (planned at 5,000 replicates).
6. Fix the carryover mismatch (M5) or treat it as an explicit
   factor, and include G9 inference.
7. Implement the path decomposition. `path` is already in the
   per-replicate data inside `run_one_rep()`; only the dropout
   indicator by path needs to be returned.
8. Add a sensitivity analysis at the analysis stage: MAR mixed
   model, delta-adjusted MI over a grid of $\delta$, and a
   tipping-point display for the interaction.
9. Generate all prose numbers with inline R from the summary file.
10. Either cover the protocol's grid ($N$, $\beta_{bm}$, carryover)
    or narrow the title and abstract to a single-scenario
    prototype.

## 8. Scope of this review

**Read in full.** `report.Rmd` (including the setup chunk; the
manuscript has no computational chunks), `README.md`,
`references.bib`, `cover-letter.md`, `title-page.md`, the project
`CLAUDE.md`, `analysis/report/STYLE.md`, the relevant parts of
`analysis/report/NOTATION.md`, the driver
`analysis/scripts/quick-sim/09-dropout-prototype-10min.R`,
`R/lme_analysis.R`, `R/censordata.R`, `R/generateData.R`
(`generateData`, `buildSigma`), and the Study 1 section of
`analysis/scripts/informative-dropout-by-design/00-ademp-protocol.md`.
I read the Hendrickson et al. (2020) text through `pdftotext`. I did
not inspect the figures.

**Executed.** Reproduction of the summary file from the replicate
file. SE ratios, $\kappa$, and loss differences with standard
errors. A positive-definiteness check of all nine path covariance
matrices under the current package with `stopOnNonPD` both FALSE
and TRUE. A fit of the OL design to identify the formula used.
Single large-$N$ fits ($N = 4000$) for three designs at
$t_{1/2}^a \in \{0, 0.5\}$. One large simulated path to characterize
the direction of the 'biased' hazard and the observation-level
missingness. `perl tools/notation-lint.pl` and the STYLE.md grep.

**Not done.** No replicated simulation was run. The bias attribution
(M5) and the hazard-direction result (M2) rest on single large
samples. The host R session ran without renv, so package versions
may differ from the container. The stored run predates the
2026-10-26 package changes. Because the matrices are PD, the repair
lock does not affect it, but I did not rerun the driver under the
current package. The stored per-replicate data lack the `warning`
and `issingular` fields, so CAR(1) fallback fits could not be
counted. I did not check Figure 4A pixel values or Hendrickson
commit 58b32a9 directly. `bullets.Rmd` was not reviewed. Bibliographic
details for added references not retrieved from PubMed are from the
reviewer's knowledge and are UNVERIFIED.
