# Phase 1 conflict-resolution log (docs/59 CR-31, CR-43, CR-44, CR-51, GP-29)
*2026-10-11 22:35 PDT*

Author: pmsimstats team

REPO = `/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng`; SP = this
scratchpad. Evidence labels: verified (run and observed), inspected
(source read), inferred, unverified. Nothing was committed.

## CR-31: anchored-change figures in docs/37 and docs/39

**Method (verified).** `docs/series/notes/phase1/cr31-compare.R` and `docs/series/notes/phase1/cr31-mc.R` merge
the pre-anchored outputs (`SP/pre-anchored-outputs/analysis/data/
quick-sim/carryover-closed-form/*.csv`) with the current anchored
outputs (`REPO/analysis/data/quick-sim/carryover-closed-form/*.csv`)
on design, construct, rule, half-life, variant and coupling, with a
`stopifnot()` that every row matches. `SP/analyze-01c.R` was rerun:
it pairs the anchored Monte Carlo with the recursive sandbox rerun
(`SP/sandbox-recursive/`, which reproduces the stored recursive CSVs
exactly; `SP/anchored-rerun-comparison.md` Section 0) on the same
per-cell seeds. `01-e9-closed-form.R` itself was not rerun: the
stored pair already isolates the carryover form.

**Closed-form power, largest absolute change, anchored against
recursive.**

| File (grid) | CO | Hybrid | OL+BDC |
|---|---|---|---|
| `e9-closed-form-curves-designs.csv` (0 to 4 weeks) | 0.0409 (CS, mean, E9, $t_{1/2} = 4$: 0.5052 to 0.4644) | 0.0016 (CS, mean, E9b, 4) | 0.0020 |
| same, restricted to $t_{1/2} \le 2$ | 0.0169 | 0.0010 | 0.0018 |
| `e9-closed-form-vs-mc.csv` (0 to 2 weeks) | 0.0171 (CS, mean, E9, 2: 0.5478 to 0.5307) | 0.0010 | 0.0018 |
| `design-power.csv` (0, 0.5, 1, 2), docs/39 | 0.0171 | 0.0001 | 0.0005 |

Other designs in `design-power.csv`: CO with open-label first period
0.0083, Hybrid timing with ABAB blocks 0.0066, ABAB 0.0041, others
0.0002 or less.

**Monte Carlo power.** Paired on the same seeds (carryover only):
CO 0.0178, Hybrid 0.0032, OL+BDC 0.0030; every $t_{1/2} = 0$ cell
identical. Stored values (old seeds against new): largest at
$t_{1/2} > 0$ is 0.0216 (Hybrid, CS, graded, E9, $t_{1/2} = 0.2$: 0.7264
to 0.7048; reseeding, its paired change is 0.0000); largest overall
0.0230 (Hybrid, CS, step, E9b, $t_{1/2} = 0$; reseeding only).

**Answer.** The 0.017 (header) and 0.041 (Appendix A) are both correct
closed-form figures for different grids: 0 to 2 weeks (the tabulated
and Monte Carlo grid) and 0 to 4 weeks (the plotted curves). The
Appendix A Monte Carlo 0.018 is the paired carryover-only change and
correct. The header's "at most 0.021" was the stored Monte Carlo change
of a reseeded Hybrid cell (0.0216, which also rounds to 0.022), not a
carryover effect, so it was misleading. The docs/39 0.017 is correct.

**Edits.** docs/37 header note rewritten with grids stated, plus a
dated "Corrected 2026-10-11" note; docs/37 Appendix A bullet gains the
grid qualifier ("checked 2026-10-11"); docs/39 header: "Values changed
by at most 0: no value" replaced by "No value" (the phrasing flagged
in docs/59 Section 7.2) and a dated "Checked 2026-10-11" note added.

## CR-43: `characterize_carryover()`

**Inspected** (`REPO/R/carryover_analysis.R`, `REPO/R/lme_analysis.R`):

- Lines 27 to 29: default grid `c(0.1, 0.25, 0.5, 1, 2, 4, 8)` weeks.
  No 0. Legacy 01 (0.1 to 8) is right; legacy 03 Section 6.2 (0.25 to
  8) is wrong. The `CLAUDE.md` example passes its own grid
  `c(0.25, 0.5, 1, 2, 4)`.
- Lines 41 to 114: for each half-life, `lme_analysis()` with
  `full_model_out = TRUE`, the interaction `bm:Dbc` in the model; the
  model-based (not CR2) `tTable` gives `drug_p` (line 76) and
  `interaction_p` (line 90); AIC, BIC, logLik recorded (lines 102 to
  104).
- Lines 116 to 118: `best_t_half` is the AIC minimum; lines 121 to 124
  refit at that value (not wrapped in `tryCatch`; if every fit fails,
  `which.min()` returns `integer(0)` and the refit errors).
- `print_carryover_summary()` lines 318 to 325 prints "Drug effect at
  best" and "Interaction at best" with their p-values as if the
  half-life were fixed.
- `lme_analysis()` lines 167 to 177: `nlme::lme` with its default REML;
  on error, refit without `corCAR1`, silently.

**Answer.** It is the F1A-type procedure (AIC choice with the
interaction in the model, then a test as if fixed), which docs/57
Sections 1 and 3.4 show is liberal (pooled 0.057 to 0.074, 0.073 at
true 1 week; derived reason: AIC_1(h) = AIC_0(h) - LR(h)). It is in
three ways weaker than docs/57's F1A: (1) REML likelihoods are compared
across half-lives although $D_{bc,it}$, and hence the fixed-effect
design, changes with the half-life, so the AIC comparison is not valid
(REML includes the $\log\lvert X^\top V^{-1} X\rvert$ term; standard
result, inspected, not run); (2) model-based instead of CR2 standard
errors; (3) the silent `corCAR1` fallback can make the AIC compare
different covariance models. Not used by any analysis driver
(`simulation-core.R` lines 557 to 562 only mention it); exported.

**Proposed code change (not implemented).**

1. Select the half-life by ML (`method = 'ML'`) and without the
   interaction term (the F1A0 rule of docs/57 Section 4.2), then refit
   by REML with the interaction at the chosen value.
2. Do not report a nominal interaction p-value at the selected
   half-life. Return the interaction test at a fixed, user-stated
   half-life (default 0.2 week, the project policy, DQ decided), and
   offer the calibrated p-value of F1A0B (parametric bootstrap,
   docs/57 Section 4.7) as an option with `n_boot`; label the
   AIC-selected p-value, if kept at all, as descriptive.
3. Use CR2 standard errors (as in the stress-test driver).
4. Record per half-life whether `corCAR1` was fitted (`cor_struct`
   column) and exclude fallback fits from the AIC comparison.
5. Guard the all-failed case before the refit; add 0 to the grid and
   align it with the grid chosen under DQ-10.
6. Update the roxygen text and `print_carryover_summary()` wording.

## CR-44: published dropout form at 58b32a9

**Inspected.** `gh api "repos/rchendrickson/pmsimstats/contents/
R/censordata.R?ref=58b32a9"` (stdout only) is identical to
`REPO/analysis/data/quick-sim/hendrickson-58b32a9/R/censordata.R`.
Lines 53 to 67:

```r
fna1<-frac_NA*(1-frac_NA_biased)        # beta0 * (1 - beta1)
fna2<-frac_NA*frac_NA_biased            # beta0 * beta1
cdt_ps2<-sapply(cdt,function(x) ((x+100)^(censorparam$eb1)))
cdt_p1<-t(t(cdt_ps1)*d$t_wk)
cdt_p2<-t(t(cdt_ps2)*d$t_wk)
cdt_r1<-runif(...,min=0,max=2*mean(cdt_p1)*(.5/fna1))
cdt_r2<-runif(...,min=0,max=2*mean(cdt_p2)*(.5/fna2))
do<-as.data.table((cdt_p1>cdt_r1)|(cdt_p2>cdt_r2))
```

Lines 79 to 86 carry each flag forward to every later visit. `cdt` are
the `D_` columns, defined in `generateData.R` line 174 at 58b32a9 as
the per-visit sum of the three component means `tv`, `pb`, `br`
(line 62); the symptom score is `BL` minus that sum (line 181), and the
Gompertz `max` parameters used by the project's Figure 4 drivers
(`analysis/scripts/figure4/02-*.R` line 120) are positive (6.5, 6.5,
11.0), so `D` is the improvement from baseline.

**Published form.** For participant $i$ at visit $k$ with interval
$\Delta t_k$ (`t_wk`) and improvement $D_{ik} = TV_{ik} + PB_{ik} +
BR_{ik}$:

$$P(F^{(1)}_{ik}) = \min\{1,\ \beta_0(1-\beta_1)\,\Delta t_k /
\overline{\Delta t}\},\qquad P(F^{(2)}_{ik}) = \min\{1,\
\beta_0\beta_1\,p_{ik}/\bar p\},\quad p_{ik} = \Delta t_k (D_{ik} +
100)^{e_{b1}},$$

independent flags, means over all participant-visits of the simulated
trial; visit $k$ and all later visits are censored if either flag is
set at any visit up to $k$.

**Verified on synthetic data** (`docs/series/notes/phase1/cr44-censor-check.R`, 20,000
participants, 8 visits of 2 weeks, half with $D = 0$ and half rising to
30): biased-only censoring ($\beta_1 \approx 1$) at the last visit 0.293
for non-improvers against 0.382 for improvers; first-visit censoring
0.043 to 0.044 at $\beta_0 = 0.05$; cumulative censoring at visit 8 is
0.335 under the "balanced" set; censoring monotone for every
participant.

**Resolution.** Legacy 13 Section 4 has the right term,
$(\text{change} + 100)^{e_{b1}}$, but the wrong direction ("higher for
participants whose symptoms have not improved", which is also what the
roxygen text says). Legacy 17 Section 7.1 item 4 has the direction
right ("responders more likely to drop out") but the form wrong
($\beta_0 + \beta_1(\text{change})^2$: there is no additive $\beta_0$,
the change enters shifted by 100 and multiplied by the interval, and
$\beta_0$ splits into the two shares). $\beta_0$ is a per-visit flag
rate, not the total dropout. The dependence on improvement is weak at
these scales: an improvement of 30 points raises the biased term by
$(130/100)^2 \approx 1.7$ relative to no change (inferred). Whether 17's
"during placebo" qualifier is meant is not supported by the code: the
term has no phase dependence beyond $D_{ik}$ (inspected).

Correction for Appendix I Section 2 (legacy 13, 17 not edited): state
the equation above with the direction.

## CR-51: failed fits and selection shares

**Verified** (`docs/series/notes/phase1/cr51-check.R`, `docs/series/notes/phase1/cr51-check2.R`) on
`REPO/analysis/data/quick-sim/carryover-closed-form/strawman-stress/
replicates-full.rds` (308,000 rows: 41 cells, 28,000 trials, 11
analyses A, C1, C2, DA, DS, EA, ES, F1, F2, R, S), `-extra.rds` (G1,
G2 plus S, 84,000 rows) and `-extra2.rds` (F1A plus S, 56,000 rows); S
identical across the three files. A failed fit is a missing estimate,
SE or p-value.

| Statement | Stored value | Verdict |
|---|---|---|
| 48 Purpose: 27 failed fits in 291,500 (39 cells, 26,500 trials, 11 analyses) | cells 1 to 39: 27 of 291,500 | correct for the 39-cell run; cells 40, 41 (`prog-BL`, added later) add 1 (F1) |
| 47 Section 9: G1, G2 12; S, C1, C2 6; F1, F1A, ES 3; DS, F2 2; A, DA, EA, R 0; over 28,000 trials | identical, each per 28,000 | correct; total 55 of 392,000 over 14 analyses, 28 of 308,000 over the 11 of the main run |
| 47 Section 9: at most 0.6% in any cell (G1, G2) | 0.006, cell 1, G1 and G2 | correct |
| 50 Section 2: DS picks 0 82% when none, 1 or 2 weeks 86% when 1 | null cells 19, 20 (1,000 each): 0.816, 0.857 | correct (DS = S with AIC) |
| 56 Section 6.1, 57 Sections 1, 3.5: F1A 87% at 0, 77% at 1 or 2 | 0.870, 0.767 | correct (F1A = F1 with AIC) |
| 57 Table 3, compound-symmetry null row: 1 or 2 weeks 0.846 | 0.382 + 0.464 = 0.846 | correct |

The docs/59 register's "47 Section 9 (sum 28)" was itself imprecise:
47's counts sum to 55 over 14 analyses; 28 is the 11-analysis main run.
No edit was made to docs/47, 48, 56 or 57 (nothing wrong); legacy 50
not edited (nothing wrong). The alternative cells (500 trials each)
give DS 0.778 / 0.856 and F1A 0.844 / 0.768 for the same quantities.

## GP-29: four suspect references

PubMed (DOI to PMID, then metadata) and Crossref (`api.crossref.org/
works/<doi>`, no mailto parameter, no email sent). Searches: PubMed
author searches for "Konings P AND Eng KH" (0 hits), "Pequignat"
(1 hit, a 1975 ascidian enzyme paper), "De Bock M AND Van den
Noortgate W" (0 hits); Crossref bibliographic search of the Chen, Li
and Yang title.

| Cited as (legacy 20, 17) | DOI resolves to (PubMed, Crossref agree) | Verdict |
|---|---|---|
| Chen Y, Li P, Yang K. Comparison of aggregated N-of-1 trials with parallel and crossover RCTs using simulation studies. *PLoS One* 2020;15(1):e0216949 | PMID 31100082: Ao Z, Wang L, Mendoza EJ, et al. Incorporation of Ebola glycoprotein into HIV particles ... *PLoS One* 2019;14(5):e0216949 | misattributed; the title is Blackston J, Chapple A, McGree J, McDonald S, Nikles J (2019; initials as given by Crossref), *Healthcare* 7(4):137, doi:10.3390/healthcare7040137 (Crossref). Not found as cited |
| Eng KH, Konings P (2018). A simulation study on estimating biomarker-treatment interaction effects in randomized trials with prognostic variables. *Trials* 19:72, doi:10.1186/s13063-018-2491-0 | PMID 29463271: Haller B, Ulm K, same title, *Trials* 2018;19(1):128 | misattributed (authors and article number); correct citation Haller and Ulm (2018), *Trials* 19:128 |
| Pequignat F, Balaam B, Sherrington R, Sherrington C, Maher C. Power analysis for idiographic (within-subject) clinical trials ... *Behav Res Methods* 2023;55:3824 to 3842, doi:10.3758/s13428-022-02012-1 | PMID 36526885: Tueller S, Ramirez D, Cance JD, Ye A, Wheeler AC, Fan Z, Hornik C, Ridenour TA, same title, *Behav Res Methods* 2023;55(8):4175 to 4199 (online 2022-12-16) | misattributed (authors, pages); `analysis/report/references.bib` key `pequignatIdiographicPower2023` already carries the corrected authors |
| De Bock M, Van den Noortgate W, Gielen M. Accurate models vs. accurate estimates: a simulation study of Bayesian single-case experimental designs. *Behav Res Methods* 2021;53:1782 to 1798, doi:10.3758/s13428-020-01522-0 | PMID 33575987: Natesan Batley P, Hedges LV, same title, volume and pages | misattributed (authors); correct citation Natesan Batley and Hedges (2021) |

The annotations in legacy 20 for these entries describe content
attributed to invented authors (for example, Chen et al. "replicate
and extend the Blackston et al. framework", while the title is
Blackston's own paper); they are not to be carried into Appendix G
without re-reading the real papers (unverified whether the annotation
content matches them).

## Files changed

- `REPO/docs/37-carryover-power-closed-form.md`: header note (lines
  15 to 35 after edit) and Appendix A final bullet.
- `REPO/docs/39-nof1-design-optimization.md`: header note.
- `REPO/docs/59-whitepaper-refactoring-plan.md`: Section 5 dated note;
  rows CR-31, CR-43, CR-44, CR-51; Section 8 dated note; rows GP-29,
  GP-34, GP-35, GP-36.
- Scratch scripts: `docs/series/notes/phase1/cr31-compare.R`, `docs/series/notes/phase1/cr31-mc.R`,
  `docs/series/notes/phase1/cr44-censor-check.R`, `docs/series/notes/phase1/cr51-peek.R`, `docs/series/notes/phase1/cr51-check.R`,
  `docs/series/notes/phase1/cr51-check2.R`.

Not edited: docs/45, docs/47, docs/48, docs/56, docs/57, anything in
`docs/legacy/`, any code.
