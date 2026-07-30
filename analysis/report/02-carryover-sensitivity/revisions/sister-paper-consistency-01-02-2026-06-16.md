# Sister-Publication Consistency Review: Papers 01 and 02
*2026-06-16 16:30 PDT*

**Papers.**

- **01** `analysis/report/01-dgp-mean-moderation-vs-mvn/` -- Two
  architectures for simulating biomarker-treatment interactions
  (mean moderation vs MVN differential correlation).
- **02** `analysis/report/02-carryover-sensitivity/` -- Robustness of
  carryover-mitigation analysis strategies for biomarker-treatment
  interaction.

**Brief.** The two are intended as sister publications (both to
Statistics in Medicine, both authored 'Temperance Persons, on behalf of
the pmsimstats team', paper 02 repeatedly cites 01 as 'the companion
manuscript'). This review checks them for mutual consistency:
shared definitions, shared numbers, shared anchors, cross-references, and
overall coherence. All numerical comparisons were recomputed from the
archived summaries (`02-grid-summary.rds`) and read directly from paper
01's results tables (`report.Rmd` lines 850-944) and abstract.

**Verdict.** The two papers' *underlying data are consistent with each
other*. The problem is that paper 02's *prose description of paper 01* is
not: paper 02 attributes to the companion a set of carryover power-loss
figures (40-60 percentage points for Architecture B, 8-13 for
Architecture A, and '0.82 at t1/2 = 1.0') that paper 01 does not report
and, in two of three cases, explicitly disclaims. This is the dominant
finding and must be fixed before the pair is submitted together. Several
secondary mismatches (different reference cells, an unacknowledged third
architecture, a missed mechanistic cross-link) should also be resolved.

---

## 1. The carryover power-loss numbers do not agree across the pair

### 1.1 What paper 02 tells the reader the companion found

Paper 02, Introduction (orig block, lines ~241-249) and Abstract bullet
(lines ~218-221):

> 'Under Architecture B, carryover at a half-life of one week reduces
> power by **40 to 60 percentage points** from the no-carryover baseline;
> under Architecture A the corresponding loss is **8 to 13 percentage
> points**.'

and (Matched-decay baseline, lines ~1203-1224):

> 'consistent with the Architecture B numbers from the companion
> manuscript (**0.82 at t1/2 = 1.0**)'.

### 1.2 What paper 01 actually reports

Paper 01's production factorial (N = 35, 1,000 replicates, exponential
decay) gives these power values (paper 01 lines 850-902):

| Architecture | Design | t1/2=0 | t1/2=0.5 | t1/2=1.0 | Relative loss |
|---|---|---:|---:|---:|---:|
| A | CO | 0.739 | 0.739 | 0.725 | 1.9% |
| A | Hybrid | 0.992 | 0.987 | 0.978 | 1.4% |
| A | OL+BDC | 0.751 | 0.748 | 0.730 | 2.8% |
| B | CO | 0.745 | 0.754 | 0.723 | 3.0% |
| B | Hybrid | 0.995 | 0.980 | 0.944 | 5.1% |
| B | OL+BDC | 0.777 | 0.694 | 0.539 | 30.6% |

Paper 01's largest Architecture-B loss is **30.6%** relative (about 24
percentage points) at OL+BDC; its Architecture-A losses are **1.4-2.8%**
(about 2 percentage points). Paper 01 then states explicitly, in its own
prose, that these fall *below* the ranges paper 02 quotes:

- Line 882: Architecture A 'fall[s] well below the 8-13% relative-loss
  range that the broader N-of-1 simulation literature reports'.
- Line 935: Architecture B 'falls below the 40-60% relative-loss range
  that the broader N-of-1 simulation literature reports'.

So in paper 01 the '40-60%' and '8-13%' figures are *literature ranges
that paper 01's results undercut*. In paper 02 the same figures are
presented as *the companion's findings*. That is a direct misattribution.

### 1.3 Paper 02's own data also contradict its introduction

At paper 02's own reference cell (N = 70, Hybrid, A2 matched analysis,
recomputed from `02-grid-summary.rds`):

| Architecture | t1/2=0 | t1/2=1.0 | Absolute loss | Relative loss |
|---|---:|---:|---:|---:|
| B (MVN) | 0.988 | 0.860 | 12.8 pp | 13.0% |
| A (mean moderation) | 0.976 | 0.960 | 1.6 pp | 1.6% |

Paper 02's own Architecture-B loss is 13% (not 40-60), and its
Architecture-A loss is 1.6% (not 8-13). The introduction is therefore
inconsistent not only with the companion but with paper 02's own results
section.

### 1.4 Three distinct errors are bundled here

1. **Misattribution.** The 40-60% / 8-13% ranges belong (in paper 01) to
   an uncited 'broader literature', not to the companion's findings.
   Paper 02 should say the companion *found* losses of roughly 1-3%
   (Architecture A) and up to ~31% (Architecture B, design-dependent),
   and that both fall below the older literature ranges.

2. **Unit conflation.** Paper 01 reports *relative* loss in percent
   (30.6%); paper 02 restates it as *percentage points* ('40 to 60
   percentage points'). Relative loss and absolute percentage-point loss
   are different quantities and must not be interchanged.

3. **The '0.82' provenance.** Paper 02 attributes '0.82 at t1/2 = 1.0' to
   the companion's Architecture B. Paper 01's Architecture-B values at
   t1/2 = 1.0 are 0.723 (CO), 0.944 (Hybrid), 0.539 (OL+BDC) -- none is
   0.82. The 0.82 is Hendrickson's (2020) original figure, which paper 02
   elsewhere correctly attributes to Hendrickson. Remove the companion
   attribution.

**Fix.** Rewrite the three passages in paper 02 (Abstract bullet,
Introduction orig block, Matched-decay baseline) to quote paper 01's
actual production numbers. Better, since the two papers share a code base,
quote a single agreed cell computed once and reused in both manuscripts.

---

## 2. The two papers do not share a reference cell

| Anchor | Paper 01 | Paper 02 |
|---|---|---|
| Reference design | OL+BDC | Hybrid |
| Reference N | 35 | 70 |
| Reference c_bm | 0.45 | 0.45 |
| Replicates / cell | 1,000 | 500 |
| Decay forms | exponential only | linear, exponential, Weibull |
| Architectures in main factorial | A, B (+ C discussed) | A, B |

Sister papers that cross-cite each other's headline numbers should either
share a reference cell or state plainly why they differ. The current
mismatch is the proximate cause of the Section 1 confusion: paper 02's
author appears to have reached for a remembered loss figure rather than
the companion's actual OL+BDC/N=35 result. A one-paragraph 'shared
calibration and reference cell' note, identical in both papers, would
prevent recurrence. Note also the cross-paper power inversion this
mismatch can produce: paper 01 reports Hybrid/Architecture B at N = 35 as
0.944 (t1/2 = 1.0), while paper 02 reports Hybrid/Architecture B at N = 70
as 0.860; paper 01 already explains (lines 1141-1149) that N = 70 is
compressed toward the ceiling at this calibration, but the two numbers
read as contradictory unless that explanation is surfaced in paper 02 as
well.

---

## 3. A mechanistic cross-link is missed, and it matters

This is the most valuable opportunity in the pair, and it connects to the
separate finding (in `readability-review/...`) that paper 02's claim of
universal A2 dominance is false.

- **Paper 01** finds the crossover (CO) design shows *no* B-versus-A
  carryover gap, because 'the CO design's within-subject AR(1) structure
  dominates the correlation-decay channel that Architecture B depends on'
  (lines 914-916; z = 0.29, not significant).

- **Paper 02** finds (recomputed) that A2 -- the exposure-weighted
  specification built specifically to exploit the carryover-correlation
  channel -- has *no advantage and is in fact dominated* under the CO
  design (A2 = 0.488 vs A1 = 0.830 at the CO/B/t1/2=1.0 cell; A2 wins in
  0% of CO cells).

These are the same mechanism seen from two sides: when the design
neutralises the carryover-correlation channel (CO), there is nothing for
the exposure-weighted predictor to recover, so A2 collapses to (or below)
the binary specification. Paper 01 supplies the explanation for paper
02's most important caveat. Yet paper 02 does not cite it, and instead
asserts A2 dominance 'across all three trial designs', contradicting both
its own data and the companion's mechanism.

**Fix.** Paper 02 should (a) correct the dominance claim (see the
readability review), and (b) cite paper 01's CO result as the mechanistic
reason A2's advantage vanishes under the crossover design. This turns a
contradiction into a mutually reinforcing pair.

---

## 4. Scope mismatch: Architecture C

Paper 01 now formalises a third architecture, Architecture C (combined
channels, parameters c_bm,A and c_bm,B), with its own derivation (Section
2.3), a combined-DGP factorial, and an Architecture-C driver (lines
~480-700, 1085-1265, 2863). Paper 02 evaluates only Architectures A and B
and states 'We retain both DGP architectures from the companion
manuscript' (Section 3.1). 'Both' is now stale: the companion has three.

This is not a numerical error, but for a coordinated pair it reads as the
two papers having drifted. Paper 02 should say it retains 'the two
single-channel architectures (A and B) from the companion' and note that
Architecture C is out of scope here, so a reader who has seen paper 01
does not wonder why C is absent.

---

## 5. Notation and definitional consistency

Largely consistent, with harmonisable drift:

| Element | Paper 01 | Paper 02 | Action |
|---|---|---|---|
| Architecture A interaction symbol | c_bm,A (and z_B, sigma_BR) | c_bm / beta_bm | Harmonise to one symbol set |
| Architecture B interaction symbol | c_bm,B | c_bm | Disambiguate; paper 02's bare c_bm is paper 01's c_bm,B |
| Decay parameterisation | exponential, lambda = ln2/t1/2 | same, plus linear and Weibull | Consistent; paper 02 is a strict extension |
| Analysis model | lme + corCAR1 AR(1), A2 (Dbc) spec | same family; A1/A2/A3 | Consistent; A2 is the shared spec |
| Calibration target | prazosin/PTSD, sigma_BR | sigma_BR = 8, theta_true = -3.6 | Confirm paper 01 uses the identical sigma_BR = 8 and cite docs/19 in both |
| Designs | CO, Hybrid, OL+BDC (OL omitted) | CO, Hybrid, OL+BDC (OL omitted) | Consistent |

The most consequential item is the c_bm symbol: in paper 01 the bare
interaction strength is split into c_bm,A and c_bm,B once Architecture C
is introduced, whereas paper 02 uses an unsubscripted c_bm that
corresponds to paper 01's c_bm,B under Architecture B and to the mean
coefficient under Architecture A. A reader moving between the papers will
trip on this. Fix the symbols centrally.

---

## 6. Shared state and shared defects

Both papers are in the same unfinished drafting state and share the same
defects, which is at least internally consistent but blocks joint
submission:

- **Both** render the three-part bullets/rgt/orig scaffold with every
  'rgt' block containing Lorem ipsum placeholder text. Neither is
  submittable until `strip-claudecode` is run on both.
- **Both** lean on the unpublished calibration companion (paper 10,
  `pmsimstats_calibration2026`) for Type I interpretation.
- **Both** attribute the 40-60% / 8-13% ranges to a 'broader N-of-1
  simulation literature' with no citation. If that range originates in
  the project's own archived exploratory runs (figure5 /
  architecture_comparison) rather than external literature, both papers
  are mis-citing an internal legacy number as published literature. Track
  down the source or remove the comparison from both.

Consistent and correct across the pair: authorship, corresponding author,
target venue, the AR(1)/corCAR1 analysis model, the prazosin-PTSD
calibration motivation, the design set, and the Type-I-conservative
finding (paper 01: [0.010, 0.074]; paper 02: mean 0.039, max 0.084).

---

## 7. Prioritised synchronisation actions

1. **Correct paper 02's description of the companion's results**
   (Section 1): replace '40-60 pp / 8-13 pp' with paper 01's actual
   production numbers (Architecture B up to ~31% relative, design
   dependent; Architecture A ~1-3%), fix the percent-vs-pp unit error, and
   remove the '0.82 = companion' attribution. This is the only item that
   rises to a factual error.
2. **Add the CO mechanism cross-link** (Section 3): cite paper 01's CO
   null result as the reason A2's advantage vanishes under the crossover
   design, and correct paper 02's universal-dominance claim in the same
   stroke.
3. **State the reference-cell and sample-size relationship** (Section 2)
   identically in both papers, including the N=70 ceiling-compression
   note, so the cross-cited numbers are reconcilable.
4. **Acknowledge Architecture C** as out of scope in paper 02 (Section 4).
5. **Harmonise the c_bm,A / c_bm,B / c_bm notation** across both
   (Section 5).
6. **Strip the scaffold and resolve the uncited '40-60% literature'**
   range in both papers (Section 6).

Items 1-3 are substantive and change what the papers claim. Items 4-6 are
editorial coherence. None requires re-running simulations; the archived
results of both papers are mutually consistent once the prose is
corrected to match them.

---
*Rendered on 2026-06-16 at 16:30 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/02-carryover-sensitivity/sister-paper-consistency-01-02-2026-06-16.md*
