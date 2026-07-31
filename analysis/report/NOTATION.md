# Canonical notation for the pmsimstats-ng manuscript compendium
*2026-07-30 11:57 PDT*

This file is the single source of truth for mathematical notation and
controlled vocabulary across the manuscripts under `analysis/report/`.
New manuscripts should adopt it rather than re-deriving symbols, and
any deviation should be stated explicitly in the deviating paper's
Methods section.

The audit that motivated this file is
`analysis/report/whitepaper-notation-audit.md`; the repairs applied on
2026-07-30 are recorded at the end of this document.

## Core symbols

| Object | Symbol | Notes |
|---|---|---|
| Participant index | $i$ | |
| Timepoint index | $t$ | weekly occasions in all papers |
| Observed outcome | $Y_{it}$ | code column `Sx` |
| Baseline level | $\mathrm{BL}_i$ | participant-specific |
| Biological-response component | $BR_{it}$ | pharmacological |
| Placebo-belief component | $PB_{it}$ | expectation-driven |
| Natural-history component | $TV_{it}$ | time-variant drift |
| Binary drug state | $D_{it}$ | 1 on drug, 0 off drug |
| Continuous exposure | $D_{bc,it}$ | 1 on drug, $e^{-\lambda t_{sd}}$ off drug; code column `Dbc` |
| Time since discontinuation | $t_{sd}$ | code `tsd` |
| Time on drug | $t_{od}$ | cumulative; Gompertz input; code `tod` |
| Biomarker | $B_i$ | participant-level, pre-treatment |
| Standardized biomarker | $b_i$ | $(B_i - \bar{B})/s_B$; code `bm_z` |
| Moderation, covariance channel | $c_{bm}$ | Architecture B; a **correlation** |
| Moderation, mean channel | $\beta_{bm}$ | Architecture A; a dimensionless **multiplier** |
| Combined-channel weights | $c_{bm,A}$, $c_{bm,B}$ | Architecture C only |
| Decay rate | $\lambda$ | $\ln 2 / t_{1/2}$; code `lambda_cor` |
| Carryover half-life | $t_{1/2}$ | code `carryover_t1half` |
| AR(1) correlation | $\rho$ | within-factor serial correlation |
| Target estimand | $\beta_{bm:D}$ | the `bm:Dbc` interaction coefficient |
| Generic interaction coefficient | $\theta$ | permitted only where the argument is deliberately generic (paper 10) |
| Sample size | $N$ | participants; **state per-path or total** |

### Rules that the symbols encode

1. **$D_{it}$ is binary; $D_{bc,it}$ is continuous.** These are
   different objects and the distinction is the subject of the
   carryover-specification manuscript. Never use one for the other.
2. **$c_{bm}$ is a correlation; $\beta_{bm}$ is a multiplier.** They
   are calibrated to a common numeric scale, because the
   Architecture-A shift is applied as $\beta_{bm}\,b_i\,\sigma_{BR}$ on
   the standardized biomarker, so the reference value $0.45$ denotes a
   matched moderation strength under either architecture. The symbols
   are nonetheless not interchangeable. A paper reporting both
   architectures at a matched value may use one label for the shared
   quantity, but must say so.
3. **Mathematics in model definitions, code identifiers in code
   contexts.** Write $Y_{ij} = \beta_0 + \beta_D D_{bc,ij} + \ldots$
   for the model; write `Sx ~ bm + t + Dbc + bm:Dbc` in `\texttt{}`
   when quoting the R formula. Do not mix the two inside one
   expression.
4. **State the sample-size convention.** $N$ per randomization path
   and $N$ in total are both in use across the series and are not
   interchangeable; every Methods section must say which it means.
5. **State the sign convention.** The three components are
   non-negative reductions in symptom severity, so an increase in any
   component lowers $Y$; treatment effects and interaction
   coefficients are consequently negative, and moderation parameters
   positive.

## Units

Carryover half-life $t_{1/2}$ is quoted in **weeks** in the
interaction-focused manuscripts (01, 02, 03, 07, 11), on the canonical
grid $\{0, 0.5, 1.0\}$, and in **days** in the main-effect manuscripts
(04, 05), on a pharmacokinetic scale with a 3-day baseline. Both
conventions are legitimate and are retained because each paper's
figures are drawn on its own scale. Every paper using days must carry
an explicit units paragraph giving the conversion (1 week = 7 days);
papers 04 and 05 now do. Timepoints $t$ are weekly in both.

## Controlled vocabulary

**Design families.** OL (open-label), OL+BDC (open-label plus blinded
discontinuation), CO (traditional two-period crossover), Hybrid (the
Hendrickson 2020 design: open-label titration, blinded
discontinuation, brief crossover). Define the abbreviation at first
use in each paper; capitalize Hybrid.

**DGP architectures.** Named by the channel that carries the
biomarker-treatment interaction. Single letters are no longer used.

| Name in prose | Short form (legends, tables) | Data value | Former label |
|---|---|---|---|
| the mean-moderation architecture | Mean | `mean_moderation` | Architecture A |
| the covariance-moderation architecture | Covariance | `mvn` | Architecture B |
| the dual-channel architecture | Dual | `combined` | Architecture C |

The mean-moderation architecture places the interaction in the
conditional mean (the biomarker scales the size of the treatment
effect); the covariance-moderation architecture places it in the
treatment-state-dependent covariance (the biomarker is a proxy for a
hidden responder state); the dual-channel architecture activates both
with independent weights $c_{bm,A}$ and $c_{bm,B}$, and recovers the
other two as its single-channel limits.

**Analysis specifications** (carryover-specification manuscript).
E1 (binary on-drug indicator), E2 (exposure-weighted continuous
predictor), E3 (binary plus lagged nuisance term). The letter E stands
for the exposure regressor that distinguishes them. These were labeled
A1/A2/A3 until 2026-07-30 and are still stored under those values in
the archived simulation output; the mapping is applied at display time
only.

**Why not A/B/C and M1/M2/M3.** The original scheme used A/B/C for
architectures and A1/A2/A3 for specifications, which collide on the
letter A within a single manuscript. Relabeling the specifications
M1/M2/M3 removed that collision but created a second one: the
calibration manuscript already uses M0-M3 for its working-covariance
model ladder, and those labels are keys into archived data. E1/E2/E3
is free across the whole corpus, and naming the architectures by
mechanism rather than by letter removes the class of collision
entirely.

**Simulation studies.** Number them (Study 1, Study 2, ...) rather
than lettering them, to keep letters free for architectures.

**Sensitivity sweeps.** Sweep labels of the form S1, S2, ... are
paper-local and are not comparable across papers. Cite them with the
paper prefix when crossing manuscripts (for example 05-S2, not S2).

**Reporting.** ADEMP follows Morris, White, and Crowther (2019);
Monte Carlo standard error is abbreviated MCSE throughout.

## Repairs applied on 2026-07-30

- Paper 03: all 14 uses of $D_{it}$ (which denoted the continuous
  decayed indicator) renamed to $D_{bc,it}$, with a definitional
  sentence added.
- Papers 01 and 04: R formulas previously typeset in math mode moved
  to `\texttt{}`; paper 04's model equation restated with $Y_{ij}$ and
  $D_{bc,ij}$.
- Paper 02: $X = \text{Dbc}_{it}$ restated as $X = D_{bc,it}$.
- Paper 09: the stray $D_b$ variant unified to $D_{bc}$.
- Papers 08 and 09 (Architecture A only): $c_{bm}$ renamed to
  $\beta_{bm}$, with a notation paragraph added to each.
- Paper 07 (both architectures): notation paragraph added defining the
  shared label and the calibration; symbols left as $c_{bm}$
  deliberately.
- Papers 04 and 05: units paragraphs added giving the day-to-week
  conversion.
- Paper 02 and its four variants: analysis specifications relabeled
  A1/A2/A3 to M1/M2/M3 in text, tables, and figures; a
  `spec_display()` helper maps the archived data values at display
  time; figure drivers `03-render-figures.R`,
  `04-sensitivity-figures.R`, and `07-type1-cr2-figure.R` updated and
  the four affected figures regenerated.

## Prose terminology

Symbol conventions are settled here; prose terminology is settled by
`docs/30-terminology-consistency-plan.md`, executed 2026-07-30 with its
D1 spelling decision locked to **US English** (overriding that
document's original frequency-based recommendation of British forms).
The two passes are complementary: this file governs mathematics, that
one governs surface forms, acronym expansion, and controlled prose
vocabulary.

Not yet addressed: the outcome symbol is still written `Sx` in some
running prose; design labels still alternate between abbreviation and
spelled-out form; simulation studies are still lettered in paper 06;
and paper 11's sample-size convention (total rather than per-path)
requires the pending boundary re-run rather than an editorial fix.
