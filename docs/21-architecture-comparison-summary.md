# Biomarker Interaction Implementation: Architecture Comparison

## Empirical Comparison

Identical simulation configuration across both architectures:
- Sample size: N = 70
- Biomarker correlation / moderation: c.bm = 0.45
- Analysis model: nlme::lme with corCAR1 and bm:Dbc interaction
- Replicates: 50 per cell

### Architecture B: Differential Correlation (MVN)

Power under carryover half-life variation:

| Design   | t_half = 0 | t_half = 0.5 | t_half = 1.0 |
|----------|:----------:|:------------:|:------------:|
| N-of-1   |    0.82    |     0.64     |     0.50     |
| CO       |    0.40    |     0.40     |     0.32     |
| OL+BDC   |    0.62    |     0.38     |     0.24     |

Power declines substantially with increasing carryover:
- N-of-1: 68% relative loss
- CO: 20% relative loss
- OL+BDC: 61% relative loss

### Architecture A: Direct Mean Moderation

Power under carryover half-life variation:

| Design   | t_half = 0 | t_half = 0.5 | t_half = 1.0 |
|----------|:----------:|:------------:|:------------:|
| N-of-1   |    0.74    |     0.72     |     0.68     |
| CO       |    0.38    |     0.36     |     0.34     |
| OL+BDC   |    0.62    |     0.58     |     0.54     |

Power is largely preserved under carryover:
- N-of-1: 8% relative loss
- CO: 11% relative loss
- OL+BDC: 13% relative loss

## Why the Architectures Diverge

### Architecture A: Biomarker scaling is proportional

Under direct mean moderation, carryover reduces the magnitude of
the drug exposure variable D during off-drug periods, but the
biomarker moderation coefficient (beta_bm) remains a fixed
population parameter. The interaction (bm * D) has reduced
variance because D has less contrast, but the signal-to-noise
ratio degrades only modestly since the noise structure is
independent of carryover.

Mathematically, the expected outcome difference between two
participants with biomarker values b1 and b2, both on drug, is:

$$E[Y | b_1, D=1] - E[Y | b_2, D=1] = \beta_1 \cdot \beta_{bm}
\cdot (b_1 - b_2) \cdot D$$

This difference is constant regardless of carryover: if D is
replaced by D * carryover_factor, the difference scales by the
same factor, but the ratio is preserved.

### Architecture B: Correlation erosion compounds mean blurring

Under the MVN approach, carryover operates on two levels
simultaneously:

1. **Mean blurring**: The BR means during off-drug periods are
   inflated by residual carryover, making them resemble on-drug
   BR means. This reduces the treatment contrast that the
   analysis model uses to distinguish drug effect from no-drug.

2. **Correlation erosion**: The BM-BR correlation during off-drug
   periods decays exponentially with time since discontinuation:
   
   $$\text{Cor}(BM, BR|t_{sd}) = c_{bm} \cdot e^{-\lambda t_{sd}}$$
   
   This is not merely a statistical modeling choice -- it is a
   structural consequence of the DGP. The correlation between
   biomarker and BR reflects the degree to which the drug effect
   is active; as the drug washes out, the correlation weakens
   because the drug-mediated component of BR variance shrinks
   relative to the drug-independent component.

The analysis model's `bm:Dbc` interaction coefficient depends on
the covariance between biomarker and Dbc-weighted outcomes. Under
Architecture B, this covariance is being attacked from both sides:
the treatment contrast in Dbc is compressed, and the biomarker-BR
correlation that generates the covariance signal is simultaneously
decaying. This dual effect produces the 40-60% power losses
observed in the simulations.

## Biological Assumptions

The choice between architectures reflects fundamentally different
assumptions about the mechanism by which the biomarker moderates
treatment response.

### When Architecture A is appropriate

Architecture A (direct mean moderation) is appropriate when the
biomarker governs the **magnitude of the drug's biological effect**
for each individual participant. Under this mechanism, the drug
effect for a participant is proportional to their biomarker value:
higher biomarker predicts larger treatment response.

Example: A pharmacogenetic biomarker that determines hepatic
metabolism capacity. Participants with high enzyme expression
metabolize the drug more rapidly, but when the drug is present
(regardless of whether it has completely washed out between
doses), they experience proportionally larger biological response.
The interaction is a property of each individual's biology.

### When Architecture B is appropriate

Architecture B (differential correlation, MVN) is appropriate when
the biomarker's predictive value emerges from **shared variance**
with a drug-responsive physiological component, and this shared
variance decays as the drug is eliminated.

Example: A dynamic biomarker that reflects current drug pathway
engagement (e.g., active drug metabolite concentration). The
biomarker and drug response are jointly determined by the same
underlying pathway state. During off-drug periods, both the
biomarker and the drug-responsive component of BR decay in
parallel, so their correlation weakens. The interaction is a
property of the joint response to the drug, not an inherent
property of the individual.

### Selecting the appropriate architecture

For simulation studies of predictive biomarker-moderated treatment
effects, the choice should be driven by the clinical hypothesis:

- **Baseline biomarker predicts intrinsic drug responsiveness:**
  Use Architecture A (direct mean moderation)

- **Biomarker reflects current drug pathway engagement:**
  Use Architecture B (differential correlation)

The power implications are substantial: the choice between
architectures determines whether carryover substantially erodes
power (Architecture B) or whether power is largely preserved
(Architecture A). This choice must be made based on clinical
reasoning, not statistical convenience.

---

*Rendered on 2026-04-08 at 09:47 PDT.*
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/21-architecture-comparison-summary.md*
