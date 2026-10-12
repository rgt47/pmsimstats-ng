# Biomarker Interaction Implementation: Architecture Comparison {.unlisted .unnumbered}
*2026-04-08 09:47 PDT (date of the original render; the creation time was not recorded separately)*

*Archived 2026-10-11; current treatment: `series/ch07` (biological rationale); see `docs/legacy/INDEX.md`.*

**Author.** pmsimstats team

**Purpose.** This early summary (April 2026) compares the power of the
biomarker-by-treatment interaction test under the two data-generating
architectures, covariance moderation and mean moderation, across
carryover half-lives, and sets out the biological assumptions that
would favor each. Its numbers come from a 50-replicate exploratory run
and are superseded by later work: see the note below.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

**Note added 2026-10-07.** The power figures below are exploratory
(50 replicates per cell, Monte Carlo standard error up to 0.07) and
predate the matched-$N$ re-run of paper 01, which reports Hybrid
carryover losses of 6.9% under mean moderation and 15.4% under
covariance moderation at $N = 70$ (`analysis/report/NOTATION.md`,
Part 2). The mechanistic account in the second section is also
superseded: `docs/37` shows, for the paired-difference statistic, that
carryover leaves the contrast variance unchanged and acts on power only
through the biomarker coupling gap, and `docs/38` argues that the
covariance construct keeps the full drug-response variance off drug.
The architecture labels A and B used in the original text have been
replaced by the canonical names (mean moderation, covariance
moderation), and the design labeled "N-of-1" is the Hybrid design.

## Notation index and glossary

| Symbol | Meaning | Status |
|----------|----------------------------------------|------|
| $N$ | total number of participants across randomization paths | canonical |
| $c_{bm}$ | covariance-moderation strength: on-drug correlation of biomarker and $BR$ | canonical |
| $\beta_{bm}$ | mean-moderation multiplier | canonical |
| $\beta_D$ | treatment main effect | canonical |
| $b$, $b_1$, $b_2$ | standardized biomarker; its values for two participants | canonical |
| $B$ | biomarker | canonical |
| $BR$ | biological (drug) response component | canonical |
| $Y$ | symptom score | canonical |
| $D$ | binary drug state | canonical |
| $D_{bc}$ | exposure-decayed drug indicator (code `Dbc`) | canonical |
| $t_{sd}$, $t_{1/2}$, $\lambda$ | time since discontinuation; carryover half-life; decay rate $\ln 2/t_{1/2}$ | canonical |
| $\sigma_{BR}$ | standard deviation of $BR$ | canonical |

Table: Notation index for the architecture comparison

**Glossary.** *Covariance moderation* (formerly Architecture B, data
value `mvn`): the interaction is a treatment-state-dependent
biomarker-$BR$ correlation. *Mean moderation* (formerly Architecture A,
`mean_moderation`): the biomarker shifts the mean drug response.
*Hybrid*: the Hendrickson et al. (2020) design (open-label titration,
blinded discontinuation, brief crossover), called "N-of-1" in the
original tables. *CO*: two-period crossover. *OL+BDC*: open-label
titration then blinded discontinuation. *Relative loss*: power at
$t_{1/2} = 0$ minus power at $t_{1/2} = 1$, divided by power at
$t_{1/2} = 0$.

## Empirical Comparison

Identical simulation configuration across both architectures (driver
`01-compare-dgp-architectures.R` in the external directory
`archive/architecture_comparison/`, outside the repository;
inspected 2026-10-07):

- Sample size: $N = 70$ in total, allocated across paths by
  `generateSimulatedResults()`
- Moderation strength: $c_{bm} = \beta_{bm} = 0.45$ (a shared value,
  per canonical rule 2)
- Data-generating correlations: AR(1) with 0.7 per week within each
  factor, cross-factor 0.1 at the same visit and 0.05 at different
  visits
- Analysis model: `nlme::lme` with `corCAR1` and the `bm:Dbc`
  interaction
- Replicates: 50 per cell

The driver ran two analysis models (half-life ignored, and $D_{bc}$
matched to the data-generating half-life); which one the tables report
was not recorded, and the output file
(`architecture_comparison_n50.rds`) was not found in the archive, so
the tables could not be re-derived. With the current package, every
covariance matrix of this configuration is positive definite (smallest
eigenvalue 0.043, OL+BDC at $t_{1/2} = 0$; verified 2026-10-07 by
`buildSigma()`), so the covariance-moderation cells did not rely on a
repaired matrix, provided the April 2026 code built the same matrices
(not verified).

### Covariance moderation (formerly Architecture B)

Power under carryover half-life variation:

| Design   | $t_{1/2}$ = 0 | $t_{1/2}$ = 0.5 | $t_{1/2}$ = 1.0 |
|----------|:----------:|:------------:|:------------:|
| Hybrid   |    0.82    |     0.64     |     0.50     |
| CO       |    0.40    |     0.40     |     0.32     |
| OL+BDC   |    0.62    |     0.38     |     0.24     |

Table: Covariance-moderation power by design and carryover half-life (50 replicates per cell)

Power declines substantially with increasing carryover (relative loss
from $t_{1/2} = 0$ to 1, computed from the table):

- Hybrid: 39% relative loss (corrected 2026-10-07; the original text
  gave 68%, which the table does not support)
- CO: 20% relative loss
- OL+BDC: 61% relative loss

### Mean moderation (formerly Architecture A)

Power under carryover half-life variation:

| Design   | $t_{1/2}$ = 0 | $t_{1/2}$ = 0.5 | $t_{1/2}$ = 1.0 |
|----------|:----------:|:------------:|:------------:|
| Hybrid   |    0.74    |     0.72     |     0.68     |
| CO       |    0.38    |     0.36     |     0.34     |
| OL+BDC   |    0.62    |     0.58     |     0.54     |

Table: Mean-moderation power by design and carryover half-life (50 replicates per cell)

Power is largely preserved under carryover:

- Hybrid: 8% relative loss
- CO: 11% relative loss
- OL+BDC: 13% relative loss

At 50 replicates the standard error of a difference between two cells
is about 0.09, so only the covariance-moderation declines in Hybrid and
OL+BDC are clearly larger than Monte Carlo error; the mean-moderation
losses and the CO loss are within it.

## Why the Architectures Diverge

### Mean moderation: the biomarker shift is confined to on-drug visits

Under mean moderation the package shifts each participant's $BR$ at
on-drug visits by $\beta_{bm}\,\sigma_{BR}\,b$ (`R/generateData.R`,
default `moderation_scaling = 'constant'`; inspected 2026-10-07). The
shift does not depend on the carryover half-life. Carryover acts only
on the off-drug $BR$ mean, which carries no biomarker information, so
the interaction signal is essentially unchanged; `docs/37` (Section
5.3) derives that, for the paired-difference statistic, power under
mean moderation does not depend on $t_{1/2}$ at all.

In the multiplicative form used by paper 01, with treatment main
effect $\beta_D$, the expected outcome difference between two
participants with standardized biomarker values $b_1$ and $b_2$ in the
same drug state $D$ is

$$E(Y \mid b_1, D) - E(Y \mid b_2, D) = \beta_D\,\beta_{bm}\,(b_1 - b_2)\,D .$$

On drug ($D = 1$) the difference is $\beta_D\beta_{bm}(b_1 - b_2)$
whatever the half-life; off drug it is zero. (Corrected 2026-10-07: the
original formula conditioned on $D = 1$ while retaining a factor $D$,
used an undefined $\beta_1$, and described the difference as both
constant and scaling with carryover.) Under the sign convention of
`NOTATION.md`, $\beta_D$ is negative, so a larger biomarker value
lowers the symptom score further on drug.

### Covariance moderation: the coupling decays off drug

Under covariance moderation the original account named two mechanisms:

1. **Mean blurring**: the $BR$ means during off-drug periods are
   inflated by residual carryover, making them resemble on-drug
   $BR$ means, which was said to reduce the treatment contrast that
   the analysis model uses to distinguish drug effect from no-drug.

2. **Correlation erosion**: the biomarker-$BR$ correlation during
   off-drug periods decays exponentially with time since
   discontinuation (`R/generateData.R`, `buildSigma()`; inspected):

   $$\mathrm{Cor}(B, BR_t) = c_{bm}\, e^{-\lambda t_{sd}} \quad \text{off drug after exposure},$$

   with $c_{bm}$ on drug and zero before first exposure.

The second rule is a specification choice of the data-generating
process, not a structural consequence of it. (Corrected 2026-10-07:
the original text called it "not merely a statistical modeling choice"
and attributed the decay to a shrinking drug-mediated share of $BR$
variance; the construct holds the $BR$ standard deviation fixed at
every visit, so no such share shrinks; `docs/38`, Section 5.)

Later work revises the two-mechanism account (note added 2026-10-07).
For the paired-difference statistic, `docs/37` derives that carryover
does not change the contrast variance (Proposition 1) and that the
mean shift is absorbed by the intercept, so the whole effect of
carryover on power runs through the coupling gap between on-drug and
off-drug visits (Propositions 2 and 3); mean blurring contributes
nothing there. For the `lme` analysis with $D_{bc}$ reported here this
was not derived. The losses in the table above are 20-61% (Hybrid
39%), not uniformly 40-60% as the original text stated.

## Biological Assumptions

The choice between architectures reflects fundamentally different
assumptions about the mechanism by which the biomarker moderates
treatment response.

### When mean moderation is appropriate

Mean moderation is appropriate when the
biomarker governs the **magnitude of the drug's biological effect**
for each individual participant. Under this mechanism, the drug
effect for a participant is proportional to their biomarker value:
higher biomarker predicts larger treatment response.

Example: a pharmacogenetic biomarker that determines hepatic
metabolism capacity, where the genotype fixes the size of each
participant's response while the drug is present. The interaction is
a property of each individual's biology. (Note added 2026-10-07: the
original example had high-expression participants metabolizing the
drug more rapidly yet responding more strongly; faster metabolism
ordinarily lowers exposure and response, so the direction of the
example depends on whether the variant affects activation or
clearance. The example is illustrative and was not checked against a
specific drug.)

### When covariance moderation is appropriate

Covariance moderation is appropriate when
the biomarker's predictive value emerges from **shared variance**
with a drug-responsive physiological component, and this shared
variance decays as the drug is eliminated.

Example: a biomarker that reflects drug pathway engagement. The
biomarker and drug response are jointly determined by the same
underlying pathway state, and their correlation weakens off drug. The
interaction is a property of the joint response to the drug, not an
inherent property of the individual. (Note added 2026-10-07: in every
current data-generating process the biomarker $B$ is a single
pre-treatment, participant-level value, so it cannot itself decay
during off-drug periods; the original example, an active drug
metabolite concentration, describes a time-varying biomarker that the
simulation does not represent. Only the coupling $c_{bm}e^{-\lambda
t_{sd}}$ decays.)

### Selecting the appropriate architecture

For simulation studies of predictive biomarker-moderated treatment
effects, the choice should be driven by the clinical hypothesis:

- **Baseline biomarker predicts intrinsic drug responsiveness:**
  use mean moderation.
- **Biomarker reflects drug pathway engagement:** use covariance
  moderation.

The power implications are substantial: the choice between
architectures determines whether carryover substantially erodes
power (covariance moderation) or whether power is largely preserved
(mean moderation). This choice must be made based on clinical
reasoning, not statistical convenience (argued).

## References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020; 2:13. doi:10.3389/fdgth.2020.00013.
2. pmsimstats team. `docs/37-carryover-power-closed-form.md`,
   `docs/38-pharmacodynamic-carryover-model.md`,
   `analysis/report/NOTATION.md`.
