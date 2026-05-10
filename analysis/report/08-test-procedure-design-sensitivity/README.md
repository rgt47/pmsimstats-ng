# 08-test-procedure-design-sensitivity

*2026-05-07 15:30 PDT*

## Working title

Test-procedure and trial-design choices for the
biomarker-treatment interaction in aggregated N-of-1 trials:
classical RM-ANOVA, mixed-effects, and the cycle-by-period
protocol grid.

## Origin

This manuscript combines two complementary methodological
threads documented separately in the project notes:

- `docs/26-rm-anova-interaction-test.md` derives the closed-form
  strict repeated-measures ANOVA test of the
  biomarker-treatment interaction in a balanced split-plot
  design and relates it to the linear-mixed Wald test.
- `docs/27-cycle-period-design-sweep-plan.md` specifies a
  parameter grid that varies the protocol-level design
  parameters (cycle count, on-drug / off-drug period length,
  on/off symmetry, optional run-in / run-out) and quantifies how
  the biomarker-treatment interaction test responds.

The two threads address distinct but coupled choices the
analyst makes (test procedure and trial-design parameters), and
treating them jointly is the contribution of this manuscript.

## Scope

1. Derive the closed-form strict RM-ANOVA $F$-statistic for the
   biomarker-stratified treatment-by-time interaction in a
   balanced split-plot design with categorical predictors and
   the sphericity assumption.
2. Relate the strict RM-ANOVA test to the Wald $t$-test on the
   biomarker-by-treatment fixed-effects coefficient in a linear
   mixed model with continuous-time AR(1) residual correlation,
   identifying the conditions on residual covariance under which
   the two tests agree and the regimes in which they diverge.
3. Conduct a Monte Carlo simulation study that crosses three
   test procedures (strict RM-ANOVA, the standard pmsimstats
   linear-mixed analysis, and a generalised-estimating-equations
   comparator) against a grid of design parameters: cycle count,
   on-drug period length, off-drug period length, on/off
   symmetry.
4. Synthesise a design-and-test recommendation table for typical
   chronic-condition N-of-1 trial scenarios.

## Simulation substrate

The simulation study uses the **`implementations/simple/`**
codebase as its substrate, rather than the full pmsimstats
package stack. The simple implementation is a single-file R
script that already implements the four Hendrickson-relevant
designs (OL, CO, N-of-1, OL+BDC) with an exposed parameter grid
covering sample size, biomarker moderation, carryover half-life,
carryover-adjustment toggle, and dropout probability. Three
features make it the appropriate substrate for this paper's
sweeps:

- The design-parameter grid is already first-class in the
  driver, which makes adding the cycle-count, period-length,
  and on/off-symmetry axes a small extension rather than a
  refactor of the full pipeline.
- The reduced-complexity DGP (mean-moderation Architecture A
  only, ANCOVA-on-phase-means analysis) isolates the design-
  parameter and test-procedure axes from the response-curve-
  family and covariance-structure axes that the full package
  exposes; these other axes are out of scope for the present
  paper.
- Pedagogical clarity: the cycle-by-period sensitivity result is
  more transparent when the substrate is simple enough that an
  external reviewer can read the simulation script end to end.

The driver scripts at `analysis/scripts/08-test-procedure-design-
sensitivity/` will be created when the simulation study is
implemented and will source `implementations/simple/simulation.R`
and extend its parameter grid with the cycle-by-period axes.

## Relationship to other papers in the series

- `01-dgp-mean-moderation-vs-mvn`: shares the
  mean-moderation DGP. The present paper's
  cycle-by-period grid uses Architecture A throughout; the
  Architecture-B comparison is held out as future work.
- `02-carryover-sensitivity`: the present paper's carryover
  axis is held at the prazosin-PTSD reference value
  ($t_{1/2} = 1$ week). Cross-paper sensitivity to varying
  carryover is reported in `02-`.
- `04-treatment-main-effect`: clinical companion that uses a
  similar simulation substrate but targets the main treatment
  effect; the present paper targets the biomarker-treatment
  interaction.
- `05-nof1-design-sensitivity`: methodologically adjacent (also
  a design-parameter sensitivity study) but focused on the main
  treatment effect rather than the interaction.

## Driver scripts

The driver directory at `analysis/scripts/08-test-procedure-
design-sensitivity/` is to be created. The simulation will
source the simple-implementation driver at
`implementations/simple/simulation.R` and extend it with:

- Cycle-count, period-length, and on/off-symmetry axes added to
  the design grid.
- Three test-procedure variants applied to each cell.
- ADEMP-formatted output following Morris, White and Crowther
  (2019).

## Author

pmsimstats team
