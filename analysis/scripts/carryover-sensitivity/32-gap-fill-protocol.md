# Gap-fill simulation for the G9 recommendation: working protocol

*2026-10-26 17:19 PDT*

**Author.** pmsimstats team

**Purpose.** Paper 02's reconstruction (commit `f2d11da14`) rests on
runs made under the corrected simulation code (jobs 21, 26, 27, 28,
30) and names five open items in its Limitations: G9's behavior
under dropout is unmeasured; the 3,000-replicate size diagnostic
predates the cross-factor constant correction; the decay-shape and
half-life misspecification checks are single cells at one effect
size; the mean-moderation and autocorrelation scenarios are a single
cell or untested under the corrected code; and the prognostic
confound has been tested only at strength 0.3 (feasible) and 0.6
(infeasible), with the prog-PB expectancy mechanism unchecked at
larger $N$. This protocol specifies one driver,
`32-run-g9-gap-fill.R`, that closes all five. It is a working
protocol for an exploratory study, not a pre-registration.

## Aims

1. Measure G9's size and power under monotone dropout, MCAR and MAR,
   at 10%, 20% and 30%, against G1, G3, G8 and the F-family on the
   same trials (the 2026-08 dropout block covered G1 to G3 only and
   predates both code corrections).
2. Replace the pre-correction size diagnostic with a 2,000-replicate
   null at two half-lives and two sample sizes under the corrected
   code, reporting the calibration ratio $\kappa$ and Type I error
   for every specification with and without CR2.
3. Extend the decay-shape, half-life, architecture and
   autocorrelation scenarios of Stage 4 from single cells to small
   grids crossing effect size and $N$.
4. Run a feasible prognostic-confound ladder, chosen from the
   closed-form ceilings of `docs/65`, and check the prog-PB
   expectancy mechanism at $N = 140$ and $280$.

## Data-generating mechanisms

The refined process throughout: covariance moderation (except the
mean-moderation block), graded coupling, separable AR(1) with
$\rho = 0.7$ unless varied, $c_1 = 0.2$ at every lag, anchored
exponential carryover at $t_{1/2} = 0.5$ week unless varied, Hybrid
design, stop on non-positive-definite matrices (the driver runs a
pre-flight eigenvalue check over every cell and refuses to start if
any cell is infeasible). Default response and baseline parameters
of `simulation-core.R`.

| Block | Factor varied | Levels | Other factors | Cells |
|---|---|---|---|---|
| A1 dropout | rate $\times$ mechanism | $\{0.1, 0.2, 0.3\} \times \{$MCAR, MAR$\}$ | $c_{bm} \in \{0, 0.45\}$, $N \in \{35, 70\}$ | 24 |
| A2 size | $t_{1/2}$ | $\{0.5, 1\}$ | $c_{bm} = 0$, $N \in \{35, 70\}$ | 4 |
| B1 decay shape | Weibull $k$ (data), exponential analysis | $\{0.5, 2\}$ | $c_{bm} \in \{0, 0.25, 0.45\}$, $N \in \{35, 70\}$ | 12 |
| B2 half-life | true $t_{1/2}$ (weeks) | $\{0.2, 2, 4\}$ | $c_{bm} \in \{0, 0.25\}$, $N \in \{35, 70\}$ | 12 |
| B3 architecture | mean moderation | one level | $c_{bm} \in \{0, 0.25, 0.45\}$, $N \in \{35, 70\}$ | 6 |
| B4 autocorrelation | $\rho$ | $\{0.5, 0.9\}$ | $c_{bm} \in \{0, 0.30\}$, $N \in \{35, 70\}$ | 8 |
| C1 prognostic ladder | pathway $\times$ strength | prog-PB, prog-TV-const at $\{0.40, 0.45\}$; prog-BL at $\{0.6, 0.8\}$ | $c_{bm} \in \{0, 0.30\}$, $N = 70$ | 12 |
| C1b prog-BL at $c_{bm} = 0.45$ | strength | $0.45$ | $c_{bm} \in \{0, 0.45\}$, $N = 70$ | 2 |
| C2 prog-PB at large $N$ | $N$ $\times$ $c_{pb}$ | $\{140, 280\} \times \{0.3, 0.45\}$ | $c_{bm} = 0$ | 4 |

Dropout is monotone: once a patient misses an occasion, all later
occasions are missing. MAR scales the per-occasion dropout
probability with baseline severity so that the sicker half of
patients is about twice as likely to leave, with the marginal rate
preserved. The B1 cells generate Weibull data and analyze with the
exponential predictor every specification assumes. In B2 the oracle
G3 uses the true half-life; G4 and G9 select over
$\{0.25, 0.5, 1, 2\}$, so $t_{1/2} = 0.2$ and $4$ lie outside the
candidate grid and $2$ is its boundary; F1 is fixed at 0.2 week. The
C1 strengths are inside the `docs/65` ceilings ($c_{pb}$, $c_{tv}$
up to 0.451 and $c_{bl}$ up to 0.815 at $c_{bm} = 0.30$; $c_{bl}$ up
to 0.496 at $c_{bm} = 0.45$); the ladder reaches to within 0.001
to 0.05 of each ceiling, which is where a loss would first appear
if one exists. C2 asks whether the small population contrast
prog-PB induces at $c_{bm} = 0$ ($-0.14$ at 0.3, about $-0.2$ at
0.45) becomes detectable Type I error inflation as $N$ grows; the
F-family, which carries the expectancy term, is the comparator. B4
uses $c_{bm} = 0.30$ rather than 0.45 because the coupling ceiling
falls with $\rho$ and the driver's pre-flight check found
$c_{bm} = 0.45$ infeasible at $\rho = 0.9$ (smallest eigenvalue
$-0.28$); the $c_{bm}$ ceiling on Hybrid at $t_{1/2} = 0.5$ week is
0.506, 0.518 and 0.391 at $\rho = 0.5$, 0.7 and 0.9 respectively.
The two $\rho$ levels therefore share one effect size.

**Pre-flight and smoke status (2026-10-26).** `--preflight` reports
all 84 cells feasible (smallest eigenvalue 0.0028, at the top of the
prognostic ladder). `--smoke` ran two replicates on two cells per
block, covering dropout, a Weibull data shape, $\rho = 0.9$ and
$N = 280$, with every specification converging.

## Estimands and methods

The interaction coefficient `bm:Dbc` (or `bm:Db` for G1, G2, G6, G7
and E0). Every specification of `simulate_cell_combined()` is fit to
the same simulated trial each replicate: G1 to G9, E0, F1, F1A, F1A0,
and F1A0B derived from F1A0's warp-speed bootstrap statistic pooled
within the cell (one resample per replicate). Reported in paper 02:
G1, G3, G4, G8, G9, E0, F1, F1A0B.

## Performance measures

Per (cell, method): rejection rate at $\alpha = 0.05$ (Type I error
at $c_{bm} = 0$, power otherwise) with its MCSE; non-convergence
rate; mean reported standard error, empirical standard error of the
estimate, and their ratio $\kappa$; mean estimate. Replicates: 500
per cell, except 2,000 in A2 and 1,000 in C2, so that the MCSE on a
rejection rate near 0.05 is about 0.005 and 0.007 respectively.
Monotone dropout removes whole tails of a patient's series, so the
denominator for convergence is the replicate, not the observation.

## Cost and scheduling

Calibrated on jobs 27 and 28 (about 6 to 7 core-seconds per
replicate for the fourteen-specification fit at $N \le 70$, six
cores on the simulation host; $N = 140$ and $280$ scale roughly
linearly):

| Block | Replicates | Estimated wall time |
|---|---|---|
| A (A1 + A2) | $24 \times 500 + 4 \times 2000 = 20{,}000$ | about 6.5 h |
| B (B1 to B4) | $38 \times 500 = 19{,}000$ | about 6 h |
| C (C1 + C1b + C2) | $14 \times 500 + 4 \times 1000$ (at 2 to 4 times cost) | about 6 h |

The driver takes `--block A|B|C|all` so each block can be submitted
as its own `simrun` job and fetched independently; checkpoints are
per (cell, replicate), so an interrupted job resumes. `--preflight`
builds every cell's covariance matrix and reports feasibility without
simulating; `--smoke` runs two replicates of two cells per block.

## Analysis plan for paper 02

- A1 replaces the dropout sentence in Limitations with a result:
  a table of G9 and G1 power by rate and mechanism, and the G9 minus
  G1 gap, to establish whether exposure weighting's advantage
  survives dropout and whether G9's size does.
- A2 replaces Table 4 (the pre-correction $\kappa$ diagnostic) with
  its corrected-code equivalent at two half-lives and two $N$.
- B1 to B4 extend Stage 4's scenario tables from one cell each to
  the small grids above, keeping the same columns.
- C1 and C1b add a "feasible ladder" row block to the Stage 5 size
  and power tables, so that the robustness claim is stated up to
  the ceiling rather than at 0.3 only. C2 is reported in the
  expectancy-terms discussion: if G9's Type I error at $N = 280$
  under $c_{pb} = 0.45$ exceeds the band while F1's does not, the
  recommendation gains an explicit caveat for large trials with a
  placebo-tracking biomarker.

## Outputs

`analysis/data/carryover-sensitivity/g9-gap-fill/` with
`all.rds` (per replicate), `summary.rds` (per cell and method),
`preflight.rds` (per cell feasibility), `meta.rds`, and
`trials/` checkpoints.

## What this protocol does not cover

OL+BDC and CO (the recommendation is for Hybrid; Stages 1 to 4
already show where it transfers); the dual-channel architecture;
informative dropout mechanisms beyond the baseline-severity MAR
rule; factor interactions (each block varies one factor against
effect size and $N$ only).
