# Plan: Expanding report.Rmd with N-of-1 Power Methods Discussion

## Objective

Add an extended section (1--2 pages, ~800--1200 words) to `report.Rmd`
covering contemporary best practices for calculating power and sample
size in N-of-1 trial designs. This section should serve dual purposes:
(a) situate the simulation study within the existing analytic
framework, and (b) justify why a simulation approach is needed beyond
closed-form solutions.

## Placement in report.Rmd

Insert as a new subsection within **Section 2 (Methods)**, between the
current "Sample Size and Power Considerations" introduction material
(lines 84--88) and the "Simulation Framework" subsection (line 112).

Proposed heading hierarchy:

```
# Methods
## Standard Approaches to Power Analysis in N-of-1 Trials   <-- NEW
### Variance decomposition in crossover designs              <-- NEW
### Closed-form power for the paired-cycles approach         <-- NEW
### Extensions to series of N-of-1 trials                    <-- NEW
### Limitations motivating simulation                        <-- NEW
## Simulation Framework  (existing)
```

## Content Outline

### 2.1 Standard Approaches to Power Analysis in N-of-1 Trials

Opening paragraph framing the problem: power in N-of-1 trials depends
on a different set of variance components than parallel-group RCTs.
The fundamental unit of replication is the treatment period within a
patient, not the patient themselves.

### 2.1.1 Variance Decomposition in Crossover Designs

- The total variance in an N-of-1 trial decomposes into:
  - Between-patient variance (sigma^2_b)
  - Within-patient, between-period variance (sigma^2_w)
  - Residual/measurement error (sigma^2_e)
- In N-of-1 designs, each patient serves as own control, eliminating
  sigma^2_b from the treatment effect estimator.
- Contrast with parallel-group RCTs where sigma^2_b inflates the
  standard error of the treatment effect.
- Key refs: Senn & Julious 2024, Schmid et al. 2014 (AHRQ Ch. 4)

### 2.1.2 Closed-Form Power for the Paired-Cycles Approach

- Senn & Julious (2024) paired-cycles framework: each AB or BA pair
  yields a within-patient treatment difference.
- For a single patient with m paired cycles, the test statistic is a
  one-sample t-test on m differences with (m-1) df.
- Power formula:
  ```
  Power = P(|T| > t_{alpha/2, m-1})
  where T ~ t_{m-1}(delta / (sigma_d / sqrt(m)))
  delta = true treatment effect
  sigma_d = SD of within-patient period differences
  ```
- Senn (2019) also gives formulae determining both the number of
  patients (N) and cycles per patient (m) jointly.
- Their key result: for a series of N-of-1 trials analyzed via
  random-effects meta-analysis, power depends on:
  ```
  N (number of patients)
  m (number of cycles per patient)
  sigma^2_w (within-patient variance of treatment differences)
  tau^2 (between-patient variance of treatment effects)
  delta (population mean treatment effect)
  ```
- Key refs: Senn & Julious 2024, Senn 2019

### 2.1.3 Extensions to Series of N-of-1 Trials

- Zucker et al. (2010): combining individual N-of-1 results via
  random-effects meta-analysis (DerSimonian-Laird).
- The effective sample size for the population-level estimate depends
  on the ratio tau^2 / sigma^2_w --- when between-patient
  heterogeneity is large relative to within-patient noise, adding
  more patients helps more than adding more cycles.
- Schmid et al. (2014, 2021): Bayesian hierarchical models provide
  an alternative that borrows strength across patients.
- Blackston et al. (2019) and Chen et al. (2020): simulation-based
  comparisons showing N-of-1 power advantages over parallel RCTs
  depend heavily on the variance ratio and number of cycles.
- Key refs: Zucker et al. 2010, Schmid et al. 2014, Schmid &
  Staudenmayer 2021, Blackston et al. 2019, Chen et al. 2020

### 2.1.4 Complications Requiring Simulation

- **Carryover effects**: closed-form solutions assume either no
  carryover or complete washout. When carryover is partial (as with
  prazosin), the effective treatment contrast is attenuated and
  variance components shift in ways not captured by standard formulas.
  (Yap et al. 2024, Tang & Landes 2020)
- **Serial correlation**: within-patient observations are typically
  autocorrelated; ignoring this inflates Type I error or distorts
  power. (Tang & Landes 2020)
- **Unbalanced designs**: practical N-of-1 trials may have unequal
  period lengths, missing periods, or adaptive stopping rules that
  violate balanced-design assumptions.
- **Non-normal outcomes**: nightmare frequency may be better modeled
  as count data; closed-form power assumes normality.
- **Joint estimation**: when the goal is both individual-level
  decisions and population inference, the power trade-off between N
  and m has no simple closed-form solution.
- Conclude: these complications collectively motivate the Monte Carlo
  simulation approach adopted in this study.

## Research Steps (New Session with Zotero MCP)

1. **Query Zotero 'nof1' collection** for all items; identify any
   additional power/sample-size references not yet in `nof1-pgt.bib`.
2. **Retrieve full text** (via PubMed PMC) for key methodological
   papers:
   - Senn 2019 (PMID: 28882093) -- core sample size paper
   - Senn & Julious 2024 (PMID: 38365817) -- paired cycles
   - Tang & Landes 2020 (PMID: 32040510) -- t-tests with correlation
   - Yap et al. 2024 -- carryover and power
   - Schmid et al. 2014 (AHRQ Ch. 4) -- design chapter
3. **Extract specific formulas and notation** from Senn (2019) and
   Senn & Julious (2024) to present in consistent notation.
4. **Check for recent (2024--2026) papers** on N-of-1 power methods
   via PubMed search that may supersede or extend Senn (2019).
5. **Write the section** in report.Rmd with LaTeX math blocks.
6. **Update nof1-pgt.bib** with any new references found.

## Style Notes

- Academic tone, no emojis (per CLAUDE.md)
- LaTeX math for all formulas (the report uses pdf_document output)
- Use `\sigma^2_w`, `\tau^2`, `\delta` notation consistently with
  existing manuscript conventions
- Wrap lines at 78 characters in the Rmd source
- Cite using `[@key]` natbib format already used in report.Rmd

## Existing References in nof1-pgt.bib Already Relevant

| Cite Key | Short Description |
|----------|-------------------|
| `sennSampleSizeConsiderations2019` | Core sample size paper |
| `sennAnalysisContinuousData2024` | Paired-cycles analysis |
| `schmidStatisticalDesignAnalytic2014` | AHRQ design chapter |
| `schmidBayesianModelsNof12021` | Bayesian N-of-1 models |
| `zuckerIndividualNof1Trials2010` | Meta-analysis of N-of-1 |
| `blackstonComparisonAggregatedNof12019` | Simulation comparison |
| `blackstonComparisonAggregatedNof12019` | Simulation comparison |
| `tangTtestsNof1Trials2020` | Serial correlation |
| `yapBehavioralCarryoverEffect2024` | Carryover and power |
| `araujoComparisonFourMethods2014` | Analysis method comparison |

## Potential New References to Find

- Duan et al. (2013) AHRQ full report on N-of-1 trials
- Percha et al. (2019) or similar -- modern power software for N-of-1
- Any 2024--2026 update to the Senn (2019) framework
- Senn (2002) "Cross-over Trials in Clinical Research" (textbook,
  foundational power formulas for crossover designs)

## Deliverables

1. Updated `report.Rmd` with new ~1200-word subsection
2. Updated `nof1-pgt.bib` with any new references
3. Consistent notation with rest of manuscript
