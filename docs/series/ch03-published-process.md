# The published process at commit 58b32a9

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** generateData() and buildtrialdesign() as published at commit
58b32a9; all parameters; Figures 4 and 5.  
**Prerequisites:** ch02

---

## Contents

1. [Notation](#notation)
2. [Summary](#summary)
3. [Construct box](#construct-box)
4. [Setting and argument](#setting-and-argument)
5. [Derivations](#derivations)
6. [Numerical evidence](#numerical-evidence)
7. [Relation to the manuscripts](#relation-to-the-manuscripts)
8. [Limitations and open items](#limitations-and-open-items)
9. [Reproduce](#reproduce)
10. [References](#references)
11. [Change log](#change-log)

---

## Notation

*Local symbols used in this chapter only. For project-wide notation see
`analysis/report/NOTATION.md`.*

<!-- notation entries -->

---

## Summary

<!-- evidence-labeled results go here -->

---

## Construct box

*No Monte Carlo simulation results in this chapter.*

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 40 Sections 3 to 6 and 8; 36 Section 3; 53 Sections 3
  to 4; 55 Sections 3 to 4; 44 Section 2; 25 Section 4; 23
  Section 3.1; CM "The published compound-symmetry construct" and
  the three "2024 revision" subsections; Appendix I Sections 1
  to 3.

- New to write: the commit attribution table (`42ac030`, `58b32a9`,
  `8609f12`, `f6ee86d`, `f70f86d`, `325314f`, `06dac83`); what
  `implementations/original/` is and is not.
-->

### Commit history (42ac030 to 06dac83)

The upstream repository (`rchendrickson/pmsimstats` on GitHub) has
seven commits relevant to this project. The local vendored copy at
`implementations/original/` holds the 2024 state of that repository
(commit `06dac83`), not the published state (`58b32a9`). All commit
dates and messages are from the GitHub commit history (inspected;
not found in the local git log, which records a single "Initial
commit" for `implementations/original/`).

| Commit | Date | Summary | Role in project |
|---|---|---|---|
| `42ac030` | 2020-02-22 | Initial publication | First commit; the DGP as submitted. Data-generating code is identical to `58b32a9`; reference for "the published process" (inspected). |
| `58b32a9` | 2020-11-19 | Last pre-publication commit | Sixth commit; the canonical reference for all "orig" analyses in this project. DGP code unchanged from `42ac030`; analysis driver may differ. |
| `8609f12` | 2024-05-06 | Ron Thomas edits | First post-publication modification. Adds `scalefactor = 2` to the carryover exponent (halves the effective half-life), replaces the step-rule biomarker coupling with proportional scaling $c_{bm}\,\mu_{BR,t}/\mu_{BR,t-1}$, and introduces a timepoint-1 guard that forces $\mathrm{Cor}(B,BR_1) = 0$ on all paths (inspected, docs/53 Section 1; docs/55 Section 4; docs/legacy/23 Section 2.4). |
| `f6ee86d` | 2024-05-08 | Autocorr typo fix | Corrects an index error in the same-visit cross-factor correlation assignment (`correlations[n2,n1]`). No change to carryover or coupling (inspected, docs/53 Section 1). |
| `f70f86d` | ca. 2024-05 | Carryover added (DGP unchanged) | Reverts the cross-factor index fix from `f6ee86d` (regression). No functional changes to the core DGP files; described as "carryover added" in the commit message but the carryover term was already present (inspected, docs/legacy/04 Section A.3.3). |
| `325314f` | 2024-05-10 | $D_{bc}$ added to analysis (broken) | Introduces the continuous drug indicator $D_{bc}$ with exponential decay into `generateSimulatedResults()`. The implementation contains two bugs: a reference to `op` instead of a scalar half-life parameter, and an off-by-one in the decay argument. Analysis results from this commit are unreliable (inspected, docs/legacy/04 Section A.3.4). |
| `06dac83` | 2024-11-27 | Current upstream HEAD | The 2024 state copied as `implementations/original/`. Retains the scale factor 2 recursion and the proportional-scaling biomarker coupling from `8609f12`. Fixes the two $D_{bc}$ bugs from `325314f`. `generateSimulatedResults()` calls `generateData()` without the `scalefactor` argument, so the default $s_f = 2$ applies (inspected, docs/53; docs/55). |

Table: Commit history of `rchendrickson/pmsimstats` from initial
publication to current upstream HEAD (inspected from GitHub; dates
for `f70f86d` inferred from sequence)

### What implementations/original/ is and is not

The directory `implementations/original/` holds a vendored copy of
the upstream `rchendrickson/pmsimstats` repository at commit
`06dac83` (2024-11-27). It is **not** the published code: it post-
dates the publication by four years and includes the 2024 revisions
described above (scale factor 2 in the carryover exponent,
proportional-scaling biomarker coupling, `Dbc` in the analysis).

The published process in this project refers to commit `58b32a9`.
A snapshot of that commit's `R/generateData.R` is archived at
`analysis/data/quick-sim/hendrickson-58b32a9/`. The `original`
collection is kept for cross-implementation comparisons against
the current upstream; it is not a reproduction of the published
simulation.

The `original` collection has no carryover sensitivity driver:
its `generateSimulatedResults()` always calls `generateData()` with
`scalefactor` defaulting to 2, so it cannot be used to reproduce
published (zero-scalefactor) results without modification. Use the
`original-extended` or `tidyverse` collections for new analyses.

### Overview of the published process

<!-- From docs/36 Section 3 "The published framework" -->

The published process is the data-generating process (DGP) of
Hendrickson et al. (2020), `pmsimstats` at commit `58b32a9` (19
November 2020, the last commit before publication). In this project
it is called "orig". The data-generating code at `58b32a9` is
identical to the initial commit `42ac030` (2020-02-22), inspected
from the GitHub commit history.

Each participant's outcome is

$$
Y_{it} = \mathrm{BL}_i - \bigl(TV_{it} + PB_{it} + BR_{it}\bigr),
$$

with the biomarker $B_i$, baseline and the three response
components drawn jointly from a multivariate normal distribution.
The parameters used for the published results were (inspected,
`Produce_Publication_Results_1_generate_data.Rmd`, and the saved
`parameterselections`):

| Parameter | Published values |
|---|---|
| Designs | OL, OL+BDC, CO, Hybrid ("N-of-1"); 8 visits over 20 weeks |
| $N$ | 35, 70 |
| $c_{bm}$ | 0, 0.3, 0.6 |
| Carryover half-life $t_{1/2}$ | 0, 0.1, 0.2 weeks |
| Within-factor correlation $\rho_{cs}$ | 0.8 for all three components, constant at every lag (CS) |
| Cross-factor correlation | $c_1 = 0.2$ same visit, $c_\times = 0.1$ otherwise |
| Biomarker-BR correlation | $c_{bm}$ wherever the carryover-adjusted BR mean is nonzero, 0 elsewhere |
| Carryover on BR mean | recursive, $\mu_{BR,t-1}(1/2)^{t_{sd}/t_{1/2}}$ with cumulative $t_{sd}$ |
| Positive definiteness | `makePositiveDefinite = TRUE` in every call: `corpcor::make.positive.definite(sigma, tol = 1e-3)` when the check fails |
| Analysis | `lmer`, random intercept, binary on-drug indicator, test of `bm:Db` |

Table: Parameters of the published Hendrickson et al. simulation framework (docs/36, Section 3)

The Hybrid schedule places visits at weeks 4, 8, 9, 10, 11, 12, 16
and 20; all four of its paths begin on drug, as do both OL+BDC
paths. CO has one path on drug first and one placebo first.

### The orig feature table

<!-- From docs/44 Section 2 "The starting point: orig" -->

The published process has seven structural features:

| Feature | orig |
|---|---|
| Joint distribution | one 26-variable multivariate normal per path |
| Within-component correlation | CS, $\rho_{cs} = 0.8$ at every lag |
| Cross-component correlation | $c_1 = 0.2$ same visit, $c_\times = 0.1$ otherwise |
| Biomarker coupling | step rule: $c_{bm}$ wherever the $BR$ mean is nonzero |
| Carryover on the $BR$ mean | recursive, cumulative $t_{sd}$, over-decays |
| Re-exposure | time on drug restarts at zero |
| Off-drug $BR$ variance | full (SD 8) at every visit |
| Invalid matrices | silently repaired (`make.positive.definite`) |

Table: Features of the published process (orig) (docs/44, Section 2)

Its two known defects share one cause, the step rule. Without
carryover the coupling has an on/off contrast, and the matrix is
invalid for any $c_{bm}$ above 0.256. With any carryover the
coupling is constant, the matrix is valid, and the interaction
disappears (docs/36, Section 5.4; docs/40, Section 5).

### The three response components and the Gompertz parameters

<!-- From docs/40 Section 3 "The three response components" and
     docs/25 Section 4 "How Hendrickson et al. (2020) use it" -->

Each component receives its own modified Gompertz curve
$G(t; m, d, r)$, with each component's own $t$ clock. The
published calibration, from the package dataset `extracted_rp`,
is (inspected):

| Component | maxr | disp | rate (per week) | sd |
|---|---|---|---|---|
| BR | 10.99 | 5 | 0.42 | 8 |
| TV | 6.51 | 5 | 0.35 | 10 |
| PB | 6.51 | 5 | 0.35 | 10 |

Table: Published (Hendrickson) modified Gompertz parameters for the
BR, TV, and PB components, from the package dataset `extracted_rp`

The TV and PB components share the same parameters -- the "tabula
rasa" assumption, because no data separately identify the
expectancy response from the time-varying natural course
(docs/25, Section 4). The BR trajectory was fitted to the
difference between the prazosin and placebo arms of the
active-duty trial; the placebo arm's trajectory was taken to be
$TV + PB$.

The inflection points are near 3.8 weeks for BR ($\ln(5)/0.42$)
and 4.6 weeks for TV and PB. At week 12 the BR curve has reached
about 97 percent of its asymptote (docs/40, Section 4). Baseline
has mean 83.07 and SD 18.48 CAPS points; the biomarker (standing
systolic blood pressure) has mean 124.33 and SD 15.36 mmHg.

Expected component means by week and Hybrid path, without
carryover (inspected, docs/40 Table 2):

| Week | 4 | 8 | 9 | 10 | 11 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|
| $TV$, all | 1.86 | 4.79 | 5.24 | 5.59 | 5.85 | 6.03 | 6.39 | 6.48 |
| $PB$, all | 1.86 | 4.79 | 2.62 | 2.79 | 2.92 | 3.02 | 3.19 | 3.24 |
| $BR$, path 1 | 4.28 | 9.22 | 9.79 | 10.19 | 0 | 0 | 4.28 | 0 |
| $BR$, path 2 | 4.28 | 9.22 | 9.79 | 10.19 | 0 | 0 | 0 | 4.28 |
| $BR$, path 3 | 4.28 | 9.22 | 9.79 | 0 | 0 | 0 | 4.28 | 0 |
| $BR$, path 4 | 4.28 | 9.22 | 9.79 | 0 | 0 | 0 | 0 | 4.28 |

Table: Expected component means by week and Hybrid path without
carryover (docs/40, Section 3)

### The biomarker coupling rule (step rule)

<!-- From docs/55 Sections 3 to 4 "The intended model" and "The
     code as published" -->

The article describes the covariance matrix in one sentence:

> The covariance matrix was assembled based on a set of modifiable
> parameters defining the correlation between the baseline biomarker
> and the BR components, the autocorrelation over time (relating
> the value of a factor for one participant at one timepoint to the
> value of that factor for that participant at subsequent
> timepoints), the correlation at a single time point between the 3
> factors, and the variance of each component.

The parameter `c.bm` is documented as "Correlation between the
biomarker and the biologic response" (`R/generateData.R`, line 11;
quoted). The article does not state what the off-drug coupling
should be. The project's reading is that the coupling should fade
with the drug response: graded coupling, $g_t = 1$ on drug,
$g_t = 2^{-t_{sd}/t_{1/2}}$ off drug after first exposure, and
$g_t = 0$ before first exposure.

The code at `58b32a9`, lines 131 to 138, takes a different route
(quoted):

```r
131      # correlation with biomarker
132      for(p in 1:nP){
133        n1<-paste(trialdesign$timeptname[p],"br",sep=".")
134        if(means[which(n1==labels)]!=0){
135          correlations[n1,'bm']<-modelparam$c.bm
136          correlations['bm',n1]<-modelparam$c.bm
137        }
138      }
```

The test on line 134 is whether the expected $BR$ mean at the
visit is nonzero. The expected mean used in the test already
includes carryover. For any $t_{1/2} > 0$ the carryover-adjusted
mean is positive at every off-drug visit after the first exposure,
so those visits receive the full on-drug correlation (inspected;
verified in all 40 path and half-life cases). This is the step
rule: $g_t = \mathbf{1}[\mu_{BR,t} \neq 0]$ (docs/55, Section 4).

### The carryover mean computation (the recursion)

<!-- From docs/53 Sections 3 to 4 "The intended model" and "The
     code as published" -->

The article states that "the value of the BR at the last timepoint
on active drug was exponentially decayed, with the half-life of
this decay being maintained as a model parameter" ([1], p. 6,
quoted). The intended model is therefore

$$
\mu_{BR,t} = \mu_{\text{last}}\; 2^{-t_{sd,t}/t_{1/2}},
$$

where $\mu_{\text{last}}$ is the mean at the last on-drug visit and
$t_{sd,t}$ is the time since discontinuation at visit $t$. This is
the anchored form.

The code at `58b32a9`, lines 83 to 92, implements a different
computation (quoted):

```r
 83      brmeans<-modgompertz(d$tod,rp$max,rp$disp,rp$rate)
 84      if(nP>1){
 85        for(p in 2:nP){
 86          if(!d[p]$onDrug){
 87            if(d[p]$tsd>0){
 88              brmeans[p]<-brmeans[p]+brmeans[p-1]*
                 (1/2)^(d$tsd[p]/modelparam$carryover_t1half)
 89            }
 90          }
 91        }
 92      }
```

Line 88 multiplies the **previous visit's expected response** by
$2^{-t_{sd}/t_{1/2}}$, where $t_{sd}$ is the **cumulative** time
since discontinuation. Because the previous visit's value has
already been decayed, elapsed time is counted again at every visit.
At the $m$-th off-drug visit the mean becomes

$$
\mu_{\text{last}}\;2^{-(t_{sd,1} + t_{sd,2} + \cdots +
t_{sd,m})/t_{1/2}},
$$

rather than $\mu_{\text{last}}\;2^{-t_{sd,m}/t_{1/2}}$ (derived;
docs/53, Section 1). The decay compounds from the second off-drug
visit on, running faster than the stated half-life. With equally
spaced visits the effective half-life at the $m$-th visit is
$2t_{1/2}/(m+1)$: two-thirds of the stated value at the second
visit, one-half at the third.

The code and the intended formula agree at every on-drug visit and
at the first off-drug visit of every run. They differ at the second
and later visits of an off-drug run.

The one-line fix is to decay by the **increment** in $t_{sd}$
since the previous visit:

```r
brmeans[p]<-brmeans[p]+brmeans[p-1]*
  (1/2)^((d$tsd[p]-d$tsd[p-1])/modelparam$carryover_t1half)
```

(docs/53, Section 1; verified to reproduce the intended means to
within $1.1 \times 10^{-14}$ in all 40 path and half-life cases.)

### Positive-definiteness failures and the repair

<!-- From docs/legacy/23 Section 3.1 "The Problem" and docs/36
     Section 4.1 -->

The DGP builds a 26 by 26 covariance matrix per trial path from
user-specified correlations and standard deviations. When the
result is not positive definite, `generateData()` projects to the
nearest PD matrix via `corpcor::make.positive.definite(sigma,
tol = 1e-3)`. In the publication code the correction was silent.

Failure rates by commit (inspected, docs/legacy/23 Section 3.1):

| Commit | Description | Failures | Rate |
|---|---|---|---|
| 42ac030 | Initial publication | 40 of 162 | 24.7 pct |
| 8609f12 | Ron Thomas edits | 102 of 162 | 63.0 pct |
| f6ee86d | Autocorr typo fix | 102 of 162 | 63.0 pct |
| f70f86d | Carryover added | 102 of 162 | 63.0 pct |

Table: Positive definiteness failure rates of covariance matrices
by commit (docs/legacy/23, Section 3.1)

All 162 matrices in the standard parameter sweep were
ill-conditioned ($\mathrm{cond}(\Sigma) > 100$) under compound
symmetry. Minimum eigenvalues of failed matrices ranged from
-33.26 (worst) to -0.52 (least severe), on a dominant eigenvalue
of approximately 700.

From docs/36, failure of the published check by design (verified):

| Design | Matrices checked | Failures | Where |
|---|---|---|---|
| OL | 18 | 0 | none |
| OL+BDC | 36 | 8 | both paths, $c_{bm} > 0$, $t_{1/2} = 0$ |
| CO | 36 | 16 | both paths at $t_{1/2} = 0$; the placebo-first path at every $t_{1/2}$ ($c_{bm} > 0$) |
| Hybrid | 72 | 16 | all four paths, $c_{bm} > 0$, $t_{1/2} = 0$ |
| Total | 162 | 40 | |

Table: Failures of the published positive definiteness check by
design (docs/36, Section 4.1)

Every matrix with $c_{bm} = 0$ passes. The response block (all
$TV$, $PB$ and $BR$ entries) is positive definite in every case.
Failure always comes from the biomarker row: the ceiling is 0.256
in the Hybrid and CO designs without carryover, and both published
effect sizes (0.3 and 0.6) exceed it.

---

## Derivations

<!-- one row per result: statement | CM label | verification
script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the
numbers agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at
docs/series/registry/open-items.csv -->

---

## Reproduce

<!-- one command per output -->

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Initial content moves from docs/36, docs/40, docs/44, docs/25, docs/53, docs/55, docs/legacy/23 | Phase 2 content move |
