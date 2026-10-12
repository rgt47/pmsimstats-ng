# pmsimstats-ng Documentation Compendium {.unlisted .unnumbered}
*2026-04-21 06:35 PDT*

*Archived 2026-10-11; current treatment: Section 3.1 in `series/ch03`, Sections 8.1 to 8.2 in `series/ch07`, Section 14 in Appendix D, Section 17 in Appendix G; see `docs/legacy/INDEX.md`.*

**Author.** pmsimstats team

**Purpose.** This compendium consolidates the novel technical content from all
source documents in `~/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/`.
It unifies mathematical derivations, empirical results, audit findings,
and implementation guidance that is currently scattered across twenty-six
numbered .md and .tex files. Each section cites the source documents that
contributed to it. Sections or claims that have been superseded by
subsequent development are flagged explicitly.

The documentation corpus was last substantially updated on 2026-04-16
(audit report) and 2026-04-15 (parity workflow, revision response to
reviewer on document 02). The current repository state reflects the
2026-04-08 merge that reorganized code into `implementations/`, added
Architecture A (mean moderation) alongside Architecture B (multivariate
normal differential correlation), and moved the pre-merge tidyverse
pipeline (`analysis/2025/`) to an external archive.

**Review note (2026-10-07).** This compendium records the state of
the project in April 2026, and several of its statements have since
been corrected or superseded. Outright errors (wrong numbers, wrong
mathematics, wrong references, notation) were fixed in place on
2026-10-07. Superseded statements are kept as written, each with a
short dated note pointing to the newer document. The main ones:
the positive-definiteness ceilings on $c_{bm}$ (Sections 2.1, 2.7
and 3.3; see `docs/34`, `docs/35` and `docs/36`); the account of the
published carryover collapse (Sections 1.5 and 5.7; see `docs/36`
and `docs/37`); the headline 40-60 percent versus 10-15 percent
power-loss contrast (Sections 1.4, 16.1 and 18; see paper 01 and
`docs/37`); and the recommended analysis model (Section 4; see
`docs/45`, `docs/46` and `docs/49`). The data-generating features
that a complete N-of-1 DGP needs are reviewed in `docs/51`.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

The notation follows the canonical table in
`analysis/report/NOTATION.md`. A symbol marked *canonical* has the
meaning given there. A symbol marked *local* is used only in this
compendium and was chosen so that it does not collide with a
canonical symbol. Code identifiers (`c.bm`, `Dbc`, `Sx`, `tsd`,
`tod`, `lambda_cor`, `carryover_t1half`) appear only in code spans.
The package passes the moderation strength as `c.bm` under both
architectures. In the mathematics it is written $c_{bm}$ under
covariance moderation and $\beta_{bm}$ under mean moderation. Where a
table reports both architectures at one value, the shared label is
$c_{bm} = \beta_{bm}$ (NOTATION rule 2). $N$ is the total sample size
across randomization paths. However, the April 2026 runs quoted here
predate the 2026-08 driver audit, so whether a quoted $N$ is a total or
a per-path count has not been re-verified (NOTATION rule 4).

| Symbol | Meaning | Status |
|---|---|---|
| $i$, $t$ | participant and timepoint (weeks) | canonical |
| $N$ | total number of participants | canonical |
| $Y_{it}$ | observed symptom score (code `Sx`) | canonical |
| $BR_{it}$, $PB_{it}$, $TV_{it}$ | biologic-response, placebo-belief and time-variant components | canonical |
| $\varepsilon_{it}$, $u_i$ | residual and participant random intercept | canonical |
| $D_{it}$ | binary drug state | canonical |
| $D_{bc,it}$ | exposure-decayed drug indicator (code `Dbc`) | canonical |
| $t_{sd}$, $t_{od}$ | time since discontinuation, time on drug (code `tsd`, `tod`) | canonical |
| $t_{1/2}$ | carryover half-life, weeks (code `carryover_t1half`) | canonical |
| $\lambda$ | decay rate $\ln 2 / t_{1/2}$ (code `lambda_cor`) | canonical |
| $B_i$ | raw biomarker (code `bm`) | canonical |
| $b_i$ | standardized biomarker $(B_i - \mu_{bm})/\sigma_{bm}$ (code `bm_z`) | canonical |
| $\mu_{bm}$, $\sigma_{bm}$ | biomarker mean and SD (124.33, 15.36) | canonical ($\sigma_{bm}$); local ($\mu_{bm}$) |
| $\sigma_{BR}$ | SD of $BR$ (8) | canonical |
| $\sigma_{PB,t}$ | SD of $PB$ at visit $t$, scaled by expectancy | local |
| $\eta_t$ | expectancy at visit $t$: 1 open label, 0.5 blinded, 0 baseline (code `e`) | canonical |
| $c_{bm}$ | covariance-moderation strength: on-drug correlation of $B_i$ with $BR_{it}$ | canonical |
| $\beta_{bm}$ | mean-moderation multiplier | canonical |
| $\beta_{bm:D}$ | biomarker-by-treatment interaction coefficient (`bm:Dbc`) | canonical |
| $\beta_D$ | treatment main effect | canonical |
| $\beta_0, \ldots, \beta_3$ | other fixed effects of the model being described | local |
| $\rho$ | AR(1) within-factor correlation per week (nlme calls the `corCAR1` parameter Phi) | canonical |
| $\phi_t$ | decay factor $2^{-t_{sd}/t_{1/2}}$ off drug, 0 before first exposure | local (as in `docs/37`) |
| $g_t$ | coupling: the fraction of $c_{bm}$ present in $\mathrm{Cor}(B_i, BR_{it})$ | local (as in `docs/37`) |
| $r_t$ | biomarker-response correlation at visit $t$, $c_{bm} g_t$ | local (as in `docs/37`) |
| $h_t$ | moderation profile of a mean-moderation shift (Section 8.1, Attempt 2) | local (as in `docs/37`) |
| $s$ | carryover scale factor in the analysis indicator (code `carryover_scalefactor`, default 1) | local |
| $\sigma_{BR,t}$ | a hypothetical time-varying SD of $BR$ (Section 5.8) | local |
| $\bar Y_{i,\text{on}}$, $\bar Y_{i,\text{off}}$, $\Delta_i$, $\sigma_\Delta$ | participant phase means, their difference, and its SD (Sections 10, 13) | local |
| $T$, $K$, $L$ | visits per participant, number and length of drug blocks (Section 13) | local |
| $\mu_{BR,t}$ | mean of $BR$ at visit $t$ | local |
| $\Delta_{\text{on}}$, $\Delta_{\text{off}}$ | spread of expected $BR$ between participants at $\pm 1$ SD of the biomarker | local |
| $\sigma_u^2$, $\sigma_\varepsilon^2$, $\sigma_T^2$ | random-intercept, residual and total within-participant variance (Section 4.3) | local |
| $\omega_0$, $\omega_1$ | intercept and biomarker slope of a change-score regression (Sections 10, 13) | local |
| $\xi_0$, $\xi_1$, $e_i$ | intercept, slope and error of the second-stage regression (Section 10) | local |
| $\pi(B_i)$ | latent-class membership probability (Section 11; paper 03 usage) | local |
| $Z_i$ | latent class label | canonical |
| $\mathrm{cond}(\Sigma)$ | condition number of a covariance matrix | local |
| $\alpha$ | nominal significance level | canonical |

Table: Notation index: symbols used in this compendium

**Glossary.**

- **Architecture A / Architecture B.** Former labels for the
  mean-moderation architecture (`dgp_architecture = 'mean_moderation'`)
  and the covariance-moderation architecture (`'mvn'`). The labels are
  kept in this historical document; NOTATION.md now names the
  architectures by mechanism.
- **N-of-1 (design).** The label these April 2026 sources use for the
  Hendrickson (2020) Hybrid design (OL, then BDC, then a brief
  crossover). It is called Hybrid elsewhere in the project.
- **OL, OL+BDC, CO.** Open-label, open-label followed by blinded
  discontinuation, and crossover designs.
- **ER.** An older name for the expectancy (placebo-belief) component,
  now $PB$; it has been renamed $PB$ throughout this document.
- **Step, proportional and graded coupling.** Rules for the off-drug
  biomarker-response correlation (Section 5). The current package uses
  graded coupling, $g_t = \phi_t$ off drug and $g_t = 0$ before first
  exposure (inspected, `R/generateData.R`, 2026-10-07).
- **PD.** Positive definite (of a covariance matrix); the *ceiling*
  is the largest $c_{bm}$ for which the matrix is PD.
- **E1/E2/E3.** The project's current labels for the three carryover
  analysis specifications (NOTATION Part 3). Section 12 predates them.

## 0. Source Document Inventory and Antiquation Map

The source files under `docs/` are listed below with their role and
their standing as of 2026-04-21. Files marked 'Current' are faithful
to the code. Files marked 'Historical' document decisions or states
that have been implemented; their technical content is preserved here
because the reasoning remains instructive. Files marked 'Superseded'
describe methodology that the revised code no longer uses; readers
should consult this compendium or the currently-labeled source rather
than those files.

| File | Role | Standing |
|---|---|---|
| 00-documentation-index.tex | Navigation index | Current |
| 01-codebase-overview.tex | Package structure and parameter table | Current |
| 02-revision-response-latent-class-expanded.md | Psychometric connection | Current |
| 02-revision-response-to-reviewer.md | Editorial revision log | Current |
| 03-audit-and-revision-report.tex | Seven DGP corrections | Current (foundational) |
| 04-revised-power-analysis.tex | Validated power estimates | Current |
| 06-ar1-residual-correlation.tex | nlme plus corCAR1 transition | Current |
| 07-positive-definiteness-failures.tex | PD failure analysis | Current |
| 08-biomarker-correlation-decay.tex | lambda_cor derivation | Current |
| 09-carryover-correlation-artifact.tex | Step-function artifact | Current (diagnostic) |
| 10-carryover-analysis-model-assessment.tex | 2x2 factorial plan | Current (proposed study) |
| 13-figure4-walkthrough.tex | Publication methodology | Superseded |
| 14-trial-analysis-implementation-guide.tex | Collaborator deliverables | Current |
| 15-tidyverse-alignment-plan.tex | Pre-merge port plan | Historical |
| 16-tidyverse-code-review.tex | Pre-merge code review | Historical |
| 17-simplification-plan.tex | Alternative simplification | Historical (not adopted) |
| 18-response-parameter-sensitivity.tex | Figure 5 validation | Current |
| 19-mean-moderation-implementation-notes.tex | Architecture A derivation | Current |
| 19-mean-moderation-implementation-notes.md | Quick-reference summary | Current |
| 20-annotated-bibliography.tex | Literature synthesis | Current reference |
| 21-architecture-comparison-summary.md | Empirical comparison | Current |
| 22-journal-target-recommendation.md | Submission strategy | Current |
| 22-tidyverse-development-parity-workflow.md | Parity workflow | Current |
| audit-2026-04-16.md | Technical audit | Current (most recent) |
| drafting-02-rgt/working.tex | Working draft of doc 02 | Working file (no material differences) |
| 07a-pd-constraints-theory.tex | PD constraints derivation (Sylvester, Gershgorin, correlation hierarchy) | Current (migrated 2026-04-23 from pmsimstats2025) |
| 09a-carryover-bias-variance.md | Non-monotonic power under model misspecification (bias-variance tradeoff) | Current (migrated 2026-04-23 from pmsimstats2025) |
| 24-component-decomposition-pedagogy.md | BR/ER/TV conceptual guide with prazosin-PTSD narrative | Current (migrated 2026-04-23 from pmsimstats2025) |
| archive/correlation-structure-alternatives.md | Alternative correlation structures (Matern, rational-quadratic, power-law) | Historical (migrated 2026-04-23; exploratory) |
| archive/hendrickson-pseudocode-comparison.md | Side-by-side pseudocode audit trail vs. Hendrickson original | Historical (migrated 2026-04-23) |

Table: Source documents with their role and current standing

*Note (2026-10-07).* The standings above are as of 2026-04-21. Several
listed files are no longer in `docs/`: the two
`02-revision-response-*.md` files,
`10-carryover-analysis-model-assessment.tex`, `audit-2026-04-16.md`
and `drafting-02-rgt/`. Documents `docs/24` to `docs/51` postdate this
compendium. Of these, `docs/34` to `docs/37`, `docs/45`, `docs/46`,
`docs/49` and `docs/51` supersede parts of Sections 1 to 5 (see the
dated notes below).

Four key antiquation points to highlight:

1. **Doc 13 (Figure 4 Walkthrough)** describes the publication
   methodology (compound symmetry, `lme4::lmer`, step-function BM-BR
   correlation, carryover half-lives of 0.1 and 0.2 weeks). The revised
   code no longer uses any of these elements. The document is retained
   as a record of the baseline from which the seven DGP corrections
   (Section 2) were derived.

2. **Doc 15 (Tidyverse Alignment Plan)** describes a pre-merge roadmap
   for porting `analysis/2025/pm_functions.R` to parity with the
   original data.table code. The merge on 2026-04-08 executed most of
   this plan, moved the pre-merge tidyverse file to an external
   archive, and reorganized the code into `implementations/`. Function
   names and paths referenced in doc 15 no longer resolve.

3. **Doc 16 (Tidyverse Code Review)** reviews the same pre-merge
   file. Its findings informed the 2026-04-08 rewrite. Line numbers
   and function names in doc 16 refer to the archived file, not the
   current `implementations/tidyverse/R/functions.R`.

4. **Doc 17 (Simplification Plan)** proposed a substantially simpler
   framework (random intercept with independent errors, ANCOVA on
   phase-mean change scores, reduced parameter space). The merge took
   a different path: it added Architecture A rather than simplifying.
   The document is retained as a record of considered alternatives
   and contains the most thorough formal specification of six
   alternative analysis methods in the project.

One execution-blocking issue was carried forward at the time of writing
(audit-2026-04-16, finding `A1`):
`analysis/carryover_factorial/02-run-factorial-2025.R:20` sources
`analysis/2025/01-pm-functions.R`, which has been moved to an external
archive. The script fails at load time until the source target is
updated to `implementations/tidyverse/R/functions.R` or a copy from the
archive is restored.

---

## 1. The Dual Data-Generating Process Architectures

Source documents: 02-dgp-mean-moderation-vs-mvn.tex;
19-mean-moderation-implementation-notes.tex; 21-architecture-comparison-summary.md;
02-revision-response-latent-class-expanded.md.

### 1.1 Architecture A: Direct Mean Moderation

The biomarker enters the mean structure of the outcome as a
multiplicative scaling of the treatment effect:

$$Y_{it} = \beta_0 + \beta_D \cdot D_{it} \cdot (1 + \beta_{bm} \cdot b_i)
  + \varepsilon_{it}$$

Expanding the product gives the additive interaction model

$$Y_{it} = \beta_0 + \beta_D D_{it} + \beta_2 b_i + \beta_{bm:D} D_{it} b_i
  + \varepsilon_{it}$$

with $\beta_2 = 0$ and $\beta_{bm:D} = \beta_D \beta_{bm}$. The
multiplicative form is therefore the special case of the additive
model that has no biomarker main effect. The identity is algebraic and
does not depend on centering. Centering ($b_i$ rather than $B_i$) only
makes $\beta_D$ the treatment effect at the mean biomarker. (Corrected
2026-10-07: the earlier text said that the two forms are equivalent
under centering.)

Properties:

- The interaction is deterministic given biomarker and treatment
  status. This is a population-level mean-structure claim, not a
  statement about individual outcomes.
- $\beta_{bm}$ is a population parameter, uniform across participants.
- The signal lives in the first moment of the outcome distribution.
- Under carryover, the ratio of drug effect between high- and
  low-biomarker participants is preserved. Power loss is modest.

Implementation in the pmsimstats-ng framework: after an MVN draw from
a covariance matrix that omits the BM-BR block, each participant's BR
factor at on-drug timepoints is shifted by a standardized biomarker
deviation:

$$BR_{it}^{\text{mod}} = BR_{it}^{\text{base}}
  + \beta_{bm}\, b_i\, \mathbf 1[D_{it} = 1]\, \sigma_{BR},
  \qquad b_i = (B_i - \mu_{bm})/\sigma_{bm}$$

The shift per one standard deviation of biomarker is
$\beta_{bm}\sigma_{BR}$. For $\beta_{bm} = 0.45$ and $\sigma_{BR} = 8$
(the value in `extracted_rp`, verified) this is 3.6 BR units, which
matches the conditional expectation under Architecture B at
$c_{bm} = 0.45$. In code the multiplier is passed as `c.bm`. *Note
(2026-10-07):* `generateData()` also offers
`moderation_scaling = 'trajectory'`. Under that option the shift follows
the on-drug BR mean curve, normalized so that its mean on-drug magnitude
is $\sigma_{BR}$ (inspected, `R/generateData.R`).

### 1.2 Architecture B: Multivariate Normal Differential Correlation

The biomarker and the biologic response factor are jointly multivariate
normal, with a treatment-state-dependent correlation:

$$\begin{pmatrix} B_i \\ BR_{it} \end{pmatrix} \sim
  \mathrm{MVN}\left( \begin{pmatrix} \mu_{bm} \\ \mu_{BR}(t, D_{it}) \end{pmatrix},
  \Sigma(D_{it}) \right)$$

$$\mathrm{Cor}(B_i, BR_{it}) = \begin{cases}
  c_{bm} & \text{if } D_{it} = 1 \\
  c_{bm} \cdot \exp(-\lambda \cdot t_{sd}) & \text{if } D_{it} = 0,\ t_{sd} > 0 \\
  0 & \text{before first exposure}
\end{cases}$$

The zero before first exposure (for example, the placebo-first CO path)
is the current package behavior (inspected, `buildSigma()`); the
earlier text omitted that case.

The conditional expectation on drug is:

$$E[BR_{it} \mid B_i, D_{it} = 1] = \mu_{BR}(t)
  + c_{bm} \cdot \frac{\sigma_{BR}}{\sigma_{bm}} \cdot (B_i - \mu_{bm})
  = \mu_{BR}(t) + c_{bm}\,\sigma_{BR}\, b_i$$

Off drug, with carryover half-life $t_{1/2}$ and
$\lambda = \ln 2 / t_{1/2}$:

$$E[BR_{it} \mid B_i, D_{it}=0, t_{sd}] = \mu_{BR,\text{off}}(t_{sd})
  + c_{bm} \cdot e^{-\lambda t_{sd}} \cdot \sigma_{BR}\, b_i$$

Properties:

- The interaction is probabilistic: a tendency, not a deterministic
  scaling.
- $c_{bm}$ is a correlation coefficient with $\lvert c_{bm}\rvert \le 1$.
- The signal lives in the second moment (the covariance structure).
- Carryover erodes the differential correlation. Power loss is
  substantial.

### 1.3 Mathematical Comparison

Architecture A yields an expected difference that scales cleanly with
drug exposure:

$$E[Y_{it} \mid b_i, D_{it}] - E[Y_{jt} \mid b_j, D_{jt}]
  = \beta_D \cdot \beta_{bm} \cdot (b_i - b_j) \cdot D$$

for two participants $i$ and $j$ in the same drug state
$D_{it} = D_{jt} = D$. The ratio of effects between high- and
low-biomarker participants is preserved as $D$ varies. If $D$ is
replaced by a carryover-decayed continuous indicator, the difference
scales by the same factor, but the high/low ratio is unchanged.

Architecture B yields an on-drug difference:

$$E[BR_{it} \mid b_i, D=1] - E[BR_{jt} \mid b_j, D=1]
  = c_{bm} \cdot \sigma_{BR} \cdot (b_i - b_j)$$

Off drug, the difference is multiplied by $e^{-\lambda t_{sd}}$. The
biomarker-specific signal erodes exponentially during washout.

### 1.4 Empirical Power Comparison

Setup (doc 21): $N = 70$, $c_{bm} = \beta_{bm} = 0.45$, `nlme::lme`
plus `corCAR1` with the `bm:Dbc` interaction, 50 replicates. Half-lives
in weeks.

Architecture B (MVN differential correlation):

| Design | $t_{1/2} = 0$ | $t_{1/2} = 0.5$ | $t_{1/2} = 1.0$ |
|---|---|---|---|
| N-of-1 | 0.82 | 0.64 | 0.50 |
| CO | 0.40 | 0.40 | 0.32 |
| OL+BDC | 0.62 | 0.38 | 0.24 |

Table: Architecture B power by design and carryover half-life, 50 replicates

Relative power loss at $t_{1/2} = 1.0$: N-of-1 39 percent, CO 20 percent,
OL+BDC 61 percent.

Architecture A (mean moderation):

| Design | $t_{1/2} = 0$ | $t_{1/2} = 0.5$ | $t_{1/2} = 1.0$ |
|---|---|---|---|
| N-of-1 | 0.74 | 0.72 | 0.68 |
| CO | 0.38 | 0.36 | 0.34 |
| OL+BDC | 0.62 | 0.58 | 0.54 |

Table: Architecture A power by design and carryover half-life, 50 replicates

Relative power loss at $t_{1/2} = 1.0$: N-of-1 8 percent, CO 11 percent,
OL+BDC 13 percent.

Monte Carlo standard error at power near 0.5 with 50 replicates is
approximately $\sqrt{0.5 \cdot 0.5 / 50} = 0.07$.

*Note (2026-10-07).* `docs/19` (both the `.md` and `.tex` versions)
reports the N-of-1 Architecture A row as 0.74 / 0.68 / 0.72, and doc 21
reports it as 0.74 / 0.72 / 0.68. The middle and last cells are
swapped, and the source run was not located, so the discrepancy is
unresolved. Either version lies within Monte Carlo error of flat power.
These are 50-replicate pilot estimates. With per-cell MCSE of about
0.07, the relative losses above are imprecise to about $\pm 15$
percentage points, and they should not be cited as the project's
result. Paper 01 later re-ran the comparison at a matched total of
$N = 70$. In that run Hybrid carryover losses were 6.9 percent under
mean moderation and 15.4 percent under covariance moderation (as
recorded in `analysis/report/NOTATION.md`, Part 2; inspected, not
recomputed). `docs/37` gives the closed-form dependence of power on
the half-life.

### 1.5 Why the Architectures Diverge

Two mechanisms operate simultaneously under Architecture B:

1. **Mean blurring** (shared with Architecture A). Off-drug BR means
   are inflated by residual carryover, compressing the on-drug versus
   off-drug treatment contrast. Both architectures experience this.

2. **Correlation erosion** (unique to Architecture B). The BM-BR
   correlation decays exponentially with $t_{sd}$. Because the
   interaction signal lives in the second moment, this directly
   weakens the signal. This mechanism has no analogue in Architecture
   A, where the interaction lives in the first moment.

Architecture A's interaction coefficient sees reduced `Dbc` contrast
under carryover, but the signal-to-noise ratio degrades only modestly
because the biomarker-by-drug signal scales proportionally with the
reduced contrast. Architecture B's interaction coefficient is attacked
from both sides: the treatment contrast in `Dbc` is compressed, and the
biomarker-BR correlation generating the covariance signal is
simultaneously decaying.

*Note (2026-10-07).* `docs/37` (derived and checked by simulation)
qualifies this two-mechanism account. For the paired-difference
statistic, carryover in the response *means* is absorbed by the
contrast's intercept. The whole effect of the half-life on power then
runs through the off-drug coupling $\bar g_{\text{off}}$, that is,
through correlation erosion alone. Under mean moderation with an
on-drug-only shift, power does not depend on the half-life
(`docs/37`, Section 1, points 1 and 4). Mean blurring therefore
explains little of the architecture contrast.

### 1.6 Biological Rationale and Architecture Selection

The choice between architectures is a claim about biology, not
statistical convenience.

Architecture A is appropriate when the biomarker governs the
**magnitude of the drug's biological effect** for each individual.
Examples:

- Renal clearance (eGFR determines drug elimination).
- Receptor density moderation (target availability determines
  pharmacodynamics).
- CYP2D6 metabolizer status.
- Dose-response with biomarker-dependent effective dose.

Architecture B is appropriate when the biomarker predicts which
participants are likely to respond, through shared variance with a
drug-responsive physiological component that decays as the drug is
eliminated. Examples:

- Blood pressure as PTSD subtype marker (elevated resting SBP indexes
  heightened central noradrenergic tone, a subtype preferentially
  responsive to prazosin as alpha-adrenergic antagonist).
- Baseline disease severity as treatment predictor.
- Inflammatory biomarkers in psychiatry.
- Genetic risk scores.

Decision table (doc 02):

| Criterion | Architecture A | Architecture B |
|---|---|---|
| Biomarker role | Causal mediator | Statistical predictor |
| Mechanism | Deterministic scaling | Probabilistic association |
| Signal location | Mean structure | Covariance structure |
| Carryover impact | Modest (10-15 pct) | Substantial (40-60 pct) |
| Appropriate when | Biomarker determines effective dose or PK | Biomarker indexes latent subtype |

Table: Decision criteria for choosing between Architecture A and Architecture B

The carryover-impact row restates the 50-replicate pilot of Section
1.4. The note there gives the later matched-$N$ figures.

Reporting recommendation: simulation studies should explicitly state
whether the biomarker-treatment interaction is mean moderation or
differential correlation, provide the biological rationale, confirm
that the carryover model operates consistently with the chosen
architecture, and include sensitivity analyses under the alternative
architecture when the biological mechanism is uncertain. When
uncertain, Architecture B produces more conservative (lower) power
estimates under carryover and is therefore the safer default for
trial design.

### 1.7 Prevalence of Each Architecture in the Literature

Architecture A is the near-universal standard in clinical trial
simulation, appearing in biomarker-stratified designs (Simon 2010),
adaptive enrichment trials (Freidlin and Korn 2014), basket and
umbrella trials (Renfro and Sargent 2017), platform trials, and
biomarker-by-treatment interaction studies (Haller and Ulm 2018).
It is standard in N-of-1 simulation (Zucker et al. 1997; Araujo et
al. 2016; Duan et al. 2013) and in N-of-1 methodology on power and
design (Wang and Schork 2019; Schork 2022).

Architecture B (MVN differential correlation) appears primarily in
latent variable and structural equation modeling, Bayesian joint
modeling, and copula-based simulation. The Hendrickson et al. (2020)
framework appears to be unique in the N-of-1 trial literature in
using differential correlation as the interaction mechanism.

Implication: published power estimates for biomarker-moderated
crossover and N-of-1 designs may be systematically optimistic for
biomarkers that operate through an Architecture B mechanism, because
those studies do not account for the additional power loss that
carryover produces when the interaction signal resides in the
covariance structure rather than the mean structure.

---

## 2. The Seven DGP Corrections

Source document: 03-audit-and-revision-report.tex (dated 2026-03-24).

The revised codebase differs from the Hendrickson et al. (2020)
publication code in seven substantive ways. Each is documented below
with the mathematical form, the reason for the change, and the effect
on inference.

### 2.1 Correction 1: AR(1) Autocorrelation

Replaces compound symmetry within each response factor.

$$\mathrm{Cor}(X_{t_i}, X_{t_j}) = \rho^{\lvert t_i - t_j\rvert}$$

Under compound symmetry the correlation is constant regardless of
time gap. Under AR(1) it decays exponentially with temporal
separation, which is clinically more realistic. The change widens
the PD-feasible range for the BM-BR correlation. *Corrected
2026-10-07* (the earlier text mixed two values of $\rho$ in one
comparison and overstated the AR(1) ceiling). Under the published
compound-symmetry specification ($\rho = 0.8$, no carryover) the
ceiling is 0.256 (`docs/36`, verified there). At the same
$\rho = 0.8$, the separable AR(1) structure raises the Hybrid ceiling
to 0.440 (`docs/36`). At the package value $\rho = 0.7$ and no
carryover, the separable AR(1) ceiling is 0.474 in OL+BDC and 0.480 in
Hybrid (`docs/35`), and the current `covar` AR(1) structure gives
0.452 and 0.456 (`docs/34`).

### 2.2 Correction 2: Exponential BM-BR Correlation Decay

Replaces the step-function rule under which the BM-BR correlation
dropped to zero the moment a participant was off drug regardless of
carryover residual.

$$\mathrm{Cor}(B_i, BR_{it}) = c_{bm} \cdot
  e^{-\lambda t_{sd}}$$

with $\lambda = \ln 2 / t_{1/2}$ by default (code `lambda_cor`). The
correlation decays in proportion to the drug's pharmacokinetic
elimination rate, and it is zero at off-drug visits before first
exposure. The derivation from first principles appears in Section 5
below. `docs/37` shows that this graded coupling removes the published
power collapse (see the note in Section 5.7).

### 2.3 Correction 3: nlme::lme with corCAR1

Replaces `lme4::lmer` with random intercept only.

```{.r}
nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc,
          random = ~1|ptID,
          correlation = nlme::corCAR1(form = ~t|ptID),
          data = datamerged,
          control = nlme::lmeControl(opt = 'optim',
                                     maxIter = 200,
                                     msMaxIter = 200))
```

Under the AR(1) DGP, `lmer` with random intercept alone inflates
Type I error to 13-17 percent for CO and OL designs. The `corCAR1`
residual structure restores nominal Type I control (Section 4
below).

### 2.4 Correction 4: Carryover Scale Factor Removal

Removes the undocumented scale factor of 2 that commit 8609f12
added to the carryover decay exponent:

- Original (commit 42ac030): `(1/2)^(t_sd / t_half)`.
- Ron Thomas edit (commit 8609f12): `(1/2)^(2 * t_sd / t_half)`,
  halving the effective half-life.
- Revised HEAD: `(1/2)^(t_sd / t_half)`. The parameter `t_half`
  means what it says.

### 2.5 Correction 5: Timepoint-1 Guard Removal

Commit 8609f12 introduced a guard `if (p > 1)` in the BM-BR
correlation assignment loop that was intended to avoid an index
error in the proportional-scaling formula. The guard inadvertently
forced the BM-BR correlation at timepoint 1 to zero even when
participants were on drug at that timepoint. The revised code
iterates from `p = 1` without the guard.

### 2.6 Correction 6: Gompertz Origin-Passing

Replaces the naive parameterization `f(t) = M * exp(-d * exp(-r*t))`,
which does not pass through the origin, with an offset-rescaled form
that guarantees `f(0) = 0` and preserves the asymptote at `M`:

```
y <- maxr * exp(-disp * exp(-rate * t))
vert_offset <- maxr * exp(-disp)
y <- y - vert_offset
y <- y * (maxr / (maxr - vert_offset))
```

Doc 15 (tidyverse alignment plan) identifies this correction as the
single largest source of divergence between the original and the
pre-merge tidyverse implementation. The naive form shifts every
component mean in the MVN distribution.

### 2.7 Correction 7: Updated Default Parameters

- Autocorrelation: 0.8 (compound symmetry) to 0.7 (AR(1)).
- Biomarker correlation: publication grid `{0, 0.3, 0.6}` reduced
  to `{0, 0.25, 0.45}`. PD-feasible maximum under AR(1) is 0.49
  for OL; `c.bm = 0.6` fails PD checks for 7 of 9 trial paths at
  zero carryover and all 9 at higher carryover under Architecture B.
  *Note (2026-10-07):* these April 2026 grid-search figures were not
  reproduced by the later closed-form ceilings. At $\rho = 0.7$ the
  current `covar` structure gives design-level ceilings of 0.619 (CO),
  0.456 (Hybrid) and 0.452 (OL+BDC) at $t_{1/2} = 0$. These rise with
  the half-life in Hybrid and OL+BDC, to 0.544 and 0.533 at
  $t_{1/2} = 1$ (`docs/34`, Table 4.1, verified there). On those
  figures CO paths admit $c_{bm} = 0.6$, so the statement that "all 9"
  paths fail at higher carryover is not supported. The reference value
  0.45 lies within 0.002 to 0.006 of the Hybrid and OL+BDC ceiling
  without carryover.
- Carryover half-life: publication grid `{0, 0.1, 0.2}` weeks
  replaced with `{0, 0.5, 1.0}` weeks. The publication values are
  too short to be clinically meaningful and were used precisely
  where the step-function correlation artifact was largest.

### 2.8 Correction Status Across Implementations

From the 2026-04-16 audit:

| Correction | original | original-extended | tidyverse | simple |
|---|---|---|---|---|
| 1. AR(1) autocorrelation | yes | yes | yes | n/a |
| 2. Exponential BM-BR decay | yes | yes | yes | n/a |
| 3. nlme::lme with corCAR1 | yes | yes | yes | no (OLS) |
| 4. Scale factor removal | yes | yes | yes | n/a |
| 5. Timepoint-1 guard removal | yes | yes | yes | n/a |
| 6. Gompertz origin-passing | yes | yes | yes | partial |
| 7. Updated defaults | yes | yes | yes | partial |

Table: Status of the seven DGP corrections across the four implementations

All seven corrections are embodied in the three production
implementations. The simple sandbox omits some by design (it uses
ANCOVA on phase-mean change scores rather than LME).

---

## 3. Positive Definiteness of the Covariance Matrix

Source documents: 07-positive-definiteness-failures.tex;
03-audit-and-revision-report.tex; 04-revised-power-analysis.tex.

### 3.1 The Problem

The DGP builds a 26-by-26 covariance matrix per trial path from
user-specified correlations and standard deviations. When the result
is not positive definite, `generateData()` projects to the nearest PD
matrix via `corpcor::make.positive.definite(sigma, tol = 1e-3)`.
In the publication code the correction was silent.

Failure rates by commit:

| Commit | Description | Failures | Rate |
|---|---|---|---|
| 42ac030 | Initial publication | 40 of 162 | 24.7 pct |
| 8609f12 | Ron Thomas edits | 102 of 162 | 63.0 pct |
| f6ee86d | Autocorr typo fix | 102 of 162 | 63.0 pct |
| f70f86d | Carryover added | 102 of 162 | 63.0 pct |

Table: Positive definiteness failure rates of covariance matrices by commit

All 162 matrices in the standard parameter sweep were ill-conditioned
($\mathrm{cond}(\Sigma) > 100$) under compound symmetry.

Condition number statistics (all 162 matrices):

| Commit | Median $\mathrm{cond}(\Sigma)$ | $> 10^2$ | $> 10^3$ | $> 10^6$ |
|---|---|---|---|---|
| 42ac030 | 229 | 162/162 | 40/162 | 40/162 |
| 8609f12 | 3.1 times 10^18 | 162/162 | 102/162 | 102/162 |

Table: Covariance matrix condition number statistics by commit across 162 matrices

Minimum eigenvalues of failed matrices ranged from -33.26 (worst) to
-0.52 (least severe), on a dominant eigenvalue of approximately 700.
The most extreme corrections displaced about 5 percent of the matrix
structure, concentrated in eigenvectors affecting the BM-BR and
cross-phase PB correlations disproportionately.

*Note (2026-10-07).* `docs/36`, working from the published commit
`58b32a9`, reports 40 of 162 failing matrices. These are every matrix
with $c_{bm} > 0$ and no carryover in CO, OL+BDC and Hybrid, plus the
placebo-first CO path at every half-life. It attributes them to the
requested biomarker-response correlation exceeding the CS ceiling of
0.256 to 0.262 rather than to the response block (`docs/36`, Section
1, point 1). This differs from the OL+BDC expectancy account in
Section 3.2.

### 3.2 Root Cause: OL+BDC Expectancy Discontinuity

Failures were concentrated in the OL+BDC design (83-92 percent
failure rate, 32 of 40 failures under the initial commit). The OL
design never fails. The mechanism:

- Expectancy (code `e`) changes abruptly: $\eta_t = 1.0$ in the
  open-label phase (4 timepoints) and $\eta_t = 0.5$ in the blinded
  phase (4 timepoints).
- The PB standard deviation scales with expectancy:
  $\sigma_{PB,t} = \sigma_{PB}\, \eta_t$, producing the pattern
  $(10, 10, 10, 10, 5, 5, 5, 5)$.
- Cross-phase covariance under compound symmetry with $\rho = 0.8$:
  $0.8 \cdot 10 \cdot 5 = 40$, yielding a pairwise correlation of
  $40 / \sqrt{100 \cdot 25} = 0.8$, consistent at the pair level.
- However, compound symmetry requires all same-factor timepoint
  pairs to have identical correlation regardless of temporal
  separation. With heterogeneous variances across phases this
  becomes self-contradictory: a single correlation coefficient
  cannot simultaneously satisfy the within-phase structure in both
  phases and the cross-phase structure.
- *Correction (2026-10-08).* The preceding bullet is incorrect.
  Heterogeneous variances only rescale the correlation matrix,
  $\Sigma = SRS$ with $S$ the diagonal matrix of SDs, and a positive
  diagonal rescaling does not affect positive definiteness: the PB
  block is positive definite whenever its compound-symmetry correlation
  matrix is (smallest eigenvalue $1 - 0.8 = 0.2$ here). The published
  compound-symmetry response block is positive definite at the
  published values (smallest eigenvalue 0.100 for 4 to 32 visits,
  independent of the schedule), and the failures arise through the
  biomarker row, where the coupling $c_{bm}$ exceeds the exact ceiling
  the response block supports (for example, 0.256 to 0.262 at
  $\rho = 0.8$ without carryover, below both published values 0.3 and
  0.6; `docs/36` attributes all 40 failing matrices of its 162-matrix
  audit to the biomarker row). The failure
  counts above are unaffected; only the mechanism stated here is
  wrong (derived; verified numerically). See
  `analysis/report/compendium-mathematics.tex`, subsections "The joint
  correlation matrix and the exact ceiling" and "The published
  compound-symmetry construct", Equations `eq:ce-ceiling` and
  `eq:ce-cs-decomp`.

OL never fails because expectancy is constant ($\eta_t = 1$) and all
variances are homogeneous within each factor block.

### 3.3 Resolution

Two mechanisms restore PD:

1. **AR(1) autocorrelation replaces compound symmetry.** Under AR(1),
   nearby timepoints are more correlated than distant ones, which
   accommodates variance heterogeneity without forcing a single
   correlation to cover all pairs. Maximum feasible $c_{bm}$
   increases. *Corrected 2026-10-07:* it is 0.256 under CS at
   $\rho = 0.8$ with no carryover (`docs/36`). At $\rho = 0.7$ it is
   0.452 to 0.456 (OL+BDC, Hybrid) under the current `covar` AR(1)
   structure and 0.474 to 0.480 under the separable AR(1) structure
   (`docs/34`, `docs/35`). The earlier "0.49 or higher" compared two
   values of $\rho$.

2. **Pre-validation via `validateParameterGrid()`.** The revised
   code tests every parameter combination for PD before simulation
   begins and reports which cells would require correction. Under
   the revised parameter grid (`c.bm <= 0.45`, `rho = 0.7`,
   clinically realistic half-lives), zero of 162 matrices require
   correction.

Maximum feasible $c_{bm}$ by design:

| Design | Original (CS, $\rho = 0.8$) | Revised (AR(1), $\rho = 0.7$) |
|---|---|---|
| OL | 0.89 | 0.49 |
| OL+BDC | 0.26 | 0.45 |
| CO | 0.25 | 0.61 |
| N-of-1 | 0.25 | 0.53 |

Table: Maximum feasible $c_{bm}$ by design under original and revised parameters (April 2026 grid search)

The publication value of $c_{bm} = 0.6$ exceeds the PD-feasible range
for every design except OL. The revised value of $c_{bm} = 0.45$ is
within range for all designs.

*Note (2026-10-07).* The original column agrees with `docs/36` (OL
0.893, OL+BDC 0.262, CO and Hybrid 0.256). The revised column does
not agree with the exact closed-form ceilings of `docs/34` at
$t_{1/2} = 0$: CO 0.619, N-of-1 (Hybrid) 0.456, OL+BDC 0.452. Under
AR(1) without graded coupling, `docs/36` gives OL as 0.576. Prefer the
closed-form values. The value 0.45 is inside the range but within
0.006 of the boundary.

### 3.4 Architecture A Removes the Constraint

Under Architecture A, the BM-BR correlation block is not entered
into the covariance matrix (it is applied post-draw as an additive
shift). The `c.bm` parameter therefore has no effect on the
eigenvalue structure. Testing `c.bm` from 0 to 0.95 in steps of 0.05
across 4 designs and 3 half-lives:

| | Architecture B | Architecture A |
|---|---|---|
| First PD failure | $c_{bm} = 0.50$ | none through $\beta_{bm} = 0.95$ |
| Feasible ceiling | approximately 0.45 | unconstrained |
| Paths failing at 0.60 | 7 of 9 | 0 of 9 |
| Paths failing at 0.95 | 9 of 9 | 0 of 9 |

Table: Positive definiteness limits on the moderation strength (`c.bm`) for Architecture B and Architecture A

### 3.5 Residual PD Instrumentation Gap (Finding C1)

The audit-2026-04-16 flagged that all three production implementations
still call `make.positive.definite()` silently when a matrix is not PD
(root `generateData.R:284-299`; corresponding lines in original,
original-extended, tidyverse). No flag is exposed to downstream
consumers. The current PD correction rate therefore cannot be verified
without instrumentation. The recommended fix is to return a
`pd_corrected` flag from `buildSigma()` and aggregate it into
per-parameter-set summaries.

*Note (2026-10-07, inspected `R/generateData.R`).* The repair is
still applied without a returned flag. `buildSigma()` now issues a
warning with the smallest eigenvalue, but only when `verbose = TRUE`.
`docs/34` (Section 1, point 6) measures the perturbation just above
the ceiling at no more than 0.004 in any correlation.

The doc 19 (mean moderation implementation notes) also flags that
`validateParameterGrid()` builds sigma matrices with the Architecture B
BM-BR correlation block regardless of the simulation's actual
`dgp_architecture` setting. Under Architecture A with high `c.bm`, the
validator may report spurious PD failures. The data generation itself
is unaffected because the simulation uses the Architecture A-specific
sigma. The warning can be safely ignored when running Architecture A.
As of 2026-10-07 `validateParameterGrid()` still takes no
`dgp_architecture` argument and calls `buildSigma()` with its default
`'mvn'` (inspected).

---

## 4. Analysis Model: nlme::lme with corCAR1

Source documents: 06-ar1-residual-correlation.tex;
03-audit-and-revision-report.tex; 04-revised-power-analysis.tex.

### 4.1 Formula Construction

The analysis model adapts to the design. For designs with
within-subject drug variation (CO, N-of-1, OL+BDC):

$$Y_{it} = \beta_0 + \beta_1 B_i + \beta_2 t
  + \beta_D D_{bc,it} + \beta_{bm:D} B_i D_{bc,it}
  + u_i + \varepsilon_{it}$$

with $D_{bc,it}$ the continuous drug indicator (code `Dbc`; 1 on
drug, exponential decay off drug) and $\beta_{bm:D}$ the interaction
of interest. The biomarker enters uncentered, as the raw `bm`
column (inspected, `R/lme_analysis.R`).

For open-label (no within-subject drug variation):

$$Y_{it} = \beta_0 + \beta_1 B_i + \beta_2 t
  + \beta_3 B_i \cdot t + u_i + \varepsilon_{it}$$

with $\beta_3$ the biomarker-by-time interaction. Random intercept
only (`random = ~1|ptID`). Residual correlation
`corCAR1(form = ~t|ptID)`.

### 4.2 Continuous-Time AR(1) Rationale

Trial timepoints are unequally spaced. The N-of-1 design uses weeks
`c(4, 8, 9, 10, 11, 12, 16, 20)` with gaps of 1 to 4 weeks. Discrete
AR(1) (`corAR1`) assumes equal spacing and is incorrect for these
designs. Continuous-time AR(1) (`corCAR1`) handles unequal spacing
via $\rho^{\lvert t_1 - t_2\rvert}$ with $t$ in weeks (nlme reports
$\rho$ as the parameter Phi). The schedule quoted is the N-of-1
example in `vignettes/01-simulate-trial-data.Rmd` (inspected).

### 4.3 Why the Random Intercept Alone Is Insufficient

Under compound symmetry the within-person correlation is constant,
and a random intercept absorbs it perfectly (residuals independent).

Under AR(1), let the within-participant deviations have total
variance $\sigma_T^2$ and correlation $\rho^{\lvert t_1 - t_2\rvert}$.
A random intercept of variance $\sigma_u^2$ leaves residuals of
variance $\sigma_\varepsilon^2 = \sigma_T^2 - \sigma_u^2$ and
covariance $\sigma_T^2 \rho^{\lvert t_1 - t_2\rvert} - \sigma_u^2$, so

$$\mathrm{Cor}(\varepsilon_{i t_1}, \varepsilon_{i t_2} \mid u_i)
  = \frac{\sigma_T^2 \rho^{\lvert t_1 - t_2\rvert} - \sigma_u^2}
         {\sigma_\varepsilon^2}
  = \rho^{\lvert t_1 - t_2\rvert}
  - \left(1 - \rho^{\lvert t_1 - t_2\rvert}\right)
    \frac{\sigma_u^2}{\sigma_\varepsilon^2}$$

(derived; corrected 2026-10-07). The earlier expression,
$\rho^{\lvert t_1 - t_2\rvert} - \sigma_u^2/(\sigma_u^2 +
\sigma_\varepsilon^2)$, omitted the rescaling by the residual
variance. The qualitative conclusion is unchanged. The random
intercept absorbs only part of the correlation. Residual correlation
remains positive for nearby timepoints (where
$\rho^{\lvert t_1 - t_2\rvert} > \sigma_u^2/\sigma_T^2$) and is
negative for distant ones.
`lme4::lmer` treats residuals as independent and therefore
underestimates standard errors.

(Qualified 2026-10-08: the direction of the misstatement depends on
the contrast and the design, so "underestimates" holds for some designs
only. For contiguous on-drug and off-drug blocks the AR(1) contrast
variance exceeds the compound-symmetry one (ratio 3.25 at four visits
per block), so a random-intercept-only model understates the standard
error and the test is anti-conservative, as in the quoted `lmer`
crossover sizes of 0.13 to 0.17; for interleaved designs the ratio is
smaller. This is argued from the contrast-variance formula, not
computed for each design. See
`analysis/report/compendium-mathematics.tex`, subsection "What a random
intercept absorbs under AR(1)", Equation `eq:inva-ri` and its note.)

### 4.4 Type I Error Correction

Under AR(1) DGP with `lme4::lmer` analysis (200 replicates, `c.bm = 0`,
no censoring, no carryover):

| Design | N=35 | N=70 |
|---|---|---|
| OL | 0.13 | 0.09 |
| OL+BDC | 0.09 | 0.05 |
| CO | 0.17 | 0.13 |
| N-of-1 | 0.06 | 0.08 |

Table: Type I error by design and N under lme4::lmer analysis

The CO design reaches 17 percent false positives at N=35, 3.4 times
nominal. CO and OL are most affected because their drug indicators
are most strongly confounded with time-dependent residual structure.

Under the revised `nlme::lme` plus `corCAR1` model:

| Design | N=35 | N=70 |
|---|---|---|
| OL | 0.03 | 0.02 |
| OL+BDC | 0.03 | 0.04 |
| CO | 0.06 | 0.08 |
| N-of-1 | 0.04 | 0.06 |

Table: Type I error by design and N under nlme::lme with corCAR1

All lie within the central 95 percent binomial range for a true size of
0.05 at 200 replicates, [0.02, 0.08] (computed), with OL at N=70 and CO
at N=70 on its edges. At 1,000 replicates (doc 04), Type I error is
2.5 to 5.8 percent across all designs.

*Note (2026-10-07).* Later work superseded the conclusion that
`nlme::lme` with `corCAR1` and model-based inference is adequate.
`docs/45` found real size excess for model-based and Kenward-Roger
inference under flexible-time and baseline-row analyses. Its own
primary recommendation (an unstructured model with Kenward-Roger) was
in turn superseded, because that model rejected at 0.057 pooled and
0.087 at N = 35 and failed every prognostic-biomarker cell. The
exploratory stress test of `docs/46` (41 cells, 28,000 trials, 14
analyses, decision rule rerun with all 14) retained ten analyses after
the size and sensitivity screens: S, C1, C2, DS, ES, EA, F1, F2, G1
and G2. Of these, ES, EA, F1, G1 and G2 are robust to a prognostic
biomarker. G2 and F1 tie on maximin shortfall (0.232 against 0.242),
so the rule selects F1, with fewer covariance parameters. F1 is the
`corCAR1` strawman with CR2 inference plus an expectancy term
(`De + bmc:De`; `docs/46` and `docs/49`). Neither reaches power 0.80 on
the grid.

### 4.5 Error Handling

`lme_analysis()` wraps the fit in a two-level fallback:

1. If `lme` with `corCAR1` fails to converge: fall back to `lme`
   without correlation structure.
2. If the fallback fails: return NA for all outputs.

Rows with NA outcomes are removed before fitting (nlme uses
`na.fail` by default). Audit finding C2 noted that the returned
`issingular` column was hard-coded FALSE
(`lme_analysis.R:178, 221` at the time); downstream consumers could
not distinguish clean fits from singular ones. The recommended fix was
to inspect `fit$apVar` for the 'Non-positive-definite approximate
variance-covariance' string and check random-effects variance for
zero eigenvalues. *Note (2026-10-07, inspected):* `lme_analysis()` now
sets `issingular` from `VarCorr(fit)`, TRUE when the random-intercept
SD is non-finite or below $10^{-6}$, and NA when both fits fail. The
`apVar` check was not adopted.

### 4.6 Coefficient Extraction

nlme returns `summary(fit)$tTable` with columns `Value`, `Std.Error`,
`p-value` (contrasted with `lme4::lmer`, which uses `Estimate`,
`Std. Error`, `Pr(>|t|)`). The extraction searches the row names for
`'bm:Dbc'` or `'Dbc:bm'` and selects the first match:

```{.r}
target <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(c))
if (length(target) == 0) {
  beta <- betaSE <- p <- as.numeric(NA)
} else {
  target <- target[1]
  p <- c[target, 'p-value']
  beta <- c[target, 'Value']
  betaSE <- c[target, 'Std.Error']
}
```

Audit finding C4 flags this as fragile; `carryover_analysis.R`
uses a more flexible pattern with a longer alias list.

---

## 5. BM-BR Correlation Decay: Derivation from First Principles

Source documents: 08-biomarker-correlation-decay.tex;
09-carryover-correlation-artifact.tex.

### 5.1 Three Historical Rules

The BM-BR correlation at an off-drug timepoint has been computed
three different ways across commits:

1. **Step function** (commit 42ac030, publication): full correlation
   $c_{bm}$ whenever the adjusted BR mean is nonzero, zero otherwise,
   $r_t = c_{bm}\, \mathbf 1[\mu_{BR,t} \neq 0]$.

2. **Proportional scaling** (commit 8609f12, Ron Thomas):
   $r_t = (\mu_{BR,t} / \mu_{BR,t-1})\, c_{bm}$. Correlation decays
   proportionally with the ratio of adjacent-timepoint BR means.

3. **Independent exponential** (current HEAD):
   $r_t = c_{bm} e^{-\lambda t_{sd}}$ (code `lambda_cor`), and
   $r_t = 0$ before first exposure. Correlation decays exponentially
   with a parameterized rate. `docs/37` calls this graded coupling,
   $g_t = \phi_t$.

### 5.2 Individual-Level Argument

Consider two participants during an on-drug period when the
population BR mean is $\mu_{BR}$:

- Participant H with biomarker one SD above mean ($b = +1$).
- Participant L with biomarker one SD below mean ($b = -1$).

From the bivariate normal conditional expectation:

$$E[BR_H] = \mu_{BR} + c_{bm} \cdot \sigma_{BR}$$
$$E[BR_L] = \mu_{BR} - c_{bm} \cdot \sigma_{BR}$$

Biomarker-specific spread:
$\Delta_{\text{on}} = 2\, c_{bm}\, \sigma_{BR}$.

### 5.3 Behavior at Drug Discontinuation

If carryover is pharmacokinetic, drug molecules clear at a rate
determined by the half-life, independent of any participant's
biomarker. Every participant's drug effect decays by the same factor
$\phi_t = (1/2)^{t_{sd}/t_{1/2}}$:

$$E[BR_H^{\mathrm{off}}] = \phi_t (\mu_{BR} + c_{bm} \sigma_{BR}),
  \quad E[BR_L^{\mathrm{off}}] = \phi_t (\mu_{BR} - c_{bm}
  \sigma_{BR})$$

Population mean off-drug: $\mu_{BR}^{\text{off}} = \phi_t \mu_{BR}$.

Biomarker-specific spread off-drug:
$\Delta_{\text{off}} = \phi_t \Delta_{\text{on}}$. The spread decays by
the same factor as the mean.

### 5.4 Implied Correlation

Setting the implied spread equal to the MVN-derived form
$\Delta_{\text{off}} = 2\, r_{\text{off}}\, \sigma_{BR}$ (which holds
when $\sigma_{BR}$ is the same on and off drug; see Section 5.8):

$$\boxed{\; r_{\mathrm{off}} = \phi_t \cdot c_{bm}
  = (1/2)^{t_{sd}/t_{1/2}} \cdot c_{bm} \;}$$

Equivalently in exponential form, since
$\phi_t = \exp(-(\ln 2 / t_{1/2})\, t_{sd})$, the correlation decay
rate (code `lambda_cor`) equals the drug elimination rate:

$$\boxed{\; \lambda = \frac{\ln 2}{t_{1/2}} \;}$$

(derived).

### 5.5 Numerical Validation

N-of-1 Path A, $c_{bm} = 0.5$, $t_{1/2} = 1.0$ week, $\sigma_{BR} = 8$:

| Metric | On drug (BD2, week 10) | Off drug (BD3, week 11) |
|---|---|---|
| $\mu_{BR,t}$ | 10.19 | 5.09 |
| $r_t$ | 0.50 | 0.25 |
| $E[BR \mid b = +1]$ | 14.19 | 7.09 |
| $E[BR \mid b = -1]$ | 6.19 | 3.09 |
| Spread $\Delta$ | 8.0 | 4.0 |
| Spread / mean | 0.785 | 0.786 |

Table: On-drug and off-drug BR mean and correlation for N-of-1 Path A

The spread halves as the mean halves (arithmetic checked 2026-10-07).
The ratio of the biomarker-specific spread to the mean is preserved at
0.79. No participant at $\pm 1$ SD has negative expected BR.

### 5.6 Comparison with Alternative Rules

Off-drug conditional expectations at BD3 ($\mu_{BR,t} = 5.09$,
$c_{bm} = 0.5$, $t_{sd} = 1$ week):

| Rule | $r_t$ | $E[BR \mid b = +1]$ | $E[BR \mid b = -1]$ | Issue |
|---|---|---|---|---|
| Zero off-drug correlation | 0 | 5.09 | 5.09 | No biomarker signal |
| Graded, $\lambda = \ln 2 / t_{1/2}$ | 0.25 | 7.09 | 3.09 | Consistent |
| Slow decay, $\lambda = 1$ | 0.18 | 6.56 | 3.62 | Under-predicts spread |
| Full correlation, $r_t = c_{bm}$ (published step rule) | 0.50 | 9.09 | 1.09 | Over-predicts spread |

Table: Off-drug conditional expectations of BR under alternative decay rules

Only graded coupling produces a spread consistent with
pharmacokinetic decay. *Corrected 2026-10-07:* the first row was
previously labeled "step function". Because $\mu_{BR,t} \neq 0$ at
BD3, the published step rule of Section 5.1 assigns full correlation
here, which is the last row. The slow-decay entries were recomputed
($5.09 \pm 0.5 e^{-1} \cdot 8$).

### 5.7 The Publication Artifact

The publication used the step-function rule at $t_{1/2} = 0.1$ and
$t_{1/2} = 0.2$ weeks. At these values the carryover residual in the
mean one week after discontinuation is at most about 0.1 percent of
peak ($2^{-10}$ at $t_{1/2} = 0.1$). The step-function rule assigned
the full $c_{bm}$ correlation to these timepoints, producing large
biomarker-specific variation around a near-zero mean. Numerical
example at BD3 with 0.1-week half-life, $c_{bm} = 0.6$,
$\sigma_{BR} = 8$, $\sigma_{bm} = 15.36$:

| Biomarker | $E[BR]$ original | $E[BR]$ revised |
|---|---|---|
| 109 (-1 SD) | -4.8 | +0.01 |
| 124 (mean) | +0.01 | +0.01 |
| 140 (+1 SD) | +4.8 | +0.01 |

Table: Expected BR by biomarker level under original and revised decay

The $\pm 1$ SD entries are $0.01 \pm 0.6 \cdot 8 = 0.01 \pm 4.8$; the
+1 SD entry was previously given as +4.9 (corrected 2026-10-07).

Under the original rule, a low-biomarker participant's conditional
expected BR is -4.8 at a timepoint where the drug effect has decayed
to about 0.1 percent of peak (previously misstated as 0.01 percent). This is clinically incoherent: it implies
that discontinuing the drug worsens symptoms beyond the natural
trajectory for low-biomarker participants, specifically because the
drug is gone. No rebound effect was hypothesized.

The analysis model sees this spurious biomarker-specific variation
as structured noise in the `bm:Dbc` term. It inflates the residual
variance of the interaction estimate, widens the standard error of
$\hat\beta_{bm:D}$, and reduces the probability of rejecting the null.
The effect is largest for designs with many off-drug timepoints
(N-of-1 has 4 of 8) and smallest for OL (always on drug).

Under the revised rule, the off-drug correlation at BD3 with a
0.1-week half-life is approximately $0.001\, c_{bm}$, effectively
zero. The Hendrickson et al. (2020) finding of a power drop from
0.95 to 0.18 at $t_{1/2} = 0.1$ weeks is therefore primarily a
modeling artifact. At clinically realistic half-lives (0.5 to 1.0
weeks), the revised correlation is substantial and carryover does
reduce power, but moderately rather than catastrophically.

*Note (2026-10-07).* The conclusion that the collapse is an artifact
of the step rule stands. The mechanism given above (added noise
inflating the standard error) has been replaced by an exact account in
`docs/36` and `docs/37`, both derived and verified by simulation.
Under the step rule the coupling is full at every off-drug visit with
any residual effect, so for every $t_{1/2} > 0$ the on/off *coupling
gap* is zero. The interaction slope is then exactly zero, and power
falls to the test size. This is a switch at $t_{1/2} = 0$, not a
decline. Graded coupling, $g_t = \phi_t = 2^{-t_{sd}/t_{1/2}}$ off
drug and 0 before first exposure, removes the collapse
(`docs/37`, Section 1, points 2 and 3). The published figure verified
in `docs/36` is Hybrid power falling from 0.74 to 0.12 at $N = 70$,
$c_{bm} = 0.3$, $t_{1/2} = 0.1$. The "0.95 to 0.18" figure quoted here
and in Section 6.3 was not re-verified.

### 5.8 The Remaining Variance Limitation

The derivation assumes $\sigma_{BR} = 8$ is constant. This is the
current implementation. The conditional variance during off-drug
periods is nearly unchanged ($\sigma_{BR}^2 (1 - r_t^2)$ with small
$r_t$), so participants have large random fluctuations around a
decaying mean. This conflates two sources of variation that should
be separated:

- **Trait variability**: between-person differences in latent drug
  responsiveness, existing before drug and persisting regardless of
  drug status. Should remain constant.
- **State-dependent variability**: variation in acute drug effect,
  scaling with effect magnitude. Should decay during carryover.

Addressing the limitation would require a time-varying $\sigma_{BR,t}$
or a variance decomposition, which would be substantial restructuring.
For the current application of comparing relative power across designs,
the fixed-variance assumption affects all designs similarly and is
acceptable as a conservative approximation.

---

## 6. Validated Results

Source documents: 04-revised-power-analysis.tex;
06-ar1-residual-correlation.tex; 14-trial-analysis-implementation-guide.tex.

### 6.1 Revised Parameter Grid

| Parameter | Value | Rationale |
|---|---|---|
| `dgp_architecture` | `'mvn'` (default); `'mean_moderation'` optional | Dual support |
| Autocorrelation | AR(1), $\rho = 0.7$ | PD-feasible, clinically realistic |
| $c_{bm}$ (Arch B) | {0, 0.25, 0.45} | Max PD-feasible approximately 0.45 for OL+BDC |
| $\beta_{bm}$ (Arch A) | up to 0.95 | PD constraint absent |
| Carryover $t_{1/2}$ | {0, 0.5, 1.0} weeks | Clinically realistic |
| $\lambda$ (`lambda_cor`) | $\ln 2 / t_{1/2}$ (auto) | Pharmacokinetic derivation |
| Cross-factor (same $t$) | `c.cf1t` = 0.2 | Unchanged from publication |
| Cross-factor (different $t$) | `c.cfct` $\times \rho^{\lvert t_i - t_j\rvert}$ | AR(1) decay |
| Analysis model | nlme::lme + corCAR1 | Required for Type I control |
| PD corrections | 0 of 162 | All matrices valid |

Table: Revised parameter grid with values and rationale

### 6.2 Type I Error at 1,000 Replicates ($c_{bm} = 0$)

| Design | N=35 | N=70 |
|---|---|---|
| OL | 0.029 | 0.025 |
| OL+BDC | 0.030 | 0.031 |
| CO | 0.053 | 0.058 |
| N-of-1 | 0.040 | 0.030 |

Table: Type I error by design and N at 1,000 replicates

### 6.3 Power at $c_{bm} = 0.45$, No Censoring, 1,000 Replicates

| Design | N | $t_{1/2} = 0$ | $t_{1/2} = 0.5$ | $t_{1/2} = 1.0$ | Interpretation |
|---|---|---|---|---|---|
| OL | 35 | 0.06 | 0.06 | 0.06 | Unaffected |
| OL | 70 | 0.09 | 0.10 | 0.09 | Unaffected |
| CO | 35 | 0.19 | 0.21 | 0.18 | Robust |
| CO | 70 | 0.36 | 0.35 | 0.30 | Robust |
| OL+BDC | 35 | 0.31 | 0.21 | 0.12 | Moderate drop |
| OL+BDC | 70 | 0.62 | 0.44 | 0.21 | Large drop |
| N-of-1 | 35 | 0.45 | 0.33 | 0.23 | Large drop |
| N-of-1 | 70 | 0.75 | 0.63 | 0.43 | Large drop |

Table: Power by design, N and half-life at $c_{bm} = 0.45$, 1,000 replicates

Without carryover the design ranking is N-of-1, then OL+BDC, then
CO, then OL, at both sample sizes. Under carryover CO draws level
with OL+BDC at $t_{1/2} = 0.5$ for N = 35 (0.21 each) and overtakes it
at $t_{1/2} = 1.0$ (0.18 against 0.12 at N = 35; 0.30 against 0.21 at
N = 70). (Corrected 2026-10-07: the earlier text ranked CO above
OL+BDC at low N, contrary to the table.) Carryover sensitivity is
moderate and gradual, not the catastrophic drop from 0.95 to 0.18 at
$t_{1/2} = 0.1$ reported in the publication (see the note in Section
5.7).

### 6.4 Runtime

Current implementation runtime on modern hardware (from doc 14):

- 200 replicates: approximately 8 minutes.
- 1,000 replicates: approximately 47 minutes.

Contributing optimizations (from doc 04 appendix): sigma caching
(pre-build 162 covariance matrices once), `furrr::future_map`
parallel execution across available cores, pre-computed Cholesky
factor, vectorized `data.table` outcome construction, pre-allocated
output collection, PD pre-validation, automatic core detection.
Together these yield roughly a 12x speedup over the publication code.

---

## 7. Repository Structure and Implementations

Source documents: 01-codebase-overview.tex; audit-2026-04-16.md;
and implementations/README.md.

### 7.1 Four Parallel Implementations

| Collection | Style | Architectures | Location | Role |
|---|---|---|---|---|
| original | data.table | B (MVN) only | `implementations/original/` | Historical baseline |
| original-extended | data.table | A + B | `implementations/original-extended/` | Production |
| tidyverse | tidyverse | A + B | `implementations/tidyverse/` | Modern alternative |
| simple | base R + tidyverse | A only | `implementations/simple/` | Pedagogical sandbox |

Table: The four parallel implementation collections with style, architectures and role

The installed package at root `R/` is effectively a copy of
`original-extended` with additional plotting, carryover analysis, and
documentation functions. The seven core R files are byte-identical
between root and original-extended except for the `dgp_architecture`
threading in `generateData.R` and `generateSimulatedResults.R`.

### 7.2 Installed Package File Inventory (root R/)

| File | Lines | Exports |
|---|---|---|
| buildtrialdesign.R | 119 | `buildtrialdesign` |
| censordata.R | 89 | `censordata` |
| datadocumentation.R | 61 | (none; roxygen stubs) |
| generateData.R | 355 | `generateData`, `buildSigma`, `validateParameterGrid` |
| generateSimulatedResults.R | 330 | `generateSimulatedResults` |
| lme_analysis.R | 237 | `lme_analysis` |
| packagedocumentation.R | 19 | (package NULL) |
| plottingfunctions.R | 341 | `PlotModelingResults`, `plotfactortrajectories` |
| utilities.R | 108 | `cumulative`, `modgompertz`, `reknitsimresults` |
| carryover_analysis.R | 377 | `characterize_carryover`, `analyze_trial_extended`, `print_carryover_summary`, `print_trial_summary` |

Table: Installed package R files with line counts and exports (April 2026)

Line counts and export lists are as of April 2026 and were not
re-counted. For example, `generateData.R` has since gained the
`'combined'` architecture, `moderation_scaling` and `components`
arguments (inspected 2026-10-07). NOTATION.md, Part 2, gives the
current list of exported functions.

### 7.3 Where Each Implementation Stands

**original**: the seven-file data.table code as corrected in the
revised repository. Architecture B only; no `dgp_architecture`
parameter. Supports OL, CO, OL+BDC, and N-of-1 designs. Has an
independent test suite (6 test files, 458 lines) using tinytest
naming conventions, driven by a custom `run_tests.R` rather than
standard package test infrastructure. Retained as the
Architecture B-only reference implementation.

**original-extended**: adds the `dgp_architecture` argument with
`match.arg` enforcement, the Architecture A additive shift block in
`generateData()`, and the BM-BR guard in `buildSigma()`. Threads the
argument through `generateSimulatedResults()`. Fifty-two additional
lines in two files relative to `original`. The `tests/` directory is
empty; coverage relies on the root `inst/tinytest/` suite and the
parity tests.

**tidyverse**: a complete reimplementation as a single 991-line
`functions.R` with snake_case naming, tidyverse selection helpers,
and three forms of carryover decay (exponential, linear, Weibull)
exposed via `carryover_decay()`. The `tests/` directory is empty.
A non-automated `test-alignment.R` lives alongside `functions.R`
and is sourced manually to validate sigma, Gompertz, and data
generation against the `original` collection.

**simple**: a single 577-line `simulation.R` that demonstrates
Architecture A mean moderation with a reduced-complexity model
(ANCOVA on phase-level means rather than LME plus corCAR1; no sigma
caching; no carryover-aware residual structure). Not package-
structured. By design, not used for publication-grade power
calculations.

### 7.4 Test Inventory

From the 2026-04-16 audit:

| Location | Files | Coverage | Assertions | Gaps |
|---|---|---|---|---|
| `inst/tinytest/` | 5 | cumulative, modgompertz, buildtrialdesign, buildSigma, generateData, lme_analysis | 57 | No Type I error check; no mean_moderation; no corCAR1 correctness; no PD recovery; no convergence/singular fit test |
| `implementations/original/tests/` | 6 | independent | 458 lines | Tinytest-style naming but sourced manually |
| `implementations/original-extended/tests/` | 0 | none | 0 | Empty |
| `implementations/tidyverse/R/test-alignment.R` | 1 | sigma, Gompertz, data parity | 5 stopifnot | Manual sourcing; no formal framework |
| `implementations/test-parity-extended-tidyverse.R` | 1 | cross-implementation | approximately 144 cells | Parity-only |

Table: Test inventory by location with coverage, assertions and gaps

The root `inst/tinytest/` suite covers core functions but has
known coverage gaps flagged as audit finding C7: no test exercises
`dgp_architecture = 'mean_moderation'`, null-effect Type I error,
`corCAR1` correctness under a known ground-truth $\rho$, or the silent
PD correction path.

Note: the project CLAUDE.md referred to the `testthat` framework at
the time. The installed package actually uses `tinytest` (configured as
`Suggests: tinytest (>= 1.4.0)` in DESCRIPTION, driven by
`tests/tinytest.R`). This discrepancy was flagged in the audit. As of
2026-10-07 the project CLAUDE.md specifies `tinytest` (inspected).

### 7.5 Bundled Data

- `results_core.rda`: publication Figure 4 results (initial commit,
  500 replicates).
- `extracted_bp.rda`: baseline parameters (BP mean/SD, CAPS mean/SD
  from pilot data).
- `extracted_rp.rda`: response parameters (Gompertz max, disp, rate,
  sd for tv, pb, br).
- `CTdata.rda`: actual clinical trial data for Vignette 3.
- `results_maxes.rda`, `results_rates.rda`, `results_trajectories.rda`:
  Figure 5 results.

### 7.6 Reproducibility Infrastructure

- **Docker**: `rocker/tidyverse:4.5.3` (ARG `R_VERSION=4.5.3`);
  non-root `analyst` user; `languageserver` for IDE support. Working
  directory `/home/analyst/project`.
- **renv.lock**: R 4.5.2 (4.4.2 in the CI workflow), approximately
  200 packages pinned from Posit PPM.
- **Makefile targets**: `help`, `r`, `rstudio`, `document`, `build`,
  `check`, `install`, `test`, `vignettes`, `manuscripts`, `deps`,
  `check-renv`, `check-system-deps`, `docker-build`, `docker-test`,
  `manuscript-01`, `manuscript-02`.
- **GitHub Actions**: `.github/workflows/r-package.yml` builds in
  `rocker/tidyverse:4.4.2`, restores from renv.lock, and runs
  `R CMD check --no-manual --ignore-vignettes` with `rcmdcheck` and
  `tinytest` installed to a separate CI tool library.
- **Pre-commit hook**: `scripts/pre-commit-parity.sh` runs
  `implementations/test-parity-original-tidyverse.R` to validate
  cross-implementation parity on each commit.

---

## 8. Architecture A: Implementation History

Source documents: 19-mean-moderation-implementation-notes.tex;
19-mean-moderation-implementation-notes.md;
02-revision-response-to-reviewer.md.

### 8.1 Four Attempts

The current Architecture A implementation was reached after three
incorrect formulations. Each is recorded because the failure modes
illuminate the calibration constraints the correct formulation must
satisfy.

**Attempt 1: Multiplicative scaling (incorrect).**

$$BR_{it}^{\text{mod}} = BR_{it}^{\text{base}} (1 + \beta_{bm} b_i)$$

Problem: during off-drug timepoints $BR_{it}^{\text{base}}$ is
approximately zero (the Gompertz evaluates to approximately 0 when
$t_{od} = 0$).
Multiplying near-zero by `(1 + small)` produces near-zero. The
signal is proportional to the magnitude of the base BR, not to the
drug exposure. Result: power dropped approximately 60 percent under
carryover, indistinguishable from Architecture B. Lesson:
multiplicative moderation conflates the interaction signal with the
response magnitude.

**Attempt 2: Additive with drug exposure decay (incorrect).**

$$BR_{it}^{\text{mod}} = BR_{it}^{\text{base}}
  + \beta_{bm}\, b_i\, h_t\, \sigma_{BR},
  \qquad h_t = \begin{cases} 1 & \text{on drug} \\
  \phi_t = (1/2)^{t_{sd}/t_{1/2}} & \text{off drug} \\
  0 & \text{never treated} \end{cases}$$

Problem: drug exposure decays during off-drug. At 2.5-week intervals
with $t_{1/2} = 1.0$, exposure values are {0.18, 0.031, 0.006, 0.001}
(computed). Moderation signal at the first off-drug timepoint is
$0.45 \cdot 1 \cdot 0.18 \cdot 8 = 0.65$ BR units, falling to 0.004.
Against per-timepoint noise $\sigma_{BR} = 8$, the signal is
undetectable.
Result: power under carryover remained low, comparable to
Architecture B. Lesson: tracking drug exposure decay in the
moderation signal reintroduces the very phenomenon Architecture A is
intended to avoid.

**Attempt 3: Additive with time-on-drug and uncentered biomarker
(incorrect scaling).**

$$BR_{it}^{\text{mod}} = BR_{it}^{\text{base}} + B_i\, t_{od,t}\, \beta_{bm},
  \qquad B_i \text{ the raw (uncentered) biomarker}$$

This was the form in the pre-merge 2025 tidyverse codebase (commit
957c8bd, `pm_functions.R` lines approximately 650-670). Problem: the
biomarker mean is approximately 124 mmHg (systolic BP) with SD
approximately 15; $t_{od}$ reaches 10 to 20 weeks. Moderation shift
at an on-drug timepoint with $t_{od} = 10$ is
$124 \cdot 10 \cdot 0.45 = 558$ BR units against a BR noise SD of 8,
a signal-to-noise ratio of about 70 to 1. Result: power = 1.00
everywhere. Lesson: dimensional analysis would have caught this
immediately. The product $B_i t_{od,t}$ has units of mmHg-weeks,
which is not commensurate with BR (a symptom-score change).

(Corrected 2026-10-08: the 558 units are the total shift at the mean
biomarker, and almost all of it is biomarker-free drift. With
$B_i = \mu_{bm} + \sigma_{bm} b_i$ the shift is
$\beta_{bm}\mu_{bm}t_{od,t} + \beta_{bm}\sigma_{bm}t_{od,t}\,b_i$; at
$\beta_{bm} = 0.45$, $\mu_{bm} = 124.33$, $\sigma_{bm} = 15.36$ and
$t_{od,t} = 10$ the drift is 559.5, which the drug and time terms
absorb and which carries no interaction information, and the moderation
that the interaction term sees is 69.1 per SD of the biomarker, about
8.6 $\sigma_{BR}$ rather than 70, growing with $t_{od,t}$. The
conclusion, a grossly overpowered signal, stands. Verified by direct
computation. See `analysis/report/compendium-mathematics.tex`,
subsection "Historical mean-moderation formulations and the calibration
check", Equation `eq:inva-attempt3`.)

**Attempt 4: Centered, on-drug-binary, calibrated (correct).**

$$BR_{it}^{\text{mod}} = BR_{it}^{\text{base}}
  + \beta_{bm}\, b_i\, \mathbf 1[D_{it} = 1]\, \sigma_{BR},
  \qquad b_i = (B_i - \mu_{bm})/\sigma_{bm}$$

Derivation: at $\beta_{bm} = c_{bm}$ this matches Architecture B's
conditional expectation from the MVN model,

$$E[BR_{it} \mid B_i, D_{it} = 1]
  = \mu_{BR} + c_{bm} \cdot \frac{\sigma_{BR}}{\sigma_{bm}}
  \cdot (B_i - \mu_{bm})$$

Why it works:

1. **Centering**: the standardized biomarker $b_i$ has mean zero.
   The population-average moderation is zero. Only deviations
   produce differential BR. Correct for an interaction, not a main
   effect.

2. **Binary indicator**: $\mathbf 1[D_{it} = 1]$ is 1 during
   treatment, 0 during off-drug. Unlike $t_{od}$ (grows with exposure)
   or $h_t$ (decays), it produces a clean constant signal during all
   on-drug timepoints and zero during off-drug. Carryover affects
   the BR population means but not the moderation term.

3. **Scaling by $\sigma_{BR}$**: a one-SD biomarker shift produces a
   $\beta_{bm}\sigma_{BR}$ shift in BR. For $\beta_{bm} = 0.45$ and
   $\sigma_{BR} = 8$ this is 3.6 BR units, calibrated to match the
   Architecture B conditional expectation.

### 8.2 Empirical Verification

Drawing 500 participants with identical parameters
($c_{bm} = \beta_{bm} = 0.45$, crossover design, no carryover):

| | Architecture A | Architecture B |
|---|---|---|
| BR shift per SD of biomarker | 3.36 | 3.48 |
| Residual SD | 5.23 | 3.89 |

Table: Empirical BR shift and residual SD for Architecture A and Architecture B

The shift magnitudes are comparable. Architecture A has slightly
higher residual variance because the BM-BR correlation in the
covariance matrix is absent (the noise components are independent),
whereas Architecture B's MVN draw produces correlated noise that
reduces the residual variance in the regression.

*Note (2026-10-07).* The two shifts lie within sampling error of the
theoretical 3.6 at $N = 500$ (SE about 0.2 to 0.3). The residual SDs
could not be checked, because the regression behind them is not
specified. A regression of a single $BR_{it}$ on $b_i$ would leave a
residual SD of about $\sigma_{BR} = 8$ under Architecture A and
$8\sqrt{1 - 0.45^2} = 7.1$ under Architecture B, not 5.23 and 3.89.

*Correction (2026-10-08).* The preceding note assumed the wrong
regression. The `.tex` version of `docs/19` (Sections 2 and 4.1) specifies a
regression of the *mean on-drug* BR on the biomarker. For the mean of
the four on-drug visits of a crossover path at 2.5-week spacing, with
AR(1) correlation $\rho = 0.7$ per week (the `buildSigma()` default),
the residual variance is $\sigma_{BR}^2\,\mathbf 1^\top A_{\mathrm{on}}
\mathbf 1/16$ under Architecture A and that quantity minus
$c_{bm}^2\sigma_{BR}^2$ under Architecture B, where $A_{\mathrm{on}}$ is
the AR(1) correlation block of the on-drug visits. This gives residual
SDs of 5.39 and 4.02, and the reported 5.23 and 3.89 lie within one
sampling SD of them (sampling SDs about 0.17 and 0.13 at $N = 500$), so
the reported values are reproduced (derived; verified). That the source
run used this schedule and $\rho$ is inferred from the match, not from
the run itself, which was not located. The explanation in the
paragraph above the note (correlated noise under Architecture B reduces
the residual variance) is consistent with this result. See
`analysis/report/compendium-mathematics.tex`, subsection "Historical
mean-moderation formulations and the calibration check", Equation
`eq:inva-calib`.

### 8.3 Hendrickson Publication Parameters Under Architecture A

With the original publication grid (`c.bm` up to 0.6, $t_{1/2}$ in
{0, 0.1, 0.2}), Architecture A makes previously infeasible
combinations evaluable. Power at $\beta_{bm} = 0.6$, 50 replicates:

| Design | N | $t_{1/2} = 0$ | $t_{1/2} = 0.1$ | $t_{1/2} = 0.2$ |
|---|---|---|---|---|
| N-of-1 | 70 | 0.98 | 0.96 | 0.96 |
| OL+BDC | 70 | 0.92 | 0.92 | 0.92 |
| CO | 70 | 0.66 | 0.62 | 0.70 |
| OL | 70 | 0.14 | 0.16 | 0.20 |
| N-of-1 | 35 | 0.68 | 0.78 | 0.70 |
| OL+BDC | 35 | 0.58 | 0.66 | 0.68 |
| CO | 35 | 0.40 | 0.28 | 0.26 |

Table: Architecture A power at $\beta_{bm} = 0.6$ by design, N and half-life

Power at $\beta_{bm} = 0.3$, 50 replicates:

| Design | N | $t_{1/2} = 0$ |
|---|---|---|
| N-of-1 | 70 | 0.48 |
| OL+BDC | 70 | 0.36 |
| CO | 70 | 0.26 |

Table: Architecture A power at $\beta_{bm} = 0.3$ by design and N without carryover

Type I error at $\beta_{bm} = 0$ is 0.00 to 0.08 across all designs
(nominal 5 percent). Observations:

1. Carryover resilience is consistent with these runs: power shows no
   systematic trend across the three half-lives for any design under
   Architecture A. The largest movement, CO at N = 35 (0.40 to 0.26),
   is about 1.5 standard errors of a difference at 50 replicates.
   Independence of power from the half-life under on-drug-only mean
   moderation is derived in `docs/37` (Section 5.3). Note that these
   half-lives are very short.
2. The design ranking is preserved:
   N-of-1 > OL+BDC > CO > OL.
3. $\beta_{bm} = 0.6$ is now evaluable. N-of-1 at $N = 70$ achieves 0.96 to
   0.98 power, previously unobtainable because the parameter failed
   PD validation under Architecture B.

---

## 9. Clinical Application and Trial Analysis

Source document: 14-trial-analysis-implementation-guide.tex.

### 9.1 Collaborator Deliverables

Three research goals mapped to implemented functions:

**Goal 1: Did prazosin improve PTSD symptoms?**

```r
result <- analyze_trial_extended(td, dat, op, threshold = 5)
print_trial_summary(result)
```

Defined in `R/carryover_analysis.R`. Outputs:

- Main drug effect (`Dbc` coefficient).
- Biomarker-treatment interaction (`bm:Dbc`).
- 95 percent confidence intervals.
- Variance components (random intercept SD, residual SD, and the
  `corCAR1` correlation $\rho$, returned as `ar1_phi`).
- Predicted response curves per participant via `predict(fit, newdata)`.
- Clinical decision threshold, returned as `threshold_bm`:
  `bm_mean + (threshold - beta_Dbc) / beta_interaction`.
  *Flag (2026-10-07, inspected):* this formula treats the `Dbc`
  coefficient as the drug effect at the mean biomarker. However,
  `lme_analysis()` fits the raw, uncentered `bm`, so that coefficient
  is the effect at $B_i = 0$. The threshold should then be
  $(\text{threshold} - \hat\beta_D)/\hat\beta_{bm:D}$ without the
  `bm_mean` offset, or the biomarker should be centered before the fit.
  The same offset appears in the `predicted_effect` function. This is a
  code defect in `R/carryover_analysis.R`, left for a code fix.

Status: complete.

**Goal 2a: Characterize carryover in actual data.**

```r
cc <- characterize_carryover(td, dat,
                             half_lives = c(0.25, 0.5, 1.0, 2.0, 4.0))
print_carryover_summary(cc)
# cc$best_t_half, cc$best_model
```

Sweeps candidate half-lives, fits `lme_analysis()` at each,
compares AIC/BIC/log-likelihood, returns the best-fitting half-life
and full model object. Validated to recover true `t_half = 1.0`
from simulated data.

Status: complete.

**Goal 2b: Identify best sequence of response measures.**

Status: not yet implemented; conceptual only.

**Goal 3a: Simulation power calculations.**

Status: complete. Validated 1,000-replicate power estimates with
Type I error 2.5 to 5.8 percent across designs. Zero PD failures
across 162 parameter combinations. Runtime 8 minutes (200 reps),
47 minutes (1,000 reps).

**Goal 3b: Actual analysis output structure.**

`lme_analysis(..., full_model_out = TRUE)` returns:

- `form`: the model formula.
- `fit`: the full nlme::lme model object.
- `data`: analysis-ready data with derived variables (`Dbc`, `tsd`,
  `De`).
- `summary`: beta, SE, p-value.

From the full model object, downstream analyses access
`summary(fit)$tTable`, `fitted(fit)`, `residuals(fit)`, `ranef(fit)`,
`intervals(fit)`, and `VarCorr(fit)`.

### 9.2 Prazosin-PTSD Context

From docs 02 and 20:

- Resting systolic blood pressure is a proxy for central noradrenergic
  tone, a PTSD subtype preferentially responsive to prazosin (an
  alpha-adrenergic antagonist).
- The biomarker is a statistical predictor of response class, not a
  causal pharmacokinetic mediator. Architecture B is the appropriate
  DGP.
- Clinical trial history: Raskind et al. 2003, 2007, 2013 (smaller
  trials, favorable results on nightmares and PTSD symptoms);
  Raskind et al. 2018 (PACT trial, N = 304, null result on primary
  endpoint); Mayer et al. 2023 (posttraumatic headache pilot;
  previously misattributed to Raskind et al.).
- The PACT null result motivates the search for a predictive
  biomarker and is the central clinical puzzle the pmsimstats
  framework addresses.

### 9.3 Implication for Trial Design

From doc 02, Section 5 (precision medicine trial design): the
architectural choice matters most in designs that deliberately
expose each participant to both treatment states to estimate
within-subject responsiveness (crossover, basket, platform, N-of-1
hybrid). Parallel-group enrichment designs escape Architecture B's
carryover sensitivity because they have no off-drug phase.
Recommendation: when the biological mechanism is uncertain, run the
power calculation under both architectures and treat the
Architecture B estimate as the conservative anchor for sample size
planning. This is consistent with the PATH Statement (Kent et al.
2020; Ann Intern Med 172(1):35-45, doi:10.7326/M18-3667).

---

## 10. Alternative Analysis Strategies Under Architecture B

Source document: 02-dgp-mean-moderation-vs-mvn.tex, Section 6.

Six candidate strategies for recovering power under Architecture B,
evaluated comparatively. The labels 6.1 to 6.7 below are the section
numbers of the source document, not of this compendium. Later
evaluations of analysis approaches are collected in `docs/43`
(catalog) and `docs/46` (stress test).

**6.1 Restrict analysis to on-drug observations only.**

Eliminates correlation erosion (all retained observations have full
`c_bm`). Costs: reduced effective sample size and loss of the
off-drug reference, which confounds biological, placebo, and
time-varying components. Best for designs with many on-drug
timepoints (OL has 8, Hybrid has 5 to 6). Least useful for CO, which
would lose half its observations.

**6.2 Weighted analysis.**

On-drug observations weight 1. Off-drug observations weight
$\exp(-\lambda t_{sd})$. Rationale: information-theoretic.
Observations with higher drug exposure carry more information about
the biomarker-treatment interaction because the correlation signal
is stronger.

**6.3 Within-subject contrast.**

For designs with both on-drug and off-drug periods (CO, Hybrid,
OL+BDC):

$$\bar{Y}_{i,\mathrm{on}} - \bar{Y}_{i,\mathrm{off}}
  = \omega_0 + \omega_1 \cdot B_i + e_i$$

(symbols renamed 2026-10-07 from $\alpha$, $\gamma$, which are
canonical for the significance level and the $BR \times PB$
interaction).

Sidesteps `Dbc` parameterization and the repeated-measures model.
Within-subject averaging reduces noise.

**6.4 Two-stage random slopes.**

Stage 1: fit a random-slopes LME.

$$Y_{it} = \beta_0 + \beta_1 t + \beta_D D_{it}
  + u_{0i} + u_{1i} D_{it} + \varepsilon_{it}$$

Stage 2: regress the estimated random deviations on the biomarker.

$$\hat{u}_{1i} = \xi_0 + \xi_1 B_i + e_i$$

(renamed from $\delta_0, \delta_1, \eta_i$, which collide with the
canonical standardized bias and expectancy covariate).

Separates within-subject treatment effect estimation (affected by
carryover) from between-subject biomarker association (which uses
cross-sectional relationships).

**6.5 Exclusion of contaminated observations.**

Exclude the first one or two off-drug timepoints, which carry the
most carryover contamination. Later off-drug timepoints provide a
cleaner contrast with the washed-out state.

**6.6 Design-level solutions.**

- Longer washout periods between drug phases.
- More on-drug timepoints relative to off-drug.
- Front-loaded drug exposure (OL+BDC design): a long initial on-drug
  block establishes the BM-BR correlation before off-drug
  observations begin.

**6.7 Comparative evaluation.**

| Strategy | Discards data | Complexity | Expected benefit |
|---|---|---|---|
| On-drug only | Yes | Low | High for OL/Hybrid; poor for CO |
| Weighted | No | Low | Moderate; depends on weight calibration |
| Within-subject contrast | Aggregates | Low | Moderate; robust to misspecification |
| Two-stage random slopes | No | Moderate | Unknown; separates estimation stages |
| Exclude early off-drug | Some | Low | Moderate; improves contrast clarity |
| Design modification | N/A | N/A | High; addresses root cause |

Table: Alternative analysis strategies by data discarded, complexity and expected benefit

---

## 11. Latent-Class and Psychometric Connections

Source document: 02-revision-response-latent-class-expanded.md.

### 11.1 The Faithful Generative Model

If the biomarker indexes an unobserved responder/non-responder
dichotomy, the mechanistically correct DGP is a finite mixture:

$$Z_i \mid B_i \sim \mathrm{Bernoulli}(\pi(B_i)),\quad
  BR_{it} \mid Z_i = z \sim f_z(t, D_{it})$$

with `Z_i` the latent class, `pi(B_i) = Pr(Z_i = 1 | B_i)` typically
a logistic function of the biomarker, and class-specific response
distributions `f_0, f_1` differing in drug-response parameters.

Key property: carryover under this DGP operates only on the
class-specific means and does not erode any correlation signal,
because class membership `Z_i` is drug-state-invariant.

### 11.2 MVN Differential Correlation as Second-Moment Approximation

A two-component mixture with class-dependent means and common
within-class covariance produces a marginal joint distribution of
`(B_i, BR_{it})` whose first two moments match, to leading order,
those of a single MVN with `Cor(B, BR) = c_bm` chosen so that
`c_bm^2` equals the between-class variance fraction. Under mild
regularity conditions (class overlap, moderate `pi`), the mixture's
second-moment structure is well-approximated by MVN with
differential correlation in treatment states where class-specific
means separate. Architecture B as implemented captures the
covariance implications of the latent-class mechanism without
committing to the discrete generative form.

(Corrected 2026-10-08: `c_bm^2` is not the between-class variance
fraction but the *product* of the two between-class variance fractions,
$c_{bm}^2 = f_B f_{BR}$ with
$f_B = \pi(1-\pi)(\mu_B^1 - \mu_B^0)^2/\mathrm{Var}(B)$ and $f_{BR}$
defined likewise for BR, under within-class independence of $B_i$ and
$BR_{it}$; equivalently $|c_{bm}|$ is their geometric mean. This
follows because the covariance is carried by class membership alone,
$\mathrm{Cov}(B, BR) = \pi(1-\pi)(\mu_B^1 - \mu_B^0)(\mu_{BR}^1 -
\mu_{BR}^0)$ (derived). See `analysis/report/compendium-mathematics.tex`,
subsection "The MVN approximation to mixture-induced moderation",
Equation `eq:latent-mvn-approx`.)

### 11.3 Where the Approximation Breaks Down

Three regimes:

1. **Bimodal responses**: when `f_0` and `f_1` are well-separated
   (non-responders show near-zero drug effect, responders show a
   large effect), the marginal BR distribution is bimodal and no
   MVN matches its tails. Power under the true mixture DGP is
   **higher** than Architecture B predicts because mixture analysis
   can exploit bimodality.

2. **Strong class-membership gradient**: when `pi(B_i)` is steep
   (approaching a step function at a threshold), the biomarker
   behaves almost deterministically as a class label and
   Architecture A becomes the better approximation.

3. **Class-varying covariance**: if autocorrelation, residual
   variance, or placebo-response parameters themselves differ
   between classes, no single-component MVN can represent the
   generative structure. Analysis models assuming constant
   covariance will be misspecified.

### 11.4 Taxonomy of Finite Mixture Models

From the psychometric tradition:

- **Latent class and latent profile analysis** (Lazarsfeld and
  Henry 1968; Goodman 1974; Clogg 1995; McLachlan and Peel 2000):
  small number of discrete subpopulations. LCA for categorical
  indicators; LPA for continuous. R: `poLCA`.
- **Growth mixture models and group-based trajectory models**
  (Muthen and Shedden 1999; Nagin 1999, 2005): class-specific
  growth curves with within-class random-effect variation. Closest
  psychometric analogue of the mixture DGP discussed above.
  Identifiability hazards (Bauer and Curran 2003; Bauer 2007)
  apply with full force at N-of-1 sample sizes.
- **Factor mixture models** (Yung 1997; Arminger, Stein, and
  Wittenberg 1999; Lubke and Muthen 2005): continuous latent
  factors combined with discrete latent classes. Subsumes both
  pure continuous (single-class factor model) and pure discrete
  (LCA with no within-class structure) as limits.
- **Regression mixtures and mixtures of experts** (DeSarbo and
  Cron 1988; Wedel and DeSarbo 1995; Jacobs et al. 1991): class-
  specific regression coefficients with covariate-dependent class
  membership probabilities.

### 11.5 The Architecture B Spectrum

Architecture B as implemented is the **most restrictive** operational
case of a broader family that ranges from mean-only moderation
(Architecture A) through combined mean-plus-covariance moderation
(Arminger et al. 1999; the canonical Psychometrika treatment) to
heterogeneous-random-effects models (Verbeke and Lesaffre 1996;
Proust-Lima et al. 2017).

Implication: the substantial power loss under carryover (40 to 60
percent) documented in Section 1.4 is a specific property of the
restrictive covariance-only case and may not generalize to
mean-plus-covariance or heterogeneous-random-effects variants.

Recommendation for a complete sensitivity analysis: report power
under at least three points on the spectrum: Architecture A
(mean-only), Architecture B as implemented (covariance-only), and a
combined mean-plus-covariance specification of the Arminger form.

### 11.6 Practical Constraint and Next Extension

Identifiability concerns for mixture-based analysis at N-of-1 trial
sample sizes (typically 30 to 100) are acute. Published FMM
applications involve hundreds to thousands of participants. The
open question is not whether mixture modeling is more faithful (it
is), but whether the fidelity gain is recoverable at realistic
sample sizes. Analysis-strategy mitigation (Section 10 above) is
more likely than mixture-model analysis to deliver recoverable
power gains at realistic trial sizes.

If the pmsimstats program pursues a richer DGP, the natural target
is the **heterogeneous-random-slopes form** rather than a full FMM.
This preserves the analysis model's compatibility with existing
N-of-1 estimation machinery (`nlme::lme` with participant-specific
random treatment effects), admits a biologically faithful latent-
class generative structure at the DGP level, preserves the audit
chain of the current framework, and permits direct comparison with
the covariance-only baseline under carryover.

R packages relevant in approximate order of fit: `lcmm` (closest
tooling for the heterogeneous-random-slopes form), `flexmix`
(regression-mixture framework), `mclust`, `gamlss`, `poLCA`,
`OpenMx` and `lavaan` with mixture extensions. Reference
implementation outside R: Mplus.

---

## 12. The 2x2 Factorial: DGP Carryover Versus Analysis-Model Carryover

Source document: 10-carryover-analysis-model-assessment.tex.

### 12.1 The Question and the Framework

Standard methodology (Senn 2002; Jones and Kenward 2014) for
assessing carryover in the analysis model is a 2x2 factorial with
DGP carryover on one axis and analysis-model carryover on the
other, reporting Type I error, power, bias, and confidence interval
coverage in each cell. Senn's modern consensus:

- The two-stage procedure (test for carryover first, then condition)
  is invalid; carryover tests have low power.
- If carryover is scientifically plausible, either design the trial
  to eliminate it via adequate washout, or always include carryover
  in the model and accept the associated power loss. Do not
  test-then-decide.
- For simulation studies: generate data under a range of carryover
  magnitudes (including zero), analyze under models that include
  and exclude the carryover term, and report operating
  characteristics for each combination.

### 12.2 Three Analysis Approaches

*Note (2026-10-07).* The approaches below predate the project's
current specification labels (NOTATION.md, Part 3). Approach A is E2
(the exposure-weighted predictor $D_{bc,it}$), and Approach C is E1
(binary indicator, carryover ignored). Approach B, a linear
time-since-discontinuation covariate, is not E3; E3 adds a lagged
just-off-drug indicator. The study proposed in Section 12.3 was carried
out, in revised form, as the carryover-sensitivity manuscript, paper 02
(`analysis/report/02-carryover-sensitivity/`).

**Approach A: Continuous drug indicator ($D_{bc,it}$, code `Dbc`).**

Formula: `Sx ~ bm + t + Dbc + bm:Dbc`. $D_{bc,it}$ equals 1 on drug
and $(1/2)^{s\, t_{sd}/t_{1/2}}$ off drug, where $s$ is
`op$carryover_scalefactor` (default 1, inspected). Conflates treatment and
carryover; captures both treatment-biomarker and carryover-
biomarker interactions simultaneously.

**Approach B: Separate carryover covariate (tsd).**

Formula: `Sx ~ bm + t + Db + bm:Db + tsd`. The linear time-since-
discontinuation term absorbs carryover. Closer to standard
crossover methodology but uses a linear carryover term, which is
misspecified relative to the exponential DGP.

**Approach C: Ignore carryover.**

Formula: `Sx ~ bm + t + Db + bm:Db`. Off-drug observations treated
as if no drug were present. Biases the treatment contrast toward
zero, reducing power, and may inflate Type I error if the bias is
systematic.

### 12.3 Proposed Simulation Study

Design:

- DGP factor: $t_{1/2}$ in {0, 0.25, 0.5, 1.0, 2.0} weeks.
- Analysis factor: three models (Ignore, Separate, Continuous).
- Fixed: N = 70, $\rho = 0.7$, all four trial designs.
- Biomarker correlation: $c_{bm}$ in {0, 0.3, 0.45}.
- 1,000 replicates per cell.
- Total cells: $5 \times 3 \times 3 \times 4 = 180$.

Metrics: Type I error (at $c_{bm} = 0$), power, bias, RMSE, coverage.

Expected findings (labels corrected 2026-10-07; they previously
numbered the approaches 1 to 3 in the reverse of the A to C order of
Section 12.2):

- Approach C (Ignore): correct Type I error at $t_{1/2} = 0$, may
  inflate at higher carryover; power declines with increasing
  carryover.
- Approach B (Separate): correct Type I error across all carryover
  values; power slightly lower at $t_{1/2} = 0$ (cost of extra
  parameter), higher at $t_{1/2} > 0$; mild bias from linear versus
  exponential misspecification.
- Approach A (Continuous): correct Type I error if decay rate
  matches DGP; power depends on whether the conflation helps or
  hurts.

Infrastructure: models supported via existing codebase through
`op$carryover_t1half` and `op$simplecarryover` parameters.

Status (April 2026): proposed simulation study; not yet complete.
See the note at the start of Section 12.2 for its later form.

---

## 13. Simplification Alternatives (Not Adopted)

Source document: 17-simplification-plan.tex.

Doc 17 proposed reducing complexity of the pmsimstats framework
substantially. The 2026-04-08 merge took a different path (adding
Architecture A rather than simplifying). The document's novel
contribution is the most thorough formal specification of six
analysis methods available in the project, including precision
equations and design implications. It is preserved here because the
analytical framework remains useful for cross-method comparisons.

### 13.1 Six Interaction Test Methods

**1. ANCOVA on phase-mean change scores (recommended in doc 17).**

Collapse to phase means, $\bar Y_{i,\text{on}}$ the mean of $Y_{it}$
over visits with $D_{it} = 1$, and similarly $\bar Y_{i,\text{off}}$.
Compute $\Delta_i = \bar Y_{i,\text{on}} - \bar Y_{i,\text{off}}$.
Regress $\Delta_i = \omega_0 + \omega_1 B_i + e_i$ (renamed from
$\gamma_0, \gamma_1$; $\gamma$ is canonical for the $BR \times PB$
interaction). Test $H_0\colon \omega_1 = 0$ by a t-test on
$\hat\omega_1$ with $N - 2$ df. Precision:
$\mathrm{Var}(\hat\omega_1) \approx \sigma_\Delta^2 / (N\,
\mathrm{Var}(B))$.

**2. Paired t-test on change scores.**

Tests `H_0: E[Delta_i] = 0` (main effect, not interaction).
Requires dichotomization of the biomarker (e.g., median split) for
the interaction. Efficiency loss from dichotomization approximately
36 percent (`2/pi` approximately 0.64 relative efficiency).

**3. Generalized Estimating Equations (GEE).**

Population-averaged model with sandwich estimator. Key precision
parameter: `Var_within(D_i) = D_bar_i (1 - D_bar_i)`, maximized at
`D_bar = 0.5` (equal on/off time). For the OL design,
`Var_within(D) = 0`, so the interaction is not estimable.

**4. Random-intercept LMM.**

Subject-specific model. Precision depends on
`SS_within(D_i) = sum_t (D_it - D_bar_i)^2`. Design-specific
values: OL has 0, OL+BDC approximately 1.5, CO 2.0, N-of-1
approximately 2.9. These values use doc 17's occasion counts: OL+BDC
with 8 visits, about 6 on drug; CO with 8 visits, 4 on drug; N-of-1
with 12 visits, about 7 on drug. With binary $D$ and $T$ visits,
$SS_{\text{within}} = T\bar D(1-\bar D) \le T/4$, so the current
8-visit Hybrid design cannot exceed 2.0 (derived 2026-10-07). In doc
17, $L$ is the block length; $K$ is not defined there and is read
from context as the number of blocks.
Short blocks (high $K$, low $L$) suffer less from
autocorrelation than long blocks because fewer consecutive
observations are highly correlated.

**5. Simple linear regression (ignoring clustering).**

Downward-biased standard errors; Type I error inflation with
ICC = 0.3 and T = 8 is approximately 15 to 25 percent at nominal
alpha = 0.05. Not recommended for inference.

(Corrected 2026-10-08: the inflation applies only to a
*between*-participant effect. With ICC 0.3 and $T = 8$ the design
effect is $1 + 7 \cdot 0.3 = 3.1$, giving a realized size of
$2\Phi(-1.96/\sqrt{3.1}) = 0.266$ (verified). For the
within-participant interaction $\beta_{bm:D}$, the project's estimand,
OLS charges the contrast the full residual variance
$\sigma_u^2 + \sigma_\varepsilon^2$ rather than $\sigma_\varepsilon^2$,
so its standard errors are overstated and the test is conservative,
not anti-conservative (derived for a pure within-participant contrast
with exchangeable errors; Paper 10 reports a standard-error ratio near
1.3). OLS remains unsuitable for inference, for the
opposite reason. See `analysis/report/compendium-mathematics.tex`,
subsections "A general decomposition of realized Type I error",
Equation `eq:kappa-ols`, and "What a random intercept absorbs under
AR(1)".)

**6. Repeated-measures ANOVA (biomarker dichotomized).**

Requires dichotomization (36 percent information loss). Sphericity
assumption violation requires Greenhouse-Geisser or Huynh-Feldt
correction, further reducing power. RM-ANOVA retains approximately
0.64 times ANCOVA efficiency. Design-specific behavior: OL uses
time as the within-subject factor; CO is a clean two-level factor
with balanced phases; N-of-1 with multiple blocks raises
sphericity concerns.

### 13.2 Cross-Design Summary (doc 17, Table 8.1)

| | OL | OL+BDC | CO | N-of-1 |
|---|---|---|---|---|
| Interaction term | bm times t | bm times D | bm times D | bm times D |
| Target | Total improvement rate | Drug-specific | Drug-specific | Drug-specific |
| Confounded with | PB, TV slopes | none | none | none |
| LMM info | sum bm^2 SS(t) | sum bm^2 (1.5) | sum bm^2 (2.0) | sum bm^2 (2.9) |
| Transitions | 0 | 1 | 1 | 3 |
| Autocorr vulnerability | Low | Medium | Medium | Low per block |
| Carryover vulnerability | None | Low | Low | High |

Table: Interaction term and target of the test across the four designs

### 13.3 The OL Design Is a Fundamentally Different Test

For the open-label design, all six methods face the same problem:
there is no within-subject drug variation. The interaction must be
formulated as biomarker-by-time. This is a weaker test than the
biomarker-by-drug interaction for two reasons: confounding (the
slope $dY/dt$ includes BR, PB, and TV; the biomarker-by-time
interaction conflates all three) and precision (the within-subject
information for the time slope depends on
`SS_within(t) = sum(t - t_bar)^2`, which is typically smaller than
`SS_within(D)`).

### 13.4 Adoption Status

Not adopted. The 2026-04-08 merge preserved the LME with corCAR1
analysis model and added Architecture A rather than simplifying to
ANCOVA on phase-mean change scores. The doc is retained because
the formal specification of analysis method precision across
designs is the most detailed in the project and remains useful for
power planning in related contexts.

---

## 14. Code Audit Findings (2026-04-16)

Source document: audit-2026-04-16.md.

The audit identified nine issues across three severity tiers.

### 14.1 P0 (Execution or Correctness Blockers)

**Finding T1**: Unqualified dplyr calls in the tidyverse Architecture
A path.

Location: `implementations/tidyverse/R/functions.R:421-423`.

```r
bm_mean <- baseline_param |> filter(cat == 'bm') |> pull(m)
bm_sd   <- baseline_param |> filter(cat == 'bm') |> pull(sd)
br_sd   <- resp_param     |> filter(cat == 'br') |> pull(sd)
```

These lines execute only when `dgp_architecture == 'mean_moderation'`.
If called without `library(dplyr)` in scope, they will either fail
at runtime or resolve `filter()` to `stats::filter`. The parity test
masks the issue because it runs in a tidyverse-loaded session. Fix:
qualify both calls with `dplyr::` and add a regression test
invoking `generate_data(..., dgp_architecture = 'mean_moderation')`
with only declared Imports loaded.

**Finding `A1`**: Broken source path.

Location: `analysis/carryover_factorial/02-run-factorial-2025.R:20`.

Sources `analysis/2025/01-pm-functions.R`, which has been moved to
an external archive. The script fails at load time. Fix: change the
source target to `implementations/tidyverse/R/functions.R` and
update downstream function-name references, or restore a copy from
the external archive if historical comparison is needed.

### 14.2 P1 (Publication and CRAN Blockers)

**Finding C1**: Silent positive-definiteness correction.

Locations: `generateData.R:284-299` (root), `orig/generateData.R:242-253`,
`ext/generateData.R:284-299`, `tid/functions.R:337-356`.

When `is.positive.definite(sigma)` returns FALSE, all
implementations call `make.positive.definite(sigma, tol = 1e-3)`
without recording the correction. Doc 04 reports 24.7 percent of
matrices in the publication required correction. The current audit
cannot verify the present rate without instrumentation. Fix: return
a `pd_corrected` flag from `buildSigma()` and aggregate into
per-parameter-set summaries.

**Finding C2**: Hard-coded `issingular = FALSE`.

Locations: `lme_analysis.R:178, 221` (copies in orig, ext);
`tid/functions.R:654`.

The `issingular` column is always FALSE even when the underlying
fit is singular. nlme exposes `fit$apVar` and
`fit$modelStruct`; singularity is inferable. Fix: detect
singularity from `fit$apVar` (character string 'Non-positive-
definite approximate variance-covariance') and check random-effects
variance for zero eigenvalues; populate `issingular` accordingly.
(Status 2026-10-07: addressed in the root package by a `VarCorr`
check; see Section 4.5.)

**Finding C5**: Unvalidated `carryover_t1half`.

Location: `generateData.R:172-177` (copies).

A negative `carryover_t1half` produces a negative `lambda_cor`,
turning `exp(-lambda_cor * tsd)` into an exponentially growing
correlation. No validation. Fix: add `stopifnot(carryover_t1half >= 0)`
at the top of `generateData()`. *Note (2026-10-07, inspected):* in
the current root `buildSigma()` a non-positive half-life sets
`lambda_cor` to 0, so the correlation cannot grow. However, the BR
carryover mean term `last_on * (1/2)^(tsd / carryover_t1half)` grows
without bound for a negative half-life, and there is still no
validation.

**Finding C6**: Stale imports.

The root DESCRIPTION Imports list includes `lme4`, `lmerTest`,
`ggpubr`, and `gridExtra`, none of which are called anywhere in
`R/`, `vignettes/`, or the three package-style implementations.
These are residues from the pre-nlme analysis model and earlier
plotting code. Fix: remove them from DESCRIPTION; if any surfaces in
a vignette, promote to Suggests. The audit notes these will generate
NOTEs on CRAN. (Status 2026-10-07: all four are still listed under
Imports in DESCRIPTION; inspected.)

**Finding C7**: Missing test coverage.

No test exercises `dgp_architecture = 'mean_moderation'`, null-effect
Type I error at nominal alpha, `corCAR1` correctness under known
`phi` ground truth, or the silent PD correction path. Fix: add four
`tinytest` cases (or migrate to `testthat` 3e) covering these gaps.

**Finding `A2`**: Missing READMEs in `analysis/` subfolders.

Locations: `analysis/02-carryover-sensitivity/`,
`analysis/architecture_comparison/`, `analysis/carryover_factorial/`,
`analysis/figure5/`. Each contains numbered driver scripts without a
README. Fix: add a one-paragraph README per subfolder describing
research question, inputs, outputs, and runtime.

### 14.3 P2 (Hygiene and Maintainability)

**Finding C3**: `eval(parse(text = ...))` for column selection.

Locations: `censordata.R:50, 84`; `plottingfunctions.R:136`;
`lme_analysis.R:93` (partial). Fix: replace with `get()`,
`dat[[name]]`, or `dplyr::all_of(name)`. No user-facing behavior
change.

**Finding C4**: Positional coefficient extraction.

Locations: `lme_analysis.R:193-217` (copies); `tid/functions.R:637-650`;
`carryover_analysis.R:71-93`. Relies on
`intersect(c('bm:Dbc', 'Dbc:bm'), rownames(...))` and `target[1]`. `carryover_analysis` uses a more
flexible pattern. Fix: harmonize around the `carryover_analysis`
pattern.

### 14.4 Parity Status

From the audit: `parity_baseline.rds` and the recent
`parity_report.rds` differ by two bytes, indicating stable parity
between `original-extended` and `tidyverse` at tolerances of 1e-10
(sigma, data) and 1e-6 (beta, SE, p-value). No expected divergences
are recorded in `parity_divergences.csv` (zero data rows).

### 14.5 Publication Readiness Summary

- CRAN blockers: T1 will fail `R CMD check --as-cran`; C6 will
  generate NOTEs; C7 weakens the submission narrative.
- Manuscript blockers: `A1` will break reviewer reproduction; C1
  (silent PD correction) should surface in methods or appendix but
  the current rate is not instrumented; C2 (singularity flag always
  FALSE) may affect reported convergence statistics.

The audit judged that, with the P0 and P1 items addressed, the
package and manuscript could proceed to submission.

---

## 15. Parity Development Workflow

Source document: 22-tidyverse-development-parity-workflow.md.

### 15.1 Baseline and Inner Loop

Before starting a tidyverse refactor:

```bash
git checkout -b tidyverse-refactor
Rscript implementations/test-parity-extended-tidyverse.R
cp implementations/parity_report.rds implementations/parity_baseline.rds
```

Commit `parity_baseline.rds` with the branch. At the end of
development, diff the current `parity_report.rds` against it to
confirm only intended differences.

Two inner-loop cadences:

- `--quick` flag: approximately 2 seconds, 16 cells. Use for
  iteration on a specific function.
- Full sweep: approximately 10 seconds. Run before every commit.

An optional pre-commit hook blocks commits on failure:

```bash
Rscript implementations/test-parity-extended-tidyverse.R --quick || \
  { echo 'parity broken; run --full for details'; exit 1; }
```

### 15.2 Edit Classification

- **Pure refactors** (rename, split function, swap dplyr for purrr,
  extract helper): parity must hold. If it breaks, a bug was
  introduced.
- **Feature additions with default values**: parity must still
  hold for default call paths.
- **Behavior changes** (different DGP parameterization, different
  MVN ordering, different sigma caching): parity will break by
  design; decide whether to apply to both implementations in
  lockstep.

Most development falls into the first two categories. The third
should be treated as a separate project with explicit scope.

### 15.3 Breakage Diagnosis

When parity breaks mid-refactor, the test prints the first failing
cell. Diagnose in this order:

1. **Sigma failure**: the change likely touched
   `build_sigma_matrix` or `build_correlation_matrix`. The BM-BR
   correlation branch is the most common culprit.
2. **Data failure**: touched MVN sampling, post-draw mean
   moderation, or participant/path assembly.
3. **Analysis failure**: touched formula construction, `Dbc`
   computation, or `nlme::lme` invocation.

The `--quick` run covers 16 cells in about 2 seconds. Re-running
after every approximately 20 lines of changes localizes bugs to
the edit that introduced them.

### 15.4 RNG Guardrails

The parity test calls `set.seed()` at the per-cell level. If a
refactor consumes random draws in a different order (for example
swapping `MASS::mvrnorm` for a Cholesky-based draw), parity will
break even though the code is semantically correct. Defenses:

- Do not change RNG-consuming loops unless `original-extended`
  is updated in the same commit.
- Keep `n_cores = 1` for parity-relevant comparison. Introduce
  parallelism only after functional parity is locked in.

### 15.5 Intentional Divergence

When a refactor intends to change tidyverse behavior in a way that
`original-extended` should not follow (for example fixing a bug
that lives only in data.table code, or using tidyverse as a
testbed for a new feature):

- Note in the commit message or `DIVERGENCES.md` which cells are
  allowed to fail and why.
- Mark those cells in the script with an `expected_divergent` flag.
- Re-baseline after the intentional divergence lands:
  `cp parity_report.rds parity_baseline.rds`.

The trap to avoid is letting an intentional divergence mask a
subsequent unintentional one.

### 15.6 Pre-Merge Checklist

- [ ] cp parity_report.rds parity_baseline.rds at start of branch.
- [ ] --quick after every focused edit (approximately 2 s).
- [ ] Full sweep before every commit (approximately 10 s).
- [ ] If parity breaks, read the first failing cell; fix before
  proceeding.
- [ ] Never change RNG consumption order in a pure refactor.
- [ ] Add edge cells when new features are added.
- [ ] Document any intentional divergence explicitly.

One-liner for diff before merging:

```r
a <- readRDS('implementations/parity_baseline.rds')
b <- readRDS('implementations/parity_report.rds')
all.equal(a$grid, b$grid)  # TRUE, or character vector of diffs
```

---

## 16. Publication Strategy

Source document: 22-journal-target-recommendation.md.

### 16.1 Manuscript Summary

Manuscript: 'Two Architectures for Simulating Biomarker-Treatment
Interactions: Implications for Statistical Power Under Carryover'.
The manuscript formalizes the distinction between direct mean
moderation (Architecture A) and differential correlation via MVN
(Architecture B), with the central finding that architecture choice
determines carryover impact on power: 10 to 15 percent loss under
Architecture A, 40 to 60 percent under Architecture B. (These
magnitudes come from the 50-replicate pilot; see the note in Section
1.4 for the matched-$N$ re-run of paper 01.)

Methodological positioning: intersection of statistical simulation
methodology, N-of-1 and crossover trial design, and biomarker-
treatment interaction modeling.

### 16.2 Ranked Journal Targets

1. **Statistics in Medicine** (primary). A leading journal for
   statistical methodology in clinical trial design. Primary
   audience (biostatisticians designing trials) would act on the
   finding. Collection precedent: Senn (2016; Stat Med 35(7):966-977).

2. **Biometrics**. Shi and Ye (2024; Biometrics 80(2):ujae023,
   previously misattributed to "Yap, Wang, and Harhay") is the closest
   comparator (carryover impact on power in crossover designs).
   The manuscript extends this work by showing that the
   mechanism determining the interaction affects the magnitude of
   power loss.

3. **Trials**. Highest concentration of N-of-1 methodology papers
   in the related literature (4 of 50 in the annotated
   bibliography). Open access, practical orientation.
   Collection precedents: Haller and Ulm (2018), Hawksworth et al.
   (2024; Trials 25:263). The latter was previously listed as two
   papers, "Chen et al. (2024)" and "Senn and Julious (2024)", neither
   of which could be found; flagged 2026-10-07.

4. **Statistical Methods in Medical Research**. Publishes
   simulation-based methodological work for trial design.
   Collection precedent: Senn (2019), doi:10.1177/0962280217726801
   (listed originally as "Arends et al. (2019)", a misattributed
   entry carrying Senn's title; verified 2026-10-07).

5. **Contemporary Clinical Trials Communications**. Practical
   framing and clinical application fit well. Collection
   precedent: Sturdevant and Lumley (2021). The originally listed
   "Liu et al. (2021)" was a misattributed entry for Lim and In
   (2021), *Korean Journal of Anesthesiology*, doi:10.4097/kja.21165,
   and is not a precedent for this journal (verified 2026-10-07).

6. **Harvard Data Science Review**. Collection precedent:
   Schork (2022). More of a perspective venue; heavy simulation
   results may not fit format.

7. **Healthcare (MDPI)**. Collection precedents: Blackston et al.
   (2019), Wang and Schork (2019). Fast review; lower stigma
   risk if strengthened before submission.

Primary recommendation: target **Statistics in Medicine**. Backup:
target **Trials** if Statistics in Medicine considers the scope
too specialized.

---

## 17. Literature Synthesis

Source document: 20-annotated-bibliography.tex.

Fifty references in the annotated bibliography organized by topic.
The summaries below focus on the claims most directly relevant to
the pmsimstats framework.

*Note (2026-10-07).* Each entry below was checked against Crossref
and PubMed. Seven were misattributed or had wrong metadata, and they
have been corrected in place, with the earlier citation recorded. The
same errors probably occur in `docs/20-annotated-bibliography.tex`,
which was not checked. DOIs are given for verified entries.

### 17.1 Core N-of-1 Simulation Methodology

- **Hendrickson, Thomas, Schork, and Raskind (2020)**, Frontiers in
  Digital Health 2:13, doi:10.3389/fdgth.2020.00013. Simulation study
  (four designs, LME analysis, 1,000 repetitions) of aggregated N-of-1
  designs for detecting predictive biomarker effects. Finding: the
  N-of-1 (Hybrid) design outperformed OL and OL+BDC but had slightly
  lower power than the traditional crossover. This is the framework
  pmsimstats-ng revises. (Corrected 2026-10-07: previously cited as
  "Thomas et al.", described as deriving analytic power formulae,
  and summarized as needing fewer participants than parallel-group
  designs, none of which matches the abstract.)
- **Blackston et al. (2019)**, Healthcare 7(4):137,
  doi:10.3390/healthcare7040137. 5,000
  simulated trials per scenario across four variance structures
  comparing aggregated N-of-1, parallel RCT, crossover.
  Finding: aggregated N-of-1 outperforms traditional designs but
  has elevated Type I error under moderate-to-strong unmodeled
  carryover.
- **Wang and Schork (2019)**, Healthcare 7(3):84,
  doi:10.3390/healthcare7030084. Analytical
  power formulae under AR(1) serial correlation. Finding: serial
  correlation substantially reduces power; increasing evaluation
  periods partially compensates.

### 17.2 Carryover and Temporal Correlation

- **Shi and Ye (2024)**, "Behavioral carry-over effect and power
  consideration in crossover trials," Biometrics 80(2):ujae023,
  doi:10.1093/biomtc/ujae023 (corrected 2026-10-07; previously
  attributed to "Yap, Wang, and Harhay"). Formal framework for
  behavioral carryover in crossover trials (psychological mechanisms,
  not pharmacokinetic). Closest comparator for the pmsimstats
  manuscript. Finding: behavioral carryover biases the basic
  estimator toward underestimating the treatment effect, costing
  power. Directly relevant to the prazosin-PTSD application.
- **Sturdevant and Lumley (2021)**, Contemporary Clinical Trials
  Communications 22:100711, doi:10.1016/j.conctc.2021.100711.
  Mixed-effects model approach for
  testing carryover. Finding: 'carryover has been largely ignored
  except as a nuisance parameter'.
- **Schork (2022)**, Harvard Data Science Review Special Issue 3,
  doi:10.1162/99608f92.f1eef6f4. Serial correlation in personalized
  and aggregated studies. Finding: ignoring serial correlation
  inflates false positives.
- **Tang and Landes (2020)**, PLoS One 15(2):e0228077,
  doi:10.1371/journal.pone.0228077. Modified
  t-tests accounting for serial correlation in N-of-1 data.
  Finding: proposed tests maintain nominal Type I where standard
  t-tests fail.

### 17.3 Biomarker-Treatment Interactions

- **Haller and Ulm (2018)**, Trials 19:128,
  doi:10.1186/s13063-018-2491-0. Simulation study of
  biomarker-treatment interaction estimation in randomized trials
  with time-to-event outcomes. Same genre as the pmsimstats
  manuscript (Architecture A paradigm).
- **Grenet et al. (2020)**, Basic and Clinical Pharmacology and
  Toxicology 126(1):59-64, doi:10.1111/bcpt.13293. Simulation comparing crossover versus
  parallel-group for binary predictive biomarker identification.
  Finding: crossover more efficient when within-subject
  correlation is high (supports N-of-1 use in pmsimstats).

### 17.4 Analysis Methods

- **Chen and Chen (2014)**, "A comparison of four methods for the
  analysis of N-of-1 trials," PLoS One 9(2):e87752,
  doi:10.1371/journal.pone.0087752 (corrected 2026-10-07;
  previously attributed to "Senn, Julious, and Araujo"; the Araujo,
  Julious and Senn paper is a different one, PLoS One 2016,
  11(12):e0167167). Four methods for N-of-1 analysis
  (fixed-effects, random-effects, Bayesian hierarchical, paired).
  Finding: all similar for population inference; differ
  substantially for individual treatment effect estimation.
- **Peterlin, Kejzar, and Blagus (2023)**, TEST 32(1):184-210,
  doi:10.1007/s11749-022-00830-1.
  Diagnostic tests for correct specification of fixed- and
  random-effects design matrices in LMM.

### 17.5 Sample Size

- **Senn (2019)**, Statistical Methods in Medical Research
  28(2):372-383, doi:10.1177/0962280217726801. Five distinct
  criteria for planning sample size in N-of-1 trials. Finding:
  criterion depends on goal (individual-level, population-level,
  heterogeneity detection).
- **Schmid and Duan (2014)**. Statistical design and analytic
  considerations for N-of-1 trials. Chapter 4 in Kravitz RL, Duan N
  (eds.), *Design and Implementation of N-of-1 Trials: A User's
  Guide*, AHRQ Publication No. 13(14)-EHC122-EF. Rockville, MD:
  Agency for Healthcare Research and Quality; 2014. The chapter title
  and its place in the guide are confirmed; the chapter's author list
  rests on the AHRQ product page as reported by a web search, because
  the AHRQ site blocks automated access (2026-10-07). The original
  entry listed Kravitz as a third author.
- **Hawksworth et al. (2024)**, "A methodological review of
  randomised n-of-1 trials," Trials 25:263,
  doi:10.1186/s13063-024-08100-1 (corrected 2026-10-07; previously
  cited as "Senn, Julious, and Chen (2024), Trials 25:264").
  Methodological review of published randomized N-of-1 trials.
  Finding: widespread inconsistencies in reporting and analytic
  practice.

### 17.6 Bayesian Approaches

- **Schmid and Yang (2022)**, "Bayesian models for N-of-1 trials,"
  Harvard Data Science Review, doi:10.1162/99608f92.3f1772ce
  (corrected 2026-10-07; previously "Schmid and Staudenmayer
  (2021)"). Review of Bayesian modeling for N-of-1 trials.
- **Daza (2018)**, "Causal analysis of self-tracked time series data
  using a counterfactual framework for N-of-1 trials," Methods of
  Information in Medicine 57(1):e10-e21, doi:10.3414/ME16-02-0044
  (corrected 2026-10-07; previously cited as Daza, Wac, and Oppezzo,
  AMIA Annual Symposium Proceedings 2018:532-541, which was not found).
  Counterfactual framework for causal inference in N-of-1 trials via
  self-tracked time series.

### 17.7 Prazosin-PTSD Clinical Context

- **Raskind et al. (2003)**, American Journal of Psychiatry
  160(2):371-373, doi:10.1176/appi.ajp.160.2.371. Initial
  proof-of-concept (N=10), double-blind crossover; significant
  reduction in nightmares.
- **Raskind et al. (2013)**, American Journal of Psychiatry
  170(9):1003-1010, doi:10.1176/appi.ajp.2013.12081133. Parallel-group RCT (N=67) in active-duty
  soldiers; significant effects; observed effect sizes and
  heterogeneity calibrate simulations.
- **Raskind et al. (2018)**, New England Journal of Medicine
  378(6):507-517, doi:10.1056/NEJMoa1507598. PACT trial (N=304); no
  significant difference versus placebo on primary outcome. The null
  result motivates the pmsimstats framework search for a predictive
  biomarker.
- **Mayer et al. (2023)**, "Randomized controlled pilot trial of
  prazosin for prophylaxis of posttraumatic headaches in active-duty
  service members and veterans," Headache 63(6):751-762,
  doi:10.1111/head.14529 (corrected 2026-10-07; previously
  attributed to Raskind et al.). Pilot for posttraumatic headache
  prophylaxis (N=48; sample size not re-checked).

---

## 18. Core Contribution Statement

The pmsimstats-ng project establishes that the choice between direct
mean moderation (Architecture A) and differential correlation
(Architecture B) as the DGP for biomarker-treatment interaction in
aggregated N-of-1 designs determines whether carryover effects
reduce statistical power modestly (approximately 10 to 15 percent,
Architecture A) or substantially (approximately 40 to 60 percent,
Architecture B). The choice is a claim about biology: Architecture
A for biomarkers that determine the magnitude of the drug's
biological effect (renal clearance, receptor density, CYP2D6 status);
Architecture B for biomarkers that predict response class through
shared variance with a drug-responsive physiological component
(blood pressure as PTSD subtype marker, baseline severity,
inflammatory markers). Most published simulation studies use
Architecture A. The Hendrickson et al. (2020) framework, which
pmsimstats-ng revises and extends, appears to be unique in the N-of-1
literature in using Architecture B (Section 1.7). Published power estimates for biomarker-
moderated crossover and N-of-1 designs may be systematically
optimistic for biomarkers operating through an Architecture B
mechanism.

The revised implementation embodies seven DGP corrections relative
to the publication (AR(1) autocorrelation, exponential BM-BR decay,
nlme with corCAR1, scale-factor removal, timepoint-1 guard removal,
origin-passing Gompertz, updated defaults), supports both
architectures via the `dgp_architecture` parameter, validates Type
I error at nominal 5 percent across all designs, and executes 1,000
replicates in 47 minutes with zero silent positive-definiteness
corrections.

*Note (2026-10-07).* This statement is as of April 2026. The power-loss
magnitudes are those of the 50-replicate pilot (Section 1.4). The
closed-form account of carryover and power is in `docs/37`. The
feasible range of $c_{bm}$ is in `docs/34` to `docs/36`, and at the
reference value 0.45 it sits within 0.006 of the boundary. The
recommended analysis is F1 from the stress test of `docs/46` and
`docs/49`, not the plain `corCAR1` model.
