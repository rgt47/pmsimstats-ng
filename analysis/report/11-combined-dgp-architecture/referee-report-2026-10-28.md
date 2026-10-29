# Referee report: paper 11, the dual-channel architecture
*2026-10-28 22:32 PDT*

**Manuscript.** "A combined data-generating architecture for
biomarker-treatment interactions: channel super-additivity and
mixing-weight sensitivity in aggregated N-of-1 trials"
(`analysis/report/11-combined-dgp-architecture/report.Rmd`, source
dated 2026-09-07).

**Target journal assumed.** *Statistics in Medicine* (methods paper on
a simulation data-generating process). The comments apply equally to
*Trials*.

**Reviewer.** pmsimstats team (internal referee pass).

**Evidence labels.** VERIFIED: I executed code and observed the
result. INSPECTED: I read the source and confirmed it. SUSPECTED:
consistent with the evidence but not confirmed. UNVERIFIED: not
checked. Scratch scripts are in the session scratchpad, prefixed
`p11-` (listed in Section 8).

## 1. Summary

The manuscript defines a data-generating process (DGP) for a
biomarker-treatment interaction with two independent weights. The
mean-channel weight $c_{bm,A}$ adds $c_{bm,A}\,b_i\,\sigma_{BR}$ to
the drug-response component at on-drug visits. The covariance-channel
weight $c_{bm,B}$ sets the on-drug biomarker-response correlation,
which decays off drug at the carryover half-life. Setting either
weight to zero recovers one of the two architectures of paper 01.

The paper makes two claims. First, the biomarker-response covariance
delivered to the analysis model is additive in the two weights, and
because power is convex in the noncentrality below about one half,
the channels combine "super-additively" in power. Second, a
$3 \times 3$ grid of weights crossed with three designs and three
half-lives (1000 replicates per cell, specification E2,
`nlme::lme` with `corCAR1`) shows monotone power, interior cells
above the probability-scale additive benchmark, comparable carryover
sensitivity of the two channels, and a preserved design ranking.
Section 4.2 states candidly that the boundary cells do not reproduce
paper 01's single-architecture baselines and attributes the gap to a
calibration mismatch. The Abstract Results are marked as provisional.

## 2. Overall assessment and recommendation

**Recommendation: reject in the present form. A substantially
reframed resubmission could be considered.**

The manuscript is not ready for external review, and the problems
are not cosmetic. The stored simulation was produced on 2026-07-31
by code that predates four corrections the program has since adopted.
One of them corrupts every crossover cell with positive half-life
(Major 1). The boundary-validation failure has a specific cause that
the manuscript misidentifies: paper 01's analysis includes the
baseline visit and this driver's does not (Major 2). The paper's
headline result, "channel super-additivity", is a property of the
probability scale. A single channel doubled shows the same excess,
and a probit-scale additive benchmark reproduces the combined cells
(Major 3). The analytical section omits the off-drug profile, which
is where the two channels actually differ, and it misstates the
direction of the variance effect (Major 4). Most fundamentally, under
a jointly normal biomarker the mean shift is linear in $b_i$. Both
channels then act on the conditional mean $\mathbb{E}[BR \mid B]$,
and the "dual-channel" model is a single Gaussian model with a
modified covariance. The mechanistic premise of two distinct channels
acting on different moments does not hold as specified (Major 5).
Several prose claims in Sections 3 and 4.3 are contradicted by the
manuscript's own table (Major 7).

Strengths, stated briefly. The boundary-reduction property is a
sound validation device, and the code implements it correctly. No
stored cell depended on repair of a non-positive-definite matrix
(verified below). Section 4.2 discloses the failed validation
openly rather than concealing it. The empirical additivity of the
fitted interaction estimate at $t_{1/2} = 0$ is real and reproduces
from the stored data.

## 3. Significance and novelty

The question is legitimate. A simulation study of power for a
predictive biomarker should state how the interaction enters the
DGP, and a family that nests the two conventions is useful as a
sensitivity framework.

The novelty claimed is narrower than stated. (a) The additivity of
the conditional slope under joint normality is elementary regression
algebra: $\mathbb{E}[BR_t \mid b] = \mu_t + (c_{bm,A} D_t +
c_{bm,B} g_t)\,\sigma_{BR}\,b$. The package source already states
this equivalence in a comment (`R/generateData.R`, lines 190-194).
(b) "Super-additivity" relative to $p(a,0) + p(0,b) - p(0,0)$ is
the generic behavior of any two additive effects on the noncentrality
scale when power lies below one half. It is not specific to the
channels (Major 3). (c) The Discussion concedes that the combined
form is introduced conceptually in paper 03. What remains is a
descriptive grid that has not passed its own validation gate.

Contemporariness. The manuscript cites no work after 1999 outside
the series. It does not engage the location-scale and
variance-function literature, the variance-ratio literature on
heterogeneity of treatment effects (HTE), or simulation-reporting
standards (Section 6). The analysis model (E2, true half-life,
model-based standard errors) is not the program's current
recommendation. Paper 02 recommends G9 or F1 with CR2 standard
errors.

## 4. Major comments

Correctness items come first.

### M1. The stored run predates four corrections, and one corrupts all crossover cells with positive half-life

*Location.* Table `tab:arch-c-power` (CO rows); Section 4.3;
Reproducibility; `analysis/data/dgp-combined/01-combined-*.rds`
(file time 2026-07-31 11:14).

*Problem.* I extracted the driver, `simulation-core.R` and
`implementations/tidyverse/R/functions.R` at commit 26b8c16a3
(2026-07-31), which matches the data timestamp. Against the current
code, the stored run differs in four ways.

1. *No pre-exposure guard in $D_{bc}$.* In July, `prepare_long_data()`
   set `Dbc = carryover_decay(tsd, h)` at every off-drug visit. The
   design sets `tsd = 0` before first exposure, and the decay at
   zero equals 1 for $h > 0$. In the CO design, path 2 is
   placebo-first. Under the July rule its participants are coded
   `Dbc = 1` at all eight visits at $t_{1/2} = 0.5$ and $1.0$. Half
   of every CO sample therefore contributes no within-participant
   on/off contrast to `bm:Dbc`. Hybrid and OL+BDC paths all begin on
   drug and are unaffected. The CO rows of the manuscript's only
   table ($t_{1/2} = 1.0$) are affected.
2. *Cross-component correlation.* `c.cfct = 0.1`, not the
   specification value $c_1 = 0.2$ (`docs/61`).
3. *Carryover rule.* The original compounding recursion
   (`mu[t] = base[t] + mu[t-1] * phi`) was used, not carryover that
   decays from the last on-drug mean.
4. *Matrix repair.* `make_positive_definite = TRUE` with no stop.
   This turned out not to matter (see "Positive definiteness" below).

The Reproducibility section implies that the current driver
reproduces the stored output. It would not, because the driver now
sources the corrected `simulation-core.R`.

*Why it matters.* The CO results at positive half-life are
uninterpretable. The other cells describe a DGP that the program no
longer recommends.

*Remedy.* Re-run the full grid under the current code, with the
extended DGP stated explicitly. Record the commit hash in the stored
output, and state in the manuscript which DGP rules were used.

*Evidence.* VERIFIED. I evaluated $D_{bc}$ for every path under both
rules (`p11-02-pd.R`). CO path 2 at $h = 1$ gave `1 1 1 1 1 1 1 1`
under July code and `0 0 0 0 1 1 1 1` under current code. The other
three differences are INSPECTED from the extracted sources.

*Positive definiteness (no finding).* I rebuilt all 81
design-by-half-life-by-weight covariance matrices (all paths) under
the July parameterization (`c.cfct = 0.1`) and under the current one
(`c.cfct = 0.2`). None was non-positive-definite. The smallest
eigenvalues were 0.144 and 2.87, both at OL+BDC with
$c_{bm,B} = 0.45$. Unlike sibling papers, no stored result depends on
silent repair. VERIFIED (`p11-02-pd.R`).

### M2. The boundary-validation gap is caused by the baseline visit, not by a calibration mismatch

*Location.* Section 4.2; Reproducibility; driver header and
`02-summarise-combined.R`.

*Problem.* Section 4.2 attributes the 10 to 19 point shortfall to a
mismatch in `resp_param`, `baseline_param`, `c.cf1t`, `c.cfct` or
`rho`. Those inputs are identical in the two drivers. Both load
`extracted_rp`/`extracted_bp` and use $\rho = 0.7$,
`c.cf1t = 0.2` and `c.cfct = 0.1`, and both allocate a total of
$N = 70$ across paths (INSPECTED: `01-dgp-prototype.R` lines
121-182; July `simulation-core.R`). The analyses differ, however.
Paper 01's baseline numbers come from the package `lme_analysis()`.
That function prepends a baseline row ($t = 0$, off drug,
$D_{bc} = 0$) to every participant (`R/lme_analysis.R` lines 82-104).
The combined driver's `prepare_long_data()` uses the eight
post-baseline visits only.

I regenerated the boundary cells at $t_{1/2} = 0$ and $N = 70$ under
the July code, 400 replicates per cell. Each dataset was analyzed
both ways (paired).

| Design | Cell | No baseline row | With baseline row | Paper 01 |
|---|---|---|---|---|
| CO | (0.45, 0) | 0.580 | 0.723 | 0.739 |
| CO | (0, 0.45) | 0.605 | 0.745 | 0.745 |
| Hybrid | (0.45, 0) | 0.694 | 0.780 | 0.842 |
| Hybrid | (0, 0.45) | 0.730 | 0.807 | 0.840 |
| OL+BDC | (0.45, 0) | 0.578 | 0.759 | 0.751 |
| OL+BDC | (0, 0.45) | 0.607 | 0.787 | 0.777 |

The Monte Carlo standard error (MCSE) is at most 0.025 per entry.
Adding the baseline row closes the CO and OL+BDC gaps to within one
MCSE. It closes about half of the Hybrid gap. A residual of 0.03 to
0.06 (1.3 to 2.5 MCSE) remains and may reflect a further difference
between the package and tidyverse generators.

Two further problems make the gate ill-posed as written. (a) The
target numbers (0.739, 0.842, ...) come from a paper 01 version
that has since been restated. The live paper 01 is restricted to the
covariance encoding under the extended DGP and no longer reports
mean-moderation baselines (`CLAUDE.md`, "Key architectural
feature"). (b) The driver header promises an "embedded" validation
gate that warns at a 3-point deviation. No such comparison exists in
either script. `02-summarise-combined.R` prints the $N = 35$
boundary cells only and compares them to nothing.

*Why it matters.* The paper's stated validation principle (Section
2.2) is sound, but the manuscript misdiagnoses its failure and points
the follow-up at the wrong inputs.

*Remedy.* Define the validation target as a run of the same driver
with `dgp_architecture = 'mean_moderation'` and `'mvn'` under common
random numbers. Boundary identity should then hold draw for draw,
which is a much stronger check than agreement within 3 points. If
agreement with a paper 01 table is also wanted, align the analysis
dataset first (baseline row in or out). Implement the warning the
header describes, or delete the claim.

*Evidence.* VERIFIED for the baseline-row effect (`p11-03-blrow.R`,
log `p11-03-blrow.log`). INSPECTED for parameter identity and for
the missing gate. SUSPECTED for the cause of the Hybrid residual.

### M3. "Channel super-additivity" is a probability-scale effect, not a property of the channels

*Location.* Title; Abstract; Section 3, paragraphs 3-4; Section 4.3,
paragraph 2; Discussion.

*Problem.* I tested the claim three ways against the stored $N = 70$
E2 grid.

1. *Probit-scale benchmark.* I mapped each single-channel power to
   an approximate noncentrality, $\hat\Lambda = \Phi^{-1}(\hat\pi) +
   1.96$. I then predicted the combined cell as
   $\Phi(\hat\Lambda_A + \hat\Lambda_B - 1.96)$. This additive-signal
   prediction matches most interior cells within Monte Carlo error.
   Examples: Hybrid $t_{1/2} = 0$, (0.22, 0.22): observed 0.711,
   predicted 0.693. CO $t_{1/2} = 1$, (0.45, 0.45): 0.813 against
   0.828. CO $t_{1/2} = 0$, (0.45, 0.22): 0.897 against 0.900. The
   largest positive departures are in OL+BDC (for example 0.576
   against 0.467 at $t_{1/2} = 0$). The deflated null rejection in
   that design (0.024 to 0.028) biases the probit inversion there.
2. *Within-channel control.* The same "super-additivity" appears when
   a single channel is doubled. Compare $p(0.45, 0)$ with
   $2p(0.22, 0) - p(0, 0)$. At Hybrid $t_{1/2} = 0$ the values are
   0.738 against 0.447 for the mean channel and 0.738 against 0.413
   for the covariance channel. The excesses (0.29 and 0.33) match the
   between-channel excess at (0.22, 0.22), which is 0.281. The
   phenomenon therefore says nothing about the interplay of the
   channels.
3. *Sign reversals.* The probability-scale benchmark exceeds 1 in six
   cells, all at (0.45, 0.45): CO, Hybrid and OL+BDC at
   $t_{1/2} = 0$, Hybrid and OL+BDC at $0.5$, and Hybrid at $1.0$. In
   those cells the combined power falls short of the benchmark by
   0.07 to 0.44. That includes the Hybrid $t_{1/2} = 1$ cell in the
   manuscript's own table (0.995 against 1.089). The claim that
   "interior cells exceed the ... benchmark" is false for these cells.

The convexity argument also has a gap. Superadditivity of a convex
$\pi$ on $[0, \Lambda_A + \Lambda_B]$ requires the sum to lie in the
convex region, not each term separately as written. The statement
that single-channel powers lie "between roughly 0.03 and 0.30" is
contradicted by the table, where they range from 0.080 to 0.656 at
$t_{1/2} = 1$ (up to 0.738 at $t_{1/2} = 0$). Several single-channel
cells therefore lie above the inflection.

*Why it matters.* The title and the stated contribution rest on this
result.

*Remedy.* Report additivity where it lives: on the scale of the
fitted interaction estimate and of $\hat\Lambda$. State that
probability-scale excess is expected for any additive effects below
$\pi \approx 0.5$. Drop "super-additivity" from the title, or define
it relative to a probit-scale benchmark and report the
within-channel control alongside it.

*Evidence.* VERIFIED (`p11-01-summary.R`, `p11-04-analytic.R`).

### M4. The analytical section omits the off-drug profile and misstates the variance effect

*Location.* Section 3, paragraphs 1-2; Section 2.2, last paragraph.

*Problem 1: variance direction.* The analysis model conditions on
the biomarker. Conditional on $B$, the mean shift is deterministic
and adds no residual variance. The $c_{bm,A}^2\sigma_{BR}^2$ term is
a marginal inflation only. The covariance channel, by contrast,
reduces the on-drug conditional variance of $BR$ to
$\sigma_{BR}^2(1 - c_{bm,B}^2)$, from 64 to 51.0 at
$c_{bm,B} = 0.45$. The manuscript has the roles reversed. At the
outcome level the reduction is small, a factor of 0.978 to 0.982 on
$\mathrm{Var}(Y_t \mid B)$ at on-drug Hybrid visits. That is
consistent with the nearly constant empirical SD of the estimate
across $c_{bm,B}$, but the stated mechanism is wrong.

*Problem 2: off-drug profile.* The derivation equates the two
channels through the on-drug covariance alone. They differ off drug.
The mean shift is applied at on-drug visits only, so the mean
channel's interaction has no carryover, while the drug-response mean
itself does carry over (`functions.R`, the `on_drug` loop). The
covariance channel's correlation decays as $c_{bm,B}\,\phi_t$. The
conditional slope at visit $t$ is therefore
$(c_{bm,A} D_t + c_{bm,B} g_t)\,\sigma_{BR}$, which is not a
function of $c_{bm,A} + c_{bm,B}$ unless $t_{1/2} = 0$. The stored
data show this. At $t_{1/2} = 1$ the mean fitted `bm:Dbc` estimate is
as follows (expected $-0.234$ per biomarker unit at a total weight of
0.45).

- Hybrid: $-0.284$ at (0.45, 0), $-0.229$ at (0, 0.45), $-0.260$ at
  (0.22, 0.22).
- OL+BDC: $-0.326$, $-0.239$, $-0.278$.

Power at the three equal-sum cells is 0.656, 0.473 and 0.577 in
Hybrid, and 0.498, 0.277 and 0.370 in OL+BDC. This contradicts the
statement that such pairs give "very similar power". The mean
channel's estimate also overshoots the nominal $-0.234$ by 21 to 39
percent. Its step-shaped interaction is analyzed with an
exposure-decayed regressor, so E2 is misspecified for that channel.

*What does hold.* At $t_{1/2} = 0$ the estimate is additive to three
decimals. For example, Hybrid gives $-0.236$, $-0.232$ and $-0.472$
at (0.45, 0), (0, 0.45) and (0.45, 0.45), against an expected
$0.45 \times 8 / 15.36 = 0.234$ per channel. The non-identifiability
statement is correct at $t_{1/2} = 0$ and only approximately correct
otherwise.

*Remedy.* Derive the conditional mean and covariance of $Y$ given
$B$ for the full visit vector, off-drug visits included. Give the
projection of the true slope profile onto $D_{bc}$ under the
analysis half-life. State the identifiability result conditional on
$t_{1/2} = 0$, or as "identifiable only through the off-drug
profile", which is itself a useful finding. Decide whether the mean
channel should carry over, as the package's "trajectory" option
partly does, and justify the choice.

*Evidence.* VERIFIED. I computed conditional moments from the
covariance matrices (`p11-04-analytic.R`). Estimates and power come
from the stored summary (`p11-01-summary.R`).

### M5. As specified, the two "channels" are one Gaussian model: the mechanistic and latent-class framing does not hold

*Location.* Sections 2.2, 2.3 and 5; Discussion; the Arminger et al.
and Verbeke and Lesaffre citations.

*Problem.* The biomarker is drawn jointly normal with the response
components, and the mean shift is linear in $b_i$. The dual-channel
DGP is therefore a single multivariate normal law. Its on-drug
$\mathrm{Cov}(B, BR_t)$ is $(c_{bm,A} + c_{bm,B})\sigma_{bm}
\sigma_{BR}$, plus a rank-one term $c_{bm,A}^2\sigma_{BR}^2$ shared
by on-drug visits. Both channels act on $\mathbb{E}[BR \mid B]$.
Neither makes $\mathrm{Var}[BR \mid B]$ depend on $B$: the
conditional variance is constant in $B$ under both. Three statements
in the manuscript fail as a result.

- Section 2.3 says the covariance channel carries information
  "through the class-conditional dispersion $\mathrm{Var}[BR \mid
  B]$". It does not.
- Sections 2.3 and 5 describe "pharmacodynamic scaling" and
  "co-regulatory coupling" as distinct mechanisms. As specified, they
  are observationally the same mechanism on the conditional mean.
  They differ only in the off-drug profile and in a small on-drug
  variance reduction (M4).
- Section 2.2 calls the combined form "the mean-and-covariance-
  structure model of Arminger et al.". That paper concerns finite
  mixtures of conditional mean- and covariance-structure models. The
  present DGP has no mixture component. Verbeke and Lesaffre (1996)
  likewise concern a mixture distribution for random effects.

*Why it matters.* The paper's reason to exist is that a drug may
"activate one, both, or neither channel". If the channels are
indistinguishable on the conditional mean at $t_{1/2} = 0$, that
sentence has no empirical content within this DGP.

*Remedy.* Either (a) reframe the paper honestly. The covariance
encoding is mean moderation with a reduced on-drug residual variance
and a decaying off-drug profile, which the program's own `CLAUDE.md`
already states. The "dual-channel" family is then a two-parameter
family of off-drug interaction profiles. Or (b) introduce a genuinely
distinct second channel, for example a heteroscedastic or
location-scale term in which the on-drug variance of $BR$ depends on
$B$ (Hedeker et al., 2008), or a responder mixture as in paper 03.
Under (b) the latent-class language and the Arminger citation become
appropriate.

*Evidence.* INSPECTED (algebra and code). VERIFIED numerically for
the conditional slope and variance (`p11-04-analytic.R`: the on-drug
slope of $Y$ on $B$ was $-0.2344$ at $c_{bm,B} = 0.45$, and the
conditional $\mathrm{Var}(BR)$ was 51.04 against a marginal 64).

### M6. The analysis model and the scenario grid are not the program's current standard, and test calibration differs by design

*Location.* Section 4.1; Discussion ("common-mode effect").

*Problem.* The simulation fits E2 at the true half-life with
model-based standard errors and no expectancy terms. Empirical
rejection under the null cell ranges from 0.017 to 0.040 in OL+BDC
and Hybrid, and from 0.046 to 0.060 in CO (MCSE 0.004 to 0.008).
That is a design-dependent conservatism. The Discussion's claim that
it applies "as a common-mode effect across the channel-weight grid"
is defensible within a design. It is not defensible for the design
comparisons the paper also draws. Paper 02 recommends G9, the E2
form with an AIC-chosen half-life and CR2 errors, with F1 as a
fallback. The half-life grid $\{0, 0.5, 1.0\}$ weeks consists, by
the program's policy, of behavioral scenarios, not the
pharmacological range of 0.1 to 0.2 weeks. Finally, the CSV
artifacts (`02-power-grid-*.csv`, `02-power-wide.csv`) and the
summary script are restricted to $N = 35$, while the manuscript
reports $N = 70$. Section 4.1 never states $N$.

*Remedy.* Re-run under the extended DGP with G9 (or F1) and CR2
errors. Report null rejection with its MCSE for every design. Add
$t_{1/2} \in \{0.1, 0.2\}$. Regenerate the CSV artifacts at the
reported $N$.

*Evidence.* VERIFIED for null rejection and for the CSV contents.
INSPECTED for the analysis specification.

### M7. Results prose contradicts the table, and the outstanding items are unaddressed

*Location.* Abstract; Introduction, paragraph 2; Section 4.3.

*Problem.*

- "The covariance and mean channels contribute comparably to
  carryover sensitivity" is contradicted by the data. From
  $t_{1/2} = 0$ to $1$ at $N = 70$, Hybrid power falls from 0.738
  to 0.656 under the mean channel but from 0.738 to 0.473 under the
  covariance channel. OL+BDC falls from 0.612 to 0.498 against 0.601
  to 0.277. The covariance channel loses two to three times as much,
  consistent with paper 01.
- "The design ranking ... is preserved throughout the combined grid"
  is false. At $t_{1/2} = 0$, CO exceeds OL+BDC at (0.22, 0.22)
  (0.620 against 0.576), while OL+BDC exceeds CO at (0, 0.45) (0.601
  against 0.581). At $t_{1/2} = 1$, CO exceeds OL+BDC at (0, 0.22)
  and (0, 0.45). Several of these contrasts are within two MCSE of
  zero, so the honest statement is that the ranking of CO and OL+BDC
  is not resolved.
- "The monotone structure is exact. This confirms that the channels
  contribute independently" is a non sequitur. Monotonicity follows
  from additivity of the slope and implies nothing about
  independence.
- The Introduction states that boundary cells "are validated
  numerically", and the Abstract says interior cells confirm "the
  analytical prediction". Section 4.2 states that validation has not
  passed.
- None of the four outstanding items in `CLAUDE.md` is complete: the
  common-driver boundary re-run, an analytical derivation, MCSE
  columns, and all three half-life panels. The table gives MCSE only
  as a per-row range and shows the $t_{1/2} = 1$ panel only.

*Remedy.* Rewrite Section 4.3 from the regenerated numbers. Report
all three half-life panels with a per-cell MCSE. Remove the Abstract
Results until validation passes.

*Evidence.* VERIFIED (`p11-01-summary.R`).

### M8. Cells use independent random streams, so the key contrasts are noisier than necessary

*Location.* Driver (`furrr_options(seed = 20260610L)` over cells);
Section 4.3.

*Problem.* `future_map_dfr` with a fixed `furrr` seed gives each grid
cell its own L'Ecuyer substream. Cells are independent rather than
paired, so every interaction contrast has a standard error near
0.02 to 0.025, the sum of four independent variances. This applies
to the benchmark excess, the channel comparisons and the boundary
identity. Common random numbers across the nine weight cells would
make the boundary reductions exact, draw for draw, and would sharpen
every contrast. The Reproducibility section also says the re-run
will "adopt the per-randomization-path sample-size convention" of
`01-dgp-prototype.R`. That driver now allocates a total, as does the
combined driver, so the sentence contradicts Section 4.2.

*Remedy.* Draw the multivariate normal deviates and the biomarker
once per replicate and apply all nine weight settings to them. Report
paired-difference MCSEs.

*Evidence.* INSPECTED.

## 5. Minor comments

1. *Style (STYLE.md section 5 grep).* There are 11 hits: "coupling"
   (lines 80, 192, 194, 264, 416, 451), "ceiling" and "sub-ceiling"
   (63, 112, 287, 289, 290), and "decoupling" (264). STYLE.md also
   maps "mean moderation" and "covariance moderation" to "the
   mean-shift encoding" and "the correlation encoding". NOTATION.md
   Part 3 retains the "-moderation architecture" names. Resolve the
   conflict once at the program level and follow it.
2. *Retired letters.* "A and B" in prose (lines 51, 102-103), "the B
   channel" and "the A channel" (419-421), and "the A-versus-B
   contrast" (464) use labels that NOTATION.md Part 3 retired. The
   letter $B$ is also the biomarker symbol ($B_i$, $\mu_B$,
   $\sigma_B$), so "the B channel" collides with it.
3. *Notation lint.* `perl tools/notation-lint.pl` reports one
   `bare-E-expect` at line 200 ($E[BR \mid B]$). Use $\mathbb{E}$.
4. *Reserved symbols.* $\eta$ is used for noncentrality but is
   reserved for design expectancy (use $\Lambda$). $\alpha$ is used
   for a mixing weight but is reserved for the significance level.
   $z_{B,i}$, $\mu_B$ and $\sigma_B$ should be $b_i$, $\mu_{bm}$ and
   $\sigma_{bm}$.
5. *Sign convention.* NOTATION.md Rule 5 is not stated. The fitted
   interaction coefficients are negative (about $-0.23$ per
   biomarker unit), and the reader needs that statement.
6. *Feasibility asymmetry.* Only $c_{bm,B}$ is bounded by positive
   definiteness. $c_{bm,A}$ is unbounded. State the largest feasible
   $c_{bm,B}$ for each design at the parameters used. At the July
   parameters the OL+BDC matrix at 0.45 had a smallest eigenvalue of
   only 0.144.
7. *Title and Section 2.2.* The title promises "mixing-weight
   sensitivity", but the mixing-weight parameterization is defined
   and then explicitly not simulated. Either simulate it, which is
   the natural design given M4 and M5, or remove it from the title.
8. *Name consistency.* "The Combined Architecture" (Section 2 title),
   "combined architecture" and "dual-channel architecture" alternate.
   Use one name.
9. *Missing Background citation.* Hendrickson et al. (2020) is in
   the bibliography but never cited. STYLE.md section 2 asks that it
   be cited once as the source of the design machinery and the
   correlation encoding.
10. *Reporting standard.* No ADEMP structure (Morris et al., 2019)
    and no Glossary (STYLE.md section 4). "E2" is used without
    definition.
11. *Table caption.* It begins in lower case ("the dual-channel
    architecture empirical power"), and its last half is revision
    history that belongs in a response letter.
12. *Section headings.* Numbers are hand-coded into `\section{}`
    while `number_sections: true` is set, so the PDF shows doubled
    numbering (INSPECTED in the source; the PDF was not opened).
    "Abstract" is also typeset as a numbered section.
13. *Inert leftovers.* `in_header: claudecode.tex` is an inert
    leftover (`CLAUDE.md`, "Paragraph structure"). `bullets.Rmd`
    exists in the directory and was not reviewed.
14. *Sentence length.* Several Abstract sentences exceed 30 words,
    against STYLE.md section 5. Semicolons join clauses at lines 303
    and 498, against STYLE.md section 1.
15. *Reproducibility statement.* It asserts R 4.6.0 in the Docker
    container. Whether the July data were produced in the container
    or locally is not recorded. BLAS differences make the two
    non-identical draw for draw (`CLAUDE.md`, "Simulation server").
    UNVERIFIED.
16. *Bibliography hygiene.*
    - The file header says "manuscript 01". It was copied from paper
      01.
    - 24 of 30 entries are uncited. Pandoc drops them, but the file
      should be pruned.
    - `arminger1999` lacks its DOI (10.1007/BF02294568), and
      `verbekelesaffre1996` lacks its DOI (10.2307/2291398). Both
      were verified by web search against publisher and index
      records.
    - The four `@unpublished` series entries carry no version or
      date, so readers cannot locate the cited claims.

## 6. Missing and questionable references

| Reference | Status | Why it is needed or questioned |
|---|---|---|
| Hendrickson et al. 2020, *Front Digit Health* 2:13 | In bib, uncited | Source of the correlation encoding and the Hybrid design |
| Arminger, Stein, Wittenberg 1999, *Psychometrika* 64:475 | Questionable | Finite mixtures; the present DGP is not a mixture (M5) |
| Verbeke and Lesaffre 1996, *JASA* 91:217 | Questionable | Mixture of random-effects laws; supports paper 03, not this DGP |
| Hedeker, Mermelstein, Demirtas 2008, *Biometrics* 64:627, doi 10.1111/j.1541-0420.2007.00924.x | Missing | Location-scale models: the standard way to put a covariate on a variance (M5 remedy b) |
| Brookes et al. 2004, *J Clin Epidemiol* 57:229, doi 10.1016/j.jclinepi.2003.08.009 | Missing | Simulation of interaction-test power; the fourfold rule |
| Kraemer et al. 2002, *Arch Gen Psychiatry* 59:877, doi 10.1001/archpsyc.59.10.877 | Missing | Moderator framework in psychiatric trials; context for the prazosin example |
| Senn 2016, *Stat Med* 35:966, doi 10.1002/sim.6739 | In bib, uncited | Variance components and HTE; what a covariance-encoded interaction implies for arm variances |
| Cortés et al. 2019, *F1000Research* (PMC6524747) | Missing | Variance-ratio evidence on constant effects. Volume and article number not verified |
| Winkelbeiner et al. 2019, *JAMA Psychiatry* 76:1063, doi 10.1001/jamapsychiatry.2019.1530 | Missing | Variance-ratio test of individual response in psychiatry |
| Plöderl and Hengartner 2019, *BMJ Open* 9:e034816, doi 10.1136/bmjopen-2019-034816 | Missing | Patient-by-treatment interaction via variance ratios; psychiatric context |
| Morris, White, Crowther 2019, *Stat Med* 38:2074, doi 10.1002/sim.8086 | In bib, uncited | ADEMP reporting for the simulation study |
| Kent et al. 2020, PATH statement, *Ann Intern Med* 172:35 | In bib, uncited | HTE framing and terminology |

PubMed records were retrieved for Hedeker et al., Brookes et al.,
Kraemer et al., Winkelbeiner et al., Burke et al. and Plöderl and
Hengartner. The DOIs above are from those records. Burke et al. 2015
(*BMJ* 351:h5651, doi 10.1136/bmj.h5651) is optional context on the
credibility of subgroup analyses. I found no published simulation
study that crosses mean- and covariance-encoded interactions in
N-of-1 designs. The search was limited (PubMed plus web search), so
absence is not established.

## 7. Suggestions for strengthening

1. Reframe the contribution around the off-drug interaction profile.
   The two weights parameterize how much of the interaction persists
   after discontinuation. That is identifiable in principle, through
   the off-drug visits and E3-type contrasts, and it matters for
   design choice.
2. Use common random numbers across weight cells (M8). The boundary
   identities then become exact unit tests that can live in
   `inst/tinytest/`.
3. Give the closed-form noncentrality for the GLS estimator under the
   true conditional covariance. This is straightforward from the
   matrices already built in `build_sigma_matrix()`. Show the
   simulated power against it for every cell. That is the analytical
   result the extraction plan asked for.
4. Simulate the mixing-weight line $c_{bm,A} + c_{bm,B} = c$, the
   equal-total parameterization. It isolates exactly the profile
   effect of M4, and it is the comparison readers will want.
5. If a second moment is wanted, add a location-scale channel and
   show whether a variance-function model or a CR2-robust test
   detects it. That would justify the "dual-channel" name.
6. Re-run under the extended DGP with G9 (or F1) and CR2 errors, with
   $t_{1/2}$ on the pharmacological and behavioral grids. Report null
   rejection by design.

## 8. Scope of this review

*Read in full:* the manuscript source with its code chunks; the
bibliography; `planning/plan-architecture-c.md`; paper 01's
`architecture-c-extraction-plan.md`; both dual-channel drivers;
`STYLE.md`; the relevant parts of `NOTATION.md` (Parts 1, 2 sample
size and 3); the combined branch of
`implementations/tidyverse/R/functions.R`; `simulate_cell()`,
`prepare_long_data()` and `fit_spec()` in `simulation-core.R`
(current and 2026-07-31 versions); paper 01's `01-dgp-prototype.R`
(design, allocation, parameters); `R/lme_analysis.R` (data assembly
and formula); and the combined branch of `R/generateData.R`. Not
reviewed: `bullets.Rmd`, the rendered PDF, paper 03, paper 06, the
full text of Arminger et al. and of Verbeke and Lesaffre (judged from
abstracts and indexing records only).

*Executed (scratchpad
`/private/tmp/claude-501/-Users-zenn-prj-res-36-pmsimstats-ng-pmsimstats-ng-analysis-report-01-dgp-mean-moderation-vs-mvn/513155aa-ec85-4a1d-889a-5d57e8607d18/scratchpad/`):*

- `p11-extract.sh`: extracted the 2026-07-31 sources (commit
  26b8c16a3).
- `p11-01-summary.R`: reproduced every number in Table
  `tab:arch-c-power` and Section 4.2 from the stored summary (exact
  match). It also produced the benchmark excesses, null rejection,
  mean estimates and SDs.
- `p11-02-pd.R`: rebuilt the covariance matrices for all 81 cells
  under the July and current parameterizations and checked the
  $D_{bc}$ coding of every path.
- `p11-03-blrow.R`: the baseline-row experiment. Six cells, 400
  paired replicates each, $t_{1/2} = 0$, $N = 70$, July code.
- `p11-04-analytic.R`: conditional moments from the covariance
  matrices, the probit-scale benchmark and the within-channel
  control.

*Not done:* no re-run of the grid under the current code (a long
simulation), so the corrected magnitudes are unknown. The Hybrid
residual in M2 was not localized. The paper 01 baseline numbers were
taken from `analysis/data/quick-sim/01-dgp-summary.txt` and not
regenerated. The probit benchmark in M3 is a normal approximation
that ignores the lower rejection tail and the conservatism of the
test. Literature searching used PubMed and general web search only.
No email or identifier was sent to any service.
