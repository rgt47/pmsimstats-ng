# Strict classical repeated-measures ANOVA for the biomarker-treatment
# interaction: closed-form derivation and limitations

*2026-05-07 06:30 PDT*

**Author.** pmsimstats team

## 1. Motivation

The founding pmsimstats question is the test of the
biomarker-treatment interaction in an aggregated N-of-1 trial
[@hendricksonOptimizing2020]. The pmsimstats-ng analysis layer
fits a continuous-time mixed-effects model
($\text{nlme::lme}(\text{Sx} \sim \text{bm} + t + \text{Dbc} +
\text{bm}{:}\text{Dbc}, \text{random} = \sim 1 \mid \text{ptID},
\text{correlation} = \text{corCAR1})$) and tests the interaction
coefficient via a Wald $t$-statistic. A reviewer trained in
classical analysis-of-variance methodology may reasonably ask
whether a strict repeated-measures ANOVA (RM-ANOVA) test of the
same interaction is available, what its closed-form expression
looks like, and how it relates to the mixed-model test the
framework currently produces. This white paper supplies that
answer for the strict classical case.

The presentation is deliberately confined to **strict** RM-ANOVA:
balanced split-plot design, categorical predictors, sphericity
assumption, method-of-moments estimation. Modern extensions (the
mixed-effects model for repeated measures, MMRM
[@mallinckrodt2008recommendations]; generalised estimating
equations, GEE [@liang1986longitudinal]; ANCOVA with random
subject effects) are tabulated for context but are not strict
RM-ANOVA in the textbook sense.

## 2. Notation and design

Consider the simplest design that supports a biomarker-treatment
interaction test in a within-subject setting:

- Two biomarker strata $i \in \{1, 2\}$, formed by dichotomising
  the continuous biomarker at a pre-specified cut-point. Strict
  RM-ANOVA requires categorical predictors; the dichotomisation
  is the price of admission to the framework
  [@royston2006dichotomizing; @altman2006prognosis;
  @sennMastering2016].
- $n$ subjects per stratum, balanced. Index $l \in \{1, \ldots, n\}$.
- Two treatment phases $j \in \{1, 2\}$ (active drug, placebo).
  In an N-of-1 design these are two phases through which every
  subject passes; in a parallel-group design they are between-
  subject and the design is not a split-plot.
- $T$ within-phase timepoints, indexed $k \in \{1, \ldots, T\}$.
  In the simplest case $T = 1$, the analyst pre-averages
  within-phase observations to a single mean per subject per
  phase.

The biomarker stratum is **between-subjects** (each subject sits
in exactly one stratum); the treatment phase and the timepoints
are **within-subjects** (each subject sees both phases at every
timepoint). This is the canonical *split-plot* repeated-measures
design [@maxwellDelaney2017designing].

## 3. The strict RM-ANOVA model

Write $Y_{ijlk}$ for the symptom outcome at biomarker stratum $i$,
treatment phase $j$, subject $l$ within stratum $i$, and timepoint
$k$. The strict classical split-plot model is

$$
\begin{aligned}
Y_{ijlk} = {} & \mu + \alpha_i + \beta_j + (\alpha\beta)_{ij}
+ \pi_{l(i)} \\
& + \tau_k + (\alpha\tau)_{ik} + (\beta\tau)_{jk}
+ (\alpha\beta\tau)_{ijk}
+ e_{ijlk} ,
\end{aligned} \qquad (1)
$$

with sum-to-zero constraints on each fixed effect, the random
subject effect $\pi_{l(i)} \sim N(0, \sigma_\pi^2)$ nested in
stratum, and within-subject error
$e_{ijlk} \sim N(0, \sigma_e^2)$. The model presupposes the
classical sphericity (compound symmetry) condition on the
within-subject covariance: every pair of within-subject
observations has the same covariance $\sigma_\pi^2$, every
observation has the same total variance
$\sigma_\pi^2 + \sigma_e^2$.

The biomarker-treatment interaction is the
$(\alpha\beta)_{ij}$ term. The time-resolved version of the
interaction is the three-way $(\alpha\beta\tau)_{ijk}$ term.

## 4. Sum-of-squares decomposition

The total sum of squares partitions into between-subjects and
within-subjects components:

$$
SS_{\text{tot}} = SS_{\text{between}} + SS_{\text{within}} . \qquad (2)
$$

The between-subjects component decomposes further into
$SS_A$ (biomarker stratum) and $SS_{S(A)}$ (subject within
stratum):

$$
SS_A = bnT \sum_{i=1}^a (\bar Y_{i\cdots} - \bar Y_{\cdots\cdot})^2 ,
\qquad (3)
$$

$$
SS_{S(A)} = bT \sum_{i=1}^a \sum_{l=1}^n
(\bar Y_{i\cdot l\cdot} - \bar Y_{i\cdots})^2 . \qquad (4)
$$

The within-subjects component decomposes into $SS_B$ (treatment
phase), $SS_{AB}$ (the **biomarker × treatment interaction**, the
target of the test), and a within-subjects error term
$SS_{B \times S(A)}$:

$$
SS_B = anT \sum_{j=1}^b
(\bar Y_{\cdot j\cdot\cdot} - \bar Y_{\cdots\cdot})^2 , \qquad (5)
$$

$$
SS_{AB} = nT \sum_{i,j}
(\bar Y_{ij\cdot\cdot} - \bar Y_{i\cdots} -
\bar Y_{\cdot j\cdot\cdot} + \bar Y_{\cdots\cdot})^2 , \qquad (6)
$$

$$
SS_{B \times S(A)} = T \sum_{i,j,l}
(\bar Y_{ijl\cdot} - \bar Y_{i\cdot l\cdot} -
\bar Y_{ij\cdot\cdot} + \bar Y_{i\cdots})^2 . \qquad (7)
$$

The corresponding degrees of freedom are
$df_A = a - 1$, $df_{S(A)} = a(n-1)$, $df_B = b - 1$,
$df_{AB} = (a-1)(b-1)$, and
$df_{B \times S(A)} = a(n-1)(b-1)$ (the latter being the
within-subjects error degrees of freedom appropriate for the
treatment effect and its interaction with the between-subjects
factor).

## 5. The test statistic: closed form

The biomarker-treatment interaction is tested by the $F$-ratio

$$
F_{AB} = \frac{MS_{AB}}{MS_{B \times S(A)}}
= \frac{SS_{AB}/df_{AB}}{SS_{B \times S(A)}/df_{B \times S(A)}} .
\qquad (8)
$$

Under the null $H_0: (\alpha\beta)_{ij} \equiv 0$ and the
sphericity assumption,
$F_{AB} \sim F(df_{AB}, df_{B \times S(A)})$.

### 5a. The $2 \times 2$ split-plot case

Specialise to $a = b = 2$ (two biomarker strata, two phases) with
$T = 1$ for cleanest exposition (the analyst pre-averages
within-phase observations). The sample interaction contrast is

$$
\hat L \;=\; (\bar Y_{11\cdot} - \bar Y_{12\cdot}) -
            (\bar Y_{21\cdot} - \bar Y_{22\cdot}) , \qquad (9)
$$

i.e.\ the difference of within-subject treatment differences
across biomarker strata. This is the natural biomarker-treatment
interaction statistic.

The interaction parameter under the sum-to-zero ANOVA
parameterisation satisfies
$\hat{(\alpha\beta)}_{11} = \hat L / 4$ (and the other three cells
take values $\pm \hat L / 4$), so the interaction sum of squares
is

$$
SS_{AB} \;=\; nT \cdot 4 \cdot (\hat L / 4)^2 \;=\; \frac{nT\,\hat L^2}{4} .
\qquad (10)
$$

With $df_{AB} = 1$, the interaction mean square equals the
interaction sum of squares.

The within-subjects error mean square is the pooled estimator of
$\sigma_e^2$ from the treatment-by-subject-within-stratum
component of variance. Under the sphericity assumption and
balanced data, $MS_{B \times S(A)} = \hat\sigma_e^2$.

Substituting into equation (8):

$$
\boxed{\,F_{AB} \;=\; \frac{nT\,\hat L^2}{4\,\hat\sigma_e^2}\,} ,
\qquad df_1 = 1, \quad df_2 = 2(n-1) . \qquad (11)
$$

This is the closed-form expression for the strict classical
RM-ANOVA test statistic for the biomarker-treatment interaction
in a $2 \times 2 \times T$ split-plot design.

### 5b. Within-subject contrast representation

Equation (11) is equivalent to a two-sample $t$-test on the
within-subject difference. Define
$\Delta_{il} = \bar Y_{i1l\cdot} - \bar Y_{i2l\cdot}$ as the
within-subject treatment difference for subject $l$ in biomarker
stratum $i$, averaged over $T$ within-phase timepoints. Then

$$
\hat L \;=\; \bar\Delta_1 - \bar\Delta_2 , \qquad
\bar\Delta_i \;=\; \frac{1}{n} \sum_{l=1}^n \Delta_{il} . \qquad (12)
$$

Under sphericity, $\Delta_{il}$ are iid within stratum with
$\text{Var}(\Delta_{il}) = 2\sigma_e^2 / T$ (the random subject
effect cancels in the within-subject difference). The variance of
$\hat L$ is therefore

$$
\text{Var}(\hat L) \;=\; \frac{2 \cdot 2\sigma_e^2}{nT}
\;=\; \frac{4\sigma_e^2}{nT} . \qquad (13)
$$

The Welch-style $t$-statistic for the contrast is

$$
t_{AB} \;=\; \frac{\hat L}{\sqrt{4\hat\sigma_e^2/(nT)}}
\;=\; \sqrt{\frac{nT}{4}} \cdot \frac{\hat L}{\hat\sigma_e} ,
\qquad (14)
$$

and $t_{AB}^2 = F_{AB}$, recovering equation (11). The strict
classical RM-ANOVA $F$-test for the biomarker-treatment
interaction is therefore identical to a two-sample $t$-test on
the within-subject treatment differences, with degrees of freedom
$2(n-1)$. This identity is well known in the split-plot
literature [@maxwellDelaney2017designing, ch.\ 12;
@kirk2013experimental] and is the analytical anchor for the
statement that the strict-RM-ANOVA test of a biomarker-treatment
interaction reduces to a simple two-sample comparison of
treatment-effect estimates between biomarker strata.

### 5c. Generalisation: $a$ strata, $b$ phases

For $a$ biomarker strata, $b$ treatment phases, $n$ subjects per
stratum, and $T$ within-phase timepoints, equation (8) generalises
to

$$
F_{AB} = \frac{SS_{AB}/[(a-1)(b-1)]}{SS_{B \times S(A)}/[a(n-1)(b-1)]} .
\qquad (15)
$$

Under sphericity,
$F_{AB} \sim F((a-1)(b-1), a(n-1)(b-1))$ under $H_0$. The closed
form of the interaction sum of squares follows directly from
equation (6).

### 5d. The time-resolved three-way case

Including the time factor $\tau_k$ explicitly, the test of
whether the biomarker modifies the treatment effect across time is
the three-way interaction $(\alpha\beta\tau)_{ijk}$. The
$F$-statistic is

$$
F_{AB\tau} = \frac{MS_{AB\tau}}{MS_{B\tau \times S(A)}} ,
\qquad (16)
$$

with $df_1 = (a-1)(b-1)(T-1)$ and $df_2 = a(n-1)(b-1)(T-1)$. The
sums of squares follow the standard split-plot pattern
[@maxwellDelaney2017designing, ch.\ 14].

## 6. Sphericity, Mauchly, and the Greenhouse-Geisser correction

Equations (11) and (15) presuppose sphericity: the within-subject
covariance matrix is compound-symmetric. When this assumption
fails, the nominal $F$-distribution is no longer exact under
$H_0$, and the test inflates Type I error. Two diagnostics and
two corrections are standard:

- **Mauchly's test** [@mauchly1940significance]. Likelihood-ratio
  test for compound symmetry against an unstructured alternative.
  Often rejects in moderate samples; not informative about the
  direction of departure.
- **Box's $\epsilon$** [@box1954someproblemsI;
  @box1954someproblemsII]. A scalar measure of departure from
  sphericity, with $\epsilon = 1$ at compound symmetry and
  $\epsilon \geq 1/(T-1)$ at maximum departure.
- **Greenhouse-Geisser correction**
  [@greenhouseGeisser1959]. Estimate $\hat\epsilon$ from the
  observed within-subject covariance and adjust degrees of
  freedom: $F_{AB} \sim F(\hat\epsilon \cdot df_1,
  \hat\epsilon \cdot df_2)$ under $H_0$. Conservative.
- **Huynh-Feldt correction** [@huynhFeldt1976]. Less
  conservative variant, recommended when $\hat\epsilon > 0.75$.

The pmsimstats DGP imposes AR(1) within-factor correlation, which
is **not** compound-symmetric. Under AR(1), Box's $\epsilon$ is
strictly below 1, and the Greenhouse-Geisser correction
substantially attenuates the available degrees of freedom in
$F$-tests on within-subject effects. In simulation, the strict-
RM-ANOVA Type I error in the present setting is documented at
13--17 percent for the CO and OL designs absent correction
(audit at `docs/06-ar1-residual-correlation.tex`); after
Greenhouse-Geisser correction the Type I error returns toward
nominal but at the cost of substantial power. The mixed-model
analysis with `corCAR1` does not pay this cost because it models
the AR(1) structure directly.

## 7. Power calculation under the closed form

For the $2 \times 2$ case, equation (11) yields a clean
closed-form power formula. Under the alternative
$H_1: \hat L = \delta$, the statistic
$F_{AB} = nT\,\hat L^2 / (4\hat\sigma_e^2)$ follows a non-central
$F$ distribution with non-centrality parameter

$$
\lambda \;=\; \frac{nT\,\delta^2}{4\sigma_e^2} . \qquad (17)
$$

Power at level $\alpha$ is

$$
1 - \beta \;=\; \Pr\!\left[ F_{AB} > F_{1, 2(n-1), 1-\alpha}
\;\Big|\; \lambda \right] , \qquad (18)
$$

which is computable via standard non-central $F$ functions
(`stats::pf` in R with the `ncp` argument; `PROC POWER` in SAS).
For sample-size calculation, equation (17) inverts to

$$
n \;=\; \frac{4\sigma_e^2 \lambda}{T \delta^2} , \qquad (19)
$$

with $\lambda$ chosen to deliver the target power at the design
$\alpha$.

This formula assumes (i) sphericity, (ii) categorical biomarker
at the dichotomised cut-point with equal stratum sizes, (iii) no
carryover, (iv) no within-treatment relapse, (v) iid within-
phase observations after within-subject differencing. All five
assumptions are violated to varying degrees in the pmsimstats
DGP; see §8.

## 8. Limitations of the strict classical formulation

The closed-form simplicity of equation (11) is bought with five
strong assumptions, each of which is empirically problematic in
the pmsimstats setting.

1. **Categorical biomarker.** The pmsimstats biomarker (resting
   SBP in the prazosin-PTSD application) is continuous and has
   no defensible cut-point. Dichotomising at the median or at a
   clinical threshold loses the within-stratum biomarker
   variation, attenuates the interaction estimator, and inflates
   Type II error
   [@royston2006dichotomizing; @altman2006prognosis;
   @sennMastering2016]. The information loss can be substantial:
   in the carryover-sensitivity production data the median split
   gives a $\hat L$ standard error roughly 1.4 to 1.7 times that
   of the continuous-biomarker mixed-model coefficient,
   translating to a 30 to 50 percent power penalty.

2. **Sphericity.** Compound-symmetric within-subject covariance
   is empirically rejected for time-series symptom data: nearby
   timepoints are more correlated than distant ones. The AR(1)
   structure that pmsimstats implements explicitly violates
   sphericity. After Greenhouse-Geisser correction, the test
   recovers nominal Type I error but loses degrees of freedom and
   power.

3. **No carryover model.** Strict RM-ANOVA has no provision for
   residual drug effect during washout; the model treats off-drug
   observations as if no drug had been administered. Under the
   pmsimstats DGP with $t_{1/2} = 1.0$ week, this misspecification
   biases $\hat L$ toward zero (attenuation by 30 to 35 percent
   in the production data, comparable to the A1 specification
   evaluated in the carryover-sensitivity manuscript).

4. **Balanced design required.** Equations (11) and (15) assume
   equal stratum sizes and equal numbers of within-phase
   timepoints. Real trials have missing data; classical RM-ANOVA
   handles missingness only through complete-case deletion. The
   MMRM framework's likelihood-based handling of MAR data
   [@mallinckrodt2008recommendations] is genuinely better-behaved.

5. **Equal time spacing.** Strict RM-ANOVA treats time as a
   categorical factor with implicitly equal-spaced levels. Real
   N-of-1 trials have unequally-spaced visits; pmsimstats
   implements continuous-time AR(1) (`corCAR1`), which the
   strict-RM-ANOVA framework cannot do.

## 9. Relationship to the pmsimstats MMRM analysis

The strict-RM-ANOVA $F$-statistic in equation (11) and the
pmsimstats `nlme::lme` Wald $t$-statistic on `bm:Dbc` are
**asymptotically equivalent** under three conditions: (i) the
biomarker is dichotomised at the cut-point that defines the
strata, (ii) the within-subject covariance is exactly compound-
symmetric, and (iii) `Dbc` is binary (the A1 specification rather
than the matched-decay A2 specification). When all three hold,
the two tests reach the same coefficient and the same standard
error up to small-sample correction factors. When any fails, the
mixed-model test extracts more information than the strict
RM-ANOVA: the continuous biomarker contributes within-stratum
variation, the AR(1) covariance is modelled directly, and the
continuous `Dbc` exposure-decay predictor (A2) leverages the
off-drug trajectory.

For the founding pmsimstats interaction question
[@hendricksonOptimizing2020], the strict-RM-ANOVA test is
therefore a **conservative reduction** of the available analysis
rather than a competing alternative. Where a reviewer requests
"RM-ANOVA results", the natural response is one of (a) point at
the existing MMRM analysis as a continuous-time variant of the
RM-ANOVA family, (b) report the strict-RM-ANOVA test from
equation (11) as a sensitivity comparator under dichotomised
biomarker and binary phase, or (c) report both with the
explicit power-loss tabulation.

## 10. Recommendation for pmsimstats

We recommend retaining the MMRM-with-continuous-time analysis as
primary and adding a strict-RM-ANOVA comparator as a planned
sensitivity:

1. **Categorise the biomarker at the median.** Conventional
   choice; the alternative pre-specified clinical cut-point is
   acceptable when one exists.
2. **Pre-average within-phase timepoints to a single mean per
   subject per phase.** Reduces the $T$-axis to $T = 1$ for the
   sensitivity test.
3. **Compute $\hat L$ per equation (9) and $F_{AB}$ per equation
   (11).** This is computable by hand from the four cell means
   and the within-subject error mean square.
4. **Report the strict-RM-ANOVA $p$-value alongside the MMRM
   $p$-value in the manuscript table.** The two should be in the
   same order of magnitude. Substantial divergence between them
   indicates that the dichotomisation, the sphericity assumption,
   or the binary phase indicator is doing meaningful work, and
   the discrepancy is informative for interpretation.

The carryover-sensitivity manuscript at
`analysis/report/02-carryover-sensitivity/` already evaluates the A1
binary-phase specification, which is the analysis-model
counterpart to the strict-RM-ANOVA closed form derived here. The
manuscript reports the A1 specification's bias and coverage under
carryover. The strict-RM-ANOVA test inherits the A1 specification's
attenuation properties.

## References

- [@gompertz1825nature] Gompertz B. On the nature of the function
  expressive of the law of human mortality. *Philosophical
  Transactions of the Royal Society of London*. 1825;115:513-583.
- [@maxwellDelaney2017designing] Maxwell SE, Delaney HD, Kelley K.
  *Designing Experiments and Analyzing Data: A Model Comparison
  Perspective*. 3rd ed. Routledge; 2017.
- [@kirk2013experimental] Kirk RE. *Experimental Design:
  Procedures for the Behavioral Sciences*. 4th ed. SAGE; 2013.
- [@mauchly1940significance] Mauchly JW. Significance test for
  sphericity of a normal $n$-variate distribution. *Annals of
  Mathematical Statistics*. 1940;11(2):204-209.
- [@box1954someproblemsI] Box GEP. Some theorems on quadratic
  forms applied in the study of analysis of variance problems, I.
  *Annals of Mathematical Statistics*. 1954;25(2):290-302.
- [@box1954someproblemsII] Box GEP. Some theorems on quadratic
  forms applied in the study of analysis of variance problems,
  II. *Annals of Mathematical Statistics*. 1954;25(3):484-498.
- [@greenhouseGeisser1959] Greenhouse SW, Geisser S. On methods
  in the analysis of profile data. *Psychometrika*.
  1959;24(2):95-112.
- [@huynhFeldt1976] Huynh H, Feldt LS. Estimation of the box
  correction for degrees of freedom from sample data in
  randomised block and split-plot designs. *Journal of
  Educational Statistics*. 1976;1(1):69-82.
- [@laird1982random] Laird NM, Ware JH. Random-effects models for
  longitudinal data. *Biometrics*. 1982;38(4):963-974.
- [@verbeke2000linear] Verbeke G, Molenberghs G. *Linear Mixed
  Models for Longitudinal Data*. Springer; 2000.
- [@mallinckrodt2008recommendations] Mallinckrodt CH, Lane PW,
  Schnell D, Peng Y, Mancuso JP. Recommendations for the
  primary analysis of continuous endpoints in longitudinal
  clinical trials. *Drug Information Journal*. 2008;42(4):
  303-319.
- [@liang1986longitudinal] Liang K-Y, Zeger SL. Longitudinal
  data analysis using generalised linear models. *Biometrika*.
  1986;73(1):13-22.
- [@royston2006dichotomizing] Royston P, Altman DG, Sauerbrei W.
  Dichotomising continuous predictors in multiple regression:
  a bad idea. *Statistics in Medicine*. 2006;25(1):127-141.
- [@altman2006prognosis] Altman DG. The cost of dichotomising
  continuous variables. *British Medical Journal*. 2006;332:1080.
- [@sennMastering2016] Senn S. Mastering variation: variance
  components and personalised medicine. *Statistics in Medicine*.
  2016;35(7):966-977.
- [@hendricksonOptimizing2020] Hendrickson RC, Thomas RG,
  Schork NJ, Raskind MA. Optimizing aggregated N-of-1 trial
  designs for predictive biomarker validation: statistical
  methods and practical considerations. *Frontiers in Digital
  Health*. 2020;2:13.

---

*Rendered on 2026-05-07 at 06:30 PDT.*<br>
*Source: ~/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/26-rm-anova-interaction-test.md*
