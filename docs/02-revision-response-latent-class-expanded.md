---
title: "Expanded response: latent-class faithfulness versus MVN approximation"
author: "pmsimstats team"
date: "2026-04-15"
---

# Expanded §4.2 response: latent-class faithfulness versus MVN approximation

This note expands the Pass 2 edit applied to §4.2 of
`02-dgp-mean-moderation-vs-mvn.tex`, which addressed the reviewer's
question in the email summary: if Architecture B is really an
"imperfect responder indicator", why not model that explicitly at the
participant level rather than via covariance structure?

The short paragraph added in Pass 2 acknowledges the tension. A fuller
treatment develops five threads.

## 1. What the faithful generative model actually looks like

If blood pressure (or any candidate predictive biomarker) genuinely
indexes an unobserved responder/non-responder dichotomy, the
mechanistically correct DGP is a finite mixture:

$$
Z_i \mid B_i \sim \mathrm{Bernoulli}\bigl(\pi(B_i)\bigr), \qquad
BR_{it} \mid Z_i = z \sim f_z(t, D_{it}),
$$

where $Z_i$ is the latent class (responder vs. non-responder),
$\pi(B_i) = \Pr(Z_i = 1 \mid B_i)$ is typically modelled as a logistic
function of the biomarker, and the class-specific response
distributions $f_0, f_1$ differ in their drug-response parameters
(peak response, onset rate, or both). Carryover under this DGP
operates only on the class-specific means and does not erode any
correlation signal, because class membership $Z_i$ is
drug-state-invariant.

## 2. Why MVN differential correlation is nonetheless a defensible second-moment summary

A two-component mixture with class-dependent means and a common
within-class covariance produces a marginal joint distribution of
$(B_i, BR_{it})$ whose first two moments match, to leading order,
those of a single MVN with $\mathrm{Cor}(B, BR) = c_{bm}$ chosen so
that $c_{bm}^2$ equals the between-class variance fraction. Under
mild regularity (class overlap, moderate $\pi$), the mixture's
second-moment structure is well-approximated by an MVN with
differential correlation in treatment states where the class-specific
means separate and shrinkage in treatment states where they do not.
Architecture B as implemented therefore captures the covariance
implications of the latent-class mechanism without committing to the
discrete generative form.

## 3. Where the approximation breaks down

The MVN approximation will be misleading in three regimes:

- **Bimodal responses.** When $f_0$ and $f_1$ are well-separated
  (e.g., non-responders show near-zero drug effect and responders show
  a large effect), the marginal distribution of $BR$ is bimodal and no
  MVN matches its tails. Power under the true mixture DGP will be
  higher than Architecture B predicts because a well-specified mixture
  analysis can exploit the bimodality.
- **Strong class-membership gradient.** When $\pi(B_i)$ is steep
  (approaching a step function at some threshold), the biomarker
  behaves almost deterministically as a class label and Architecture A
  with a thresholded or dose-response moderation becomes the better
  approximation.
- **Class-varying covariance.** If the autocorrelation, residual
  variance, or placebo-response parameters themselves differ between
  classes, no single-component MVN can represent the generative
  structure and analysis models assuming constant covariance will be
  misspecified.

## 4. Implications for the power estimates in Section 3

The substantial Architecture B power loss under carryover (40-60%)
reflects erosion of a second-moment signal. A mixture-model analysis
fitted to mixture-DGP data would, in principle, recover additional
power by estimating class membership directly, because the class label
is carryover-invariant. The gap between Architecture B's
carryover-sensitive estimate and the mixture-faithful estimate is an
upper bound on what a mixture analysis could recover; it is not
currently known whether this gap is large. A focused simulation
contrasting (i) mixture-DGP with MVN analysis, (ii) mixture-DGP with
mixture analysis, and (iii) MVN-DGP with both analyses would
characterise the loss attributable to second-moment approximation
versus the loss attributable to analysis-model misspecification.

## 5. Practical recommendation and scope of the open question

For the prazosin-PTSD application, mixture modelling is attractive in
principle but carries identifiability costs that are acute in N-of-1
settings with modest participant counts: the EM algorithm for finite
mixtures of mixed-effects models requires either informative priors on
class-membership probabilities or strong separation between
class-specific response curves, neither of which is assured in a trial
designed to detect an interaction whose existence is itself uncertain.

The open question, therefore, is not whether mixture modelling is more
faithful (it is), but whether the fidelity gain is recoverable at
realistic sample sizes, or whether analysis-strategy mitigations
(Section 6: restricted analysis, weighting, within-subject contrasts)
deliver a better power-per-complexity trade. Resolving this question
empirically is out of scope for the present white paper but is a
natural next step in the pmsimstats programme.

## 6. Connection to the psychometric latent variable literature

The reviewer's question, translated into psychometric vocabulary,
asks whether biomarker-moderated treatment response is better modelled
as a continuous-trait phenomenon (classical factor-analytic / latent
regression tradition) or as a discrete-class phenomenon (latent class
/ mixture tradition). The psychometric literature has worked this
distinction over for more than fifty years, originally in the context
of individual differences in education, personality, and
psychopathology, and the taxonomy it developed maps directly onto the
Architecture A versus Architecture B distinction drawn in this paper.

The core generative question is whether unobserved heterogeneity in
outcome distributions is best represented by (i) a single population
with continuous individual-differences variation along one or more
latent dimensions (factor-analytic models; Spearman 1904;
Lawley & Maxwell 1971), (ii) a small number of qualitatively distinct
subpopulations with within-class homogeneity (latent class models;
Lazarsfeld & Henry 1968; Goodman 1974), or (iii) a hybrid in which
classes differ along continuous dimensions within each class (factor
mixture models; Lubke & Muthén 2005; Clark et al. 2013). The same
trichotomy appears throughout the clinical trial simulation
literature, usually implicitly: Architecture A encodes the continuous
individual-differences view (treatment response is a deterministic
function of a continuous biomarker), while Architecture B as
implemented in Hendrickson et al. (2020) is most naturally interpreted
as a second-moment approximation to the latent-class view.

## 7. Finite mixture taxonomy relevant to biomarker moderation

Four families of finite mixture model from the psychometric and
econometric literatures bear directly on the biomarker-treatment
interaction problem.

### 7.1 Latent class and latent profile analysis

Latent class analysis (LCA; categorical indicators) and latent profile
analysis (LPA; continuous indicators) assume a small number of
discrete subpopulations, each with its own marginal distribution over
observed variables. Class membership is itself unobserved and is
inferred from the pattern of covariation among indicators. In the
biomarker-trial context, an LPA-style model would posit two (or more)
responder classes with class-specific marginal distributions of
on-drug and off-drug response, and estimate class-membership
probabilities from the observed response trajectories. Classical
references include Lazarsfeld & Henry (1968), Goodman (1974), Clogg
(1995), and McLachlan & Peel (2000); the Vermunt & Magidson Latent
GOLD software and the `poLCA` R package implement the standard EM
estimation routines.

### 7.2 Growth mixture models and group-based trajectory models

Growth mixture models (GMM; Muthén & Shedden 1999; Muthén 2001, 2004)
extend LPA to longitudinal data by specifying class-specific growth
(or, in the trial setting, drug-response) curves with within-class
random-effects variation. Group-based trajectory models (Nagin 1999,
2005) are the degenerate case in which within-class random-effect
variances are set to zero. GMM is the closest psychometric analogue
of the mixture DGP discussed in §4.2: each participant belongs to a
latent class with its own drug-response trajectory, and the class
label is drug-state-invariant. Bauer & Curran (2003) and Bauer
(2007) document the identifiability and overextraction hazards that
accompany GMM estimation at realistic sample sizes; these cautions
apply with full force to N-of-1 trials, where the number of
participants is typically well below the sample sizes at which GMM
class-count recovery has been validated.

### 7.3 Factor mixture models

Factor mixture models (FMM; Yung 1997; Arminger, Stein & Wittenberg
1999; Lubke & Muthén 2005) combine continuous latent factors with
discrete latent classes, yielding a generative structure in which
within-class variation is governed by a factor model and between-class
variation is governed by class-specific means and (optionally)
class-specific factor loadings or residual covariances. The FMM
framework subsumes both the pure continuous (single-class factor
model) and the pure discrete (LCA with no within-class structure)
cases as limits. For the biomarker-moderation problem, FMM provides
a principled way to represent a biomarker as both a noisy class
indicator (via class-membership probability) and a continuous
modulator of within-class response (via factor loadings on the
biological-response factor).

### 7.4 Regression mixtures and mixtures of experts

Regression mixtures (DeSarbo & Cron 1988; Wedel & DeSarbo 1995;
Jedidi, Jagpal & DeSarbo 1997) and the closely related mixture of
experts architecture (Jacobs, Jordan, Nowlan & Hinton 1991) allow
class-specific regression coefficients with class-membership
probabilities that depend on covariates (the "gating" function in
mixture-of-experts terminology). In the trial-simulation idiom, this
corresponds to a DGP in which the biomarker governs both the
probability of being a responder and the magnitude of the
within-responder drug effect, which is precisely the generative form
the reviewer's question invites.

## 8. The Architecture B spectrum: covariance, mean, and combined biomarker moderation

Architecture B as implemented in Hendrickson et al. (2020) and
pmsimstats specifies biomarker-dependent covariance structure with
population means held constant across biomarker values. This is one
point on a broader spectrum of biomarker-moderated DGPs recognised in
the psychometric literature; the spectrum is usefully organised by
which moments of the joint distribution depend on the biomarker.

### 8.1 Covariance-only moderation (Architecture B as implemented)

The Hendrickson-Schork specification holds $E[BR_{it}]$ fixed across
biomarker strata and encodes the interaction entirely in
$\mathrm{Cov}(B_i, BR_{it} \mid D_{it})$. This is an unusual
restriction within the psychometric mixture tradition: most standard
FMM and GMM parameterisations permit class-specific means at a
minimum, and the means-constant-covariance-varies case is typically
encountered only as a constrained submodel used for testing invariance
hypotheses (Meredith 1993; Widaman & Reise 1997). Its appeal in the
Hendrickson context is tractability for power simulation under
multivariate normality; its cost, as the reviewer observes, is that
it does not correspond cleanly to any natural generative mechanism.

### 8.2 Combined mean and covariance moderation

Arminger, Stein & Wittenberg (1999), in the canonical Psychometrika
treatment of mixtures with covariate-dependent structure, develop
finite mixtures of conditional mean- and covariance-structure models
in which both the expected-response vector and the residual covariance
matrix may differ between classes and may depend on observed
covariates within class. This is the general psychometric form of
Architecture B and reduces to Architecture A in the limit of a single
class with covariate-dependent mean only, and to the
Hendrickson-implemented submodel in the limit of class-invariant means
with covariate-dependent covariance only. Verbeke & Lesaffre (1996)
provide the analogous treatment for linear mixed-effects models with
heterogeneity in the random-effects population, which is the more
direct generalisation of the repeated-measures structure used in
pmsimstats.

### 8.3 Location-scale moderation (GAMLSS and related)

Generalised additive models for location, scale, and shape (GAMLSS;
Rigby & Stasinopoulos 2005) and the related distributional-regression
tradition (Klein et al. 2015) allow the biomarker to enter multiple
parameters of the outcome distribution simultaneously through separate
link functions. In the biomarker-trial context, this yields a DGP in
which high-biomarker participants may have both a shifted mean drug
response (location effect) and altered response variability (scale
effect), without requiring discrete class structure. The GAMLSS
framework is continuous-heterogeneity analogue of FMM and occupies an
intermediate position between pure Architecture A (mean only) and
pure Architecture B (covariance only).

### 8.4 Heterogeneous random-effects models

Verbeke & Lesaffre (1996), Zhang & Davidian (2001), and
Proust-Lima et al. (2017) develop linear mixed-effects models in which
the random-effects distribution is itself a finite mixture, allowing
class-specific random-intercept and random-slope distributions.
Applied to aggregated N-of-1 trial data, this framework represents
responder heterogeneity directly at the participant-specific random
treatment effect, which is arguably the most clinically transparent
representation: responders and non-responders have distinct
random-slope distributions on the treatment indicator, and the
biomarker predicts class membership rather than entering the
analysis-model mean structure. The `lcmm` R package (Proust-Lima
et al. 2017) implements this family.

## 9. Implications for the pmsimstats framework

Three implications follow from locating Architecture B within the
broader psychometric spectrum.

**First**, the covariance-only parameterisation adopted in
Hendrickson et al. (2020) is the most restrictive operational case of
a much richer family. The substantial power loss under carryover
documented in Section 3 is a specific property of this restrictive
case and may not generalise to the mean-plus-covariance or
heterogeneous-random-effects variants that are standard in
psychometrics. A complete sensitivity analysis for the prazosin-PTSD
application should accordingly report power under at least three
points on the spectrum: Architecture A (mean-only), Architecture B as
implemented (covariance-only), and a combined mean-plus-covariance
specification of the Arminger et al. (1999) form.

**Second**, the identifiability concerns documented in the
psychometric GMM and FMM literatures (Bauer & Curran 2003; Bauer 2007;
Lubke & Muthén 2005) bound what can be asked of mixture-based analysis
models at N-of-1 trial sample sizes. The published FMM applications
typically involve several hundred to several thousand participants;
biomarker-trial applications with $N \in [30, 100]$ are well outside
this range. Analysis-strategy mitigation (Section 6 of the main
document) is therefore more likely than mixture-model analysis to
deliver recoverable power gains at realistic trial sizes.

**Third**, if the pmsimstats programme pursues a richer DGP in
subsequent work, the natural target is the heterogeneous-random-slopes
form (§8.4) rather than a full FMM. The random-slopes
parameterisation preserves the analysis model's compatibility with
existing N-of-1 estimation machinery (`nlme::lme` with
participant-specific random treatment effects) while admitting the
biologically faithful latent-class generative structure at the DGP
level. This would preserve the audit chain of the current framework
and permit direct comparison with the covariance-only baseline under
carryover.

## 10. Software and implementation references

R packages relevant to estimating the models surveyed above, in
approximate order of fit to the pmsimstats use case:

- `lcmm` (Proust-Lima et al. 2017). Joint latent class mixed models
  for longitudinal and survival data; the closest available tooling
  for heterogeneous random-effects specification (§8.4).
- `flexmix` (Leisch 2004; Grün & Leisch 2008). General regression
  mixture framework with user-definable component models; suitable
  for regression-mixture and mixture-of-experts specifications (§7.4).
- `mclust` (Scrucca et al. 2016). Gaussian finite mixture models with
  a range of class-specific covariance parameterisations; supports the
  Arminger et al. (1999) mean-plus-covariance specifications
  non-longitudinally.
- `gamlss` (Rigby & Stasinopoulos 2005). Distributional regression
  with covariate-dependent location, scale, and shape parameters
  (§8.3).
- `poLCA` (Linzer & Lewis 2011). Polytomous latent class analysis for
  the categorical-indicator case.
- `OpenMx` and `lavaan` (with mixture extensions). Structural
  equation mixture modelling for FMM-style specifications; more
  general but with steeper configuration cost than `lcmm` or
  `flexmix`.

The psychometric reference implementation outside R is Mplus (Muthén &
Muthén 1998-2017), which remains the standard for GMM, FMM, and
related mixture-SEM specifications.

## Relation to the existing §4.2 paragraph

The five threads in Sections 1-5 above could extend the §4.2 Pass 2
paragraph in place, or be promoted to a new §4.3 ("Mixture-modelling
alternative to MVN differential correlation") so that §4.2 retains
its current clinical-examples focus. Either approach preserves the
document's existing architecture-comparison narrative while giving
the reviewer's concern the technical treatment it warrants.

The psychometric-connection material in Sections 6-10 is more
naturally placed either as a new §4.4 ("Architecture B in the
psychometric latent-variable tradition") or as an appendix. It
answers a distinct question from the §4.2 edit: whereas §4.2 concerns
the biological faithfulness of Architecture B, Sections 6-10 concern
its methodological provenance and the spectrum of related
specifications already developed in the psychometric literature.
Treating the two as separate document additions preserves the logical
structure of Section 4 (biological assumptions) while giving the
methodological-provenance discussion the space it requires.

## References (additions beyond the main document bibliography)

Arminger G, Stein P, Wittenberg J. Mixtures of conditional mean- and
covariance-structure models. *Psychometrika*. 1999;64(4):475-494.

Bauer DJ. Observations on the use of growth mixture models in
psychological research. *Multivariate Behav Res*. 2007;42(4):757-786.

Bauer DJ, Curran PJ. Distributional assumptions of growth mixture
models: implications for overextraction of latent trajectory classes.
*Psychol Methods*. 2003;8(3):338-363.

Clark SL, Muthén B, Kaprio J, et al. Models and strategies for factor
mixture analysis: an example concerning the structure underlying
psychological disorders. *Struct Equ Modeling*. 2013;20(4):681-703.

Clogg CC. Latent class models. In: Arminger G, Clogg CC, Sobel ME,
eds. *Handbook of Statistical Modeling for the Social and Behavioral
Sciences*. Plenum; 1995:311-359.

DeSarbo WS, Cron WL. A maximum likelihood methodology for
clusterwise linear regression. *J Classif*. 1988;5(2):249-282.

Goodman LA. Exploratory latent structure analysis using both
identifiable and unidentifiable models. *Biometrika*.
1974;61(2):215-231.

Grün B, Leisch F. FlexMix version 2: finite mixtures with concomitant
variables and varying and constant parameters. *J Stat Softw*.
2008;28(4):1-35.

Jacobs RA, Jordan MI, Nowlan SJ, Hinton GE. Adaptive mixtures of
local experts. *Neural Comput*. 1991;3(1):79-87.

Jedidi K, Jagpal HS, DeSarbo WS. Finite-mixture structural equation
models for response-based segmentation and unobserved heterogeneity.
*Mark Sci*. 1997;16(1):39-59.

Klein N, Kneib T, Klasen S, Lang S. Bayesian structured additive
distributional regression for multivariate responses. *J R Stat Soc
Ser C*. 2015;64(4):569-591.

Lawley DN, Maxwell AE. *Factor Analysis as a Statistical Method*.
2nd ed. Butterworth; 1971.

Lazarsfeld PF, Henry NW. *Latent Structure Analysis*. Houghton
Mifflin; 1968.

Leisch F. FlexMix: a general framework for finite mixture models and
latent class regression in R. *J Stat Softw*. 2004;11(8):1-18.

Linzer DA, Lewis JB. poLCA: an R package for polytomous variable
latent class analysis. *J Stat Softw*. 2011;42(10):1-29.

Lubke GH, Muthén B. Investigating population heterogeneity with
factor mixture models. *Psychol Methods*. 2005;10(1):21-39.

McLachlan GJ, Peel D. *Finite Mixture Models*. Wiley; 2000.

Meredith W. Measurement invariance, factor analysis and factorial
invariance. *Psychometrika*. 1993;58(4):525-543.

Muthén B. Latent variable mixture modeling. In: Marcoulides GA,
Schumacker RE, eds. *New Developments and Techniques in Structural
Equation Modeling*. Erlbaum; 2001:1-33.

Muthén B. Latent variable analysis: growth mixture modeling and
related techniques for longitudinal data. In: Kaplan D, ed.
*Handbook of Quantitative Methodology for the Social Sciences*.
Sage; 2004:345-368.

Muthén B, Shedden K. Finite mixture modeling with mixture outcomes
using the EM algorithm. *Biometrics*. 1999;55(2):463-469.

Muthén LK, Muthén BO. *Mplus User's Guide*. 8th ed. Muthén &
Muthén; 1998-2017.

Nagin DS. Analyzing developmental trajectories: a semiparametric,
group-based approach. *Psychol Methods*. 1999;4(2):139-157.

Nagin DS. *Group-Based Modeling of Development*. Harvard University
Press; 2005.

Proust-Lima C, Philipps V, Liquet B. Estimation of extended mixed
models using latent classes and latent processes: the R package
lcmm. *J Stat Softw*. 2017;78(2):1-56.

Rigby RA, Stasinopoulos DM. Generalized additive models for
location, scale and shape. *J R Stat Soc Ser C*. 2005;54(3):507-554.

Scrucca L, Fop M, Murphy TB, Raftery AE. mclust 5: clustering,
classification and density estimation using Gaussian finite mixture
models. *R J*. 2016;8(1):289-317.

Spearman C. "General intelligence", objectively determined and
measured. *Am J Psychol*. 1904;15(2):201-292.

Verbeke G, Lesaffre E. A linear mixed-effects model with
heterogeneity in the random-effects population. *J Am Stat Assoc*.
1996;91(433):217-221.

Vermunt JK, Magidson J. Latent class cluster analysis. In: Hagenaars
JA, McCutcheon AL, eds. *Applied Latent Class Analysis*. Cambridge
University Press; 2002:89-106.

Wedel M, DeSarbo WS. A mixture likelihood approach for generalized
linear models. *J Classif*. 1995;12(1):21-55.

Widaman KF, Reise SP. Exploring the measurement invariance of
psychological instruments: applications in the substance use domain.
In: Bryant KJ, Windle M, West SG, eds. *The Science of Prevention*.
APA; 1997:281-324.

Yung YF. Finite mixtures in confirmatory factor-analysis models.
*Psychometrika*. 1997;62(3):297-330.

Zhang D, Davidian M. Linear mixed models with flexible distributions
of random effects for longitudinal data. *Biometrics*.
2001;57(3):795-802.

---

*Rendered on 2026-04-15 at 12:35 PDT.*
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/02-revision-response-latent-class-expanded.md*
