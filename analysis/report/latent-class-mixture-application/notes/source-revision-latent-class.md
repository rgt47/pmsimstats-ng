---
title: "Expanded response: latent-class faithfulness versus MVN approximation"
author: "pmsimstats team"
date: "2026-04-22"
---

# Expanded §4.2 response: latent-class faithfulness versus MVN approximation

This note expands the Pass 2 edit applied to §4.2 of
`02-dgp-mean-moderation-vs-mvn.tex`, which addressed the reviewer's
question in the email summary: if Architecture B is really an
'imperfect responder indicator', why not model that explicitly at the
participant level rather than via covariance structure?

The short paragraph added in Pass 2 acknowledged the tension. This
expanded treatment develops ten threads, the first five concerning
statistical substance and the last five concerning the psychometric
and econometric literatures from which the relevant machinery is
drawn. The exposition is deliberately paced for readers whose
primary training is in biostatistics and who may have had limited
contact with latent-variable methods; where psychometric terminology
is unavoidable, it is glossed with a biostatistical analogue on first
use.

## 0. Preamble for biostatistician readers

The reviewer's question lives at the interface between two
traditions. One is the regression-and-mixed-effects tradition
familiar from clinical trial biostatistics, in which participant
heterogeneity is represented by random effects and individual
response to treatment is estimated from repeated measurements within
person. The other is the latent-variable tradition developed within
psychometrics and econometrics, in which heterogeneity that is
structural (i.e., believed to reflect distinct underlying
subpopulations or unobserved traits) is represented by explicit
latent classes or latent factors that are estimated jointly with the
outcome model. The two traditions study overlapping problems but
adopt different conceptual defaults and different notation.

For a biostatistician reading this note, the key translations are
the following. A 'latent class' is an unobserved categorical
subgroup membership: the closest familiar analogue is the
responder/non-responder dichotomy that in ordinary trial practice is
only available post hoc through responder analysis with an arbitrary
cutoff. A 'latent factor' is an unobserved continuous trait that
governs covariation among observed variables: the closest familiar
analogue is a frailty in a survival model, or a random intercept
that explains correlation across repeated measurements on the same
participant. A 'mixture model' is a density that is a weighted sum
of simpler densities, where the weights are the unknown
class-membership probabilities: the closest familiar analogue is a
two-component gaussian where one component represents
treatment-responsive participants and the other represents
non-responsive participants, both sharing the observed outcome
scale. The EM algorithm, which is the standard estimator for these
models, is the same algorithm that biostatisticians encounter in
multiple imputation and in missing-data maximum likelihood for mixed
models; the only conceptual difference is that the 'missing data'
being imputed in a mixture fit is the unobserved class label
$Z_i$ for each participant rather than a missing covariate or
outcome value.

Against this background, the question Architecture B answers is
narrower than it may appear. Architecture B does not ask whether
discrete responder classes exist; it asks whether the
second-moment consequences of unobserved responder classes, in
particular the biomarker-by-response covariance that a mixture model
would induce, can be captured by a single multivariate normal with
treatment-state-dependent correlation. Sections 1-5 below argue
that the answer is generally yes for power calculation purposes,
with specific caveats. Sections 6-10 place this answer within the
broader psychometric spectrum and clarify what a more faithful
generative model would look like, drawing on latent-class, factor
mixture, and heterogeneous-random-effects literatures that have
already worked through the identifiability and estimation problems
in detail.

## 0.1 Notation and term crosswalk

The following shorthand is used throughout. $B_i$ denotes the
biomarker value for participant $i$; in the prazosin-PTSD
application, this is a baseline blood-pressure summary. $BR_{it}$
denotes the biological-response component of outcome for participant
$i$ at time $t$; this is one of three latent response factors in the
Hendrickson decomposition alongside a time-variant component
$TV_{it}$ and a placebo-belief component $PB_{it}$. $D_{it}$ denotes
drug state at time $t$, with the continuous-decay analysis form
$D_{it} = \exp(-\lambda t_{sd,it})$ during off-drug periods and
$D_{it} = 1$ during on-drug periods. $Z_i$ denotes an unobserved
class label.

The principal psychometric-to-biostatistical term correspondences,
given as a running glossary rather than a table to accommodate the
longer biostatistical paraphrases:

- **Latent class.** An unobserved categorical subgroup
  (biostatistical analogue: a true responder versus non-responder
  stratum that is not directly observed).
- **Latent factor.** An unobserved continuous trait that governs
  covariation among observed variables (analogue: a frailty in a
  survival model, or a random effect in a linear mixed model).
- **Class-conditional distribution $f_z$.** The outcome
  distribution within a single subgroup (analogue: a
  subgroup-specific sampling distribution).
- **Mixing proportion $\pi$.** The prevalence of a given
  subgroup in the population (analogue: a stratum proportion).
- **Class indicator function.** Membership in an unobserved
  stratum, treated as a random variable.
- **Measurement invariance.** Whether a model applies identically
  across subpopulations (analogue: the assumption underlying
  pooled versus stratum-specific analyses).
- **Growth mixture model (GMM).** A longitudinal mixed model in
  which the population of random effects is itself a finite mixture
  (analogue: a longitudinal LME with subgroup-specific trajectory
  classes).
- **Factor mixture model (FMM).** A hybrid model with discrete
  subgroups at the upper level and within-subgroup continuous
  heterogeneity at the lower level (analogue: a two-level
  hierarchical model where the upper level is categorical and the
  lower level is a standard LME).
- **Mixture of experts, regression mixture.** Covariate-dependent
  class probabilities combined with class-specific regression
  slopes (analogue: a propensity-score-like gating model paired
  with within-stratum outcome regressions).
- **EM algorithm.** The same expectation-maximisation familiar
  from multiple imputation, applied here to the missing class
  labels $Z_i$ rather than to missing covariates or outcomes.

With this vocabulary in hand, the technical development follows.

## 1. What the faithful generative model actually looks like

If blood pressure (or any candidate predictive biomarker) genuinely
indexes an unobserved responder/non-responder dichotomy, the
mechanistically correct data-generating process (DGP) is a finite
mixture:

$$
Z_i \mid B_i \sim \mathrm{Bernoulli}\bigl(\pi(B_i)\bigr), \qquad
BR_{it} \mid Z_i = z \sim f_z(t, D_{it}),
$$

where $Z_i$ is the latent class indicator (responder vs.
non-responder), $\pi(B_i) = \Pr(Z_i = 1 \mid B_i)$ is typically
modelled as a logistic function of the biomarker, and the
class-specific response distributions $f_0$ and $f_1$ differ in
their drug-response parameters (peak response, onset rate, or both).

Unpacking this for readers unfamiliar with mixture notation, the
generative recipe is the following. For each participant $i$, first
draw a latent label $Z_i$ that equals 1 with probability
$\pi(B_i)$ and 0 otherwise. Participants with higher biomarker
values have a higher chance of being drawn into the responder class
but not a deterministic assignment, which is why the biomarker is
described in the main text as an 'imperfect indicator' of the
underlying responder status. Once $Z_i$ is drawn, the within-person
trajectory $BR_{it}$ is generated from the response curve
appropriate to that class: responders receive a curve with a
substantial drug effect; non-responders receive a curve that is
essentially flat in $D_{it}$. Carryover, modelled in this
framework through the continuous drug indicator $D_{it}$, operates
only on the class-specific mean trajectories. It does not erode the
biomarker-outcome association at the level of the marginal
distribution, because the biomarker-mediated mechanism is class
membership $Z_i$, which is fixed over time and drug-state-invariant.

This feature (class membership being carryover-invariant) is the
intuitive kernel of the reviewer's concern. Under a faithful
mixture DGP, carryover washes out the within-class drug effect but
leaves the between-class mean separation intact. Any analysis model
that can see class membership, even imperfectly, retains a signal
proportional to the between-class difference. By contrast,
Architecture B encodes the biomarker-response association purely as
a second-moment object (a covariance), and carryover erodes that
covariance directly. The two DGPs therefore predict different
patterns of power loss as a function of carryover intensity, and
this is what Section 3 of the main document documents.

## 2. Why MVN differential correlation is nonetheless a defensible second-moment summary

A two-component mixture with class-dependent means and a common
within-class covariance produces a marginal joint distribution of
$(B_i, BR_{it})$ whose first two moments match, to leading order,
those of a single multivariate normal with
$\mathrm{Cor}(B, BR) = c_{bm}$ chosen so that $c_{bm}^2$ equals the
between-class variance fraction (the proportion of total variance in
$BR$ that is attributable to class separation rather than to
within-class residual variation).

The biostatistical intuition is the following. If responders and
non-responders differ by some amount $\Delta$ in mean drug response
and the biomarker has some accuracy $\pi(B_i)$ in predicting class
membership, then marginally (i.e., averaging over the unobserved
$Z_i$) there is a linear association between $B_i$ and $BR_{it}$
whose strength is determined by two factors: how far apart the two
classes are in drug response ($\Delta$) and how well the biomarker
discriminates them ($\pi$). A single multivariate normal with an
appropriately chosen correlation coefficient $c_{bm}$ can reproduce
this marginal linear association even though no such gaussian
actually generated the data. This is the sense in which Architecture
B 'approximates' the mixture DGP: it reproduces the marginal
biomarker-response association that a mixture would produce, without
committing to the discrete generative mechanism.

Under mild regularity conditions (adequate class overlap and
moderate class-prevalence probabilities $\pi$ away from 0 and 1),
the mixture's second-moment structure is well-approximated by an MVN
with differential correlation: high correlation in treatment states
where class-specific means separate (i.e., on-drug) and shrinkage
towards zero in treatment states where class-specific means
coincide (i.e., off-drug after full carryover washout). Architecture
B as implemented therefore captures the covariance implications of
the latent-class mechanism without committing to the discrete
generative form. For power calculations driven by linear
association statistics (t-tests for $\beta_{bm:D}$ in a linear mixed
model), this approximation is often adequate.

## 3. Where the approximation breaks down

The MVN approximation will be misleading in three identifiable
regimes, each of which has a concrete biostatistical signature.

- **Bimodal responses.** When $f_0$ and $f_1$ are well-separated
  (non-responders show near-zero drug effect and responders show a
  large effect), the marginal distribution of $BR$ is bimodal: a
  histogram of on-drug response would show two peaks rather than a
  single peak with gaussian scatter around it. No MVN matches the
  tails of a bimodal distribution, and power under the true mixture
  DGP will be higher than Architecture B predicts because a
  well-specified mixture analysis can exploit the bimodality.
  Biostatisticians working with clinical response data will
  recognise this pattern from settings in which a subpopulation
  genuinely does not respond to a drug (e.g., certain targeted
  oncology treatments), where response histograms are
  characteristically non-gaussian.
- **Strong class-membership gradient.** When $\pi(B_i)$ is steep
  (approaching a step function at some biomarker threshold), the
  biomarker behaves almost deterministically as a class label.
  Architecture A (mean moderation) with a thresholded or
  dose-response moderation then becomes the better approximation,
  because the biomarker-response relationship is dominated by its
  mean-structure component. The familiar biostatistical scenario
  here is a diagnostic biomarker with near-perfect sensitivity and
  specificity for the responsive phenotype, where dichotomising the
  biomarker at its clinical threshold would lose little information.
- **Class-varying covariance.** If the autocorrelation, residual
  variance, or placebo-response parameters themselves differ between
  classes (i.e., responders and non-responders differ not only in
  mean drug effect but in how variable their day-to-day symptom
  scores are, or in how autocorrelated those scores are), no
  single-component MVN can represent the generative structure.
  Analysis models assuming constant covariance will be misspecified
  regardless of whether they model the biomarker interaction in the
  mean or in the covariance. This is the latent-class analogue of
  the heteroscedastic residual problem that biostatisticians address
  with `varIdent` or `weights` arguments in `nlme::lme`, except that
  the stratifying variable is itself unobserved.

## 4. Implications for the power estimates in Section 3

The substantial Architecture B power loss under carryover (40-60%
across the designs reported in Section 3 of the main document)
reflects erosion of a second-moment signal. A mixture-model analysis
fitted to mixture-DGP data would, in principle, recover additional
power by estimating class membership directly, because the class
label is carryover-invariant. The gap between Architecture B's
carryover-sensitive power estimate and the mixture-faithful power
estimate is an upper bound on what a mixture analysis could recover
at the same carryover level; whether this gap is large in the
prazosin-PTSD parameter regime is an empirical question not yet
answered by the existing pmsimstats infrastructure.

A focused simulation contrasting three cells would characterise the
loss attributable to second-moment approximation versus the loss
attributable to analysis-model misspecification:

1. Mixture DGP analysed with the current MVN linear mixed model
   (i.e., the analysis model makes no attempt to recover class
   structure).
2. Mixture DGP analysed with a matched mixture analysis model (e.g.,
   a `lcmm`-style joint latent class mixed model).
3. MVN DGP (Architecture B as currently implemented) analysed with
   both the current MVN linear mixed model and the `lcmm`-style
   mixture model, to verify that the mixture model does not gain
   power on data for which no latent class structure exists.

The difference between cells 1 and 2 quantifies the power recovery
attributable to using a faithful analysis model; the difference
between cells 1 and 3 (with MVN analysis) quantifies the
approximation error in Architecture B itself. This factorial design
maps cleanly onto the existing pmsimstats simulation harness and
would be a natural extension of the present work.

## 5. Practical recommendation and scope of the open question

For the prazosin-PTSD application, mixture modelling is attractive
in principle but carries identifiability costs that are acute in
N-of-1 settings with modest participant counts. The EM algorithm
for finite mixtures of mixed-effects models requires either
informative priors on class-membership probabilities or strong
separation between class-specific response curves, neither of which
is assured in a trial designed to detect an interaction whose
existence is itself uncertain. The standard diagnostic for mixture
identifiability, the entropy of posterior class-membership
probabilities, typically requires sample sizes in the hundreds
before it stabilises at values that support confident class
extraction (Bauer and Curran 2003, Bauer 2007). Aggregated N-of-1
trials of the size contemplated in the prazosin-PTSD application
(typically $N \in [30, 100]$) are an order of magnitude below this
threshold.

The open question, therefore, is not whether mixture modelling is
more faithful (it is), but whether the fidelity gain is recoverable
at realistic sample sizes, or whether analysis-strategy mitigations
(Section 6 of the main document: restricted analysis, weighting,
within-subject contrasts) deliver a better power-per-complexity
trade. These analysis-strategy mitigations are closer to standard
biostatistical practice and have better-understood inferential
properties at the relevant sample sizes. Resolving the
mixture-versus-mitigation question empirically is out of scope for
the present white paper but is a natural next step in the
pmsimstats programme, and the factorial simulation outlined in
Section 4 above would provide the evidence base required.

## 6. Connection to the psychometric latent variable literature

The reviewer's question, translated into psychometric vocabulary,
asks whether biomarker-moderated treatment response is better
modelled as a continuous-trait phenomenon (classical factor-analytic
or latent regression tradition) or as a discrete-class phenomenon
(latent class or mixture tradition). The psychometric literature
has worked this distinction over for more than fifty years,
originally in the context of individual differences in education,
personality, and psychopathology, and the taxonomy it developed maps
directly onto the Architecture A versus Architecture B distinction
drawn in this paper.

The core generative question is whether unobserved heterogeneity in
outcome distributions is best represented by (i) a single population
with continuous individual-differences variation along one or more
latent dimensions (factor-analytic models; Spearman 1904;
Lawley and Maxwell 1971), (ii) a small number of qualitatively
distinct subpopulations with within-class homogeneity (latent class
models; Lazarsfeld and Henry 1968; Goodman 1974), or (iii) a hybrid
in which classes differ along continuous dimensions within each
class (factor mixture models; Lubke and Muthén 2005;
Clark et al. 2013). For biostatistician readers, the three
alternatives correspond respectively to: a standard linear mixed
model with a continuous random effect (option i), a stratified
analysis with unobserved strata (option ii), and a stratified
analysis with within-stratum random effects (option iii). The same
trichotomy appears throughout the clinical trial simulation
literature, usually implicitly: Architecture A encodes the
continuous individual-differences view (treatment response is a
deterministic function of a continuous biomarker), while
Architecture B as implemented in Hendrickson et al. (2020) is most
naturally interpreted as a second-moment approximation to the
latent-class view.

Recognising this connection has practical consequences for the
pmsimstats programme. The decision between Architecture A and
Architecture B is not merely a parameterisation convenience; it
reflects a substantive biological hypothesis about whether
biomarker-outcome coupling arises from graded continuous variation
(Architecture A, consistent with a mechanism in which every
participant has some degree of drug response, scaled by biomarker)
or from discrete subpopulation structure (Architecture B in its
latent-class interpretation, consistent with a mechanism in which
some participants respond substantially and others not at all,
with the biomarker predicting which are which). The two hypotheses
yield comparable marginal biomarker-outcome associations at a
single timepoint but diverge sharply under carryover, which is the
core finding of Section 3.

## 7. Finite mixture taxonomy relevant to biomarker moderation

Four families of finite mixture model from the psychometric and
econometric literatures bear directly on the biomarker-treatment
interaction problem. The exposition below introduces each family
with its biostatistical analogue, its psychometric origin, and its
specific relevance to the pmsimstats framework.

### 7.1 Latent class and latent profile analysis

Latent class analysis (LCA) and the closely related latent profile
analysis (LPA) are mixture models for cross-sectional data in which
each participant belongs to one of a small number of unobserved
classes, and the observed variables have class-specific marginal
distributions. LCA is the version with categorical indicators (e.g.,
symptom-present / symptom-absent items in a psychiatric assessment),
and LPA is the version with continuous indicators (e.g., biomarker
measurements or symptom severity scores).

The biostatistical reader should picture this as follows. Suppose
one wished to identify a responder subgroup using baseline
characteristics alone, without access to any treatment-response
information. If the responder subgroup has a characteristic
baseline biomarker profile that differs from the non-responder
profile, a mixture model fitted to the baseline biomarker set can
estimate (a) the proportion of participants in each class, (b) the
class-specific mean biomarker profiles, and (c) each participant's
posterior probability of belonging to each class. These posterior
probabilities function as soft subgroup-membership predictors that
can subsequently be used as covariates in outcome analysis. This is
the two-stage analytic strategy standard in psychometric
applications. In the biomarker-trial context, an LPA-style model
would posit two or more responder classes with class-specific
marginal distributions of on-drug and off-drug response, and
estimate class-membership probabilities from the observed response
trajectories rather than from baseline covariates alone.

Classical references include Lazarsfeld and Henry (1968), Goodman
(1974), Clogg (1995), and McLachlan and Peel (2000); the Vermunt
and Magidson Latent GOLD software and the `poLCA` R package
implement the standard EM estimation routines. For a biostatistician
coming to this literature, McLachlan and Peel (2000) is the most
accessible entry point, as it is structured as a statistical
monograph rather than a substantive-psychometric one.

### 7.2 Growth mixture models and group-based trajectory models

Growth mixture models (GMM; Muthén and Shedden 1999; Muthén 2001,
2004) extend LPA to longitudinal data by specifying class-specific
growth (or, in the trial setting, drug-response) curves with
within-class random-effects variation. In biostatistical vocabulary,
a GMM is a longitudinal linear mixed model in which the population
distribution of random effects is itself a mixture: rather than
assuming a single multivariate normal for random intercepts and
slopes, GMM assumes a mixture of multivariate normals, each
corresponding to a distinct trajectory class. Group-based trajectory
models (Nagin 1999, 2005) are the degenerate case in which
within-class random-effect variances are set to zero, so that all
participants within a class share an identical trajectory up to
observation-level residual noise.

GMM is the closest psychometric analogue of the mixture DGP
discussed in §4.2: each participant belongs to a latent class with
its own drug-response trajectory, and the class label is
drug-state-invariant. For the pmsimstats biomarker moderation
setting, a GMM that uses the biomarker as a class-membership
predictor is essentially the faithful generative model described in
Section 1 of this note. Bauer and Curran (2003) and Bauer (2007)
document the identifiability and overextraction hazards that
accompany GMM estimation at realistic sample sizes, showing in
particular that GMM can extract apparent class structure from data
generated by single-class non-gaussian processes. These cautions
apply with full force to N-of-1 trials, where the number of
participants is typically well below the sample sizes at which GMM
class-count recovery has been validated.

### 7.3 Factor mixture models

Factor mixture models (FMM; Yung 1997; Arminger, Stein and
Wittenberg 1999; Lubke and Muthén 2005) combine continuous latent
factors with discrete latent classes, yielding a generative
structure in which within-class variation is governed by a factor
model and between-class variation is governed by class-specific
means and (optionally) class-specific factor loadings or residual
covariances. The FMM framework subsumes both the pure continuous
(single-class factor model) and the pure discrete (LCA with no
within-class structure) cases as limits.

For biostatisticians, the practical intuition is that FMM is a
two-level hierarchy: at the upper level, participants are assigned
to one of a small number of classes (as in LCA); at the lower
level, each class has its own within-class random-effects structure
(as in a standard linear mixed model). This corresponds to the
biological hypothesis that there are genuinely distinct responder
subtypes, and within each subtype there is additional continuous
heterogeneity driven by covariates or unmeasured individual
differences. For the biomarker-moderation problem, FMM provides a
principled way to represent a biomarker as both a noisy class
indicator (via class-membership probability) and a continuous
modulator of within-class response (via factor loadings on the
biological-response factor). This hybrid representation is a better
fit to the biology that motivates the prazosin-PTSD application
than either pure LCA or pure factor analysis, but estimation
requirements are correspondingly higher.

### 7.4 Regression mixtures and mixtures of experts

Regression mixtures (DeSarbo and Cron 1988; Wedel and DeSarbo 1995;
Jedidi, Jagpal and DeSarbo 1997) and the closely related mixture of
experts architecture (Jacobs, Jordan, Nowlan and Hinton 1991) allow
class-specific regression coefficients with class-membership
probabilities that depend on covariates (the 'gating' function in
mixture-of-experts terminology). For a biostatistician, the clearest
way to see this family is through its relationship to standard
interaction modelling. An ordinary linear model with a
treatment-by-biomarker interaction posits a single regression
surface in which the slope of response on treatment varies smoothly
with the biomarker. A regression mixture instead posits two or
more distinct regression surfaces (one per class), with a separate
covariate-dependent model governing which surface applies to which
participant. The two-stage structure (classification by the gating
model, then prediction by the within-class regression) is directly
analogous to the two-stage structure of causal inference methods
that use propensity scores to stratify before estimating
within-stratum treatment effects.

In the trial-simulation idiom, regression mixtures correspond to a
DGP in which the biomarker governs both the probability of being a
responder and the magnitude of the within-responder drug effect,
which is precisely the generative form the reviewer's question
invites. `flexmix` in R (Leisch 2004; Grün and Leisch 2008)
implements this family with user-definable component models,
including mixed-effects components suitable for repeated-measures
data.

## 8. The Architecture B spectrum: covariance, mean, and combined biomarker moderation

Architecture B as implemented in Hendrickson et al. (2020) and
pmsimstats specifies biomarker-dependent covariance structure with
population means held constant across biomarker values. This is one
point on a broader spectrum of biomarker-moderated DGPs recognised
in the psychometric literature; the spectrum is usefully organised
by which moments of the joint distribution depend on the biomarker.
For biostatistician readers, the moment-by-moment organisation is
the natural way in: any multivariate normal distribution is
characterised by its mean vector and covariance matrix, and
biomarker moderation can in principle enter either object or both.

### 8.1 Covariance-only moderation (Architecture B as implemented)

The Hendrickson-Schork specification holds $E[BR_{it}]$ fixed across
biomarker strata (i.e., high- and low-biomarker participants have
the same expected drug response at any given time) and encodes the
interaction entirely in $\mathrm{Cov}(B_i, BR_{it} \mid D_{it})$
(i.e., high- and low-biomarker participants differ in the
covariance between biomarker and drug response, specifically in how
strongly the biomarker co-varies with response during on-drug
periods versus off-drug periods).

This is an unusual restriction within the psychometric mixture
tradition. Most standard FMM and GMM parameterisations permit
class-specific means at a minimum, and the
means-constant-covariance-varies case is typically encountered only
as a constrained submodel used for testing invariance hypotheses
(Meredith 1993; Widaman and Reise 1997). The biostatistical analogue
is a model in which two subgroups have identical expected outcomes
but differ in residual variance or in within-subject
autocorrelation: such models are occasionally fitted in sensitivity
analyses but rarely assumed as the primary DGP. Architecture B's
appeal in the Hendrickson context is its tractability for power
simulation under multivariate normality (the entire DGP can be
reduced to drawing from a single MVN with a cleverly constructed
covariance). Its cost, as the reviewer observes, is that it does
not correspond cleanly to any natural generative mechanism.

### 8.2 Combined mean and covariance moderation

Arminger, Stein and Wittenberg (1999), in the canonical Psychometrika
treatment of mixtures with covariate-dependent structure, develop
finite mixtures of conditional mean- and covariance-structure models
in which both the expected-response vector and the residual
covariance matrix may differ between classes and may depend on
observed covariates within class. This is the general psychometric
form of Architecture B and reduces to Architecture A in the limit
of a single class with covariate-dependent mean only, and to the
Hendrickson-implemented submodel in the limit of class-invariant
means with covariate-dependent covariance only.

The biostatistical interpretation of the combined form is that
responders and non-responders differ both in their typical drug
response (mean moderation, Architecture A) and in how noisily that
response is expressed within class (covariance moderation,
Architecture B). Both channels carry information about class
membership: the mean-moderation channel contributes via the
conditional expectation $E[BR \mid B]$ and the covariance-moderation
channel contributes via the conditional variance
$\mathrm{Var}[BR \mid B]$. Under carryover, the mean-moderation
channel is comparatively robust because class means need not
converge off-drug; the covariance-moderation channel is eroded
because off-drug covariance is by construction attenuated. A
combined-specification simulation would therefore be expected to
show power loss intermediate between the Architecture A and
Architecture B extremes, and this intermediate case is arguably the
most biologically plausible regime for the prazosin-PTSD
application. Verbeke and Lesaffre (1996) provide the analogous
treatment for linear mixed-effects models with heterogeneity in the
random-effects population, which is the more direct generalisation
of the repeated-measures structure used in pmsimstats and is the
entry point most likely to feel familiar to biostatistician readers.

### 8.3 Location-scale moderation (GAMLSS and related)

Generalised additive models for location, scale, and shape (GAMLSS;
Rigby and Stasinopoulos 2005) and the related distributional
regression tradition (Klein et al. 2015) allow the biomarker to
enter multiple parameters of the outcome distribution simultaneously
through separate link functions. In plain terms, a standard
regression models the conditional mean of an outcome as a function
of covariates; a GAMLSS-style distributional regression models the
conditional mean, the conditional variance, and if appropriate the
conditional skewness and kurtosis, each as a separate function of
covariates. For biostatisticians, the familiar entry point is the
heteroscedastic linear regression in which the residual variance is
allowed to depend on a covariate; GAMLSS generalises this to
arbitrary distributional shape.

In the biomarker-trial context, a GAMLSS-based DGP would allow
high-biomarker participants to have both a shifted mean drug
response (a location effect, which is Architecture A) and altered
response variability (a scale effect, which is analogous to
Architecture B's covariance moderation but operating on marginal
rather than joint moments), without requiring discrete class
structure. The GAMLSS framework is thus the continuous-heterogeneity
analogue of FMM and occupies an intermediate position between pure
Architecture A (mean only) and pure Architecture B (covariance
only). This intermediate position is attractive for sensitivity
analysis because it does not require committing to a specific
number of latent classes.

### 8.4 Heterogeneous random-effects models

Verbeke and Lesaffre (1996), Zhang and Davidian (2001), and
Proust-Lima et al. (2017) develop linear mixed-effects models in
which the random-effects distribution is itself a finite mixture,
allowing class-specific random-intercept and random-slope
distributions. The biostatistical setting this extends is the
standard linear mixed model with random intercept and random
treatment slope, familiar from any longitudinal trial analysis
using `nlme::lme` or `lme4::lmer`. In the standard model, the
random treatment slope is drawn from a single normal distribution
with mean equal to the average treatment effect and variance
capturing individual heterogeneity. In the heterogeneous
random-effects extension, that single distribution is replaced by
a mixture: a fraction of participants are drawn from a 'responder'
distribution with a large mean treatment slope and another fraction
from a 'non-responder' distribution with a near-zero mean treatment
slope. The biomarker predicts which mixture component a participant
comes from.

Applied to aggregated N-of-1 trial data, this framework represents
responder heterogeneity directly at the participant-specific random
treatment effect, which is arguably the most clinically transparent
representation available: responders and non-responders have
distinct random-slope distributions on the treatment indicator, and
the biomarker predicts class membership rather than entering the
analysis-model mean structure at the fixed-effect level. The `lcmm`
R package (Proust-Lima et al. 2017) implements this family and is
the closest available tooling to the pmsimstats use case: it
accepts standard longitudinal mixed-model syntax, fits a mixture on
the random-effects distribution, and returns posterior class
probabilities for each participant. For the biostatistician
planning an extension of the current pmsimstats simulation, `lcmm`
is the recommended entry point.

## 9. Implications for the pmsimstats framework

Three implications follow from locating Architecture B within the
broader psychometric spectrum; each is stated first in summary form
and then unpacked in biostatistical terms.

**First**, the covariance-only parameterisation adopted in
Hendrickson et al. (2020) is the most restrictive operational case
of a much richer family. The substantial power loss under carryover
documented in Section 3 is a specific property of this restrictive
case and may not generalise to the mean-plus-covariance or
heterogeneous-random-effects variants that are standard in
psychometrics. Unpacked, this means that a sensitivity analysis
based solely on Architecture B risks overstating the
carryover-induced power loss that a biomarker-guided aggregated
N-of-1 trial should expect, because the DGP on which Architecture B
is built puts all of the biomarker-outcome signal into a channel
(the covariance) that carryover erodes. If, in reality, some of
the biomarker-outcome signal flows through the mean structure (as
in Architecture A or in the combined specification of §8.2), that
component of the signal survives carryover, and overall power loss
is smaller than Architecture B predicts. A complete sensitivity
analysis for the prazosin-PTSD application should accordingly
report power under at least three points on the spectrum:
Architecture A (mean-only), Architecture B as implemented
(covariance-only), and a combined mean-plus-covariance specification
of the Arminger et al. (1999) form. The existing pmsimstats
infrastructure already supports the first two; adding the third
requires a modest extension to `generateData` that permits
simultaneous mean scaling and covariance moderation.

**Second**, the identifiability concerns documented in the
psychometric GMM and FMM literatures (Bauer and Curran 2003; Bauer
2007; Lubke and Muthén 2005) bound what can be asked of
mixture-based analysis models at N-of-1 trial sample sizes. The
published FMM applications typically involve several hundred to
several thousand participants; biomarker-trial applications with
$N \in [30, 100]$ are well outside this range. For a
biostatistician, this means that analysis plans built around
fitting an `lcmm` or `flexmix` mixture model directly to trial data
are unlikely to yield reliable class-membership inference and
should not be relied upon as the primary analysis. Analysis-strategy
mitigation (Section 6 of the main document) is therefore more
likely than mixture-model analysis to deliver recoverable power
gains at realistic trial sizes. Mixture modelling remains useful as
a sensitivity check (to see whether class-structured residuals are
consistent with the pre-specified analysis) and as a DGP for
simulation (to verify robustness of the primary analysis to
class-structured departures from gaussianity).

**Third**, if the pmsimstats programme pursues a richer DGP in
subsequent work, the natural target is the heterogeneous-random-
slopes form (§8.4) rather than a full FMM. The random-slopes
parameterisation preserves the analysis model's compatibility with
existing N-of-1 estimation machinery (`nlme::lme` with
participant-specific random treatment effects) while admitting the
biologically faithful latent-class generative structure at the DGP
level. For a biostatistician, this is the most conservative path
forward: the simulation DGP is extended in a familiar direction (a
mixture over random-slope distributions, which is a natural
relaxation of the homogeneous-population assumption built into
standard mixed-effects models) while the analysis model remains the
pre-specified `nlme::lme` fit. This preserves the audit chain of
the current framework and permits direct comparison with the
covariance-only baseline under carryover.

## 10. Software and implementation references

R packages relevant to estimating the models surveyed above, in
approximate order of fit to the pmsimstats use case:

- `lcmm` (Proust-Lima et al. 2017). Joint latent class mixed models
  for longitudinal and survival data; the closest available tooling
  for heterogeneous random-effects specification (§8.4). For the
  biostatistician coming to this literature, `lcmm`'s syntax is
  essentially `nlme::lme` with a latent-class extension, making it
  the lowest-barrier entry point.
- `flexmix` (Leisch 2004; Grün and Leisch 2008). General regression
  mixture framework with user-definable component models; suitable
  for regression-mixture and mixture-of-experts specifications
  (§7.4). `flexmix` is more general than `lcmm` but requires more
  configuration.
- `mclust` (Scrucca et al. 2016). Gaussian finite mixture models
  with a range of class-specific covariance parameterisations;
  supports the Arminger et al. (1999) mean-plus-covariance
  specifications non-longitudinally. `mclust` is useful for
  cross-sectional LPA-style analyses of baseline biomarker panels
  but does not natively handle repeated-measures structure.
- `gamlss` (Rigby and Stasinopoulos 2005). Distributional regression
  with covariate-dependent location, scale, and shape parameters
  (§8.3). For biostatisticians familiar with `mgcv`, the syntax and
  model-fitting workflow of `gamlss` will be largely familiar.
- `poLCA` (Linzer and Lewis 2011). Polytomous latent class analysis
  for the categorical-indicator case. Less directly relevant to the
  pmsimstats continuous-outcome setting, but useful if baseline
  class indicators (e.g., symptom-checklist items) are incorporated
  into class assignment.
- `OpenMx` and `lavaan` (with mixture extensions). Structural
  equation mixture modelling for FMM-style specifications; more
  general but with steeper configuration cost than `lcmm` or
  `flexmix`.

The psychometric reference implementation outside R is Mplus
(Muthén and Muthén 1998-2017), which remains the standard for
GMM, FMM, and related mixture-SEM specifications. Most of the
psychometric literature on mixture-model identifiability and
overextraction hazards (Bauer and Curran 2003; Lubke and Muthén
2005) uses Mplus for its primary simulations, and investigators
cross-validating pmsimstats simulation results against published
psychometric benchmarks may need occasional access to Mplus for
direct comparison.

## Relation to the existing §4.2 paragraph

The five threads in Sections 1-5 above could extend the §4.2 Pass 2
paragraph in place, or be promoted to a new §4.3
('Mixture-modelling alternative to MVN differential correlation')
so that §4.2 retains its current clinical-examples focus. Either
approach preserves the document's existing architecture-comparison
narrative while giving the reviewer's concern the technical
treatment it warrants.

The psychometric-connection material in Sections 6-10 is more
naturally placed either as a new §4.4 ('Architecture B in the
psychometric latent-variable tradition') or as an appendix. It
answers a distinct question from the §4.2 edit: whereas §4.2
concerns the biological faithfulness of Architecture B,
Sections 6-10 concern its methodological provenance and the
spectrum of related specifications already developed in the
psychometric literature. Treating the two as separate document
additions preserves the logical structure of Section 4 (biological
assumptions) while giving the methodological-provenance discussion
the space it requires.

For a biostatistician preparing the manuscript revision, the
recommended sequencing is: incorporate Sections 1-5 into the main
response to reviewer (as §4.3), defer Sections 6-10 to a
methodological appendix that can be cited but need not be read for
the primary argument to go through, and commit to the factorial
simulation described in Section 4 as a follow-on paper. This
preserves the white paper's focus on the prazosin-PTSD application
while placing the broader methodological context on record for
future reviewers and collaborators.

## References (additions beyond the main document bibliography)

Arminger G, Stein P, Wittenberg J. Mixtures of conditional mean-
and covariance-structure models. *Psychometrika*. 1999;64(4):
475-494.

Bauer DJ. Observations on the use of growth mixture models in
psychological research. *Multivariate Behav Res*. 2007;42(4):
757-786.

Bauer DJ, Curran PJ. Distributional assumptions of growth mixture
models: implications for overextraction of latent trajectory classes.
*Psychol Methods*. 2003;8(3):338-363.

Clark SL, Muthén B, Kaprio J, et al. Models and strategies for
factor mixture analysis: an example concerning the structure
underlying psychological disorders. *Struct Equ Modeling*. 2013;
20(4):681-703.

Clogg CC. Latent class models. In: Arminger G, Clogg CC, Sobel ME,
eds. *Handbook of Statistical Modeling for the Social and
Behavioral Sciences*. Plenum; 1995:311-359.

DeSarbo WS, Cron WL. A maximum likelihood methodology for
clusterwise linear regression. *J Classif*. 1988;5(2):249-282.

Goodman LA. Exploratory latent structure analysis using both
identifiable and unidentifiable models. *Biometrika*. 1974;61(2):
215-231.

Grün B, Leisch F. FlexMix version 2: finite mixtures with
concomitant variables and varying and constant parameters.
*J Stat Softw*. 2008;28(4):1-35.

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

Leisch F. FlexMix: a general framework for finite mixture models
and latent class regression in R. *J Stat Softw*. 2004;11(8):1-18.

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

Muthén LK, Muthén BO. *Mplus User's Guide*. 8th ed. Muthén and
Muthén; 1998-2017.

Nagin DS. Analyzing developmental trajectories: a semiparametric,
group-based approach. *Psychol Methods*. 1999;4(2):139-157.

Nagin DS. *Group-Based Modeling of Development*. Harvard University
Press; 2005.

Proust-Lima C, Philipps V, Liquet B. Estimation of extended mixed
models using latent classes and latent processes: the R package
lcmm. *J Stat Softw*. 2017;78(2):1-56.

Rigby RA, Stasinopoulos DM. Generalized additive models for
location, scale and shape. *J R Stat Soc Ser C*. 2005;54(3):
507-554.

Scrucca L, Fop M, Murphy TB, Raftery AE. mclust 5: clustering,
classification and density estimation using Gaussian finite mixture
models. *R J*. 2016;8(1):289-317.

Spearman C. 'General intelligence', objectively determined and
measured. *Am J Psychol*. 1904;15(2):201-292.

Verbeke G, Lesaffre E. A linear mixed-effects model with
heterogeneity in the random-effects population. *J Am Stat Assoc*.
1996;91(433):217-221.

Vermunt JK, Magidson J. Latent class cluster analysis. In:
Hagenaars JA, McCutcheon AL, eds. *Applied Latent Class Analysis*.
Cambridge University Press; 2002:89-106.

Wedel M, DeSarbo WS. A mixture likelihood approach for generalized
linear models. *J Classif*. 1995;12(1):21-55.

Widaman KF, Reise SP. Exploring the measurement invariance of
psychological instruments: applications in the substance use
domain. In: Bryant KJ, Windle M, West SG, eds. *The Science of
Prevention*. APA; 1997:281-324.

Yung YF. Finite mixtures in confirmatory factor-analysis models.
*Psychometrika*. 1997;62(3):297-330.

Zhang D, Davidian M. Linear mixed models with flexible distributions
of random effects for longitudinal data. *Biometrics*. 2001;57(3):
795-802.

---

*Rendered on 2026-04-22 at 11:20 PDT.*<br>
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/02-revision-response-latent-class-expanded.md*
