# N-of-1 Trials and Precision Medicine: A Detailed Lexicon
*2026-06-17 17:42 PDT*

A working glossary of the terminology used across the N-of-1 / single-case
trial and precision-medicine literature. Each entry pairs a technical term
with a plain-language definition pitched at an undergraduate statistics
major, and lists the alternative names the same concept carries in
different fields (biostatistics, psychometrics / single-case research, and
precision / personalized medicine).

**Provenance.** Terms were mined from the full text of N-of-1 and
precision-medicine articles in the project Zotero library
(`nof1` collection, 72 curated articles; full text obtained and
machine-read for 42). The statistical and psychometric supplement was
compiled from standard methodology references. See the closing
*Provenance and epistemic status* section for scope and caveats.

**How to read the synonym column.** Where biostatistics, psychometrics,
and precision medicine use different names for the *same* underlying
concept, all variants are listed together so the same idea can be
recognised whatever vocabulary a given paper adopts. A condensed
same-concept / different-name map appears in Section 12.

---

## 1. Trial designs and architectures

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| N-of-1 trial | A controlled experiment run inside a single patient, who is switched between treatments multiple times so the best treatment *for that one person* can be estimated. | single-patient trial; single-subject trial; single-case experimental design (SCED, psychology); individualized / idiographic clinical trial; one-person trial; n=1 RCT |
| Aggregated N-of-1 trial | A set of separate single-patient trials pooled together to estimate both a population-average effect and sharpened individual effects. | series / combined / pooled N-of-1 trials; aggregated personalized studies; multiple N-of-1 study |
| Single-case experimental design (SCED) | The behavioural-science family of designs that study one subject repeatedly across baseline and treatment phases; N-of-1 trials are the randomized-crossover member of this family. | single-subject design; single-case design; single-subject research |
| Crossover trial | A trial in which each participant receives the treatments in sequence over separate periods, serving as their own control. | cross-over / change-over trial; within-subject / within-patient design; repeated-measures comparison |
| AB/BA (2x2) crossover | The simplest crossover: half the subjects get A then B, half get B then A, in two periods. | two-period two-treatment design; classical crossover |
| Multiple-crossover design | A crossover in which the patient switches back and forth between treatments several times rather than once. | multi-crossover; repeated-period crossover; double crossover |
| Parallel-group trial | A trial in which each participant receives only one treatment and group averages are compared between people. | parallel-group RCT (PG-RCT); between-participant / between-patient design |
| Randomized controlled trial (RCT) | A study that assigns participants to treatments by chance to give an unbiased comparison. | controlled trial; randomized trial |
| AB / ABAB design | Single-case notation for alternating a baseline/treatment (A) and an alternate treatment (B); ABAB introduces, withdraws, and reintroduces treatment. | reversal design; withdrawal design; challenge-withdrawal design |
| Multiple-baseline design | A single-case design that staggers when treatment starts across people, settings, or behaviours to rule out coincidental change. | staggered-baseline design |
| Changing-criterion design | A single-case design in which the success threshold is shifted in steps and the outcome is shown to track each new target. | -- |
| Alternating-treatments design | A single-case design that rapidly switches among two or more treatments within one subject to compare them. | multi-element design |
| Open-label (OL) design | A study or phase in which both patient and staff know which (active) treatment is given. | unblinded design |
| Blinded discontinuation (BDC) | A phase in which active treatment is secretly replaced by placebo to see whether benefit fades. | open-label + blinded discontinuation (OL+BDC) |
| Micro-randomized trial (MRT) | A within-person trial (often via phone or wearable) that randomizes a brief intervention many times to build just-in-time adaptive interventions. | micro-randomised trial |
| Interrupted time-series design | Repeated measurements on one unit compared in level/slope before and after a defined intervention point. | ITS; single-case time series |
| Adaptive design | A trial allowing prospectively planned changes (e.g., dropping arms, reallocation) based on accumulating data. | response-adaptive / sequential design; group-sequential design |
| Play-the-winner / response-adaptive randomization | An adaptive scheme that shifts the randomization odds toward the treatment performing better. | response-adaptive randomization |
| Biomarker-stratified trial | A trial that splits randomization by a baseline biomarker and tests the treatment within each biomarker group and overall. | biomarker-driven / stratified design |
| Enrichment / adaptive-enrichment design | A trial that preferentially enrols (or, after an interim look, narrows to) the patients most likely to respond. | population enrichment; adaptive enrichment |
| Basket trial | One treatment tested across multiple diseases that share a common molecular target. | -- |
| Umbrella trial | Multiple treatments matched to different biomarkers within a single disease. | -- |
| Platform trial | A standing trial framework evaluating multiple treatments against a shared control, adding or dropping arms over time. | -- |
| Two-stage / multi-stage design | A trial run in phases with an interim analysis between them that can change later conduct. | two-stage adaptive design |

## 2. Periods, sequence, carryover, and washout

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Treatment period | A block of time during which a patient stays on one assigned treatment before switching. | period; treatment block; phase; intervention administration period |
| Cycle | One full pass through the treatments being compared (e.g., one A-then-B pair). | pair; treatment pair; block |
| Sequence | The full ordered list of treatment periods a participant goes through. | treatment sequence; allocation sequence |
| Washout period | A treatment-free gap inserted so leftover effects of the prior treatment disappear before the next is assessed. | wash-out; clearance / rest interval |
| Active (analytical) washout | Switching directly to the next treatment but not *measuring* until the previous drug's effect has gone, used when stopping treatment would harm the patient. | analytical washout |
| Run-in / titration phase | A preliminary phase before randomization used to stabilize dose, confirm eligibility, clear prior drugs, or assess adherence. | lead-in; stabilization phase; dose titration / escalation |
| Baseline phase | The pre-treatment stretch of repeated measurements used as the comparison reference. | pre-intervention phase; phase A |
| Carryover effect | When a treatment's effect persists into and contaminates the following period. | carry-over; residual / lingering effect |
| Simple carryover | A carryover whose size is the same regardless of which treatment follows; assumed to last only one period. | first-order carry-over |
| Complex carryover | A carryover whose size or direction depends on the specific order or pairing of treatments. | treatment-sequence-dependent carryover |
| Behavioural / psychometric carryover | A non-pharmacological carryover where earlier treatment changes habits, adherence, or how a patient remembers rating themselves, mimicking a drug carryover. | behavioural carry-over; recall carryover |
| Period effect | A systematic change in outcome tied to which time period it is, independent of treatment. | time / secular / nominal trend; cycle effect |
| Order / sequence effect | A bias arising because the order in which treatments are given affects the response. | treatment-ordering effect |
| Treatment-by-period interaction | When the treatment difference is not the same across periods, often a sign of carryover or instability. | period-by-treatment interaction |
| Counterbalancing | Building a balanced (non-random) treatment order such as ABBA to offset known time trends. | counterbalanced design |
| Condition- vs time-dependent crossover rule | Whether to switch treatments once the patient returns to baseline (condition-dependent) or after a fixed wait (time-dependent). | crossover rule |
| On/off kinetics | How quickly a drug starts working and wears off, which sets how long a washout must be. | rapid onset / offset; drug half-life |

## 3. Randomization, blinding, and conduct

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Randomization | Assigning treatment or treatment order by chance to remove systematic bias. | random allocation / assignment |
| Block (permuted-block) randomization | Randomizing within small balanced groups so each treatment appears equally often and bad runs are avoided. | blocking; restricted randomization; counterbalancing |
| Stratified randomization | Randomizing separately within subgroups (site, sex) to keep those factors balanced. | stratification; minimization; covariate-adaptive allocation |
| Allocation concealment | Hiding the upcoming assignment from those enrolling participants so it cannot be foreseen or manipulated. | concealed allocation |
| Blinding | Concealing which treatment is given from patients and/or assessors to prevent expectation bias. | masking; single / double / triple blind |
| Placebo / placebo effect | An inert comparator made to look identical to the active treatment; the placebo effect is improvement wrongly credited to that inert treatment. | sham / dummy treatment; expectancy response |
| Active comparator / control | The reference treatment (placebo, usual care, or an established drug) against which the test treatment is judged. | comparator; control condition; reference arm |
| Equipoise | Genuine uncertainty about which compared treatment is better, which ethically justifies randomizing. | clinical equipoise |

## 4. Analysis models

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Linear mixed-effects model | A regression combining population-level fixed effects with subject-specific random effects, suited to repeated measurements per person. | mixed / mixed-effects model; multilevel model (MLM); hierarchical linear model (HLM); random-effects model; LME; MMRM (mixed model for repeated measures) |
| Multilevel / hierarchical model | A model with parameters nested at several levels (measurements within patients within studies) that share information across levels. | nested model; HLM; intensive hierarchical model |
| Fixed effect | A model term whose coefficient is a single constant assumed to apply to everyone. | population-average / marginal effect |
| Random effect | A model term whose value varies across subjects, treated as a draw from a distribution. | subject-specific effect; random coefficient |
| Random intercept | A per-subject baseline offset letting each individual start at their own level. | subject-specific intercept |
| Random slope | A per-subject deviation in how strongly the outcome changes with a predictor (e.g., over time). | subject-specific slope; random trend |
| Design matrix | The table of predictor values mapping coefficients to observations; separate ones for the fixed (X) and random (Z) parts. | model matrix; X / Z matrix |
| Marginal model | A model describing the *population-average* outcome rather than a specific individual's. | population-average model |
| Conditional model | A model describing the outcome for a specific individual given their own random effects. | subject-specific model |
| Generalized estimating equations (GEE) | A method for correlated/repeated outcomes that estimates the population-average effect using a chosen working correlation, staying valid even if that correlation is mis-specified. | GEE; Zeger-Liang / marginal approach; PGEE (penalized) |
| Working correlation structure | The assumed (possibly approximate) within-subject correlation pattern GEE uses. | independence / exchangeable / autoregressive working structure |
| Generalized linear (mixed) model | A regression for non-normal outcomes (counts, proportions, binary), optionally with random effects. | GLM; GLMM |
| Ordinary least squares (OLS) | Fitting a regression by minimizing the sum of squared residuals. | least-squares regression |
| Generalized least squares (GLS) | Regression that accounts for correlated or unequal-variance errors, giving the best linear unbiased estimate. | weighted least squares (related); BLUE estimation |
| Analysis of covariance (ANCOVA) | A regression that adjusts the outcome for a baseline covariate to improve precision. | baseline adjustment |
| ANHECOVA | Covariate adjustment that includes all treatment-by-covariate interactions, guaranteeing an efficiency gain. | analysis of heterogeneous covariance |
| Change score | The outcome minus its baseline value, a simple (often inefficient) way to adjust for baseline. | difference score; baseline subtraction |
| Proportional hazards (Cox) model | A survival regression in which covariates multiply a baseline rate of an event over time. | Cox regression; coxph |
| Proportional odds model | An ordinal-outcome regression assuming covariate effects are constant across category thresholds. | cumulative logit model |
| Dynamic regression model | A regression that explains the present outcome using lagged (past) values of itself and of predictors, absorbing autocorrelation. | autoregressive distributed-lag model; dynamic modelling |
| Paired t-test | A test comparing two related measurements (e.g., one patient on A vs B) by analyzing their differences. | matched-pairs / Student's t-test |
| Serial t-test | A t-test for one person's repeated data that corrects the variance and degrees of freedom for serial correlation. | -- |
| Randomization (permutation) test | A significance test that compares the observed effect with effects from many random reassignments of the treatment schedule. | permutation test; randomization-test reference set |
| Visual analysis | Judging a single-case treatment effect by inspecting the graphed time series (level, trend, variability) rather than by formal statistics. | visual inspection; graphical analysis |
| Nonoverlap index | A single-case effect size based on how little the treatment-phase data overlap the baseline data. | overlap statistic (see Tau-U, NAP, PND in Section 9) |
| Meta-analysis | Statistically pooling results across studies (or N-of-1 trials) into one overall estimate. | research / quantitative synthesis; mega-analysis |
| Individual participant data (IPD) meta-analysis | Pooling the raw per-patient data across studies, usually via mixed models, rather than just their summaries. | one-stage meta-analysis |
| Network meta-analysis | Combining trials of different treatment pairs to compare treatments never tested head-to-head. | mixed-treatment / indirect comparison |
| Intention-to-treat (ITT) | Analyzing participants by their assigned treatment regardless of adherence, preserving randomization. | ITT (modified ITT excludes those with no outcome data) |
| Per-protocol analysis | Analyzing only participants who followed the protocol, which can introduce bias. | completers / per-protocol population |

## 5. Serial correlation and time series

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Serial correlation | The tendency for measurements taken close together in time on the same person to be related rather than independent. | autocorrelation; serial dependence; within-subject / intra-subject correlation |
| Autocorrelation function (ACF) | A function or plot giving the correlation of a series with itself at each successive lag. | -- |
| Lag | The number of time steps separating two observations being compared. | -- |
| First-order autoregressive process, AR(1) | A model where each measurement correlates with the previous one and correlation decays as `rho` raised to the time gap. | AR(1); lag-1 autocorrelation; corAR1 / corCAR1 (continuous-time) |
| ARMA / ARIMA | Time-series models combining autoregressive and moving-average (and, for ARIMA, differencing) terms. | autoregressive integrated moving average |
| Stationarity | The assumption that a series' mean, variance, and correlation structure stay constant over time. | stationary time series |
| Detrending | Removing a systematic upward/downward trend so fluctuations can be studied. | -- |
| Pre-whitening | Removing autocorrelation from a series before analysis, which can accidentally remove a real effect. | whitening |
| Cochrane-Orcutt / Prais-Winsten | Older regression procedures that adjust estimates for first-order autocorrelated errors. | -- |
| Effective sample size | The smaller number of *independent* observations that correlated data are statistically worth; shrinks under positive correlation. | effective number of independent observations |
| Slope-autocorrelation indeterminacy | The difficulty, in short series, of telling apart a genuine trend from trend produced by autocorrelated errors. | confounding of slope and autocorrelation |

## 6. Variance, power, and estimation

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Components of variation | The separate sources of variability in an outcome (between-treatment, between-patient, patient-by-treatment, within-patient). | variance components; sources of variation |
| Between-patient variance | Variability in average outcomes from one patient to another. | between-subject variance; tau-squared (meta-analysis) |
| Within-patient variance | Variability across repeated occasions in the same patient on the same treatment. | within-subject / residual / error variance |
| Intraclass correlation coefficient (ICC) | The fraction of total variability due to differences between subjects; equivalently, how alike a subject's repeated measures are. | ICC; rho; intra-subject correlation |
| Covariance / variance-covariance matrix | A table of the variances and pairwise correlations of repeated measures or parameter estimates. | Sigma / Omega matrix; dispersion matrix |
| Compound symmetry (CS) | A covariance pattern assuming every pair of repeated measures is equally correlated. | exchangeable correlation |
| Unstructured covariance | A covariance pattern estimating every variance and correlation freely with no imposed shape. | general covariance |
| Toeplitz covariance | A covariance pattern where the correlation depends only on how many timepoints apart two measures are. | banded covariance |
| Kronecker correlation | Building a large correlation matrix as the structured product of smaller ones (e.g., within-period x between-period). | Kronecker-product structure |
| Positive definite | A required property of a valid covariance matrix (no impossible negative variances). | PD |
| Heteroscedasticity | Unequal variance of the outcome across conditions or over time. | heteroskedasticity; unequal variances |
| Variogram | A plot of how the variance of differences between measurements grows with the time gap, used to detect autocorrelation. | semivariogram |
| Maximum likelihood estimation (MLE) | Choosing parameter values that make the observed data most probable. | ML; likelihood-based estimation |
| Restricted maximum likelihood (REML) | A maximum-likelihood variant that estimates variance components with less bias by adjusting for the fixed effects. | restricted / residual maximum likelihood |
| Likelihood ratio test (LRT) | A test comparing two nested models by the difference in their fit. | deviance chi-squared test |
| Wald test | A test judging whether a coefficient differs from zero using its estimate divided by its standard error. | Wald statistic |
| Kenward-Roger approximation | A small-sample degrees-of-freedom correction improving fixed-effect tests in mixed models. | KR correction |
| Empirical Bayes / BLUP | A subject-specific estimate that borrows strength from the group mean, shrinking noisy individual estimates toward it. | best linear unbiased predictor; shrinkage estimate |
| Shrinkage | Pulling individual or extreme estimates toward a common value to reduce noise and overfitting. | regularization; penalization; regression to the mean; James-Stein shrinkage |
| Borrowing strength | Using data from other individuals or studies to sharpen a sparse unit's estimate. | partial pooling; shrinkage |
| Penalized splines / B-splines | Flexible curve-fitting building blocks with a roughness penalty, used to model smooth nonlinear effects without overfitting. | B-spline basis; smoothing / roughness penalty |
| Smoothing parameter (lambda) | A tuning knob controlling how wiggly versus smooth a fitted curve is. | penalization / tuning parameter |
| Non-centrality parameter | A quantity in a test statistic's distribution under the alternative hypothesis; larger values mean higher power. | NCP |
| Type I error | Falsely declaring an effect that does not exist (a false positive); its rate is the significance level alpha. | alpha; size of the test; false-positive rate |
| Type II error | Failing to detect a real effect (a false negative); its complement is power. | beta; false negative |
| Statistical power | The probability of correctly detecting a true effect of a given size. | 1 - beta; sensitivity of a test |
| Sample-size / power calculation | Computing how many patients (or periods/measurements) are needed for adequate power or precision. | power analysis; study planning |
| Effect size | A standardized measure of how large a treatment difference is, often in standard-deviation units. | standardized mean difference (SMD); Cohen's d; tau |
| Bias (of an estimate) | A systematic tendency for an estimate to be too high or too low on average. | systematic error |
| Mean squared error (MSE) | An accuracy summary combining an estimator's variance and squared bias. | RMSE (its square root); bias-variance trade-off |
| Confidence interval | A range that, over many repeats, would contain the true value a stated percentage of the time. | CI; interval estimate |
| Coverage | How often computed intervals actually contain the true value across repeated studies/simulations. | interval coverage |
| Monte Carlo simulation | Repeatedly generating random datasets with known truth to study how a method or design behaves. | simulation study |
| Multivariate normal distribution | A bell-shaped joint distribution for several correlated variables, set by a mean vector and covariance matrix, used to generate repeated-measures data. | MVN; joint normal |
| Bootstrap | Resampling the observed data (with replacement) to approximate the sampling distribution of a statistic. | parametric / wild / multiplier bootstrap |
| Pitman asymptotic relative efficiency | A ratio comparing how many subjects two designs need to reach the same power. | relative efficiency |
| Information criterion (AIC / BIC / QIC) | Scores ranking models by fit while penalizing extra parameters; lower is better (QIC is the GEE analogue of AIC). | Akaike / Bayesian / quasi-likelihood information criterion |

## 7. Bayesian methods

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Bayesian inference | Updating prior beliefs with data via Bayes' rule to obtain a posterior probability distribution over parameters. | Bayesian analysis |
| Prior distribution | The probability distribution encoding belief about a parameter before seeing the data. | prior |
| Posterior distribution | The updated distribution for a parameter after combining prior and data. | posterior |
| Posterior predictive distribution | The distribution of new/future data implied by the fitted model, used for prediction and model checking. | predictive distribution |
| Noninformative prior | A deliberately vague prior that lets the data dominate. | flat / diffuse / weak / reference prior |
| Informative prior | A prior carrying real external knowledge or expert opinion. | -- |
| Weakly informative prior | A prior that gently constrains estimates to plausible ranges without strongly steering results. | -- |
| Conjugate prior | A prior whose form yields a posterior in the same family, simplifying computation. | (e.g., inverse-gamma prior for a variance) |
| Credible interval | A Bayesian interval containing the parameter with a stated probability given the data. | posterior interval; Bayesian confidence interval |
| Posterior probability of benefit | The Bayesian probability that a treatment's true effect is favourable, read off the posterior. | -- |
| Hierarchical / partial pooling | Estimating group-specific values as a compromise between each group's own data and the overall average. | partial pooling; Bayesian shrinkage |
| Exchangeability | Treating individuals' effects as interchangeable draws from a common distribution. | -- |
| Markov chain Monte Carlo (MCMC) | A simulation algorithm that draws samples from a complex posterior to enable Bayesian computation. | MCMC; Gibbs / Metropolis sampling; JAGS / WinBUGS / Stan estimation |
| Posterior predictive check | A model check comparing data simulated from the fitted model against the observed data. | Bayesian p-value |
| Prior elicitation | Formally turning experts' beliefs into a prior distribution, especially influential when data are scarce. | expert-opinion elicitation |
| Prior effective sample size | A way to express how much information a prior contributes by equating it to a number of observations. | prior ESS; ECSS / EHSS (effective current / historical sample size) |
| Prior-data conflict | Disagreement between the prior and the observed data, which can distort inference unless the prior discounts it. | prior-likelihood conflict; commensurability |
| Robust / mixture / power / commensurate prior | Priors designed to down-weight historical or external information when it conflicts with current data. | meta-analytic-predictive (MAP) prior; discounting prior |
| Sequential probability ratio test (SPRT) | A method that tests a hypothesis after each new observation and stops once evidence is conclusive. | SPRT; bSPRT (bootstrapped); sequential analysis |

## 8. Missing data and dropout

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Dropout / attrition | Participants leaving the study before completion, leaving incomplete data. | withdrawal; loss to follow-up; discontinuation; noncompletion |
| Differential dropout | Unequal (or systematically different) dropout between treatment arms. | differential / relative attrition |
| Biased / informative dropout | When those who leave differ systematically (e.g., responders quitting when switched to placebo), distorting results. | informative dropout |
| Missing completely at random (MCAR) | Missingness unrelated to any observed or unobserved data. | MCAR |
| Missing at random (MAR) | Missingness explainable by observed data, so it can be validly modelled / imputed. | ignorable missingness |
| Missing not at random (MNAR) | Missingness that depends on the unseen value itself and cannot be fully corrected from observed data. | nonignorable / informative missingness (IM) |
| Last observation carried forward (LOCF) | Filling a missing value with the patient's last recorded value; simple but often biased. | -- |
| Multiple imputation (MI) | Filling missing values several times from their predicted distribution and combining results to reflect uncertainty. | imputation |
| Complete-case analysis | Analyzing only subjects with no missing data, which can bias results. | listwise deletion; completers analysis |
| Inverse probability weighting (IPW) | Reweighting observed cases to stand in for similar dropouts. | inverse-probability weights |
| Selection model | A missing-data model that predicts the chance of dropping out from the (possibly unseen) outcome. | -- |
| Pattern-mixture model | A missing-data model letting completers and dropouts have different outcome distributions. | -- |
| Sensitivity analysis | Re-running an analysis under alternative assumptions to see how robust the conclusions are. | robustness analysis |

## 9. Single-case / psychometric effect measures

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Level change | An abrupt shift in the average outcome from one phase to the next. | mean shift; phase jump; intercept shift |
| Trend / rate change | A change in the slope (rate of improvement or decline) between phases. | slope change; time-by-phase interaction; drift |
| Tau-U | A nonoverlap effect size for single-case data that quantifies improvement, optionally correcting for baseline trend. | -- |
| Nonoverlap of all pairs (NAP) | The proportion of treatment-phase points exceeding baseline points across all pairwise comparisons. | -- |
| Percentage of nonoverlapping data (PND) | The percentage of treatment points beyond the most extreme baseline point. | -- |
| Percentage of all nonoverlapping data (PAND) | The proportion of data that would need removing to eliminate overlap between phases. | -- |
| Standardized mean difference (SMD) | The mean difference between phases or treatments divided by a standard deviation, allowing comparison across scales. | between-case SMD; standardized effect size |
| Multilevel meta-analysis of single cases | Combining many single-case datasets in a hierarchical model to estimate average effects and their variation across cases. | MultiSCED-type analysis |
| Simulation Modelling Analysis (SMA) | A bootstrap-like resampling method generating null replicates to get empirical p-values for single-case data. | -- |
| Piecewise / segmented regression | Fitting separate linear segments to different phases (baseline vs treatment) of a series. | spline / discontinuity growth model |
| Idiographic | Focused on patterns within a single individual. | within-person approach |
| Nomothetic | Focused on general laws or averages across many individuals. | population approach |
| Ecological momentary assessment (EMA) | Repeatedly sampling a person's experiences in real time in their natural environment. | experience sampling |
| Reliability | The degree to which a measurement yields consistent, reproducible results. | (test-retest, internal-consistency reliability) |
| Standard error of measurement | The expected spread of measurement errors around a person's true score. | -- |
| Measurement error | The discrepancy between an observed measurement and the true underlying value. | observation error; noise |
| Floor / ceiling effect | Scores clustering at the lowest or highest possible value, hiding true differences. | -- |

## 10. Biomarkers, precision medicine, and treatment-effect heterogeneity

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Precision medicine | Tailoring treatment to an individual's biological, genetic, physiological, and behavioural profile. | personalized / personalised / individualized / stratified / tailored medicine |
| Patient-by-treatment interaction | Genuine variation among patients in how strongly they respond to a treatment; the core target of personalized medicine. | treatment-effect heterogeneity (HTE); subject-by-treatment interaction; differential / individual response; effect modification |
| Heterogeneity of treatment effects (HTE) | The fact that the same treatment helps different patients by different amounts. | treatment-effect heterogeneity; differential treatment effect |
| Effect modifier | A patient characteristic that changes the size or direction of a treatment's effect. | moderator (psychology); interaction variable; predictive marker |
| Biomarker-treatment interaction | The statistical effect showing a treatment's benefit depends on a patient's biomarker value; the key signal these designs test. | treatment-by-biomarker interaction; predictive interaction; effect modification |
| Predictive biomarker | A measurable patient trait that forecasts who will benefit (differentially) from a specific treatment and guides its selection. | predictive marker / signature / covariate; treatment-effect modifier; theranostic biomarker |
| Prognostic biomarker | A trait that predicts the likely course of disease regardless of treatment. | prognostic factor / variable / covariate; baseline covariate |
| Qualitative (biologic) interaction | An interaction where treatment helps some patients but harms others (the effect reverses direction); cannot be removed by rescaling. | crossover / non-removable interaction |
| Quantitative (statistical) interaction | An interaction in degree only (benefit varies in size but keeps direction); often removable by transforming the scale. | removable / scale-dependent interaction |
| Additivity | When the treatment effect is constant across patients on a chosen scale, so there is no interaction. | additive model |
| Individual treatment effect (ITE) | The difference in outcome a specific patient would have under treatment versus control. | individual-level effect |
| Conditional average treatment effect (CATE) | The average treatment effect within a subgroup defined by particular characteristics. | -- |
| Average treatment effect (ATE) | The mean treatment effect across the whole population. | population / main effect |
| Subgroup analysis | Examining treatment effects within predefined patient subsets; prone to spurious findings if unplanned. | subgroup comparison |
| Risk modeling (predictive HTE) | Stratifying patients by their predicted baseline risk to see who benefits most. | risk-stratified / risk-based subgrouping |
| Effect modeling (predictive HTE) | Building a model with treatment-by-covariate interactions to predict individual benefit. | interaction-based modeling |
| Responder / non-responder | A patient classified as benefiting or not, often by a threshold. | response classification |
| Biomarker-positive / -negative subgroup (B+ / B-) | The patients whose biomarker indicates they are more (B+) or less (B-) likely to benefit. | targeted / non-targeted subgroup |
| Biomarker prevalence | The fraction of the population that is biomarker-positive. | biomarker frequency |
| Companion diagnostic | A laboratory test required to identify patients eligible for a specific targeted therapy. | -- |
| Potential outcomes / counterfactual | The pair of outcomes a unit would show under each treatment (only one is observed); the counterfactual is the unobserved one. | Neyman-Rubin model; potential outcome |
| SUTVA | The assumption that one unit's treatment does not affect another's outcome and there is only one version of each treatment. | stable unit treatment value assumption |
| Pharmacogenomics | Using a patient's genetics to predict and tailor treatment response. | individualized / precision therapy |
| Learning health system | An infrastructure that continually uses routine clinical data to learn which patients benefit from which treatments. | rapid-learning health system |

## 11. Clinical, pharmacology, and outcome terms

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| Pharmacokinetic / pharmacodynamic (PK/PD) model | Models of how a drug's concentration changes over time (PK) and how concentration produces an effect (PD). | PK/PD model; Emax / dose-response model |
| Half-life | The time for a drug's concentration (or effect) to fall by half; governs carryover length. | t-half; elimination rate; clearance |
| Cumulative dose effect | A drug effect that depends on the whole dosing history, weighted toward recent doses. | dose-history dependence |
| Modified Gompertz function | A three-parameter S-shaped curve (maximum response, rate, displacement) modelling a response rising to a plateau. | three-parameter Gompertz |
| Exponential decay | A model where an effect shrinks by half over a fixed time, used for carryover persistence (`exp(-lambda * t)`). | half-life / lambda decay |
| Regression to the mean | The tendency of extreme initial values to drift toward average on re-measurement, mimicking improvement. | -- |
| Bioequivalence | Showing two drug formulations deliver the active ingredient similarly enough to be interchangeable. | average / population / individual bioequivalence |
| Prescribability / switchability | Whether a generic can be given to new patients (prescribability) or substituted in patients already on the brand (switchability). | population vs individual equivalence |
| Minimal clinically important difference (MCID) | The smallest change in an outcome that patients perceive as meaningful. | clinical-importance cut-off |
| Clinical vs statistical significance | Whether an effect is large enough to matter to a patient versus merely unlikely to be chance. | clinical significance / importance |
| Patient-reported outcome measure (PROM) | An outcome reported directly by the patient (symptoms, quality of life). | PRO; self-report outcome |
| Surrogate outcome | An early or proxy measurement standing in for the true outcome of interest. | surrogate endpoint |
| Time-to-event endpoint | An outcome measured as the time until an event occurs, analyzed with survival methods. | survival endpoint; PFS / OS |
| Hazard rate / hazard ratio | The instantaneous risk of an event over time; the ratio comparing that risk between groups. | HR; baseline hazard |
| Absolute risk difference | The plain difference in outcome probability between treated and control patients. | risk difference; absolute treatment effect |
| Relative risk / risk ratio | The ratio of event probability in one group to another. | RR; relative risk reduction (its complement) |
| Odds ratio | A measure comparing the odds of an outcome between two groups. | OR |
| Censoring | When the event of interest is not observed within follow-up (e.g., the study ends first). | administrative / right censoring |
| Comparative effectiveness research (CER) | Research comparing real-world treatments head-to-head to see which works better for whom. | head-to-head comparison; patient-centered outcomes research (PCOR) |
| Evidence-based medicine (EBM) | Practicing medicine by applying the best current research evidence to patient care. | -- |
| Therapeutic window | The dose range that is effective without causing harmful side effects. | therapeutic range |
| Rare / ultra-rare disease | A condition affecting so few people that conventional group trials are infeasible, favouring N-of-1 designs. | low-volume treatment context |
| Actigraphy | Using a wrist-worn motion sensor to objectively estimate sleep and activity. | actigraph; accelerometry |
| Prazosin | An alpha-1 adrenoceptor antagonist studied for PTSD trauma nightmares and sleep disturbance (the project's running clinical example). | alpha-1 blocker |
| PTSD / CAPS | Post-traumatic stress disorder; the Clinician-Administered PTSD Scale (CAPS) is a structured interview rating its symptom frequency and intensity. | post-traumatic stress disorder; CAPS-IV |

## 12. Cross-field synonym map (same concept, different name)

The single concept on the left is named differently depending on whether
a paper is written from a biostatistics, single-case / psychometric, or
precision-medicine perspective.

| Underlying concept | Biostatistics term | Single-case / psychometric term | Precision-medicine term |
|---|---|---|---|
| One-patient repeated-crossover study | N-of-1 trial; single-patient trial | single-case experimental design (SCED); single-subject design | idiographic / individualized clinical trial |
| Genuine person-to-person variation in response | patient-by-treatment interaction; variance component | subject-by-treatment interaction | heterogeneity of treatment effects (HTE) |
| A variable that changes a treatment's effect | treatment-by-covariate interaction; effect modification | moderator | predictive biomarker / effect modifier |
| Subject-specific deviations in a model | random effects; variance components | between-case variance | (individual effects) |
| Pulling individual estimates toward the mean | shrinkage; empirical Bayes / BLUP | -- | borrowing strength / partial pooling |
| Correlation of nearby-in-time measures | serial correlation; AR(1); ICC | autocorrelation | -- |
| Combining many small studies/patients | mixed-model / IPD meta-analysis | multilevel meta-analysis of single cases | aggregated N-of-1 analysis |
| The unobserved alternative outcome | counterfactual; potential outcome | -- | individual treatment effect (ITE) |
| Multiple measurements per unit over time | repeated measures; longitudinal data | single-case time series | -- |
| Standardized magnitude of an effect | standardized mean difference; Cohen's d | Tau-U / NAP / PND (nonoverlap indices) | -- |
| Real-time repeated outcome capture | -- | ecological momentary assessment (EMA) | micro-randomized trial sampling |

## 13. Reporting guidelines and ethics

| Term | Plain-language definition | Synonyms / aliases (by field) |
|---|---|---|
| CONSORT | The standard checklist for transparently reporting randomized controlled trials. | Consolidated Standards of Reporting Trials |
| CENT 2015 | The CONSORT extension giving reporting standards specifically for N-of-1 trials. | CONSORT Extension for N-of-1 Trials |
| SPENT / SPIRIT | SPIRIT is the trial-protocol reporting guideline; SPENT is its N-of-1 protocol extension. | SPIRIT extension for N-of-1 |
| RoBiNT scale | A risk-of-bias / methodological-quality rating scale for single-case and N-of-1 studies. | Risk of Bias in N-of-1 Trials scale |
| PATH statement | The Predictive Approaches to Treatment effect Heterogeneity reporting guidance for HTE analyses. | -- |
| TRIPOD / PROGRESS / ICEMAN | Reporting/appraisal frameworks for prediction models (TRIPOD), prognosis research (PROGRESS), and credibility of subgroup/interaction claims (ICEMAN). | -- |
| EQUATOR Network | An organization promoting accurate, transparent reporting of health research via guidelines. | Enhancing Quality and Transparency Of health Research |
| Operating characteristics | A design's statistical performance profile (type I error, power, bias) across scenarios. | frequentist properties |
| Pre-specification (a priori analysis) | Defining hypotheses and analyses before seeing data, to avoid data-driven false findings. | pre-registration |
| Selective reporting / data dredging | Cherry-picking favourable analyses or subgroups, inflating false positives. | multiplicity; subgroup-fishing |
| Family-wise error rate (FWER) | The probability of at least one false positive across multiple tests. | multiplicity control |
| Alpha allocation / spending | Dividing the overall significance level across hypotheses or interim looks to control multiplicity. | alpha splitting / spending; gatekeeping; fixed-sequence testing |
| Internal vs external validity | Whether conclusions are valid within the studied sample (internal) versus generalize beyond it (external). | generalizability; applicability; representativeness |
| Equipoise | Genuine uncertainty about which treatment is better, ethically justifying randomization. | clinical equipoise |
| Informed consent / IRB approval | Ethical safeguards requiring participant agreement and institutional ethics-board review before research. | ethics-committee approval; written consent |
| Patient and public involvement (PPI) | Involving patients and the public in designing, running, and interpreting research. | patient representatives |

---

## Provenance and epistemic status

- **Corpus.** The `nof1` collection of the project Zotero library holds
  72 curated N-of-1 / precision-medicine articles. Local full text was
  obtained for 46; 42 were successfully machine-read and mined for terms
  (4 were skipped owing to a content-filter false positive on one
  processing batch, but their concepts are redundantly covered by other
  articles). Three of the read files held mismatched or banner-only text
  (an unrelated BMJ news page, an arbovirus review, and a publisher
  download banner) and contributed no terms; one file was mislabelled in
  the catalogue (a Bayesian effective-sample-size paper) and was mined on
  its actual content.
- **Definitions** were extracted by reading the source texts; plain-language
  wording is the lexicon authors' paraphrase, *verified* against the
  mined passages. The statistical / psychometric supplement (Sections 5-9
  in part) draws on standard methodology references (Pinheiro & Bates;
  Fitzmaurice, Laird & Ware; Box, Jenkins & Reinsel; Gelman et al.;
  Imbens & Rubin; de Vet et al.; CENT; Senn) and is *inferred from
  established domain knowledge* rather than verified against a specific
  retrieved source.
- **Synonym groupings** reflect usage observed across the corpus plus
  conventional cross-field equivalences; they are aids to recognition,
  not claims of exact mathematical identity (e.g., a quantitative
  interaction and an additive-scale effect coincide only on a particular
  scale).

---
*Rendered on 2026-06-17 at 17:42 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/docs/29-nof1-precision-medicine-lexicon.md*
