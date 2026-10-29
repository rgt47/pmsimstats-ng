# Referee report: Critical evaluation of the modified Gompertz response function as a parametric template for treatment-response trajectories in N-of-1 trial simulation
*2026-10-28 17:47 PDT*

Target venue (per cover letter): *Statistics in Medicine*, Methods Note
or Tutorial track.
Manuscript source: `analysis/report/07-gompertz-evaluation/report.Rmd`
(Rmd dated 2026-10-08; stored results dated 2026-09-06).
Prior report: `referee-report-2026-06-17.md`. A `bullets.Rmd` draft
scaffold is present in the directory and was not reviewed.

Evidence labels used throughout: VERIFIED (I executed code and observed
the result), INSPECTED (I read the source or artifact), SUSPECTED
(consistent with what I saw but not tested), UNVERIFIED (not checked).

---

## 1. Summary

The manuscript asks whether the parametric family used for the latent
drug-response trajectory in the pmsimstats Monte Carlo engine affects
inference on the biomarker-by-treatment interaction `bm:Dbc`. Four
families (modified Gompertz, symmetric logistic, hyperbolic tangent,
piecewise-linear breakpoint) are matched at the asymptote and at a
week-5 anchor and crossed with two interaction encodings (mean
moderation with trajectory-scaled moderation; covariance moderation)
and two effect sizes ($c_{bm} \in \{0, 0.45\}$) in a 16-cell factorial
at $N = 70$ (35 per path), a two-path open-label plus blinded
discontinuation design, zero carryover, 1,000 replicates per cell. The
headline claim, revised since the previous round, is that the family
effect is architecture-dependent: exactly nil under covariance
moderation (verified, the authors say, by a matched-seed check) and a
real 0.039 power spread under mean moderation, with Gompertz the
lowest-power family. Type I error is reported as "slightly
conservative" (0.019 to 0.038) and the interaction estimate as
family-invariant.

## 2. Overall assessment and recommendation

**Recommendation: Reject in current form; encourage resubmission after
re-running the comparison with common random numbers and rewriting the
findings.**

The previous referee's load-bearing objection (the family manipulation
was inert by construction) has been addressed: the family now drives
the BR mean inside `buildSigma`, and the mean-moderation shift rides
the trajectory. The stored tables reproduce exactly from the committed
replicate file, all 16 cells are positive definite under the project's
current repair-refusal default, and the driver is reproducible from a
stated seed. Those are real improvements.

The central empirical claim nevertheless does not survive
verification. The 16 production cells use a different seed stream for
every family, so each between-family contrast carries a standard error
of about 0.019. I re-ran the four families under common random numbers
(1,000 matched seeds, every cell at $c_{bm} = 0.45$) with the shipped
package code. Under mean moderation the Gompertz and logistic families
returned identical rejection decisions in all 1,000 seeds (paired
difference 0.000), the piecewise-linear family differed in one seed,
and the hyperbolic tangent was lower by 0.010 (paired $z = 3.2$). Under
covariance moderation the pattern was the same: Gompertz and logistic
identical, piecewise-linear within 0.001, hyperbolic tangent lower by
0.008 (paired $z = 2.8$), with eight seeds in which the families
disagreed on the decision. The reported 0.039 Gompertz-logistic gap is
therefore Monte Carlo noise from independent seeding, the claim that
Gompertz is the lowest-power family is wrong, and the claim that the
covariance architecture is "exactly inert" with "bit-for-bit identical"
decisions is also wrong. The architecture-dependence that the title,
abstract, Discussion and Conclusions are built around is not supported;
what the data support is a small (about 0.01) and architecture-
independent sensitivity to the one family whose calibrated curve
actually differs from Gompertz at the sampled visits.

A corrected paper is publishable and arguably more useful than the
present one, because the corrected finding (the test is insensitive to
family at this design, in both encodings, to within one percentage
point) is clean and the method for establishing it (paired common
random numbers) is the one the authors already used for half their
argument. But the rewrite touches every section, several of the
previous referee's comments remain unaddressed, and the Type I error
and bias reporting also need correction, so I cannot recommend
revision of the existing text.

## 3. Significance and novelty

The gap the paper identifies, that no N-of-1 simulation study compares
response-curve families, survived my searches (PubMed and web; see
Section 8). No prior work compares trajectory families for an N-of-1
power simulation. The surrounding claim is overstated, however. The
manuscript says the N-of-1 methodology corpus is "silent on parametric
trajectory shape" and that papers "specify a linear or AR(1) form that
is too coarse to capture saturation." According to PubMed, Percha et
al. (2019, [DOI](https://doi.org/10.2196/12641)) built an N-of-1
simulation model with an explicit saturating response, exponential
approach to a treatment-driven target with separate wash-in and
wash-out time constants, and released software; Gärtner et al. (2023,
[DOI](https://doi.org/10.1186/s12874-023-02012-5)) released a synthetic
N-of-1 data generator with carryover; Liao et al. (2023,
[DOI](https://doi.org/10.1002/sim.9676)) simulate several lag shapes.
None of these compares families, so the novelty claim stands in its
narrow form, but the Introduction's characterization of the field must
be corrected and these papers cited.

The significance as delivered is low because the headline finding is
an artifact (Section 4, M1). The significance of the corrected finding
is modest but genuine for a Methods Note: it tells users of this and
similar simulators that at a saturating design the family can be left
at its default, and it supplies a worked example of why between-DGP
comparisons must be paired. The identifiability half of the
contribution promised in the first submission is now explicitly
deferred, which is honest but leaves a short paper with one result.

## 4. Major comments

Correctness first.

**M1. The headline 0.039 mean-moderation spread, and the ordering that
places Gompertz lowest, are Monte Carlo noise from unpaired seeds.
Under common random numbers the four families are within 0.010 of one
another in both architectures.** VERIFIED (executed).
*Location:* Abstract (Results); Section 5.1 table and text; Section
6.1 to 6.4; Conclusions; cover letter.
*Problem:* The driver assigns each of the 16 cells its own seed offset
(`seed_off` keyed by family, architecture and $c_{bm}$), so the
Gompertz-versus-logistic contrast at $c_{bm} = 0.45$ is a difference of
two independent binomial proportions with SE 0.019. The manuscript
reports the gap as 0.039, "approximately two standard errors," and
builds the architecture-dependence narrative on it, while
acknowledging that it is the largest of six post hoc contrasts. I ran
the shipped `generateData`/`lme_analysis` pipeline with the driver's
exact design, calibration and parameters, for 1,000 matched seeds per
family:

| Architecture | Gompertz | Logistic | Piecewise | Tanh | Discordant seeds |
|---|---|---|---|---|---|
| Mean moderation (trajectory) | 0.782 | 0.782 | 0.781 | 0.772 | 10 |
| Covariance moderation | 0.797 | 0.797 | 0.796 | 0.789 | 8 |

Gompertz versus logistic: zero discordant decisions in 1,000 seeds in
both architectures. Gompertz versus piecewise-linear: one discordant
seed in each. Gompertz versus hyperbolic tangent: 10 and 8 discordant
seeds, all in the same direction (paired $z = 3.2$ and $2.8$). The
mean estimates were $-0.2355$, $-0.2353$, $-0.2347$, $-0.2332$ under
mean moderation. The marginal powers in this run (0.78 to 0.80) differ
from the production cells (0.74 to 0.80) by amounts of the same order
as the "spread" the paper interprets, which is what independent seeding
at $n = 1{,}000$ produces.
*Why it matters:* every interpretive paragraph in the paper (the
"architecture-dependent finding," the "mildly conservative default,"
the "precision effect," the two "mechanisms," the "conditions that
would change the answer," the prazosin reference-cell advice) explains
a quantity that does not exist. The true family effect is about one
percentage point, is confined to the hyperbolic tangent, and is present
in both architectures.
*Remedy:* re-run the production study with common random numbers
across families (the same seed for the same replicate index in every
family, as the authors already did for the covariance check and as
Morris et al. 2019 recommend for method comparisons), report paired
differences with McNemar-type standard errors, and rewrite Abstract,
Results, Discussion and Conclusions around the corrected finding. The
re-run costs about seven minutes on the authors' hardware. Delete the
claim that Gompertz is the lowest-power family.

**M2. The covariance architecture is not "exactly inert," the
rejection decisions are not "bit-for-bit identical," and no shipped
code performs the matched-seed check the paper describes.** VERIFIED
(executed) and INSPECTED (driver).
*Location:* Abstract; Section 4 ("as we verify by a matched-seed
check, cannot move the interaction test"); Section 5.1 ("bit-for-bit
identical rejection decisions"); Section 6.2 ("This is not a
small-effect statement but an exact one"); deviations log,
"Verification."
*Problem:* (a) In 1,000 matched seeds under covariance moderation the
families disagreed on the decision in 8 seeds and the largest p-value
difference against Gompertz was 0.056; in a 40-seed check the
decisions happened to agree while the p-values differed by up to
0.003. The "exactness" argument in Section 6.2 (the family enters only
the BR mean, "to which the interaction test is invariant") is wrong as
stated: the family-specific BR mean curve is a function of visit and
path that is not in the column space of the fixed effects (intercept,
`bm`, linear `t`, `Dbc`, `bm:Dbc`), so it enters the residuals, the
REML variance components, the working CAR(1) parameter, and hence the
GLS weights and the test statistic. Exact invariance would require the
mean shift to be spanned by the design, which it is not. (b) The
driver's `validate_block()` runs 300 seeds under `moderation_scaling =
'constant'` and prints only cell powers; it does not compare p-values
or decisions across families, saves nothing, and does not run the
trajectory-scaled arm. The "bit-for-bit" claim therefore rests on an
unrecorded check.
*Why it matters:* a mathematical claim of exact invariance is made and
is false; a verification is claimed that the repository cannot
reproduce.
*Remedy:* replace "exact" with the empirical statement (family
differences of at most 0.01 under common random numbers, with the
mechanism above), add the paired check to the driver as a saved
artifact, and correct the deviations log.

**M3. At the sampled visits the calibrated logistic is the Gompertz
curve, and the hyperbolic tangent is the logistic at twice the rate;
the design has at most two distinct shapes where the data are
observed.** VERIFIED (computed) for the first part; the second was
VERIFIED algebraically by the previous referee (M2) and is unaddressed.
*Location:* Section 3; Figure 1 and its caption; Section 6.1 ("the
four families span the qualitative range").
*Problem:* With the logistic rate fixed at 0.5 and the anchor at week
5, the logistic-minus-Gompertz deviation at the six on-drug durations
sampled by the design (weeks 4, 8, 12, 16, 17, 18) is at most 0.13
CAPS-IV points against a residual SD of 8; the two curves reach 50, 90
and 99 percent of the asymptote at 4.7/4.7, 9.2/8.9 and 14.8/13.7
weeks. This is why they return identical decisions in 1,000 seeds (M1).
The hyperbolic tangent, $(m/2)(1 + \tanh(r(t - t_0))) = m/(1 +
\exp(-2r(t - t_0)))$, is the logistic with rate $2r$; the manuscript
still calls it a family that "neither the logistic nor the Gompertz
family can span with their standard parameterizations," which is false
(set the logistic rate to 1.0). The von Bertalanffy curve is still
named as a comparator in Section 1.1 ("The comparators are the
symmetric logistic, the von Bertalanffy growth curve, and a
piecewise-linear breakpoint model") and is not used; the previous
referee flagged this (m3).
*Why it matters:* the paper's claim to span "asymmetric, symmetric
smooth, and non-smooth" shapes is not realized at the visits where it
would matter; the only family that differs from Gompertz by more than
one point at any sampled visit is the tanh (by about $-1.0$ at week 4
and $+1.2$ at week 8), and that is the only family that moves power.
*Remedy:* either (i) reframe as a rate-and-shape sensitivity with the
calibration stated as the limiting factor, or (ii) replace tanh with a
structurally different curve (Richards with a free shape parameter,
Hill/Emax in time, or a delayed-onset form) and choose calibration
constraints that leave visible differences at weeks 4 and 8. Remove von
Bertalanffy from the comparator list or use it.

**M4. Type I error is roughly half nominal in every null cell, the
mechanism offered is unsupported, and the null cells are duplicated
across architectures without the paper noticing.** VERIFIED
(computed from the stored replicates) and INSPECTED (`lme_analysis.R`).
*Location:* Abstract; Section 5.2 table and text; Section 6.
*Problem:* The eight null rejection rates are 0.019 to 0.038, pooled
0.027. Against the nominal 0.05 with binomial SE 0.0069 the cells sit
1.7 to 4.5 standard errors low (Gompertz under mean moderation: $z =
-4.5$). "Slightly conservative" understates a size of about half
nominal. The explanation, "consistent with the corCAR1 small-sample
degree-of-freedom accounting," was flagged as SUSPECTED by the previous
referee (m9) and is repeated without a diagnostic; `nlme::lme` uses
containment degrees of freedom with no small-sample correction, so the
conservatism must come from elsewhere (for instance, mismatch between
the random-intercept-plus-CAR(1) working model and the three-component
AR(1) DGP with expectancy-scaled placebo variance). Separately, at
$c_{bm} = 0$ the mean-moderation and covariance-moderation DGPs are
identical (no shift is applied and no coupling is built), so the eight
null cells are four DGPs run twice with different seeds. The paper
treats them as eight; the same-DGP pairs differ by up to 0.016
(logistic, 0.038 versus 0.022, $z = 1.9$), which shows the noise level
and shows that the mean-moderation null ordering (Gompertz 0.019 lowest,
logistic 0.038 highest), which coincidentally tracks the reported power
ordering, is noise. Finally, because the test does not hold its size,
the power comparisons are made at an effective $\alpha$ near 0.027,
which the paper should state.
*Remedy:* pool the duplicated null cells (2,000 replicates per family),
report the size deficit as a finding with its magnitude, either test
the stated mechanism (fit with a fixed correlation or at larger $N$, or
compare against a Kenward-Roger or Satterthwaite reference) or demote
it to a conjecture, and say what the size deficit means for the power
numbers.

**M5. A closed-form reference value for the interaction coefficient
exists, so bias at the alternative can and should be reported.**
VERIFIED (computed) for the reference value; INFERRED for its status
as the LME estimand under trajectory scaling.
*Location:* Abstract (Results, last sentence); Section 4, first
paragraph; Section 5.3.
*Problem:* The paper states that the calibration "yields no
closed-form true $\beta_{bm:D}$ on the outcome scale to difference
against." With the analysis on the raw biomarker and $D_{bc}$ binary at
$t_{1/2} = 0$, the population interaction slope is $-c_{bm}\,
\sigma_{BR}/\sigma_{bm} = -0.45 \times 8 / 15.36 = -0.2344$ points per
biomarker unit under covariance moderation (the conditional mean of
$BR$ on drug given $B_i$) and under constant-scaled mean moderation;
NOTATION.md already records $\theta_{\text{true}} = -c_{bm}\sigma_{BR}$
on the standardized scale for paper 10. Under trajectory scaling the
on-drug mean multiplier is normalized to $\sigma_{BR}$, so the same
value is the visit-averaged target, and the LME's weighting across
visits is the only reason it is approximate. The stored cell means
($-0.2296$ to $-0.2393$) and the common-random-number means
($-0.2332$ to $-0.2355$) sit on this value.
*Remedy:* report bias and its MCSE against $-0.2344$ in every
alternative cell, state the approximation under trajectory scaling,
and delete the claim that no reference exists.

**M6. Several mathematical characterizations of the curves are
reversed or wrong at the calibrated parameters.** VERIFIED (computed).
*Location:* Section 2 (equation 1 and the sentence after it); Section
3.1; Section 3.2.
*Problem:* (a) Section 2 says the Gompertz "rise [is] relatively slower
before the inflection and faster after." The opposite holds: the
inflection is at $t^\ast = \ln d / r = 3.83$ weeks ($d = 5$, $r =
0.42$), where the curve has reached $1/e = 37$ percent of its
asymptote; it then needs 11 more weeks to reach 99 percent (14.8
weeks). The Gompertz asymmetry is a short acceleration phase and a long
deceleration tail. (b) Section 3.1 says that at matched asymptote and
anchor "the logistic family rises faster than Gompertz before the
inflection and slower after." At the calibrated parameters the
logistic lies above the Gompertz on essentially the whole range (below
only between weeks 5.5 and 6, by 0.01 points) and reaches 90 and 99
percent earlier. (c) Section 3.2 presents the tanh's "qualitatively
faster approach to the asymptote" as a property of the family; it is a
consequence of fixing both rates at 0.5 (M3). (d) Equation (1) is
written with the code names `maxr`, `disp`, `rate` as mathematical
symbols; NOTATION.md (decision 5, 2026-10-09) requires $m$, $d$, $r$ in
equations and the code names only in code. The inflection formula and
the three implementations (offset subtraction then rescaling so the
asymptote is restored, which the manuscript's equations (2) and (3)
omit) are otherwise correct.
*Remedy:* correct the three sentences with the numbers above, write
equation (1) in $m$, $d$, $r$, and show the rescaling in equations (2)
and (3) so the "asymptote held constant" statement is true of what the
code does.

**M7. Reproducibility and provenance hygiene: tables remain
hand-typed, the stated runtime does not match the artifact, and the
title page, README, driver header and deviations log describe
superseded studies.** INSPECTED.
*Location:* Section 5 (tables; "597 s"); Reproducibility; `README.md`;
`title-page.md`; driver header; `02-deviations-log.md`.
*Problem:* The previous referee's M5 (results hard-coded rather than
read from the artifact) is not addressed: both result tables are
static markdown, only the introductory figure is computed. The
manuscript reports 597 s wall-clock; the committed summary file
reports 421.5 s (generated 2026-09-06), so the prose was not updated
when the run was repeated. The driver header still says "N=35" and
"Last updated 2026-06-17"; the deviations log says "at N = 35" twice
(NOTATION rule 4 makes $N$ the total, 70). `README.md` says the driver
directory "does not yet exist." `title-page.md` describes "Prototype
Monte Carlo outputs (200 reps x 2 trajectory families under
Architecture B)," says the Results section is brief because the
comparators "have not yet been run," gives a placeholder author and
e-mail, and a word count from 2026-05-08. The protocol footer still
points at `~/prj/alz/...` (previous m6). Two bibliography entries
(`corbo2005shelflife`, `achiron2025piecewise`) are still uncited
(previous m5). The self-citation relationship to Hendrickson et al.
(2020) is still unstated (previous m8).
*Remedy:* wire the tables to `07-gompertz-faithful-summary.txt` or the
replicate file; regenerate all run metadata from the artifact; bring
README, title page, driver header, deviations log and protocol footer
into line with the executed study; remove the orphan entries; add one
sentence on authorship overlap with the evaluated framework.

**M8. The manuscript violates the project's notation and style rules
in ways a cross-manuscript reader will notice.** INSPECTED; the STYLE
grep and `tools/notation-lint.pl` were run (the linter is clean).
*Location:* throughout.
*Problem:* (a) "Architecture-A" and "Architecture-B" appear about
fifteen times; NOTATION.md Part 3 retired the letters. (b) The design
is called "Hybrid OL+BDC" (Sections 4 and 6.1); in NOTATION.md "Hybrid"
is the four-path design with a crossover block, while the driver builds
the two-path OL+BDC design (visits at weeks 4, 8, 12, 16, 17, 18, 19,
20; discontinuation at week 18 or 19; no crossover). The protocol
describes yet another design ("16-week total, 8-week open label +
4-week blinded discontinuation + 4-week crossover"). Three documents,
three designs. (c) STYLE.md section 4 requires a Glossary section;
there is none, and $c_{bm}$, $\beta_{bm}$, the anchor, and the analysis
formula are explained in running prose. (d) The STYLE grep returns one
prose hit ("coupling," line 374). (e) The Abstract contains at least
five sentences over 30 words. (f) Voice: "we shall describe," "the
reader will note," "it is worth stating plainly," "Two features of the
design conspire," a rhetorical question closing Section 3.3, "lives
in," "rides," "loads into," parentheses holding clauses. (g) The
analysis formula is written once with lme4 syntax `(1|ptID)` alongside
`corCAR1`; the fit is `nlme::lme` with `random = ~1|ptID`. (h) The
estimand is written $\hat\beta_{bm:D}$ in most places and
$\hat\beta_{bm:D_{bc}}$ in Section 5.3. (i) The sign convention (rule
5) is not stated, so a reader meets a negative interaction estimate
with no explanation.
*Remedy:* rename the architectures per NOTATION.md, name the design
correctly and reconcile the protocol, add the Glossary, and edit the
prose to STYLE.md section 1.

**M9. The literature framing overstates the field's silence and omits
the closest prior simulators.** INSPECTED; searches in Section 8.
See Section 3 and Section 6. Cite Percha et al. (2019), Gärtner et al.
(2023), Liao et al. (2023), and Wang and Schork (2019,
[DOI](https://doi.org/10.3390/healthcare7030084)), and replace "silent
on parametric trajectory shape" with the accurate statement that these
simulators each fix one response form (exponential approach, step with
lag, AR(1)) and none compares forms.

### Status of the 2026-06-17 comments

| Item | Status |
|---|---|
| M1 inert family manipulation | Addressed (family drives BR mean; trajectory-scaled moderation); but see M1 and M2 above |
| M2 tanh is a rescaled logistic | Not addressed; the contradictory sentence remains |
| M3 sloppiness study promised, not delivered | Addressed by explicit deferral (Section 1.6) |
| M4 undocumented protocol deviations | Addressed; `02-deviations-log.md` exists, but states $N = 35$ |
| M5 hard-coded tables | Not addressed; the inaccurate "read back by chunks" sentence was removed |
| m1 MCSE range | Addressed (0.013 to 0.014) |
| m2 SE of the difference | Addressed (0.019), though the gap itself is noise (M1) |
| m3 von Bertalanffy named, unused | Not addressed |
| m4 tanh "used in PK modeling" unsupported | Not addressed; no source found in my search either |
| m5 orphan bib entries | Not addressed |
| m6 stale protocol path | Not addressed |
| m7 draft scaffolding | Addressed (`report.Rmd` is plain prose; `bullets.Rmd` remains as a separate file) |
| m8 self-citation statement | Not addressed |
| m9 corCAR1 d.f. mechanism asserted | Not addressed |

## 5. Minor comments

1. Abstract: "mean $\hat\beta_{bm:D}$ within 0.004 of zero in every
   null cell." The largest null mean is $-0.00403$ (Gompertz,
   covariance moderation). Say "within 0.005" or report the value.
   VERIFIED.
2. Section 1.2, Figure 1 caption and the computed sentence after it
   are correct: the six on-drug durations are 4, 8, 12, 16, 17, 18
   weeks; the Gompertz reaches 39, 84 and 97 percent of its asymptote
   at the first three on-drug visits; the breakpoint is at week 9.28.
   VERIFIED. The sentence "any answer we report should be read as
   conditional on a visit schedule that samples the discriminating
   region twice" is the most important limitation in the paper and
   belongs in the Abstract and Conclusions.
3. Section 1.3: the statement that the analysis routine returns a
   `bm:t` slope "in the degenerate case where a design leaves no
   within-participant variation in drug state" is consistent with
   `lme_analysis.R` (INSPECTED) but is irrelevant to this design and
   can be cut.
4. Section 4, "8-core `parallel::mclapply`": the summary artifact
   confirms 8 workers. VERIFIED. Note that the stored run and a run on
   another platform are not draw-for-draw identical (CLAUDE.md,
   "Simulation server"); the reproducibility section should say so.
5. Section 4, "Notation for the moderation parameter": the shared
   $c_{bm}$ label is declared, as NOTATION rule 2 requires. Good. The
   paragraph then uses "Architecture-A shift," which rule 2's own
   example avoids.
6. Section 5, "Cell MCSE on power runs from 0.013 to 0.014": correct
   for the alternative cells (0.0126 to 0.0138). VERIFIED.
7. Section 5.3: "SD in the range 0.073 to 0.083" matches the artifact
   (0.0733 to 0.0831). VERIFIED. "Roughly equal between architectures"
   is fair (0.080 to 0.083 versus 0.073 to 0.078), though the
   covariance cells are uniformly tighter; say so with the numbers.
8. Section 6.1, the week-4 and week-8 orderings of the families
   (tanh lowest then highest; Gompertz lowest at week 8) are correct.
   VERIFIED from the trajectory-scaling multipliers (week 4: 0.34,
   0.45, 0.46, 0.49 for tanh, Gompertz, logistic, piecewise; week 8:
   0.97, 0.98, 0.98, 1.10). The paragraph that follows, declining to
   name a mechanism, becomes moot after M1.
9. Section 6.3, "Second, larger or smaller $c_{bm}$ ... Pushing
   $c_{bm}$ to 0.30 (where power lies near 0.5)": no evidence is given
   for the 0.5. UNVERIFIED; cite a source cell or drop the number.
10. Data availability says the data and code "are openly available";
    the title page says the public organization is pending. Reconcile
    before submission.
11. `claudecode.tex` is still included via `in_header`; it defines
    environments the manuscript no longer uses (CLAUDE.md, "Paragraph
    structure"). Remove the include.
12. The Introduction's description of Hendrickson et al.'s parameter
    set ("tabula rasa," "the authors having declined to assume any
    difference") is attributed to the package documentation; cite the
    documentation object or the paper for that characterization.
13. "1{,}000" renders correctly in LaTeX math but appears inside prose
    in the Abstract and Results; use plain "1,000" outside math.
14. `laird1964dynamics`: the bibliography's volume 13, issue 3, pages
    490 to 502 agrees with the PubMed record (PMID 14219541,
    [DOI](https://doi.org/10.1038/bjc.1964.55)). All seven entries I
    checked by citation lookup resolve to a PubMed record (Laird 1964;
    Ebling 1996; Triggs 1996; Jagadeesan 2023; Hockenberry 2017;
    Achiron 2025; Schork 2023). VERIFIED. The remaining entries were
    checked against my own knowledge only (INSPECTED).

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| Sec. 1.1, 1.4: N-of-1 corpus "silent on parametric trajectory shape"; papers use "linear or AR(1)" forms | Overstated. Prior N-of-1 simulators specify explicit saturating or lagged response forms; none compares families | Percha et al. 2019, J Med Internet Res, [DOI](https://doi.org/10.2196/12641); Gärtner et al. 2023, BMC Med Res Methodol, [DOI](https://doi.org/10.1186/s12874-023-02012-5); Liao et al. 2023, Stat Med, [DOI](https://doi.org/10.1002/sim.9676); Wang and Schork 2019, Healthcare, [DOI](https://doi.org/10.3390/healthcare7030084); Blackston et al. 2019, Healthcare, [DOI](https://doi.org/10.3390/healthcare7040137) |
| Sec. 3.2: tanh "has been used in pharmacokinetic modeling for drugs with rapid receptor saturation" | Still uncited (previous m4). My web search found hyperbolic (Emax) saturation models and no tanh time-course model | Cite a specific source or delete |
| Sec. 1.1, contribution statement: von Bertalanffy as comparator | Cited (`vonbertalanffy1957`) but not used (previous m3) | Use or delete |
| Sec. 5.2, 6: conservatism attributed to "corCAR1 small-sample degree-of-freedom accounting" | No source; `nlme` uses containment d.f. | Kenward and Roger 1997, Biometrics 53:983-997 for small-sample d.f.; or demote to conjecture |
| Sec. 1.3, 3: sigmoidal families in longitudinal mixed models, robustness to misspecification | Directly relevant methodological precedent, not cited | Capuano et al. 2016, Stat Methods Med Res 27:863-875, [DOI](https://doi.org/10.1177/0962280216645632) |
| Sec. 1.3: model-selection and identifiability consequences in growth-curve simulation | Optional, recent | Mohsin et al. 2023, J Theor Biol 576:111656, [DOI](https://doi.org/10.1016/j.jtbi.2023.111656) |
| Sec. 4, 5: paired comparison across DGPs with common random numbers | The design recommendation that M1 turns on | Morris, White and Crowther 2019 (already cited); cite the relevant section on using the same datasets across methods |
| `corbo2005shelflife`, `achiron2025piecewise` | In `references.bib`, never cited | Remove |

According to PubMed, I found no citation in the bibliography that is
mischaracterized in the claim it supports; the growth-curve, tumor,
microbial and pharmacology citations are used accurately. The
`morris2019simulation` entry lacks a DOI (10.1002/sim.8086).

## 7. Suggestions for strengthening

1. **Re-run under common random numbers and report paired
   differences.** This is the only change that is strictly necessary,
   and it converts a noisy and wrong headline into a precise and
   correct one. Record the paired check as a saved artifact.
2. **Make the family contrast real at the sampled visits.** Either
   choose a calibration that leaves the families visibly different at
   weeks 4 and 8 (for instance, match the asymptote and the time to
   half-maximum but let the early slope differ), or add a delayed-onset
   or Hill-in-time family whose early visits actually differ from the
   Gompertz, and say in the Abstract that the design samples the
   discriminating region at two visits.
3. **Report size and bias properly.** Pool the duplicated null cells,
   report the size deficit with its magnitude and a tested or
   conjectured cause, and report bias against $-0.2344$.
4. **Add one slice of the deferred grid.** A matched-seed run at
   $t_{1/2} = 0.5$ weeks and at $N = 35$ total would cost about 15
   minutes and would answer the two "conditions that would change the
   answer" that the Discussion can only speculate about.
5. **State the design conditionality.** The finding is a property of
   a schedule with three on-drug visits at or beyond saturation; a
   design with weekly visits in weeks 1 to 8 could give a different
   answer. Say this in the Abstract.
6. **Housekeeping before resubmission.** Wire the tables to the
   artifact, update the runtime, title page, README, driver header,
   deviations log and protocol; remove orphan references and the
   inert header include; add the Glossary and the sign-convention
   sentence; rename architectures and the design per NOTATION.md.

## 8. Scope of this review

- **Verified (executed):**
  - Re-summarized `analysis/data/quick-sim/07-gompertz-faithful-replicates.rds`
    (16,000 rows, `N = 70` throughout) and reproduced every power,
    MCSE, mean-estimate, SD and convergence entry in both manuscript
    tables and in `07-gompertz-faithful-summary.txt`.
  - Rebuilt the covariance matrix for all 16 cells and both paths with
    `buildSigma(..., stopOnNonPD = TRUE)` under the current package:
    all positive definite (minimum eigenvalue 0.14 at $c_{bm} = 0.45$,
    covariance moderation, path A). The stored run did not depend on
    repaired matrices.
  - Recomputed the calibration (anchor 5.92 points at week 5; $m =
    10.99$; logistic $t_0 = 4.295$; tanh $t_0 = 4.829$; breakpoint
    9.28 weeks), the Gompertz inflection (3.83 weeks), the family
    values and deviations at the six sampled durations, the 50/90/99
    percent times, and the trajectory-scaling multipliers.
  - Ran the shipped pipeline under common random numbers: 40 seeds for
    both architectures (first pass), then 1,000 seeds for both
    architectures at $c_{bm} = 0.45$, all four families, with the
    driver's design, parameters, caches and `moderation_scaling`
    settings. Results in M1 and M2.
  - Computed the implied interaction reference value and the pairwise
    and null-cell $z$ statistics from the stored replicates.
  - Ran the STYLE.md terminology grep (one prose hit) and
    `perl tools/notation-lint.pl` (clean).
  - Resolved seven bibliography entries to PubMed records by citation
    lookup.
- **Inspected (read and judged):** `report.Rmd` in full;
  `02-faithful-trajectory-sweep.R`; `R/generateData.R` (`buildSigma`,
  `generateData`, PD handling, the trajectory-scaling and coupling
  branches); `R/utilities.R` (`modgompertz` and the three shape
  functions); the `Dbc` construction and p-value extraction in
  `R/lme_analysis.R`; `00-ademp-protocol.md`; `02-deviations-log.md`;
  `README.md`; `title-page.md`; `cover-letter.md`; `references.bib`;
  the 2026-06-17 referee report; the series 07 entries of
  `docs/reviews/pub_review_whitepaper_2026-08-16.md`; CLAUDE.md,
  STYLE.md and NOTATION.md.
- **Not checked:** I did not re-run the production driver end to end
  or the 300-seed validation block; I did not read
  `01-architecture-a-trajectory-sweep.R` or `bullets.Rmd`; I did not
  verify the `lme` degrees-of-freedom claim beyond reading how
  `lme_analysis` extracts the p-value; I did not test the conjectured
  cause of the size deficit; I did not independently confirm the
  prazosin reference parameters against Hendrickson et al. (2020); the
  Gärtner, Liao, Percha and Capuano papers were judged from their
  PubMed abstracts, not full texts.
- **Searches run (no user identifiers were sent to any service):**
  PubMed: (a) N-of-1 trial simulation with trajectory, response-curve,
  time-course or onset terms, 2015 onward (zero hits on the strict
  query); (b) N-of-1 or single-case simulation with carryover, washout
  or data-generating terms, 2017 onward (12 hits, of which Percha 2019,
  Gärtner 2023, Liao 2023, Wang and Schork 2019, Blackston 2019 and
  Hendrickson 2020 are relevant); (c) Gompertz or sigmoid
  misspecification or model selection in longitudinal mixed models
  (Capuano 2016; Mohsin 2023); (d) onset-of-action and Emax or
  indirect-response models in trial simulation (nothing on point); (e)
  Gompertz or growth-curve comparisons with logistic, Richards or von
  Bertalanffy, 2015 onward (eleven agricultural growth-curve
  comparisons; none in a trial-simulation setting). Web: (f)
  hyperbolic-tangent pharmacodynamic time-course models (no tanh
  source found); (g) N-of-1 simulation response-trajectory family
  sensitivity (Percha 2019, Gärtner 2023 and Liao 2023 again; nothing
  comparing families).

---
*Rendered on 2026-10-28 at 17:47 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/07-gompertz-evaluation/referee-report-2026-10-28.md*
