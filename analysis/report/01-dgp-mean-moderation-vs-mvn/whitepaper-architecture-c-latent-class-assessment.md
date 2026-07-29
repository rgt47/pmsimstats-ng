---
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - \linespread{0.97}
  - \setlength{\parskip}{0.45em}
---

# White paper: is Architecture C biologically motivated? A latent-class assessment, with application to the prazosin/PTSD setting

*2026-07-29 10:38 PDT*

**pmsimstats team.** A cross-paper note relating the combined DGP
architecture of paper 01 to the latent-class formulation of paper 03,
and assessing which specification best fits the prazosin/PTSD
reference application. Written to be readable without prior exposure
to either paper; terms are defined at first use and collected in a
glossary at the end.

## 1. Purpose and scope

Paper 01 defines a third data-generating architecture, Architecture C,
which switches on a mean-moderation channel (parameter $c_{bm,A}$) and
a differential-correlation channel (parameter $c_{bm,B}$) at the same
time. Its stated justification is mechanistic: the two channels are
'mechanistically distinct and need not scale in proportion', so a drug
could activate one, both, or neither.

Paper 03 develops a latent-class account of the same phenomenon. This
note asks what paper 03's theory implies for paper 01's Architecture
C. The conclusion qualifies paper 01's justification: combining the
two channels is defensible, but not for the reason paper 01 gives, and
not in the form paper 01 implements.

Two results in this note are new and were computed for it. Section 5
works out, by simulation, which combinations of the two parameters a
two-class biological mechanism can actually produce, and finds that
the answer is a restricted wedge that excludes pure Architecture B
entirely. Section 9.1 resolves a calibration question about the
existing code and finds that paper 01's two architectures are not
matched on effect size as closely as the manuscript implies.

A note on terminology. A *data-generating process* (DGP) is the recipe
the simulation uses to manufacture fake trial data. The *analysis
model* is the statistical model we then fit to that fake data, exactly
as we would fit it to a real trial. Power is estimated by generating
many trials, fitting the analysis model to each, and counting how
often the test comes out significant. The subject of this note is the
recipe, not the analysis.

## 2. What we are estimating, and the three things that can change it

Everything here follows from one question: what is the trial actually
testing?

The analysis model is fixed across all of paper 01:

```
Sx ~ bm + t + Dbc + bm:Dbc
```

where `Sx` is the symptom score, `bm` is the patient's baseline
biomarker, `t` is time, and `Dbc` is drug exposure, equal to 1 while
the patient is on drug and decaying smoothly toward 0 after they stop.
The term of interest is `bm:Dbc`, the interaction. Its coefficient
answers the clinical question: **does the drug work better in patients
with a higher biomarker?** That coefficient is the *estimand*, the
quantity we are trying to estimate, and the trial succeeds when we can
distinguish it from zero.

Power for that test is governed, to a good approximation, by

$$\text{power} \approx f\left(\frac{\text{size of the true
coefficient}}{\text{standard error of its estimate}}\right)$$

so anything the recipe does must show up in one of three places:

1. **Signal: the size of the coefficient.** How much does a one-unit
   increase in the biomarker actually change the drug's effect on the
   outcome? This is the numerator.
2. **Noise: the scatter around it.** How variable are patients'
   responses once you know their biomarker and drug status? More
   scatter means a larger standard error, which is the denominator.
3. **Information: how much the off-drug timepoints contribute.** The
   interaction is identified by contrasting on-drug against off-drug
   observations. If carryover blurs that contrast, or if the biomarker
   stops predicting anything once the drug is withdrawn, the off-drug
   observations stop pulling their weight.

These three buckets organize the whole note. Two recipes that agree in
all three give the same power; two that differ in any one may not.

### 2.1 The quantities that populate the three buckets

**Conditional mean.** The average outcome among patients who share a
given biomarker value, written $E[BR \mid B = b]$. Read it as: 'among
patients whose biomarker equals $b$, what is the mean biological
response?' The word 'conditional' means only 'restricting attention to
patients with that biomarker value'.

**Conditional mean gradient.** How fast that average changes as the
biomarker increases. If patients with a biomarker of 2 respond one
unit better than patients with a biomarker of 1, the gradient is 1 per
biomarker unit. *This is why we care:* the gradient is essentially
what the `bm:Dbc` coefficient estimates. It is bucket 1, the signal.
Two recipes with the same gradient feed the same signal into the test,
whatever else differs between them.

**Conditional variance.** How spread out responses are among patients
who share a biomarker value, written $\mathrm{Var}(BR \mid B = b)$.
*This is why we care:* it is the residual scatter the model must see
through, so it drives the standard error. It is bucket 2, the noise.

Buckets 1 and 2 are what statisticians call the *first moment* (the
mean) and the *second moment* (the variance and covariance) of a
distribution. When paper 01 says Architecture A puts the signal 'in
the first moment' and Architecture B puts it 'in the second moment',
this is the distinction being drawn. Section 3 shows that the
distinction, stated that way, is not accurate.

**Carryover and washout.** *Washout* is the period after a patient
stops the drug. *Carryover* is residual drug effect persisting into
it. The *carryover half-life* is how long that residual takes to fall
by half. Bucket 3 concerns what happens to the biomarker's predictive
value during washout.

## 3. What Architectures A and B actually differ in

The usual summary is that Architecture A encodes the biomarker effect
in the mean and Architecture B in the variance. That is not what the
two recipes do.

**On drug, the two produce the same conditional mean gradient.** Under
Architecture A the recipe draws the response independently of the
biomarker and then adds a shift of $c_{bm} z_B \sigma_{BR}$, where
$z_B$ is the patient's biomarker in standard-deviation units and
$\sigma_{BR}$ is the response standard deviation. Under Architecture B
the recipe instead builds a correlation between biomarker and response
into the joint draw, giving a conditional mean of
$\mu_{BR} + c_{bm}(\sigma_{BR}/\sigma_B)(b - \mu_B)$. Substituting
$z_B = (b - \mu_B)/\sigma_B$ shows these are the same quantity. The
calibration comment at `R/generateData.R:129-131` states the
correspondence explicitly, and it is the reason the two architectures
are comparable at matched $c_{bm}$ at all.

**So both are identical in bucket 1.** They differ in buckets 2 and 3:

- **Bucket 2, the noise.** Under Architecture B, knowing the biomarker
  genuinely tells you something about the response, which narrows the
  range of plausible responses:
  $\mathrm{Var}(BR \mid B) = \sigma_{BR}^2 (1 - c_{bm}^2)$, shrunk.
  Under Architecture A the biomarker and response are drawn
  independently and the effect is bolted on as a fixed shift, sliding
  each patient's response up or down without making them any more
  predictable: $\mathrm{Var}(BR \mid B) = \sigma_{BR}^2$, unchanged.
- **Bucket 3, the off-drug information.** Architecture A's shift is
  controlled by a hard on/off switch, `onDrug := (tod > 0)`
  (`R/generateData.R:141`, `:166`), so it is exactly zero at every
  off-drug timepoint. Architecture B's correlation fades gradually, as
  $c_{bm} e^{-\lambda t_{sd}}$, where $t_{sd}$ is time since
  discontinuation.

The accurate description is therefore: **Architectures A and B are two
different noise-and-washout completions of one shared signal.** This
reframing decides the biological question, because it means
Architecture C's two parameters are not independent knobs. They both
turn the same signal dial.

## 4. What a latent-class mechanism produces

### 4.1 The vocabulary

A *latent class* is an unobserved subgroup. 'Class' means a discrete
category, such as responder versus non-responder; 'latent' means we
never observe which category a patient is in. It is the formal version
of the clinical intuition that a drug works in some patients and not
others.

A *gating function*, written $\pi(B)$, is the probability that a
patient with biomarker value $B$ belongs to the responder class. It
'gates' the drug effect in the sense of a valve: it controls how
likely the effect is to be switched on for that patient. Its shape is
the whole story:

- A *flat* gating function means the biomarker carries no information
  about who responds.
- A *step* gating function means the biomarker is a perfect classifier
  with a threshold: everyone above a cutoff responds, nobody below.
- A *graded* function sits in between: higher biomarker values make
  response more likely without guaranteeing it.

The *gating slope* measures how sharply the function rises, that is,
how close it is to a step. This single quantity turns out to determine
which architecture is the better approximation.

### 4.2 The three signatures

Suppose the truth is: each patient belongs to a hidden class
$C_i \in \{0, 1\}$; the chance of being a responder is $\pi(B_i)$; and
on drug the two classes have mean responses $\mu_0$ and $\mu_1$,
separated by $\Delta = \mu_1 - \mu_0$.

We never observe the class, so what the trial sees is the average over
both possibilities, weighted by how likely each is. That produces one
signature in each bucket:

1. **A mean gradient (bucket 1).** The observed average response at
   biomarker value $b$ is $\mu_0 + \Delta\,\pi(b)$: the non-responder
   mean, plus the extra benefit, weighted by the chance of being a
   responder. Because $\pi$ rises with $b$, so does average response.
   The steepness of that rise is roughly $\Delta$ times the gating
   slope, so a bigger class difference and a sharper gate both
   strengthen the signal.
2. **Variance that depends on the biomarker (bucket 2).** The spread
   of responses is $\Delta^2 \pi(b)(1 - \pi(b))$ plus ordinary
   within-class variability. The extra term is largest where
   $\pi(b) \approx 0.5$. The intuition is direct: a patient whose
   biomarker leaves them on a coin flip between the two classes is the
   least predictable, while patients at either extreme are nearly
   certain to be in one class and vary less. Variance that changes
   with a predictor this way is called *heteroscedasticity*.
3. **A class label that persists (bucket 3).** A patient's class is a
   property of the patient, not of their current drug status. It does
   not wash out. Paper 03 identifies this as the reason some
   information about who responds survives discontinuation and stays
   recoverable from the patient's full trajectory.

Both architectures capture signature 1. **Neither captures signature
2:** Architecture B shrinks variance by the same factor for every
patient, Architecture A does not shrink it at all, and the truth
shrinks it by an amount depending on where the patient sits on the
gating curve. **Neither captures signature 3:** both are built
timepoint by timepoint from correlations and shifts, with no
persistent label attached to the patient.

### 4.3 Summary comparison

| Recipe | Bucket 1: signal | Bucket 2: noise | Bucket 3: off drug |
|---|---|---|---|
| Architecture A | gradient $c\,z_B\sigma_{BR}$ | unchanged | exactly zero |
| Architecture B | same gradient | shrunk by $1-c^2$ | fades as $e^{-\lambda t_{sd}}$ |
| Architecture C | sum of both gradients | partial shrinkage | mixed |
| Two-class mixture | gradient $\Delta\pi'(b)$ | varies with $b$ | label persists |

## 5. Worked example: which parameter pairs can a mixture actually produce?

Section 3 argued that Architecture C's two parameters are coupled
rather than independent. That argument can be checked numerically.

The experiment is simple. Simulate a genuine two-class mixture with a
known gating slope, then ask what parameter value each architecture
would need in order to reproduce what the mixture produced. For
Architecture A the relevant summary is the regression slope of
response on standardized biomarker, divided by the within-class
response SD. For Architecture B it is the observable correlation
between biomarker and response. Call these $c_A$ and $c_B$. Two
million patients per cell were simulated, so the values below are
precise to about three decimal places.

**Table 1.** Gating slope against the implied parameter pair, at class
separation $\Delta = 1$ within-class SD, with half of patients in each
class. The last two columns give the conditional SD of response at the
centre of the biomarker distribution and two SDs out.

| Gating slope | $c_A$ | $c_B$ | $c_B / c_A$ | SD at $b=0$ | SD at $b=2$ |
|---|---|---|---|---|---|
| 0.25 | 0.061 | 0.055 | 0.895 | 1.119 | 1.107 |
| 0.5  | 0.118 | 0.106 | 0.895 | 1.120 | 1.091 |
| 1.0  | 0.206 | 0.184 | 0.895 | 1.119 | 1.049 |
| 2.0  | 0.302 | 0.270 | 0.896 | 1.117 | 1.012 |
| 4.0  | 0.366 | 0.327 | 0.895 | 1.116 | 0.995 |
| 8.0  | 0.390 | 0.349 | 0.894 | 1.117 | 0.998 |

Three things are visible.

**First, the two parameters rise together and their ratio never
moves.** Sharpening the gate from nearly flat to nearly a step raises
$c_A$ six-fold and $c_B$ six-fold, but the ratio stays pinned at
0.895. The mixture traces a straight ray out from the origin. It never
travels along either axis.

**Second, heteroscedasticity appears exactly as predicted, and it is
modest.** At a nearly flat gate the conditional SD is the same at the
centre and in the tail. As the gate sharpens, the tail SD falls to
about 1.00 while the centre stays at about 1.12, so patients near the
50/50 point end up roughly 12 percent more variable than patients at
the extremes. That is a real departure from both architectures, both
of which assume constant conditional variance, but at this class
separation it is not a dramatic one. This is a point in favor of the
correlation approximation in the graded regime, and it is consistent
with paper 03's claim that a single multivariate normal reproduces a
mixture to second order.

**Third, paper 03's bound is exact.** Paper 03 proves that the
achievable correlation satisfies $c_{bm}^2 = f_B \times f_{BR}$, the
product of the between-class variance fractions in biomarker and in
response.
Computing that bound for each row reproduces the simulated $c_B$ to
three decimals in every case (0.055, 0.106, 0.185, 0.271, 0.326,
0.348 against the observed values above). The proposition holds
numerically.

### 5.1 What sets the direction of the ray

Table 1 fixes the class separation. Varying it shows what controls the
ray's angle.

**Table 2.** Implied parameter pairs at two gating slopes, across
class separations.

| $\Delta$ | $c_A$ (slope 1) | $c_B$ (slope 1) | $c_A$ (slope 8) | $c_B$ (slope 8) | ratio |
|---|---|---|---|---|---|
| 0.5 | 0.103 | 0.100 | 0.195 | 0.190 | 0.970 |
| 1.0 | 0.206 | 0.185 | 0.390 | 0.348 | 0.895 |
| 1.5 | 0.309 | 0.247 | 0.583 | 0.467 | 0.801 |
| 2.0 | 0.412 | 0.292 | 0.778 | 0.550 | 0.708 |
| 3.0 | 0.620 | 0.345 | 1.167 | 0.648 | 0.555 |

The ratio depends only on the class separation, and is identical at
both gating slopes. So the structure is clean: **class separation sets
the direction of the ray, gating slope sets how far along it you
travel.** The closed form behind the ratio column is

$$\frac{c_B}{c_A} \;=\; \frac{\sigma_w}{\sqrt{\sigma_w^2 +
\Delta^2 p(1-p)}}$$

where $\sigma_w$ is within-class SD and $p$ the responder fraction.
The ratio is 1 when there is no class separation and falls toward 0 as
separation grows.

### 5.2 The consequence for Architecture C's grid

This ratio is always at most 1, so every reachable point satisfies
$c_B \le c_A$. The set of parameter pairs a two-class mechanism can
produce is the **wedge between the 45-degree line and the $c_A$
axis**. Within that wedge:

- **Pure Architecture B is unreachable.** The $c_B$ axis, where the
  mean channel is zero and only the correlation channel is active,
  corresponds to no two-class mixture at all. Producing correlation
  without a mean gradient would require the gating function to be
  uncorrelated with the biomarker, which defeats the purpose of having
  a predictive biomarker.
- **Pure Architecture A is the limiting case** of very large class
  separation, where the mean gradient dominates and the correlation
  contribution becomes negligible by comparison.
- **The mean channel always leads.** For any two-class mixture,
  $c_A > c_B$. A mixture is therefore always closer to Architecture A
  than to Architecture B in these coordinates.

Table 2 also puts paper 01's headline parameter value in context.
Reaching $c_B = 0.45$ requires a nearly step-like gate together with a
class separation of about 1.5 within-class SDs (giving $c_B = 0.467$).
At a moderate gate, even a separation of 3 SDs only reaches
$c_B = 0.345$. So $c_{bm} = 0.45$ is a demanding assumption: it needs
both an unusually clean classifier and a large responder advantage.
And the mixture that delivers $c_B = 0.467$ delivers $c_A = 0.583$ at
the same time. The pure-B cell that carries paper 01's most quoted
result corresponds to no two-class biology whatever.

## 6. Assessment: is combining the two channels biologically motivated?

### 6.1 The unconstrained grid is not produced by any single mechanism

Since both parameters feed the same conditional mean gradient (Section
3), turning both up adds two signals together. Paper 01 observes the
consequence without drawing the inference: at
$(c_{bm,A}, c_{bm,B}) = (0.45, 0.45)$, 'the total interaction effect
on the BR mean is approximately double that of either pure
architecture'.

Across most of the $3 \times 3$ grid, moving toward the interior is
therefore mostly **making the effect bigger**, with noise and washout
structure changing alongside. It is not switching on a second
biological channel while holding the first fixed. Section 5 sharpens
this: the reachable set is a wedge, so the grid's corner cells and
much of its interior correspond to no two-class mechanism at all.
Paper 01's stated reason for preferring the unconstrained grid, that
it 'spans the full $(c_{bm,A}, c_{bm,B})$ space', is the wrong
criterion if much of that space corresponds to no biology.

### 6.2 The mixing-weight version is coherent

Paper 01 defines, then declines to use, a constrained alternative: set
$c_{bm,A} = (1 - \alpha)c$ and $c_{bm,B} = \alpha c$, where $c$ is a
fixed total effect size and $\alpha$ runs from 0 to 1. At $\alpha = 0$
this is pure Architecture A, at $\alpha = 1$ pure Architecture B, and
in between a genuine blend.

This version holds bucket 1 fixed and varies only buckets 2 and 3,
which is exactly what separates the two architectures, so sweeping
$\alpha$ answers a clean question rather than a confounded one. It is
also the sweep paper 03's theory calls for, since paper 03 shows that
a graded gate is well approximated by the correlation recipe whereas a
near-step gate is better approximated by a thresholded mean-moderation
recipe (`03-latent-class-mixture-application/report.Rmd:1202`). The
weight $\alpha$ can thus be read as a stand-in for gating sharpness, a
property of the biology one could investigate, rather than as a
mixture of two separate mechanisms.

One caveat follows from Section 5.2. Because reachable pairs satisfy
$c_B \le c_A$, the biologically supported part of the mixing line is
its lower half, $\alpha \le 0.5$. Values above that correspond to
correlation exceeding mean gradient, which no two-class mechanism
produces.

### 6.3 What no Gaussian recipe captures

Neither pure architecture, and so no point on the mixing line,
reproduces the biomarker-dependent variance or the persistent class
label of Section 4.2. Architecture C should accordingly not be
described as covering the latent-class possibility. If the worry is
that the biomarker tags a responder subtype, the right response is to
simulate the subtype structure directly. Paper 03's
heterogeneous-random-slopes form is the appropriate target, since it
admits latent-class structure in the recipe while leaving the analysis
model an ordinary `nlme::lme` fit.

## 7. When does the correlation approximation fail?

Paper 03 names three regimes in which a single multivariate normal
cannot stand in for a mixture. Each has a distinct empirical
signature, which makes the taxonomy usable rather than merely
cautionary.

**Bimodal responses.** When the two classes are well separated, so
non-responders show near-zero effect and responders a large one, the
distribution of on-drug response has two peaks rather than one. No
single normal distribution matches the tails of a two-peaked
distribution. Paper 03 notes that in this regime power under the true
mixture *exceeds* the normal approximation's prediction, because a
correctly specified mixture analysis can exploit the two peaks. The
clinical analogue is a subpopulation that genuinely does not respond,
as in some targeted oncology indications.

**A strong class-membership gradient.** When the gate approaches a
step, the biomarker behaves almost like a class label. A thresholded
mean-moderation specification then becomes the better approximation,
because the relationship is dominated by its mean component. The
clinical analogue is a biomarker with near-perfect sensitivity and
specificity for the responsive phenotype. Section 5 supports this from
the other direction: as the gate sharpens, both implied parameters
saturate, and the mean channel remains the larger throughout.

**Class-varying covariance.** If autocorrelation, residual variance,
or placebo-response parameters differ between classes, no
single-component normal can represent the structure, and an analysis
assuming constant covariance is misspecified regardless of where it
puts the biomarker interaction. This is the latent-class analogue of
heteroscedastic residual variance, ordinarily handled with `varIdent`
or `weights` in mixed-model practice.

Only the second of these regimes is addressed by moving along the
mixing line. The first and third require the mixture recipe itself.

## 8. Assessment for the prazosin/PTSD setting

Paper 01 concludes that the prazosin application 'is an Architecture B
story' because resting systolic blood pressure (SBP) is 'a proxy for
an underlying neurobiological state, not a direct determinant of drug
pharmacokinetics'. Four considerations qualify that conclusion.

### 8.1 The classification is drawn on the wrong axis

*Pharmacokinetics* (PK) is what the body does to the drug: absorption,
distribution, clearance, and hence the concentration achieved.
*Pharmacodynamics* (PD) is what the drug does to the body: how much
effect a given concentration produces, which depends on how much
target is available to act on.

Paper 01's argument rules out Architecture A on PK grounds. But its
own Architecture A criteria are not confined to PK: the worked
examples include receptor density, 'if a biomarker measures how many
target receptors a patient has, more receptors means proportionally
more drug effect', which is a PD quantity. Noradrenergic tone is a PD
substrate of exactly this kind. Prazosin blocks alpha-1 adrenergic
receptors, so a patient with more alpha-1 signaling to block
plausibly gets a proportionally larger effect. On paper 01's own
taxonomy, that is an Architecture A mechanism.

The section is internally inconsistent here. It lists 'alpha-1
receptor distribution' among the sources of slippage that supposedly
make SBP a mere proxy, but alpha-1 receptor density is the paper's own
canonical Architecture A example. The sound argument for Architecture
B in this setting is not that SBP fails to determine drug levels; it
is that SBP measured in the arm is a noisy readout of noradrenergic
activity in the brain, and prazosin's benefit in PTSD is thought to be
centrally mediated. That is an argument about measurement slippage,
and it coexists with a real PD component rather than replacing it.
Both being true argues for an intermediate mixing weight, not for pure
Architecture B. Section 5.2 reinforces this from theory: pure
Architecture B is not something a responder-subtype mechanism can
produce in the first place.

### 8.2 The clinical evidence points past both architectures to a mixture

The evidence base paper 01 cites is itself the signature of a
responder subgroup diluted by non-responders. The earlier single-site
trials reported benefit; the multisite PACT trial did not replicate,
and paper 01 notes this 'has been interpreted as consistent with
effect heterogeneity, namely that prazosin works in some patients but
not others'. That is a latent-class hypothesis stated in plain
language.

This points the modeling error in an unexpected direction. Under
paper 03's bimodal-response regime (Section 7), power under the true
mixture exceeds what the normal approximation predicts. If the
prazosin responder subgroup is real and reasonably distinct,
Architecture B may be **understating** achievable power, the opposite
of the direction paper 01's carryover warning points.

This is checkable cheaply. Plotting the distribution of on-drug
response in the existing Raskind datasets would distinguish a single
broad peak, which supports a correlation-leaning specification, from
two separate peaks, which supports neither pure architecture and calls
for the mixture recipe. Section 5 gives the quantitative version of
the same check: an estimate of the responder advantage in within-class
SD units and of the responder fraction is enough to locate the
application on the wedge and read off the implied parameter pair.

### 8.3 Carryover here is behavioral, and neither washout model fits

Prazosin is cleared from the blood in a few hours, yet the simulations
use carryover half-lives of 0.5 to 1.0 weeks. Whatever persists across
those weeks is not the drug. It is the behavioral and physiological
state the drug produced: sleep architecture and nightmare frequency
re-equilibrating after discontinuation.

That undermines both washout models, in opposite ways. Architecture B
fades the biomarker-response *association* exponentially, which is
right when the drug effect is gone and the association goes with it.
But under behavioral carryover the benefit itself persists while the
drug is absent, so a responder keeps responding for some weeks and the
association should persist too, plausibly at near full strength
through the first off-drug week. Architecture B's fade is thus too
fast here. Architecture A is worse: its hard on/off switch sets the
biomarker-scaled benefit to exactly zero at every off-drug timepoint
while the average benefit lingers, which is not a coherent account of
behavioral carryover under any reading.

The specification that would fit is a mean-moderation channel fading
with the same exposure profile the analysis model already uses for
`Dbc`, giving
$\Delta BR_{it} = c_{bm,A} z_{B,i} \sigma_{BR} e^{-\lambda t_{sd}}$
off drug instead of zero. **This variant is not implemented.** The
loop at `R/generateData.R:165-170` applies the shift only where
`onDrug` is true, and neither the `constant` nor the `trajectory`
setting of `moderation_scaling` changes that gate.

The consequence deserves emphasis. Because Architecture A's signal
does not decay at all in the recipe, the 1 to 3 percent power losses
paper 01 reports for Architecture A come entirely from the analysis
model's exposure variable shrinking, not from any erosion of the
generated signal. That is a coding choice being reported as a finding
about mean-moderation biology, and it should be labeled as such.

### 8.4 Recommended specification

Three independent considerations point the same way, and all three
suggest the pure-Architecture-B OL+BDC figure of 0.539, a 30.6 percent
loss, is a pessimistic bound rather than a planning estimate:

1. SBP plausibly carries a real pharmacodynamic component alongside
   its proxy role (Section 8.1), so the mean channel is not empty, and
   Section 5.2 shows it is always the larger channel under a
   responder-subtype mechanism.
2. Behavioral carryover implies the association persists longer than
   Architecture B allows (Section 8.3).
3. If the responder subgroup is real and distinct, the true mixture
   supports higher power than the correlation approximation predicts
   (Section 8.2).

Sizing a trial on 0.539 will over-power it. Sizing on the pure
Architecture A figure of 0.730 will under-power it if the gate is
smooth and arm-to-brain slippage is large. The defensible procedure is
to report power along the constrained mixing line at fixed total
effect size, restrict attention to $\alpha \le 0.5$ per Section 6.2,
treat the resulting interval as the planning range, and adopt an
intermediate weight as the working estimate pending evidence on gating
sharpness.

### 8.5 The design decision does not have to wait on any of this

The choice of trial design is robust to the architecture question even
though the power estimate is not.

In paper 01's Architecture C grid, at full signal and one-week
carryover, the Hybrid design attains 0.870 against 0.688 for OL+BDC
and 0.528 for the crossover. In the Section 3.1 tables the Hybrid
design is likewise first under both pure architectures at every
carryover level. Hybrid wins under every specification examined.

The second-versus-third ordering is *not* robust: under Architecture A
at one-week carryover, OL+BDC (0.730) marginally exceeds the crossover
(0.725), whereas under Architecture B the crossover (0.723) exceeds
OL+BDC (0.539) by a wide margin. The architecture question governs the
power estimate and the runner-up ranking, but not the primary
recommendation. For the prazosin application the recommendation is
Hybrid, and it can be made now.

## 9. How to tell which architecture is right

### 9.1 A calibration check on the existing code

An open question from an earlier draft has now been resolved by direct
simulation against the package. Architecture A adds its shift on top
of a draw that already has full variance, so it should inflate total
response variance by a factor of $(1 + c_{bm}^2)$, whereas
Architecture B builds the association into the draw and leaves total
variance alone. If so, the two architectures are not matched on
observable effect size at nominal parity.

Running `generateData` at $N = 200{,}000$ on an eight-timepoint
open-label design confirms it.

**Table 3.** Observable correlation and response SD at an on-drug
timepoint, by architecture and nominal parameter value.

| Architecture | nominal $c_{bm}$ | induced correlation | SD of response |
|---|---|---|---|
| B (MVN) | 0.45 | 0.446 | 4.99 |
| A (mean moderation) | 0.45 | 0.405 | 5.47 |
| B (MVN) | 0.30 | 0.294 | 4.99 |
| A (mean moderation) | 0.30 | 0.281 | 5.21 |

Architecture B lands on its nominal value. Architecture A at nominal
0.45 delivers an observable correlation of 0.405, about 9 percent low,
and inflates the response SD from 4.99 to 5.47. **Paper 01's
head-to-head comparison is therefore not exactly effect-size matched,
and Architecture A is the handicapped arm in every cell**, carrying
both a smaller induced association and more response variance.

The direction of the headline conclusion survives, since Architecture
A is disadvantaged and still shows the greater carryover robustness.
But two things follow. The baseline power gap at no carryover
(0.751 for A against 0.777 for B at OL+BDC) is at least partly this
artifact rather than a property of the architectures. And because
power loss is not scale-free, comparing losses measured from different
baselines contaminates the difference-of-differences contrast to some
degree. The comparison should be re-run with A's parameter adjusted so
the two induce matching observable correlations.

### 9.2 Diagnostics available from real trial data

Three checks discriminate among the candidates without needing a new
trial:

1. **Shape of the on-drug response distribution.** Two peaks indicate
   a separated mixture, for which neither architecture is adequate.
   One broad peak is consistent with graded gating.
2. **Conditional variance across biomarker strata.** Split patients
   into biomarker tertiles and compare response SD. A mixture predicts
   the middle stratum is the most variable (Section 5, Table 1). Both
   architectures predict no difference. The effect is modest, roughly
   12 percent at a class separation of 1 SD, so this check needs a
   reasonable sample.
3. **Does the biomarker still predict response during washout?** If
   the association survives discontinuation at close to full strength,
   the persistent-label picture is right and both architectures'
   washout models are wrong. This is the most decisive check and the
   most directly relevant to the prazosin setting.

None of these has been run here.

## 10. Recommendations

1. Keep Architecture C, but move the production grid onto the
   constrained mixing line at fixed total effect size, restricted to
   $\alpha \le 0.5$. Reserve the unconstrained grid for the two
   boundary cells, which are useful as implementation checks because
   they reproduce the pure architectures exactly.
2. Re-run the Section 3.1 head-to-head with Architecture A's parameter
   adjusted for the calibration gap in Table 3, so the two
   architectures induce matching observable correlations.
3. Implement an exposure-decayed version of the mean-moderation
   channel and re-run Architecture A under it. Until then, describe
   the 1 to 3 percent Architecture A losses as conditional on a hard
   on/off gate.
4. State explicitly that Architecture C does not span the latent-class
   possibility, and cross-reference paper 03's
   heterogeneous-random-slopes recipe as the construct that does.
5. Revise paper 01 Section 4.3 to argue the prazosin classification on
   the central-versus-peripheral measurement axis rather than the
   pharmacokinetic axis, and to acknowledge the pharmacodynamic
   component.
6. Reconcile paper 01 and paper 03 on carryover before either is
   submitted. Paper 03 warns that a covariance-only recipe overstates
   carryover-induced power loss; paper 01's headline result is that
   loss. The two are reconcilable, since paper 01 conditions on the
   biomarker genuinely behaving like Architecture B, but the
   reconciliation should be written down rather than left to the
   reader.

## 11. Epistemic status

**Verified by direct execution.** The calibration results in Table 3
were produced by calling `generateData` from the installed package.
The mixture results in Tables 1 and 2 were produced by simulation at
two million patients per cell. Paper 03's variance-fraction bound was
checked numerically against Table 1 and reproduces the simulated
correlation to three decimals.

**Verified by inspection of source.** The shared on-drug conditional
mean gradient, the binary `onDrug` gate, and the absence of an
exposure-decayed mean channel were confirmed against
`R/generateData.R` at the line references given.

**Analytic but not independently tested.** The closed form for the
ratio $c_B/c_A$ in Section 5.1 is derived, and matches the simulated
ratio column, but has not been checked against a separate derivation.
The three signatures in Section 4.2 are standard results for finite
mixtures.

**Argued, not tested.** The prazosin assessment in Section 8 argues
from mechanism and from the cited clinical literature. No re-analysis
of the Raskind trials was performed, and none of the three diagnostics
in Section 9.2 has been run. All power figures quoted from paper 01
are read as published and were not recomputed.

**Scope limits.** The mixture simulations use a two-class logistic
gate, a normally distributed biomarker, equal class sizes, and
constant within-class variance. Section 7's third regime,
class-varying covariance, is therefore outside the scope of Tables 1
and 2, and the wedge result of Section 5.2 should not be assumed to
hold when class sizes are strongly unequal. The gating intercept was
held at zero throughout, which is what pins the responder fraction at
one half; an offset gate would move the ratio.

## 12. Glossary

**Carryover.** Residual drug effect persisting after discontinuation.

**Conditional mean.** Average outcome among patients sharing a given
biomarker value.

**Conditional variance.** Spread of outcomes among patients sharing a
given biomarker value.

**Data-generating process (DGP).** The recipe used to manufacture
simulated trial data.

**Estimand.** The quantity a study is trying to estimate; here, the
`bm:Dbc` interaction coefficient.

**Gating function.** The probability that a patient with a given
biomarker value belongs to the responder class.

**Gating slope.** How sharply the gating function rises; a measure of
how nearly the biomarker acts as a perfect classifier.

**Heteroscedasticity.** Variance that changes with the value of a
predictor.

**Latent class.** An unobserved patient subgroup, such as responders
versus non-responders.

**Pharmacodynamics (PD).** What the drug does to the body; how much
effect a given concentration produces.

**Pharmacokinetics (PK).** What the body does to the drug; absorption,
distribution, and clearance.

**Washout.** The period after a patient stops the drug.

## Appendix: reproducing Tables 1 and 2

```r
set.seed(20260729)
n <- 2e6; sigma_w <- 1.0
b <- rnorm(n)

row_for <- function(delta, s) {
  cls <- rbinom(n, 1, plogis(s * b))
  br  <- delta * cls + rnorm(n, 0, sigma_w)
  slope <- coef(lm(br ~ b))[2]
  c(delta = delta, s = s,
    cA = slope / sigma_w,
    cB = cor(b, br),
    ratio = cor(b, br) / (slope / sigma_w))
}

grid <- expand.grid(delta = c(0.5, 1.0, 1.5, 2.0, 3.0),
                    s = c(1.0, 8.0))
round(t(mapply(row_for, grid$delta, grid$s)), 3)
```

---
*Rendered on 2026-07-29 at 10:38 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/01-dgp-mean-moderation-vs-mvn/whitepaper-architecture-c-latent-class-assessment.md*
