# Why Architecture B Departs from the Published Implementation {.unlisted .unnumbered}

*2026-09-09 10:57 PDT*

**Author.** pmsimstats team

**Purpose.** This paper takes the differences between the
covariance-moderation process of this package (`covar`) and the
published Hendrickson et al. [1] implementation (`orig`) one at a time,
and asks of each whether it repairs a defect, is a defensible modeling
choice, or follows from another change. It supports the comparison of
the two processes in paper 01.

**Review note (2026-10-07).** Three later results bear on this paper.
(a) The ceilings quoted in Sections 3.3 and 4.3 were computed at
$\rho = 0.7$; the published code uses $\rho = 0.8$, at which the
published construct's ceiling is 0.256 (Hybrid and CO) and 0.262
(OL+BDC) without carryover (`docs/36`, Section 4.2; recomputed in this
review). Under the separable AR(1) form of `docs/35` the ceiling is
0.440 at $\rho = 0.8$, and 0.474 (Hybrid 0.480) at $\rho = 0.7$.
(b) `docs/37` shows that under the published step coupling the
on-drug versus off-drug coupling gap is zero, so the interaction
vanishes at *any* positive half-life; graded coupling,
$g_t = 2^{-t_{sd}/t_{1/2}}$ off drug and 0 before first exposure,
removes the collapse. (c) Paper 01 has since been restructured: the
section numbers cited here for it (2.2.4, B.5 to B.7) refer to its
September 2026 draft. Its current Appendix B lists the code changes to
`58b32a9` (graded coupling, stop instead of repair, anchored
carryover). A copy of `58b32a9` is now held at
`analysis/data/quick-sim/hendrickson-58b32a9/`.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined here and are not part of the canonical set.

| Symbol | Meaning | Status |
|----------|----------------------------------------|------|
| $t$, $s$ | measurement-occasion indices ($s$ local) | canonical |
| $n$ | number of measurement occasions on a path | local |
| $w_t$ | cumulative calendar week of occasion $t$ | local |
| $B$, $BR$, $PB$, $TV$ | biomarker and the three response components | canonical |
| $c_{bm}$ | covariance-moderation correlation: on-drug correlation of $B$ with $BR$ | canonical |
| $\rho$, $\rho_c$ | AR(1) or compound-symmetry within-factor correlation, common or for component $c$ (subscript local) | canonical |
| $A_c$ | within-factor correlation block of component $c$ | local |
| $C_{cc'}$ | cross-factor correlation block between components $c$ and $c'$ | local |
| $c_1$ | same-occasion cross-factor correlation (code `c.cf1t`) | local |
| $c_\times$ | cross-factor correlation at distinct occasions (code `c.cfct`) | local |
| $\tilde r$, $\tilde r_t$ | biomarker coupling vector: correlations of $B$ with $BR_t$ | local |
| $v$ | coupling pattern, $\tilde r = c_{bm} v$ | local |
| $g_t$ | graded coupling factor off drug, $2^{-t_{sd}/t_{1/2}} = e^{-\lambda t_{sd}}$ (`docs/37`) | local |
| $\mathbf{1}_n$ | vector of $n$ ones | local |
| $M$ | correlation block of the three response factors | local |
| $c_{bm}^{\ast}$ | ceiling: largest $c_{bm}$ with a positive definite matrix | local |
| $t_{sd}$, $t_{sd,t}$ | time since discontinuation (at occasion $t$) | canonical |
| $t_{1/2}$ | carryover half-life, weeks | canonical |
| $\lambda$ | carryover decay rate, $\ln 2 / t_{1/2}$ | canonical |
| $\mu_{BR}$ | mean of the drug-response component | local |

Table: Notation index

**Glossary.** *`orig`*: the published implementation at commit
`58b32a9`. *`covar`*: the package's covariance-moderation process
(`dgp_architecture = 'mvn'`, formerly Architecture B). *Step coupling*:
$\tilde r_t = c_{bm}$ wherever $\mu_{BR} \neq 0$, else 0. *Graded
coupling*: $\tilde r_t = c_{bm}$ on drug, $c_{bm} e^{-\lambda t_{sd,t}}$
off drug after exposure, 0 before first exposure. *Compound symmetry
(CS)*: one constant correlation for every pair of occasions. *Ceiling*:
$c_{bm}^{\ast}$ above. *Recursive* versus *anchored* carryover: decay of
the previous, already decayed mean versus decay from the mean at
discontinuation.

## 1. Purpose

Architecture B, which we call `covar` in code and in the appendices of
paper 01, is not the covariance-moderation process that Hendrickson et
al. [1] published. It differs from the published implementation in three
places, and a fourth difference has since been introduced by a
correction of our own. A reader who notices this is entitled to ask
whether each departure was necessary or merely preferred.

We take that question seriously enough to answer it one difference at
a time. Our conclusion is that the three differ in kind. One repairs a
defect that would otherwise make the interaction unidentifiable under
carryover. One is a modeling choice which we believe is clearly
better, though the published alternative is defensible. One is a
consequence of the second, and would be difficult to justify on its
own.

Throughout, `orig` denotes the published implementation at commit
`58b32a9` of 19 November 2020, the last commit before the 2024
revisions. (An earlier version of this sentence called it the last
commit preceding publication. The article appeared on 28 August 2020
[1], so `58b32a9` postdates it; its data-generating code is identical
to that of the initial commit `42ac030`, per `docs/36`, Section 2.)
Section 7 records what happened to that repository afterward, and why
its later history does not settle anything.

## 2. The four differences

| Component | `orig` | `covar` |
|---|---|---|
| Within-factor $A_c$ | $\rho_c$ | $\rho_c^{\lvert w_t - w_s \rvert}$ |
| Cross-factor $C_{cc'}$, $t \neq s$ | $c_\times$ | $c_\times \rho^{\lvert w_t - w_s \rvert}$ |
| Coupling $\tilde r_t$, off drug | $0$ or $c_{bm}$ | $c_{bm} e^{-\lambda t_{sd,t}}$ |
| Carryover on $\mu_{BR}$ | recursive | anchored |

Table: Four differences between `orig` and `covar` by covariance component

Everything else agrees. The same-occasion cross-factor entry, the
on-drug value of the coupling vector, the baseline row and column, and
the biomarker's zero correlation with the two non-response components
are common to both.

Sections 3 through 6 take the four in order of how strongly we would
defend them.

## 3. The coupling vector: a defect, not a preference

This is the one substantive repair, and we would make it even if
nothing else were at stake.

### 3.1 What the published gate does

The published code sets the biomarker-response correlation by asking
whether the response mean at that occasion is nonzero. Where it is,
the correlation takes the full value $c_{bm}$; where it is not, zero.
That is a sensible-looking rule. Its difficulty appears only once
carryover is switched on.

Carryover keeps the response mean nonzero after discontinuation. That
is what carryover is. So the gate, which asks whether any drug effect
is present rather than how much, returns yes at off-drug occasions
too, and assigns them the same $c_{bm}$ it assigns on-drug occasions.

In the Hybrid design at any positive half-life (one week is the
example used below), the result is

$$\tilde r^{\text{orig}} = c_{bm}\mathbf{1}_n,$$

a constant vector, because every Hybrid path begins on drug, so every
off-drug occasion follows exposure and carries a nonzero residual mean
(inspected, `58b32a9` `R/generateData.R`, lines 83 to 92 and 131 to
138). Every occasion, on drug and off, carries the same
biomarker-response correlation.

### 3.2 Why this is fatal rather than merely inaccurate

The biomarker-treatment interaction is identified by $\tilde r$
varying between on-drug and off-drug occasions. That variation is the
whole signal. A constant $\tilde r$ describes a biomarker that
correlates with the response equally whether or not the patient is
taking the drug, which is a biomarker main effect and not an
interaction at all.

The consequence runs the wrong way round. Under the published gate,
any carryover at all flattens $\tilde r$, so at every positive
half-life the interaction is not identifiable through the coupling
channel. The loss is a jump at $t_{1/2} = 0$, not a gradual decline
(derived, `docs/37`; observed in the published results, `docs/36`,
Section 5). But carryover ought to blur the on-drug and off-drug
contrast, not equalize it at full strength. A patient two weeks off
drug should look partly like a treated patient and partly like an
untreated one. The published gate says they look exactly like a
treated patient.

### 3.3 The symptom in the feasible range

The defect leaves a numerical fingerprint which is worth recording,
because it is initially misleading. Evaluating the exact ceiling on
the moderation parameter (section 8) on an evenly spaced schedule of
eight occasions two weeks apart, on drug at the first four and off at
the last four, with $\rho_c = 0.7$, $c_1 = 0.2$ and $c_\times = 0.1$
(verified: recomputed from the closed form of section 8 in the
2026-10-07 review):

| Configuration | Ceiling on $c_{bm}$ |
|---|---|
| `orig`, $t_{1/2} = 0$ | 0.346 |
| `orig`, $t_{1/2} > 0$ | 0.841 |
| `covar`, $t_{1/2} = 1$ | 0.670 |

Table: Ceiling on $c_{bm}$ for `orig` and `covar` on an evenly spaced eight-occasion schedule, $\rho_c = 0.7$

At the published $\rho_c = 0.8$ the `orig` values on this schedule are
0.256 and 0.893, and on the paper 01 Hybrid design they are the same,
0.256 and 0.893 (verified; `docs/36`, Section 4.2). The values in the
table are not those of the published parameters.

At positive carryover the published process admits a *larger*
nominal moderation than Architecture B does. That looks like an
advantage and is the opposite. A constant coupling vector aligns with
the direction in which the covariance is best conditioned, so the
matrix will tolerate a large $c_{bm}$ precisely because that $c_{bm}$
is no longer buying an interaction. A process that can carry a bigger
number while carrying no signal is not thereby preferable.

### 3.4 The replacement

Architecture B gates on the drug indicator itself and decays with time
since discontinuation,

$$\tilde r_t = \begin{cases} c_{bm}, & \text{on drug} \\
c_{bm}e^{-\lambda t_{sd,t}}, & \text{off drug after exposure,} \\
0, & \text{before first exposure,} \end{cases}
\qquad \lambda = \ln 2 / t_{1/2}.$$

In the same Hybrid cell ($c_{bm} = 0.45$, $t_{1/2} = 1$) this gives
$0.45$ on drug against $0.225$ and $0.113$ at one and two weeks off
(verified). The contrast survives, and it
narrows as residual exposure decays, which is what the phenomenon
being modeled actually does.

We should be clear about the standing of this change. It is not a
recovery of what the published code was trying to do; it is a
departure from what the published code does. Our warrant is
Hendrickson et al.'s own stated rationale [1], that the biomarker is
informative because it shares variance with drug-responsive
physiology. A coupling that persists undiminished through a washout
does not express that rationale.

## 4. The within-factor form: a modeling choice we would defend

Compound symmetry assigns one constant to every pair of occasions
within a factor, whatever their separation. Architecture B uses AR(1)
decay in cumulative calendar time.

### 4.1 What compound symmetry asserts

Take the Hybrid schedule, whose occasions fall at weeks
$4, 8, 9, 10, 11, 12, 16, 20$. With $\rho_c = 0.7$:

| Pair | Gap | AR(1) | Compound symmetry |
|---|---|---|---|
| week 9 vs 10 | 1 wk | 0.700 | 0.70 |
| week 4 vs 8 | 4 wk | 0.240 | 0.70 |
| week 4 vs 20 | 16 wk | 0.003 | 0.70 |

Table: AR(1) and compound-symmetry correlations for selected Hybrid visit pairs

(Derived: $0.7^{1}$, $0.7^{4}$ and $0.7^{16}$.)

Compound symmetry asserts that a symptom measurement in week 4
resembles one taken sixteen weeks later exactly as strongly as it
resembles one taken the following week. Under AR(1) the near pair is
some two hundred times more strongly correlated than the far pair;
under compound symmetry the ratio is one.

For a fluctuating condition measured over five months, we do not find
the compound-symmetric claim credible. Symptoms persist and then
drift. A patient's state next week is largely their state this week;
their state next season is not.

### 4.2 What each form implies about the patient

The two forms encode different accounts of why a patient's
measurements resemble one another.

Compound symmetry is the correlation a participant random intercept
induces. It says the only thing making a patient's occasions alike is
a fixed patient-level offset, present equally at every visit and
unchanging. That is a real phenomenon, and in a trial with two or
three widely spaced visits it may be the dominant one.

AR(1) adds the thing that serial measurement is for. States persist
and decay. A bad week tends to be followed by a bad week and not by a
bad month five months later.

An N-of-1 trial measures the same patient eight times over twenty
weeks, with some occasions a week apart and others four. It exists in
order to exploit within-patient temporal structure. A correlation
model that flattens that structure discards the feature the design was
built to use.

### 4.3 The numerical consequence, compared like with like

Holding the coupling vector fixed (graded, $t_{1/2} = 1$) on the
eight-occasion schedule of section 3.3, the ceiling on $c_{bm}$ is
$0.370$ under compound symmetry with a constant cross-factor block. It
is $0.646$ when only the within-factor form is changed to AR(1), and
$0.670$ when the cross-factor block is also made to decay, as in
`covar` (verified, $\rho_c = 0.7$). An earlier version of this
paragraph set 0.370 against 0.670 and attributed the whole gain to the
within-factor form; the cross-factor change contributes about 0.02 of
it. On this schedule the AR(1) form nearly doubles the moderation the
covariance can carry. The gain is schedule- and $\rho$-dependent: on
the paper 01 designs without carryover at $\rho = 0.7$ it is 0.09 to
0.25 with a constant cross-factor block, and at $\rho = 0.3$ compound
symmetry admits the higher ceiling (`docs/34`, Sections 4.3 and 4.4).

We stress the qualification. This is a like-for-like comparison, both
arms using the decaying coupling. It is not the `orig` against `covar`
comparison of section 3.3, which mixes the two changes and reverses
direction. Reporting either figure without saying which comparison it
belongs to would be misleading, and an earlier draft of paper 01's
appendix did exactly that.

### 4.4 What we do not claim

Compound symmetry is not an error. It is the standard implication of a
random intercept, it is exchangeable and therefore indifferent to
occasion ordering, and it has one parameter rather than a parameter
plus a time metric. In a design with few and evenly spaced occasions
the two forms nearly coincide. Our claim is that the Hybrid schedule
is not such a design, and that its one-week and four-week gaps should
not be treated alike.

## 5. The cross-factor form: a consequence, not an independent choice

Architecture B applies the same decay to the cross-factor
off-diagonal. We would not defend this change on its own, and section
2.2.4 of paper 01 does not list it among the four differences it
enumerates. It follows from section 4.

### 5.1 The incoherence of mixing the two

Suppose the within-factor block decays and the cross-factor block does
not. At the Hybrid schedule's widest gap, with $\rho = 0.7$ and
$c_\times = 0.1$:

- BR at week 4 with BR at week 20: $0.7^{16} = 0.003$
- BR at week 4 with PB at week 20: $0.1$

The pharmacological response would be some thirty times more strongly
related to a *different* latent factor sixteen weeks away than to
*itself* sixteen weeks away. We can construct no account of a
patient under which that is true.

### 5.2 The numerical signature

Incoherence of that kind shows up in the covariance. Applying the
cross-factor decay without the AR(1) change lowers the ceiling on
$c_{bm}$ from $0.34$ to $0.16$, the worst of the sixteen
configurations examined in the crossed decomposition of
`docs/54`, Section 5 (worst cell over the
paper 01 designs at $\rho = 0.7$, under the 2024 vendored coupling,
which leaves the first occasion uncoupled; with every occasion coupled
the worst cell is 0.151, in OL+BDC, `docs/34`, Table 4.1). Applying it alongside AR(1) it contributes $+0.004$
to $+0.042$ across the graded-coupling cells (verified, `docs/34`,
Table 4.1: $+0.004$ in Hybrid and $+0.028$ in CO (0.591 to 0.619)
without carryover, up to $+0.042$ in OL+BDC at $t_{1/2} = 1$).
(Corrected 2026-10-08: this sentence first said "roughly $+0.02$ to
$+0.03$", which covers only part of the range in `docs/34`, Table 4.1.
The reading below is unchanged. See
`analysis/report/compendium-mathematics.tex`, Section "Ceilings at the
package value $\rho = 0.7$: the crossed structures".)

Read together, the two numbers say the cross-factor form is not a
lever in its own right. It is nearly free when it agrees with the
within-factor form and expensive when it does not.

### 5.3 Our position

We take the two rows together for consistency and claim nothing more
for the second than that. A reader who rejects the argument of section
4 should reject this one too, and would then have compound symmetry
throughout, which is at least coherent.

## 6. The carryover recursion: our correction, not architecture

The fourth row of the table is not a departure from the published
implementation at all. It records a correction we made in September
2026 to `covar` and to the mean-moderation process, and did not make
to `orig`.

All three processes originally adjusted the response mean at an
off-drug occasion by multiplying the previous occasion's
*already adjusted* mean by a decay factor keyed to
*cumulative* time since discontinuation. Elapsed time is
therefore counted twice, and every off-drug occasion after the first
in a run is over-decayed. At a one-week half-life the second off-drug
occasion receives half the residual effect it should (derived: with
occasions one and two weeks after discontinuation, the recursion
applies $2^{-2}$ to a mean already decayed by $2^{-1}$, giving
$2^{-3}$ where $2^{-2}$ is intended).

The corrected form anchors the decay to the mean at discontinuation,
so the cumulative factor is applied once. We left `orig` alone because
the purpose of that arm is to describe the published code rather than
to improve it.

Two consequences follow for anyone comparing the two processes. The
first off-drug occasion of any run agrees under both forms, so
divergence begins at the second, and designs with long uninterrupted
washouts show it most. And a comparison of `covar` against `orig` now
carries one difference of our making. We recommend restoring the
recursion in `covar` for the headline comparison, so that the four
architectural differences are isolated, and reporting the corrected
version alongside.

*Note added 2026-10-07.* Paper 01 has since taken the other course: it
applies the anchored form to the published code itself, as one of its
listed code changes to `58b32a9` (paper 01, Appendix B.3), rather than
restoring the recursion in `covar`. In all current processes the time
on drug restarts at zero on re-exposure.

## 7. What the repository's later history does and does not settle

The repository was revised in 2024. Commit `8609f12` of 6 May 2024
replaced the step gate with a graded form, added a scale-factor
constant to the carryover adjustment, and left the first measurement
occasion unassigned.

The first of those changes points in the same direction as section 3.
We do not offer it as corroboration. That commit is authored by
Hendrickson but its message describes it as the present author's edits
curated into working code, so the two lines of reasoning are not
independent and it would be circular to cite one in support of the
other. The argument of section 3 rests on identification and on
nothing else.

We record the later history for a different reason. The vendored copy
of the published source in this repository is taken from the 2024
head, not from `58b32a9`, and any re-run of the comparison arm must
use the latter. (Since 2026-10, a copy of the `58b32a9` sources is held
at `analysis/data/quick-sim/hendrickson-58b32a9/`.) A reader who reproduces the appendix worked examples
against the current head will not obtain the published vectors, and
will additionally find the first occasion silently unassigned.

## 8. The exact ceiling used above

The bound quoted in sections 3.3 and 4.3 is not a sweep. Dropping the
baseline row and column, which are those of the identity and factor
out, write $M$ for the block of the response factors and
$\tilde r = c_{bm} v$ for the biomarker coupling extended by zeros.
Positive definiteness requires $M$ positive definite and
$1 - \tilde r^\top M^{-1} \tilde r > 0$, giving

$$c_{bm}^{\ast} = \bigl(v^\top M^{-1} v\bigr)^{-1/2}.$$

The ceiling depends on the correlation structure only through how the
coupling pattern sits relative to $M^{-1}$, which is why the coupling
change of section 3 moves it as much as the correlation change of
section 4 does.

We note one further result, since it corrects a natural but mistaken
intuition. Under compound symmetry the smallest eigenvalue of $M$ is
$(1 - \rho_c) - (c_1 - c_\times)$, which contains neither the number
of occasions nor the number of factors (derived; 0.200 at the
production parameters for $n = 4$ to 32, verified). The response block
itself therefore does not become harder to sustain as occasions are
added.

*Correction (2026-10-07).* An earlier version of this paragraph went
on to say that compound symmetry does not become harder to sustain as
occasions are added, and that paper 01's draft appendix was wrong to
attribute the narrower feasible range to the number of occasions. That
does not follow. The ceiling depends on $v^\top M^{-1} v$, not on the
smallest eigenvalue of $M$, and under compound symmetry it falls with
the number of occasions at a fixed on-drug share, roughly as $n^{-1/2}$:
0.480, 0.346, 0.247 and 0.176 for $n = 4$, 8, 16 and 32 with half the
occasions on drug (verified; closed form in `docs/34`, Section 5.1).
What the eigenvalue result does show is that the loss comes through
the coupling vector, not through the response block.

## 9. Conclusions

In conclusion, four points are to be emphasized.

First, the coupling vector change is a repair and not a preference.
The published gate flattens the interaction channel under carryover,
so that at any positive half-life the Hybrid design carries no on-drug
versus off-drug contrast at all. An interaction that cannot vary with
treatment state is not an interaction.

Second, the within-factor change is a modeling choice, and we would
defend it while acknowledging that compound symmetry is defensible.
The Hybrid schedule mixes one-week and four-week gaps, and a form that
treats them alike discards the temporal structure the design exists to
exploit. Holding the coupling fixed, the change nearly doubles the
feasible range of the moderation parameter on the evenly spaced
schedule of section 3.3 at $\rho = 0.7$; the gain is smaller on the
paper 01 designs and reverses at small $\rho$ (section 4.3).

Third, the cross-factor change is a consequence of the second and
should be presented as such. On its own it lowers the feasible range
substantially, which is the numerical signature of asserting that a
factor is more closely related to a different factor than to itself at
the same remove.

Fourth, the carryover recursion is our correction rather than an
architectural difference, and a comparison intended to isolate the
architectural differences should hold it fixed across arms.

## References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020;2:13. doi:10.3389/fdgth.2020.00013. (Verified via
   Crossref, 2026-10-07.) Code: `github.com/rchendrickson/pmsimstats`,
   commit `58b32a9`.
