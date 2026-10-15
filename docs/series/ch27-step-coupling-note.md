# Note: step coupling in the published code

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** Standalone note for the original authors: step coupling and the coupling gap in the published code.  
**Prerequisites:** None (standalone)

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

<!-- table: process, carryover, coupling, allocation, design, N, analysis ID,
t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 55 (whole); 54 Section 3.6.2.

- New to write: the 54 Section 3.6.2 item.
-->

<!-- Spine: docs/55-step-coupling-in-published-code.md (copied verbatim; Phase 2) -->

---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# Step Coupling in the Published pmsimstats Code: The Biomarker Stays Fully Coupled to the Drug Response After the Drug Is Stopped {.unlisted .unnumbered}
*2026-10-09 09:15 PDT*

**Author.** pmsimstats team

**Purpose.** This note describes a second specific feature of the code
that accompanies Hendrickson et al. (2020) [1], the `pmsimstats`
package at commit `58b32a9` [2]. In the simulation the biomarker
predicts the drug response through a correlation between the two.
The code sets that correlation to its full value at every visit where
the expected drug response is nonzero. Once carryover is switched on,
the expected drug response is nonzero at every off-drug visit after
the first exposure, so the biomarker remains as strongly associated
with the drug response when the drug has been stopped as when it is
being taken. The on-drug versus off-drug contrast then carries no
biomarker signal, and the power to detect the biomarker-by-treatment
interaction falls to the test size at any positive half-life. The note
quotes the code verbatim, derives the size of the effect, works two
examples with every number shown, states what the feature does and
does not change downstream, and gives a short correction with a check
that it reproduces the intended coupling exactly. It is a companion to
`docs/53`, which treats the other defect in the same file, the
compounding carryover recursion of the drug-response mean. Like
`docs/53`, it is written to be read on its own, and in particular by
the authors of the article and the original code, who are best placed
to confirm or correct the reading given here.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\listoffigures
\clearpage
```

## Notation index and glossary

Notation follows the project's canonical table
(`analysis/report/NOTATION.md`) and, for the coupling quantities, the
closed-form paper `docs/37`. Symbols marked *local* are defined here
(or in `docs/37`) and are not in the canonical table. (Notation,
2026-10-09: the implied biomarker slope is now $s^{BR}_t$, since
$\phi_t$ is the canonical decay weight; the coupling pattern $v$ is
now $g$, padded $\tilde g$ in the ceiling formula; the
compound-symmetry correlation is $\rho_{cs}$, and expectation is
$\mathbb E$. The coupling quantities
of `docs/37` are now in `NOTATION.md` and are marked canonical.)

| Symbol | Meaning | Status |
|------------|--------------------------------------------------|----------|
| $t$ | visit index within a randomization path | canonical |
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug | canonical |
| $t_{od}$ | cumulative time on drug at a visit (code `tod`); the input to the Gompertz curve | canonical |
| $t_{sd}$ | time since discontinuation (code `tsd`), cumulative within an off-drug run, zero on drug and before first exposure | canonical |
| $t_{1/2}$ | carryover half-life, in weeks (code `carryover_t1half`) | canonical |
| $\lambda$ | decay rate, $\ln 2 / t_{1/2}$ (code `lambda_cor` in pmsimstats-ng) | canonical |
| $B$ | participant-level biomarker (code `bm`; $B_i$ in `NOTATION.md`); $\bar B$ its sample mean | canonical |
| $\mu_{bm}$, $\sigma_{bm}$ | biomarker mean and standard deviation, 124.33 and 15.36 | canonical |
| $BR_t$ | the biological (drug) response component at visit $t$ | canonical |
| $\sigma_{BR}$ | standard deviation of $BR_t$, 8 | canonical |
| $\mu_{BR,t}$ | expected $BR$ at visit $t$ (code `brmeans[p]`) | canonical |
| $\mu_{\text{last}}$ | $\mu_{BR}$ at the last on-drug visit before an off-drug run | canonical |
| $m$, $t_{sd,m}$ | position of an off-drug visit within its run, and its cumulative $t_{sd}$ | local (`docs/53`) |
| $c_{bm}$ | interaction strength of the covariance construct: the on-drug correlation of $B$ with $BR_t$ (code `c.bm`) | canonical |
| $g_t$ | coupling: the fraction of $c_{bm}$ present at visit $t$, so that $\mathrm{Cor}(B, BR_t) = c_{bm}\,g_t$ | canonical |
| $\bar g_{\text{on}}$, $\bar g_{\text{off}}$ | mean coupling over a path's on-drug and off-drug post-baseline visits | canonical |
| $\bar g_{\text{on}} - \bar g_{\text{off}}$ | the coupling gap ($\mathcal G$ in `NOTATION.md`) | canonical |
| $n_{\text{on}}$, $n_{\text{off}}$ | numbers of on-drug and off-drug post-baseline visits of a path | canonical |
| $\beta_{bm:D}$ | the biomarker-by-treatment interaction: the slope of a participant's on-minus-off outcome contrast on $B$ | canonical |
| $s^{BR}_t$ | implied slope of $BR_t$ on $B$, $c_{bm}\,g_t\,\sigma_{BR}/\sigma_{bm}$ (points of $BR$ per biomarker unit; formerly $\phi_t$, which is the canonical decay weight) | local |
| $Y_{it}$, $\mathrm{BL}_i$, $TV_{it}$, $PB_{it}$ | outcome, baseline level and the other two response components | canonical |
| $\Delta_i$ | participant $i$'s on-minus-off outcome contrast, $\bar Y_i^{\text{on}} - \bar Y_i^{\text{off}} = \sum_t a_t Y_{it}$ | canonical |
| $a_t$ | contrast weights, $1/n_{\text{on}}$ on drug and $-1/n_{\text{off}}$ off drug | canonical |
| $\beta_D$ | treatment main effect (the `Db` coefficient of the published model) | canonical |
| $g$, $\tilde g$ | coupling pattern over visits, the vector of $g_t$, and the same vector padded with zeros at the $TV$ and $PB$ entries (formerly $v$) | canonical |
| $M$ | correlation matrix of the response components | canonical |
| $\rho_{cs}$ | common within-component correlation of the published compound-symmetry construct | canonical |
| $s_f$ | the scale factor added in 2024 (code `scalefactor`, default 2) | canonical |
| $c_{bm}^{\ast}$ | ceiling: the largest $c_{bm}$ for which the correlation matrix is positive definite | canonical |

Table: Notation index for this note

**Glossary.**

- **Coupling.** How strongly the biomarker is tied to the drug response
  at a visit, measured as the fraction $g_t$ of the on-drug
  correlation $c_{bm}$ that is present there. The term is ours, not the
  article's.
- **Step coupling.** The rule of the `58b32a9` code:
  $g_t = 1$ wherever $\mu_{BR,t} \neq 0$, and $g_t = 0$ elsewhere.
- **Graded coupling.** The rule recommended here: $g_t = 1$ on drug,
  $g_t = 2^{-t_{sd}/t_{1/2}}$ off drug after first exposure, and
  $g_t = 0$ before first exposure.
- **Coupling gap.** $\bar g_{\text{on}} - \bar g_{\text{off}}$, the
  difference between how strongly the biomarker is tied to the drug
  response on drug and off drug, averaged over a path's visits.
- **E9, E9b.** The paired-difference statistic of `docs/37`: each
  participant's mean outcome on drug minus mean outcome off drug,
  regressed on the biomarker by ordinary least squares. E9b counts
  the baseline as one more off-drug visit, as the published analysis
  does. Both have exact closed-form moments, and E9b reproduces the
  mean estimates and the Hybrid power of the published mixed-model
  analysis (`docs/37`, Section 6).
- **Hybrid design.** The design that the article calls the "N-of-1
  design": eight visits at weeks 4, 8, 9, 10, 11, 12, 16 and 20, with
  open-label titration, a blinded discontinuation and a brief
  crossover, and four randomization paths (`pathA` to `pathD`, here
  paths 1 to 4), with 18, 18, 17 and 17 of $N = 70$ participants.
- **OL+BDC.** Open-label treatment followed by blinded discontinuation;
  two paths. **CO.** Two-period crossover, eight visits 2.5 weeks
  apart; two paths, one of which starts on placebo.
- **Ceiling.** The largest $c_{bm}$ at which the correlation matrix the
  code builds is positive definite; above it the code repairs the
  matrix silently (`docs/36`).
- **The published process.** The data-generating process of the
  `58b32a9` code (`NOTATION.md`). **The refined process.** The
  pmsimstats-ng data-generating process, which among other changes uses
  graded coupling (`NOTATION.md`). Both are used only where a name is
  needed.

**Evidence labels.** Each claim carries one of four labels: *quoted*
(text reproduced from a source), *inspected* (read from source code),
*derived* (follows algebraically) or *verified* (computed by running
code, with the script named in Section 13).

## 1. Summary for the busy reader

**What the code does.** In `R/generateData.R` at `58b32a9`, lines 131
to 138 set the correlation between the biomarker and the drug response
at visit $t$ to $c_{bm}$ if the expected drug response at that visit is
nonzero, and leave it at 0 otherwise (inspected). The expected drug
response used in the test already includes carryover. For any
$t_{1/2} > 0$ it is positive at every off-drug visit after the first
exposure, so those visits receive the full on-drug correlation
(derived; verified in all 40 path and half-life cases of Section 9.3).

**What it should do.** The biomarker is meant to predict the size of
the drug effect. A drug effect that has faded cannot be predicted at
full strength, so the coupling should fade with it:
$\mathrm{Cor}(B, BR_t) = c_{bm}$ on drug,
$c_{bm}\,2^{-t_{sd}/t_{1/2}}$ off drug after exposure, and 0 before
first exposure. The article does not state an off-drug rule (Section
3), so this graded rule is our reading of the pharmacological intent.

**This is a defect in the coupling, triggered by the mean.** The
expected drug response itself is not at fault here; its carryover is
the subject of `docs/53`. What goes wrong is the second moment: the
code uses "is the mean nonzero?" as a switch for the correlation, and
carryover turns the switch on everywhere after exposure.

**How large.** With full coupling on and off drug, the coupling gap is
0 and the slope of the on-minus-off contrast on the biomarker is
exactly 0 (derived). In the closed form for the Hybrid design at
$N = 70$, $c_{bm} = 0.25$, power is 0.721 without carryover and the test
size, 0.050, at every half-life from 0.0125 weeks upward; for E9b, the
proxy for the published analysis, 0.654 falls to 0.081 (verified). On
Hendrickson's own code with the published analysis, power falls from
0.628 at $t_{1/2} = 0$ to 0.078 at 0.1 week (`docs/53`, Section 8.3;
quoted; reproduced in Section 7.3). With graded coupling it is 0.658
at 0.1 week and declines gradually, to 0.486 at 1 week and 0.344 at 2
(quoted, same source). The collapse is a switch, not a decline: it is
discontinuous at $t_{1/2} = 0$, and any positive half-life, however
small, triggers it, until the half-life is so small that the expected
response underflows to exactly zero in floating point (Section 7.2).

**What it does not touch.** The means, the variances and every
correlation not involving the biomarker are unchanged, so any summary
that does not use the biomarker, including the population-average
drug effect, has the same distribution under either rule (derived;
verified, Section 8.3). The rule does change the positive-definiteness
ceiling: under step coupling the Hybrid ceiling jumps from 0.256 to
0.893 at any positive half-life, precisely because the coupling no
longer varies between visits (verified, Section 8.2).

**The fix.** Replace the test on the mean with a coupling computed from
the drug state and $t_{sd}$, which `buildtrialdesign()` already supplies:

```r
if(d[p]$tod>0){
  g<-1
}else if(d[p]$tsd>0 && modelparam$carryover_t1half>0){
  g<-(1/2)^(d[p]$tsd/modelparam$carryover_t1half)
}else{
  g<-0
}
correlations[n1,'bm']<-modelparam$c.bm*g
correlations['bm',n1]<-modelparam$c.bm*g
```

The patched file reproduces the intended coupling to within
$1.1 \times 10^{-16}$ in all 40 cases checked (verified, Section 9.3).
The 2024 upstream revision replaces the step with a ratio of
consecutive means, which removes the collapse but fades at half the
stated half-life, never couples visit 1, and works only because the
means compound (Section 4.5).

## 2. Scope

This note concerns one statement: the off-drug correlation between the
biomarker and the drug response. It does not address the carryover of
the drug-response mean (`docs/53`), the compound-symmetry correlation
structure, the positive-definiteness repair or the analysis model,
except where they interact with the coupling rule; those interactions
are stated in Sections 8 and 10.

## 3. The intended model

The article describes the covariance matrix in one sentence (quoted;
[1], Methods, p. 6):

> The covariance matrix was assembled based on a set of modifiable
> parameters defining the correlation between the baseline biomarker
> and the BR components, the autocorrelation over time (relating the
> value of a factor for one participant at one timepoint to the value
> of that factor for that participant at subsequent timepoints), the
> correlation at a single time point between the 3 factors, and the
> variance of each component.

and describes the carryover of the drug response, but not of its
correlation with the biomarker (quoted; [1], Methods, p. 6):

> The BR factor was set to zero at times when participants had never
> been on active drug; however, a carryover effect was built in such
> that when a participant moved from active drug to being off active
> drug, the value of the BR at the last timepoint on active drug was
> exponentially decayed, with the half-life of this decay being
> maintained as a model parameter.

The parameter is called "the correlation coefficient relating the
biomarker to the biologic response to drug" in the caption of Figure 4
([1], p. 8; quoted), and the package documents `c.bm` as "Correlation
between the biomarker and the biologic response" (`R/generateData.R`,
line 11; quoted). Neither source says what the correlation should be
at an off-drug visit. The rule below is therefore our reading of the
intent, not a quotation.

In words: the biomarker predicts how large a participant's drug
response is. While the drug is taken, the correlation is $c_{bm}$. When
the drug is stopped, the drug response decays with half-life $t_{1/2}$,
and a biomarker that predicts the size of the drug effect should
predict a fading effect with fading strength; "the biomarker predicts
the size of the drug effect, so it cannot predict a drug effect that is
no longer present" (report 01, Methods, "Pharmacological grounding";
quoted). Before the first exposure there is no drug effect and no
coupling. Writing the coupling as $g_t$,

$$
\mathrm{Cor}(B, BR_t) = c_{bm}\, g_t, \qquad
g_t = \begin{cases}
1, & D_{it} = 1, \\
2^{-t_{sd,t}/t_{1/2}} = e^{-\lambda t_{sd,t}}, & D_{it} = 0 \text{ after first exposure}, \\
0, & \text{before first exposure, or } t_{1/2} = 0 .
\end{cases}
$$

The decay of the coupling is the decay the article states for the drug
response, so the ratio of the biomarker-related spread of $BR_t$ to its
mean stays at its on-drug value off drug (derived;
`analysis/report/compendium-mathematics.tex`, section "The
pharmacokinetic argument for graded coupling"). The interaction the
analysis is meant to detect is the difference between on-drug and
off-drug coupling: if the biomarker predicts the response as strongly
off drug as on drug, there is no biomarker-by-treatment interaction to
find (`docs/37`, Section 2; quoted in substance). Carryover in
within-subject designs has a long methodological literature [6]; what
is specific here is that carryover must be specified twice, for the
mean of the drug response and for its association with the biomarker.

The article's discussion of the result is consistent with this reading:
"Although this decrease in power is consistent with expectation, the
magnitude of this drop with even short half-life carryover effects
underscored the critical nature of this parameter" ([1], pp. 11 to 12;
quoted). A decrease in power that grows with the half-life is what
graded coupling produces (Section 7); a drop to the test size at 0.1
week is what step coupling produces.

## 4. The code as published

### 4.1 The coupling lines

`R/generateData.R` at `58b32a9`, lines 131 to 138 (quoted;
[permalink](https://github.com/rchendrickson/pmsimstats/blob/58b32a92b99a8beb145de00c7bd8a0880f06f0f0/R/generateData.R#L131-L138)):

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

The block sits inside the loop `for(c in cl)` that opens at line 101,
so it runs three times with the same effect. The matrix it writes into
starts as the identity (lines 98 to 100), so a visit that fails the
test keeps a biomarker correlation of exactly 0, and the covariance
matrix is formed from it at line 141
([permalink](https://github.com/rchendrickson/pmsimstats/blob/58b32a92b99a8beb145de00c7bd8a0880f06f0f0/R/generateData.R#L97-L101)):

```r
 97    # Turn this into a matrix of correlations
 98    correlations<-diag(length(labels))
 99    rownames(correlations)<-labels
100    colnames(correlations)<-labels
101    for(c in cl){
...
141    sigma<-outer(sds,sds)*correlations
```

### 4.2 Where the means in the test come from

The vector `means` tested on line 134 is assembled at lines 73 to 95;
its drug-response part is built at lines 83 to 93 (quoted;
[permalink](https://github.com/rchendrickson/pmsimstats/blob/58b32a92b99a8beb145de00c7bd8a0880f06f0f0/R/generateData.R#L83-L93)):

```r
 83        brmeans<-modgompertz(d$tod,rp$max,rp$disp,rp$rate)
 84        if(nP>1){
 85          for(p in 2:nP){
 86            if(!d[p]$onDrug){
 87              if(d[p]$tsd>0){
 88                brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^(d$tsd[p]/modelparam$carryover_t1half)
 89              }
 90            }
 91          }
 92        }
 93        means<-c(means,brmeans)
```

with the drug state set at line 55, `d[,onDrug:=(tod>0)]`, and $t_{od}$
and $t_{sd}$ supplied by `buildtrialdesign()` (lines 86 to 106 of that
file; quoted and explained in `docs/53`, Section 4.2). Line 88 is the
compounding recursion of `docs/53`: it decays the previous visit's
already-decayed mean by the cumulative $t_{sd}$. For the present note
only one property of it matters: for $t_{1/2} > 0$ it is strictly
positive at every off-drug visit after an exposure, whatever its exact
value (derived: the Gompertz mean at the last on-drug visit is
positive, and $2^{-x}$ is positive for every finite $x$). The same is
true of the anchored mean that `docs/53` recommends. The step rule's
trigger is therefore pulled by any carryover, compounding or not.

### 4.3 Line-by-line reading

All statements in this subsection are inspected, except where marked.

- **Line 83.** `brmeans` is the Gompertz mean $G(t_{od})$ at every
  visit. $G(0) = 0$, so it is 0 at every off-drug visit before the
  carryover lines run.
- **Lines 85 to 88.** At an off-drug visit after exposure
  ($t_{sd} > 0$), the mean becomes $\texttt{brmeans[p-1]} \times
  2^{-t_{sd}/t_{1/2}}$. With $t_{1/2} > 0$ this is positive. With
  $t_{1/2} = 0$, R evaluates $2^{-t_{sd}/0} = 2^{-\infty} = 0$, so the
  mean stays 0.
- **Line 134.** The test reads the mean of the drug-response variable
  at visit $p$ (`means[which(n1==labels)]`) and compares it with 0.
  Without carryover it is nonzero exactly at the on-drug visits, so
  the test separates on-drug from off-drug visits, as intended. With
  carryover it is nonzero at every visit from the first exposure
  onward.
- **Lines 135 and 136.** Every visit that passes the test receives the
  same correlation, `modelparam$c.bm`; there is no intermediate value.
  The coupling is two-valued: 0 or 1 in units of $c_{bm}$.
- **Consequence (derived).** For $t_{1/2} > 0$ the code's coupling is
  $g_t = 1$ at every on-drug visit and at every off-drug visit after
  first exposure, and $g_t = 0$ only before first exposure. In the
  Hybrid and OL+BDC designs every path starts on drug, so after the
  first visit the coupling is 1 everywhere.

### 4.4 A coupling defect with a mean trigger

It may be natural to read this as a problem with the means, since the
test is written on the means. It is not. The means at the off-drug
visits are small and positive, as carryover intends: at
$t_{1/2} = 0.1$ week the drug-response mean one week after stopping is
0.0099 points, 0.1% of its last on-drug value of 10.19 (verified,
Table 3). The defect is that a quantity of that size switches on a
correlation of full strength. The mean decides *whether* the visit is
coupled; nothing decides *how much*. Correcting the means, as `docs/53`
does, leaves the switch exactly where it was (verified: the pattern of
nonzero means is identical under the published and the anchored means
in all 40 cases of `docs/53`, Section 8.2).

The coupling enters the data only through the covariance matrix. Under
joint normality,

$$
\mathbb E(BR_t \mid B) = \mu_{BR,t} + s^{BR}_t\,(B - \mu_{bm}), \qquad
s^{BR}_t = c_{bm}\, g_t\, \frac{\sigma_{BR}}{\sigma_{bm}},
$$

so $s^{BR}_t$ is the number of points by which a participant's expected
drug response at visit $t$ rises per unit of biomarker (derived). Under
step coupling, at $c_{bm} = 0.25$, $s^{BR}_t = 0.25 \times 8 / 15.36 =
0.1302$ at every visit after first exposure, whether the drug is on or
off, and whether the mean drug response there is 10 points or 0.01
(verified, Table 3). The off-drug mean is the carryover of the drug
effect, and it fades; the step rule makes the biomarker's prediction of
that effect carry over undiminished.

### 4.5 What the 2024 revision changed

Three later commits changed `R/generateData.R` (inspected, GitHub
commit history and per-commit patches; `docs/53`, Section 4.4). Only
one, `8609f12` (2024-05-06, committed by RC Hendrickson with the
message "Ron Thomas' edits, currated into working code";
[commit](https://github.com/rchendrickson/pmsimstats/commit/8609f124eef1770763ecbb4a8978a474c382f615)),
changes the coupling. The current upstream file (`master` at
`06dac83`, 2024-11-27), lines 145 to 154 (quoted;
[permalink](https://github.com/rchendrickson/pmsimstats/blob/06dac83f31ca85c994a6a2c1a0b98bc18ef233bd/R/generateData.R#L145-L154)):

```r
145      # correlation with biomarker
146      for(p in 1:nP){
147        n1<-paste(trialdesign$timeptname[p],"br",sep=".")
148        ## RON THOMAS VERSION:
149        if (p > 1) {
150          n0 <- paste(trialdesign$timeptname[p - 1], "br", sep = ".")
151          mm1 <- means[which(n1 == labels)]
152          mm0 <- means[which(n0 == labels)]
153          correlations["bm", n1] <- correlations[n1, "bm"] <- ifelse(brtest[p], ifelse(brmeans[p] == 0, 0, (mm1 / mm0) * modelparam$c.bm), modelparam$c.bm)
154        }
```

The original lines 134 to 137 follow as comments (lines 156 to 161).
`brtest` is set before the carryover lines, at line 87,
`brtest <- brmeans == 0`, so it is `TRUE` exactly at the off-drug
visits; line 89 names it with `labels[19:26]`, which assumes eight
visits (inspected). Reading line 153 (derived, and verified in 40 cases
to $1.1 \times 10^{-16}$, Section 9.3):

- **On drug, from visit 2:** $g_t = 1$.
- **Off drug, after exposure:** $g_t = \mu_{BR,t}/\mu_{BR,t-1}$, the
  ratio of consecutive carryover-adjusted means. Under the current
  recursion with $s_f = 2$ (line 96) that ratio is
  $2^{-s_f t_{sd,t}/t_{1/2}}$ with the cumulative $t_{sd}$, so
  $g_t = 2^{-2t_{sd,t}/t_{1/2}}$: graded coupling at half the stated
  half-life.
- **Off drug before exposure, or with $t_{1/2} = 0$:** $g_t = 0$.
- **Visit 1:** never assigned, because the block is guarded by
  `if (p > 1)`; its correlation stays at the identity's 0, even though
  visit 1 is on drug in every path of the Hybrid and OL+BDC designs.
  `docs/54`, Section 3.4, analyzes this side effect ("the first
  occasion is never assigned"; that note calls the 2024 vendored copy
  `orig` and uses its line numbers).

So the current upstream code no longer collapses: its coupling gap is
positive at every half-life (Table 6). It is not the intended rule in
three respects (derived; verified, Sections 7 and 10).

1. **It fades twice as fast as the stated half-life**, through the
   default $s_f = 2$, which `generateSimulatedResults()` does not
   override (`docs/53`, Section 4.4; inspected).
2. **It never couples visit 1.** The on-drug mean coupling is
   $(n_{\text{on}} - 1)/n_{\text{on}}$, 0.8 on Hybrid paths 1 and 2 and
   0.75 on paths 3 and 4, so the interaction slope without carryover is
   $-0.101$ instead of $-0.130$, and closed-form E9 power is 0.495
   instead of 0.721 (verified, Table 6). The direct graded rule assigns
   $c_{bm}$ at visit 1.
3. **It is right only because the means compound.** The ratio of
   consecutive means equals $2^{-s_f t_{sd,t}/t_{1/2}}$ only under the
   compounding recursion of `docs/53`. If the means are corrected to
   the anchored form and this rule is kept, the ratio becomes a
   constant factor per visit and the coupling stops fading after the
   first off-drug visit (Section 10, Table 11; `docs/53`, Section 9.5).

| Visit | $D_{it}$ | $t_{sd}$ | 58b32a9 (step) | 2024 upstream (ratio, $s_f = 2$) | Intended (graded) |
|------|------|------|------|------|------|
| OL1 | 1 | 0 | 1 | 0 | 1 |
| OL2 to BD2 | 1 | 0 | 1 | 1 | 1 |
| BD3 | 0 | 1 | 1 | 0.250 | 0.500 |
| BD4 | 0 | 2 | 1 | 0.0625 | 0.250 |
| COd | 1 | 0 | 1 | 1 | 1 |
| COp | 0 | 4 | 1 | 0.0039 | 0.0625 |

Table: Coupling $g_t$ on Hybrid path 1 at $t_{1/2} = 1$ week under the 58b32a9 code, the current upstream code and the intended rule (verified)

## 5. Why it removes the interaction signal

The algebra is short, and it is exact for the paired-difference
statistic E9 of `docs/37` (all statements derived; the source is
`docs/37`, Sections 4.2 and 5.1).

For participant $i$, let $\Delta_i = \bar Y_i^{\text{on}} -
\bar Y_i^{\text{off}} = \sum_t a_t Y_{it}$, with weights
$a_t = 1/n_{\text{on}}$ on drug and $-1/n_{\text{off}}$ off drug. Since
$Y_{it} = \mathrm{BL}_i - (TV_{it} + PB_{it} + BR_{it})$ and only
$BR_{it}$ is correlated with the biomarker, with covariance
$c_{bm} g_t \sigma_{bm}\sigma_{BR}$,

$$
\beta_{bm:D} = \frac{\mathrm{Cov}(\Delta, B)}{\mathrm{Var}(B)}
= -\sum_t a_t\, s^{BR}_t
= -c_{bm}\frac{\sigma_{BR}}{\sigma_{bm}}
  \Big(\bar g_{\text{on}} - \bar g_{\text{off}}\Big).
$$

The interaction slope is the implied biomarker slope on drug minus the
implied biomarker slope off drug, and it is proportional to the
coupling gap. Without carryover $\bar g_{\text{on}} = 1$,
$\bar g_{\text{off}} = 0$, and $\beta_{bm:D} = -c_{bm}\sigma_{BR}/\sigma_{bm}
= -0.1302$ at $c_{bm} = 0.25$.

**Step coupling.** For $t_{1/2} > 0$ every visit after first exposure
has $g_t = 1$ (Section 4.3), so in the Hybrid and OL+BDC designs
$\bar g_{\text{on}} = \bar g_{\text{off}} = 1$, the gap is 0, and

$$
\beta_{bm:D}^{\text{E9}} = 0 \quad \text{for every } t_{1/2} > 0 .
$$

The biomarker slope is the same on and off drug, and the on-minus-off
contrast removes it. E9b, which counts the baseline (where $g_0 = 0$)
as one more off-drug visit, keeps a remnant:

$$
\beta_{bm:D}^{\text{E9b}} = -c_{bm}\frac{\sigma_{BR}}{\sigma_{bm}}\,
\Big(1 - \frac{n_{\text{off}}}{n_{\text{off}} + 1}\,\bar g_{\text{off}}\Big)
= -c_{bm}\frac{\sigma_{BR}}{\sigma_{bm}}\cdot\frac{1}{n_{\text{off}} + 1},
$$

one quarter of the full slope on Hybrid paths 1 and 2
($n_{\text{off}} = 3$) and one fifth on paths 3 and 4. That remnant,
pooled $-0.029$, is what the published mixed model estimates under
carryover: $-0.026$ to $-0.032$ in simulation (`docs/37`, Section 5.1;
quoted).

**Graded coupling.**
$\bar g_{\text{off}} = n_{\text{off}}^{-1}\sum_{\text{off}} 2^{-t_{sd,t}/t_{1/2}}$
is strictly increasing in $t_{1/2}$ and tends to 0 as $t_{1/2} \to 0$,
so the gap starts at 1 and declines smoothly. Because every off-drug
visit in Hybrid has $t_{sd} \ge 1$ week, at $t_{1/2} = 0.1$ week
$\bar g_{\text{off}} \le 2^{-10} \approx 0.001$, and the gap is still
above 0.999. As $t_{1/2} \to \infty$, $\bar g_{\text{off}} \to 1$: the
step rule is the infinite-half-life limit of graded coupling, applied
at every positive half-life (`docs/37`, Section 5.2).

**Power depends on the coupling only through the gap.** The variance of
$\Delta$ involves only the response variances and correlations, which
neither rule touches, so it does not depend on $t_{1/2}$ or on the
coupling rule (`docs/37`, Proposition 1). The residual variance is
$\mathrm{Var}(\Delta) - \beta_{bm:D}^2\sigma_{bm}^2$, and power is a
strictly increasing function of $|\beta_{bm:D}|$, hence of the gap
(`docs/37`, Section 5.4). Neither which off-drug visits are coupled
nor how the coupling is spread over them matters for E9; only the
average does. Under step coupling the gap is 0 at every positive
half-life, the bottom of that curve.

## 6. Worked example 1: Hybrid path 1

### 6.1 Set-up

The design is the article's Hybrid design as defined in the package's
publication vignette (`vignettes/Produce_Publication_Results_1_generate_data.Rmd`,
lines 48 to 61 at `58b32a9`): path 1 (`pathA`) is on drug at visits 1
to 4 and 7 (weeks 4, 8, 9, 10 and 16) and off at visits 5, 6 and 8
(weeks 11, 12 and 20), one, two and four weeks after a discontinuation
(inspected). The parameters are those of the article: $c_{bm} = 0.25$
(below the ceiling, so no matrix is repaired), $\sigma_{BR} = 8$,
$\sigma_{bm} = 15.36$, and the shipped `extracted_rp` drug-response
curve (inspected).

Every coupling below was read from the correlation matrix that the
code itself builds: the published `generateData()` and the patched
version of Section 9 were run with the multivariate normal sampler
replaced by a function that records the mean vector and covariance
matrix it is given, and $g_t$ was computed as
$\mathrm{Cor}(B, BR_t)/c_{bm}$ from that matrix (verified). The
graded column equals $2^{-t_{sd}/t_{1/2}}$ off drug to within
$10^{-16}$.

### 6.2 Step by step at $t_{1/2} = 0.1$ week

- **OL1 to BD2 (weeks 4 to 10), on drug.** The Gompertz means are
  4.281, 9.223, 9.793 and 10.187 points. Both rules: $g_t = 1$,
  $s^{BR}_t = 0.1302$.
- **BD3 (week 11), off drug, $t_{sd} = 1$.** Published mean
  $10.187 \times 2^{-1/0.1} = 10.187 \times 2^{-10} = 0.00995$. It is
  not zero, so the step rule sets $g_t = 1$ and $s^{BR}_t = 0.1302$: one
  biomarker standard deviation (15.36 units) raises the expected drug
  response by 2.0 points, at a visit where the expected drug response
  itself is 0.01 points. Graded: $g_t = 2^{-10} = 0.000977$,
  $s^{BR}_t = 0.000127$.
- **BD4 (week 12), off drug, $t_{sd} = 2$.** Published mean
  $0.00995 \times 2^{-20} = 9.5 \times 10^{-9}$, still not zero: step
  $g_t = 1$. Graded: $2^{-20} = 9.5 \times 10^{-7}$.
- **COd (week 16), on drug again.** $t_{od}$ restarts at 4 weeks; both
  rules $g_t = 1$.
- **COp (week 20), off drug, $t_{sd} = 4$.** Published mean
  $4.281 \times 2^{-40} = 3.9 \times 10^{-12}$: step $g_t = 1$.
  Graded: $2^{-40} = 9.1 \times 10^{-13}$.

| Visit | Week | $D_{it}$ | $t_{sd}$ | $\mu_{BR,t}$ (published) | $\mu_{BR,t} \neq 0$ | $g_t$, step | $g_t$, graded | $s^{BR}_t$, step | $s^{BR}_t$, graded |
|---|---|---|---|---|---|---|---|---|---|
| OL1 | 4 | 1 | 0 | 4.281 | yes | 1 | 1 | 0.1302 | 0.1302 |
| OL2 | 8 | 1 | 0 | 9.223 | yes | 1 | 1 | 0.1302 | 0.1302 |
| BD1 | 9 | 1 | 0 | 9.793 | yes | 1 | 1 | 0.1302 | 0.1302 |
| BD2 | 10 | 1 | 0 | 10.187 | yes | 1 | 1 | 0.1302 | 0.1302 |
| BD3 | 11 | 0 | 1 | 0.00995 | yes | 1 | $9.8 \times 10^{-4}$ | 0.1302 | $1.3 \times 10^{-4}$ |
| BD4 | 12 | 0 | 2 | $9.5 \times 10^{-9}$ | yes | 1 | $9.5 \times 10^{-7}$ | 0.1302 | $1.2 \times 10^{-7}$ |
| COd | 16 | 1 | 0 | 4.281 | yes | 1 | 1 | 0.1302 | 0.1302 |
| COp | 20 | 0 | 4 | $3.9 \times 10^{-12}$ | yes | 1 | $9.1 \times 10^{-13}$ | 0.1302 | $1.2 \times 10^{-13}$ |

Table: Hybrid path 1 at $t_{1/2} = 0.1$ week: published drug-response mean, the step test, coupling and implied biomarker slope $s^{BR}_t$ (points of $BR$ per biomarker unit) under the step and graded rules, $c_{bm} = 0.25$ (verified)

### 6.3 The same path at $t_{1/2} = 1$ week

At one week the off-drug means are large enough to see, 5.094, 1.273
and 0.268 points (the compounding values of `docs/53`; the intended
anchored values are 5.094, 2.547 and 0.268). The step rule gives the
same coupling as at 0.1 week, $g_t = 1$ everywhere; the graded rule
gives 0.5, 0.25 and 0.0625, the fraction of the drug effect that
remains.

| Visit | Week | $D_{it}$ | $t_{sd}$ | $\mu_{BR,t}$ (published) | $\mu_{BR,t} \neq 0$ | $g_t$, step | $g_t$, graded | $s^{BR}_t$, step | $s^{BR}_t$, graded |
|---|---|---|---|---|---|---|---|---|---|
| OL1 | 4 | 1 | 0 | 4.281 | yes | 1 | 1 | 0.1302 | 0.1302 |
| OL2 | 8 | 1 | 0 | 9.223 | yes | 1 | 1 | 0.1302 | 0.1302 |
| BD1 | 9 | 1 | 0 | 9.793 | yes | 1 | 1 | 0.1302 | 0.1302 |
| BD2 | 10 | 1 | 0 | 10.187 | yes | 1 | 1 | 0.1302 | 0.1302 |
| BD3 | 11 | 0 | 1 | 5.094 | yes | 1 | 0.500 | 0.1302 | 0.0651 |
| BD4 | 12 | 0 | 2 | 1.273 | yes | 1 | 0.250 | 0.1302 | 0.0325 |
| COd | 16 | 1 | 0 | 4.281 | yes | 1 | 1 | 0.1302 | 0.1302 |
| COp | 20 | 0 | 4 | 0.268 | yes | 1 | 0.0625 | 0.1302 | 0.0081 |

Table: Hybrid path 1 at $t_{1/2} = 1$ week: published drug-response mean, the step test, coupling and implied biomarker slope under the step and graded rules, $c_{bm} = 0.25$ (verified)

### 6.4 The contrast slope

Path 1 has $n_{\text{on}} = 5$ and $n_{\text{off}} = 3$. By Section 5,
$\beta_{bm:D}^{\text{E9}} = -0.1302\,(\bar g_{\text{on}} - \bar g_{\text{off}})$:

- **Step, any $t_{1/2} > 0$.** $\bar g_{\text{on}} = \bar g_{\text{off}} = 1$,
  gap 0, slope $-0.1302 \times 0 = 0$. E9b:
  $-0.1302 \times (1 - \tfrac34 \times 1) = -0.0325$.
- **Graded, $t_{1/2} = 0.1$.** $\bar g_{\text{off}} = (2^{-10} + 2^{-20} +
  2^{-40})/3 = 0.000326$, gap 0.99967, slope $-0.13015$. E9b $-0.13016$.
- **Graded, $t_{1/2} = 1$.** $\bar g_{\text{off}} = (0.5 + 0.25 +
  0.0625)/3 = 0.2708$, gap 0.7292, slope $-0.1302 \times 0.7292 =
  -0.0949$. E9b: $-0.1302 \times (1 - \tfrac34 \times 0.2708) = -0.1037$.

The same slopes were computed a second way, directly from the code's
covariance matrix as $\mathrm{Cov}(\Delta, B)/\mathrm{Var}(B)$; the two
agree to within $10^{-16}$ (verified). Under the step rule the direct
value is $-2 \times 10^{-17}$, zero to rounding.

| Rule | $t_{1/2}$ (weeks) | $\bar g_{\text{on}}$ | $\bar g_{\text{off}}$ | Gap | E9 slope | E9b slope |
|------|------|------|------|------|------|------|
| Either | 0 | 1 | 0 | 1 | $-0.1302$ | $-0.1302$ |
| Step | 0.1 | 1 | 1 | 0 | 0 | $-0.0325$ |
| Graded | 0.1 | 1 | 0.00033 | 0.99967 | $-0.1302$ | $-0.1302$ |
| Step | 1 | 1 | 1 | 0 | 0 | $-0.0325$ |
| Graded | 1 | 1 | 0.2708 | 0.7292 | $-0.0949$ | $-0.1037$ |

Table: Hybrid path 1: mean coupling, coupling gap and slope of the on-minus-off contrast on the biomarker, step and graded rules (verified)

![Hybrid path 1: coupling $g_t$, the fraction of $c_{bm}$ present in the correlation between the biomarker and the drug response, under the published step rule (dashed, triangles) and graded coupling (solid, circles), at $t_{1/2} = 0.1$ week (left) and 1 week (right). Points are visits; baseline is plotted at 0; shaded intervals are on drug. The step rule is 1 at every visit after the first exposure, whatever the half-life, so there is no on-drug against off-drug contrast for the analysis to detect.](figures/55-fig1-coupling-by-visit.png)

## 7. Worked example 2: power

### 7.1 Closed form

Table 6 gives closed-form power for the Hybrid design at $N = 70$ and
$c_{bm} = 0.25$, computed from the code's own mean vector and covariance
matrix for each path and pooled over the four paths with the published
allocation (18, 18, 17, 17), following `docs/37`, Section 4.4
(verified). It reproduces the `docs/37` values (its Section 5) to the
three decimals printed, from an independent
implementation that reads the matrices from Hendrickson's code rather
than from a reconstruction.

| Rule | Statistic | $t_{1/2}$ = 0 | 0.1 | 0.2 | 0.5 | 1 | 2 | 4 |
|------|------|------|------|------|------|------|------|------|
| Step (58b32a9) | E9 | 0.721 | 0.050 | 0.050 | 0.050 | 0.050 | 0.050 | 0.050 |
| | E9b | 0.654 | 0.081 | 0.081 | 0.081 | 0.081 | 0.081 | 0.081 |
| Graded (intended) | E9 | 0.721 | 0.721 | 0.712 | 0.631 | 0.473 | 0.278 | 0.140 |
| | E9b | 0.654 | 0.653 | 0.647 | 0.586 | 0.473 | 0.328 | 0.210 |
| 2024 upstream (ratio) | E9 | 0.495 | 0.495 | 0.495 | 0.475 | 0.399 | 0.260 | 0.127 |
| | E9b | 0.438 | 0.438 | 0.438 | 0.424 | 0.372 | 0.273 | 0.166 |

Table: Closed-form power of the paired-difference statistic by coupling rule and $t_{1/2}$, Hybrid design, $N = 70$, $c_{bm} = 0.25$, compound symmetry (verified)

Two features stand out. Under the step rule power is flat at every
positive half-life, at the test size for E9 and at the baseline
remnant for E9b; the magnitude of the half-life plays no role. Under
graded coupling power is unchanged at the article's half-lives (0.721
at 0.1 week, 0.712 at 0.2) and declines only when carryover is long
compared with the one-week spacing of the off-drug visits. The coupling
gap, patient-weighted over paths, is 0.9997 at 0.1 week, 0.991 at 0.2,
0.906 at 0.5, 0.754 at 1 and 0.551 at 2 under graded coupling, and 0
at every positive half-life under step coupling (verified; the graded
values agree with `docs/37`, Section 5.4).

### 7.2 The boundary: $t_{1/2} = 0$ and very small half-lives

At $t_{1/2} = 0$ the two rules agree exactly: there is no carryover,
the off-drug means are exactly 0 (Section 4.3), and the graded rule's
last branch also gives 0 (verified; identical matrices). The step
rule's collapse is therefore a discontinuity at $t_{1/2} = 0$: the
limit of its power as $t_{1/2} \to 0^+$ is 0.081 (E9b), not the 0.654
it takes at 0 (derived).

At very small half-lives the floating-point arithmetic intervenes. The
smallest positive double-precision number is about
$4.9 \times 10^{-324}$, or $2^{-1074}$ [5]; a decayed mean smaller than
that is stored as exactly 0, and the step test then reads it as "no
carryover". The published mean at the $m$-th off-drug visit of a run is
$\mu_{\text{last}}\,2^{-(t_{sd,1} + \cdots + t_{sd,m})/t_{1/2}}$
(`docs/53`, Section 5), so it underflows when the summed exponent
exceeds about 1077 (derived). The half-life below which each visit's
mean is exactly 0 was found by bisection on the restated recursion and
confirmed on the code itself just above and below each threshold
(verified): in the Hybrid design 0.00093 weeks for the first off-drug
visit of every run with $t_{sd} = 1$, up to 0.0121 weeks for path 4's
fourth off-drug visit; in OL+BDC up to 0.0056; in CO up to 0.0232. Above
0.0121 weeks no Hybrid visit underflows, and the step rule's power is
exactly its plateau. Table 7 shows the transition.

| $t_{1/2}$ (weeks) | Step coupling at the off-drug visits, paths 1 / 2 / 3 / 4 | Gap, step | E9b, step | E9b, graded |
|------|------|------|------|------|
| 0 | 0 0 0 / 0 0 0 / 0 0 0 0 / 0 0 0 0 | 1 | 0.654 | 0.654 |
| 0.0005 | 0 0 0 / 0 0 0 / 0 0 0 0 / 0 0 0 0 | 1 | 0.654 | 0.654 |
| 0.001 | 1 0 0 / 1 0 0 / 1 0 0 0 / 1 0 0 0 | 0.707 | 0.437 | 0.654 |
| 0.005 | 1 1 1 / 1 1 0 / 1 1 0 1 / 1 1 0 0 | 0.268 | 0.170 | 0.654 |
| 0.01 | 1 1 1 / 1 1 1 / 1 1 1 1 / 1 1 1 0 | 0.061 | 0.096 | 0.654 |
| 0.0125 and above | 1 1 1 / 1 1 1 / 1 1 1 1 / 1 1 1 1 | 0 | 0.081 | 0.654 at 0.0125 |

Table: The boundary of the step rule: step coupling at each path's off-drug visits (1 coupled, 0 not), patient-weighted coupling gap and closed-form E9b power at very small half-lives, Hybrid, $N = 70$, $c_{bm} = 0.25$ (verified)

At 0.0005 weeks (about 5 minutes) every off-drug mean underflows, and
the step rule happens to reproduce the no-carryover result; at 0.001
weeks only the first visit of each run survives; by 0.0125 weeks
nothing underflows and the collapse is complete. Graded coupling is
unaffected throughout: its off-drug coupling is at most
$2^{-1/0.01} = 7.9 \times 10^{-31}$ at 0.01 week, and its own underflow,
when $2^{-t_{sd}/t_{1/2}}$ falls below $2^{-1074}$, replaces a
negligible number by zero (verified). The step rule thus has three
regimes in $t_{1/2}$, set by the visit schedule and the limits of
double precision rather than by the model: no carryover (exactly 0, or
small enough to underflow), a ragged transition below about 0.012
weeks, and the collapse above it. Figure 1 of `docs/37` starts at 0.02
weeks for this reason (`docs/37`, Section 5.1; quoted in substance).

The thresholds depend on the recursion. With the anchored means of
`docs/53` the mean at a visit underflows when $t_{sd}/t_{1/2}$ alone
exceeds about 1077, so path 4's fourth off-drug visit ($t_{sd} = 7$)
would underflow below 0.0065 weeks rather than 0.0121 (derived). This
is the only way in which the recursion touches the step rule, and it
lies entirely below the half-lives of any practical interest.

### 7.3 Hendrickson's own code

The closed form is for E9 and E9b, not for the published mixed model.
Three simulations ran the published `generateData()` and the published
`lme_analysis()` (the `lmer` model of `lme4` [4], with `useDE = FALSE`)
on the Hybrid
design at $N = 70$, $c_{bm} = 0.25$, without censoring, 500 replicates
per cell, with the same seed for each version at a half-life (Table 8).

| $t_{1/2}$ (weeks) | Published, run A | Graded, run A | Published, run B | Graded, run B | Published, run C | Graded, run C | 2024 upstream, run C |
|------|------|------|------|------|------|------|------|
| 0 | 0.668 | 0.668 | 0.628 | 0.628 | 0.670 | 0.670 | 0.456 |
| 0.001 | | | | | 0.464 | 0.660 | 0.482 |
| 0.01 | | | | | 0.098 | 0.676 | 0.458 |
| 0.1 | 0.126 | 0.638 | 0.078 | 0.658 | 0.100 | 0.656 | 0.464 |
| 0.2 | 0.092 | 0.656 | 0.082 | 0.670 | | | |
| 1 | 0.094 | 0.462 | 0.088 | 0.486 | 0.076 | 0.448 | 0.368 |
| 2 | | | 0.078 | 0.344 | | | |

Table: Power of the published analysis for the interaction on data from Hendrickson's generateData() as published (step), with graded coupling only, and in the 2024 upstream version; Hybrid, $N = 70$, $c_{bm} = 0.25$, no censoring, 500 replicates per cell. Run A: `14-graded-patch-58b32a9.R` (verified, stored output). Run B: `docs/53`, Section 8.3 (verified, stored output). Run C: `55-power-check.R` (verified). Monte Carlo standard error at most 0.023

Four things follow (verified unless stated).

- **At $t_{1/2} = 0$ the published and graded versions coincide** in
  every run, as they must: identical matrices and a shared seed give
  identical data and identical decisions.
- **At any practical half-life the published code collapses and the
  patch does not.** At 0.1 week published power is 0.078 to 0.126 in
  the three runs, against 0.638 to 0.658 with graded coupling; at
  1 week, 0.076 to 0.094 against 0.448 to 0.486. The spread between
  runs is within about two combined Monte Carlo standard errors except
  run A at 0.1 week, which `docs/37` (Appendix A) discusses.
- **The boundary of Section 7.2 appears in the published code
  itself.** At 0.001 week, where only the first off-drug visit of each
  run escapes underflow, published power is 0.464, between the
  no-carryover value and the plateau, as the closed form predicts
  (E9b 0.437); at 0.01 week it is already 0.098 (closed form 0.096).
- **The 2024 data generator does not collapse, but loses power at
  every half-life.** Its power is 0.456 to 0.482 at half-lives up to
  0.1 week, against 0.670 for the published and graded versions
  without carryover, because visit 1 is never coupled (Section 4.5);
  at 1 week it is 0.368 against 0.448 for graded coupling. The mean
  interaction estimate is $-0.102$ at $t_{1/2} = 0$, against $-0.131$
  for the other two versions, matching the closed-form $-0.101$.

No matrix was repaired in any cell of run C.

![Power for the biomarker-by-treatment interaction against the carryover half-life on a logarithmic axis (the point at the left edge is $t_{1/2} = 0$). Lines: closed form for E9b, the proxy for the published analysis (left), and for E9 (right), under the published step rule, graded coupling and the 2024 upstream ratio rule. Open circles (left panel): power of the published `lmer` analysis on data from Hendrickson's own `generateData()` in each version, runs A and C of Table 8. Hybrid design, $N = 70$, $c_{bm} = 0.25$. The step rule's power is a switch between its value at 0 and its plateau, with a ragged transition below about 0.012 weeks caused by floating-point underflow; graded coupling declines smoothly.](figures/55-fig2-power-vs-halflife.png)

## 8. Consequences

### 8.1 What changes: interaction power under carryover

Every carryover scenario of the article's Figure 4B uses $t_{1/2} = 0.1$
or 0.2 week (the `coremodelparams` grid, vignette line 318; inspected).
In each of them the step rule leaves no on-drug against off-drug
contrast in the Hybrid and OL+BDC designs, and half of it in CO, where
the placebo-first path's off-drug visits come before any exposure and
keep $g_t = 0$ (`docs/37`, Section 5.7; quoted). This is the pattern the
article reports: "a precipitous decline in power in both the N-of-1
and, to a lesser but still very significant degree, the open label +
blinded discontinuation design ... A decrease in power in the cross
over design was also seen, but this was significantly less severe"
([1], p. 7; quoted). Rerun on the published code at $c_{bm} = 0.25$ for
all four designs, both sample sizes and every censoring pattern,
graded coupling removes the collapse in each design (report 01,
Results, "The correction in the published code", and its Figure 4B
replication; quoted).

At the article's own $c_{bm} = 0.6$ there is a further asymmetry. In
the Hybrid design the correlation matrix is not positive definite at
$t_{1/2} = 0$ on any path, and it is positive definite at
$t_{1/2} = 0.1$ on every path (verified). The code repairs the matrix
silently when it is not (lines 145 to 150), so the no-carryover column
of Figure 4B was simulated from repaired matrices, whose on-drug
correlation is 0.54 to 0.57 rather than 0.6 (`docs/36`, Section 4.3;
quoted), and the carryover columns from unrepaired matrices with full
coupling at every visit (derived from Section 8.2). The comparison
across those columns is therefore a comparison between two different
models, neither of which is the stated one.

### 8.2 What changes: the positive-definiteness ceiling

A single biomarker can be correlated with a set of partly independent
visits only up to a limit, the ceiling $c_{bm}^{\ast} =
(\tilde g^\top M^{-1}\tilde g)^{-1/2}$, with $\tilde g$ the coupling
pattern padded with zeros at the $TV$ and $PB$ entries and $M$ the
correlation matrix of the response components (`docs/36`, Section 4;
derived there). The coupling rule changes $g$, so it changes the
ceiling (Table 9; verified, with the code's own matrices).

| Design | Rule | $t_{1/2}$ = 0 | 0.1 | 0.2 | 0.5 | 1 | 2 |
|------|------|------|------|------|------|------|------|
| Hybrid | Step | 0.256 | 0.893 | 0.893 | 0.893 | 0.893 | 0.893 |
| | Graded | 0.256 | 0.256 | 0.257 | 0.274 | 0.307 | 0.369 |
| | 2024 upstream | 0.256 | 0.256 | 0.256 | 0.260 | 0.273 | 0.295 |
| OL+BDC | Step | 0.262 | 0.893 | 0.893 | 0.893 | 0.893 | 0.893 |
| | Graded | 0.262 | 0.262 | 0.265 | 0.289 | 0.348 | 0.469 |
| | 2024 upstream | 0.256 | 0.256 | 0.256 | 0.260 | 0.274 | 0.307 |
| CO | All three | 0.256 | 0.256 | 0.256 | 0.256 | 0.256 | 0.256 |

Table: Positive-definiteness ceiling of $c_{bm}$ (minimum over paths) by design, coupling rule and $t_{1/2}$, compound symmetry, $\rho_{cs} = 0.8$ (verified)

Under the step rule the ceiling jumps from about 0.26 to 0.893 at any
positive half-life, because a coupling that is constant over visits
aligns with the large person-level component of the compound-symmetric
structure and is cheap to accommodate (`docs/36`, Section 4.2;
`analysis/report/compendium-mathematics.tex`, Table "Design-level
ceilings by configuration"; both quoted, and reproduced here). The
cells that stop needing repair are exactly the cells that have lost
their signal. Under graded coupling the ceiling rises only slowly with
the half-life; under compound symmetry it stays below 0.3 up to
$t_{1/2} = 0.5$, so the article's 0.3 and 0.6 remain infeasible there
under any coupling rule. That is the separate positive-definiteness
problem of `docs/36` (problem 1), which graded coupling does not solve
and is not meant to. In CO the binding path is the placebo-first one,
whose pre-exposure visits have zero coupling under every rule, so its
ceiling does not move (`docs/34`, Section 4.6; quoted). The 2024 rule's
ceilings at $\rho_{cs} = 0.7$, 0.346, 0.351 and 0.368 in Hybrid at
$t_{1/2} = 0$, 0.5 and 1, reproduce those of `docs/54` and the
compendium exactly (verified).

### 8.3 What does not change

- **The means and the response block.** The coupling rule writes only
  the biomarker row and column of the correlation matrix (lines 135
  and 136). The mean vector, the variances and every correlation among
  the baseline and the three response components are identical under
  the step and graded rules (inspected; verified, largest difference 0
  in all 40 cases). Below the ceiling no repair intervenes, so the
  joint distribution of the outcomes $Y_{it}$, marginal over the
  biomarker, is the same under both rules (derived).
- **The population-average drug effect.** It follows that any summary
  of the outcomes that does not use the biomarker has the same
  distribution under both rules: the expected on-minus-off contrast
  $\mathbb E(\Delta)$ is identical in every path (verified), and so is the drug
  effect averaged over participants (derived). In the published mixed
  model this is the `Db` effect evaluated at the mean biomarker,
  $\hat\beta_D + \hat\beta_{bm:D}\bar B$; its mean over replicates was
  $-8.24$ under all three versions at $t_{1/2} = 0$ and $-8.16$ and
  $-8.14$ under the graded and published versions at 0.1 week, a
  difference of 0.02 against a Monte Carlo standard error of about
  0.04 (Table 10; verified). The interaction itself carries none of
  the main drug effect: the coupling moves where the biomarker places
  each participant around the mean response, not the mean response.
- **One qualification.** The published model enters the biomarker
  uncentered (`Sx~bm+Db+t+bm*Db`; lines 117 and 137 of
  `R/lme_analysis.R`, inspected), so its `Db` coefficient is the drug
  effect at a biomarker value of 0, far outside the data (mean 124.33).
  That coefficient equals the average effect minus
  $\beta_{bm:D}\,\mu_{bm}$, and so it moves with the interaction: at
  0.1 week the mean interaction estimates differ by
  $-0.128 - (-0.031) = -0.097$, so the `Db` coefficient should differ
  by about $0.097 \times 124.3 = 12$ points (derived); the simulated
  means are $+7.73$ (graded) and $-4.26$ (published), a difference of
  12.0 (verified, Table 10). Results that report the `Db` coefficient
  itself, rather than the effect at the biomarker mean, would be
  affected.

| $t_{1/2}$ (weeks) | Version | Mean `bm:Db` estimate | Mean `Db` coefficient (at $B = 0$) | Mean drug effect at mean biomarker | Its SD over replicates |
|------|------|------|------|------|------|
| 0 | Published | $-0.131$ | 8.10 | $-8.23$ | 0.86 |
| 0 | Graded | $-0.131$ | 8.10 | $-8.23$ | 0.86 |
| 0 | 2024 upstream | $-0.102$ | 4.48 | $-8.23$ | 0.88 |
| 0.1 | Published | $-0.031$ | $-4.26$ | $-8.14$ | 0.85 |
| 0.1 | Graded | $-0.128$ | 7.73 | $-8.16$ | 0.86 |
| 0.1 | 2024 upstream | $-0.099$ | 4.21 | $-8.15$ | 0.86 |
| 1 | Published | $-0.032$ | $-2.81$ | $-6.81$ | 0.89 |
| 1 | Graded | $-0.105$ | 6.30 | $-6.72$ | 0.91 |
| 1 | 2024 upstream | $-0.091$ | 3.71 | $-7.63$ | 0.92 |

Table: Drug coefficients of the published analysis by data-generator version, run C: the uncentered `Db` coefficient moves with the interaction, the drug effect at the mean biomarker does not (500 replicates per cell; verified)

The 2024 version's drug effect at 1 week ($-7.63$ against $-6.81$ and
$-6.72$) differs for a different reason: its mean recursion uses
$s_f = 2$ (`docs/53`, Section 4.4), which changes the off-drug means,
not the coupling (inferred from the code; not separated here).

### 8.4 Where the step rule could matter more

The following are inferred, not tested here.

- **Analyses with a decayed drug indicator.** An analysis that weights
  off-drug visits by $2^{-t_{sd}/t_{1/2}}$ (the continuous indicator
  $D_{bc,it}$ that the article's discussion anticipates, "analysis
  methods that incorporate an expectation of carryover", [1], p. 12)
  is built for graded coupling. Under step coupling the true
  interaction off drug is at full strength, so such an analysis would
  be misspecified in the opposite direction from the one it guards
  against.
- **Designs with long off-drug runs.** Under graded coupling a design
  gains contrast from off-drug visits late in a run, where the
  coupling has faded. Under step coupling those visits add no contrast
  at all, so design comparisons made with the published code
  understate the value of washout periods.
- **Estimating the half-life from the biomarker channel.** Under step
  coupling the correlation structure carries no information about
  $t_{1/2}$ except whether it is zero.

## 9. The fix

### 9.1 The corrected code

The replacement couples in proportion to the drug effect that remains,
in the style of the surrounding code (report 01, Appendix B.1; quoted):

```r
  # correlation with biomarker: full on drug, decaying with the drug
  # effect off drug (graded coupling), zero before first exposure
  for(p in 1:nP){
    n1<-paste(trialdesign$timeptname[p],"br",sep=".")
    if(d[p]$tod>0){
      g<-1
    }else if(d[p]$tsd>0 && modelparam$carryover_t1half>0){
      g<-(1/2)^(d[p]$tsd/modelparam$carryover_t1half)
    }else{
      g<-0
    }
    correlations[n1,'bm']<-modelparam$c.bm*g
    correlations['bm',n1]<-modelparam$c.bm*g
  }
```

As a minimal diff against `58b32a9`, which replaces lines 134 to 137
and nothing else:

```diff
--- a/R/generateData.R
+++ b/R/generateData.R
@@ -131,9 +131,14 @@
     # correlation with biomarker
     for(p in 1:nP){
       n1<-paste(trialdesign$timeptname[p],"br",sep=".")
-      if(means[which(n1==labels)]!=0){
-        correlations[n1,'bm']<-modelparam$c.bm
-        correlations['bm',n1]<-modelparam$c.bm
-      }
+      if(d[p]$tod>0){
+        g<-1
+      }else if(d[p]$tsd>0 && modelparam$carryover_t1half>0){
+        g<-(1/2)^(d[p]$tsd/modelparam$carryover_t1half)
+      }else{
+        g<-0
+      }
+      correlations[n1,'bm']<-modelparam$c.bm*g
+      correlations['bm',n1]<-modelparam$c.bm*g
     }
   }
```

Three details of the replacement are deliberate.

- **The coupling is computed from the trial design, not from the
  means.** `d` is the trial-design table of line 53, which already
  carries `tod` and `tsd` from `buildtrialdesign()` (`docs/53`, Section
  4.2). No new input is needed, and the rule is independent of how the
  means are computed.
- **The pre-exposure branch.** `tsd` is 0 both on drug and before first
  exposure (`buildtrialdesign.R`, lines 104 to 105). The `tod>0` test
  comes first, so on-drug visits get 1, and the remaining `tsd = 0`
  visits, the baseline and the placebo-first crossover visits, get 0.
  A rule written as "coupling equals $2^{-t_{sd}/t_{1/2}}$" would give
  them $2^0 = 1$, full coupling with no drug ever taken (report 01,
  Appendix B.4, describes the same defect in an analysis indicator).
- **The $t_{1/2} = 0$ branch.** With no carryover, off-drug visits get
  0, so the corrected code reproduces the published matrices exactly
  at $t_{1/2} = 0$ (verified, Section 7.2).

### 9.2 The package's version

pmsimstats-ng (version 0.3.0) implements the same rule in `buildSigma()`,
`R/generateData.R`, lines 341 to 355 (inspected):

```r
    if(cc=="br" && dgp_architecture %in% c("mvn", "combined")){
      c_bm_cov<-if(dgp_architecture == "combined") modelparam$c.bm_b else modelparam$c.bm
      for(p in 1:nP){
        n1<-paste(trialdesign$timeptnames[p],"br",sep=".")
        if(d[p]$onDrug){
          correlations[n1,'bm']<-c_bm_cov
          correlations['bm',n1]<-c_bm_cov
        } else if(d[p]$tsd>0 && lambda_cor>0){
          decay<-exp(-lambda_cor * d$tsd[p])
          correlations[n1,'bm']<-c_bm_cov * decay
          correlations['bm',n1]<-c_bm_cov * decay
        }
      }
    }
```

with `lambda_cor` $= \ln 2 / t_{1/2}$ by default (lines 224 to 229), so
`exp(-lambda_cor * tsd)` $= 2^{-t_{sd}/t_{1/2}}$, and pre-exposure
visits left at the identity's 0. The rest of that package's process
differs from `58b32a9` (it uses AR(1) correlation, among other
changes), so it is not a clean test of the coupling change alone; the
patched `58b32a9` file is.

### 9.3 Check

The patched `58b32a9` file was run on every path of the Hybrid, OL+BDC
and CO designs (eight paths) at $t_{1/2} \in \{0.1, 0.2, 0.5, 1, 2\}$,
40 cases, and the coupling read from its correlation matrix was
compared with the intended rule of Section 3 computed directly. The
largest absolute difference was $1.1 \times 10^{-16}$ (verified). In the
same 40 cases the package's `buildSigma()` matched the intended rule to
$1.1 \times 10^{-16}$ (verified), the published file's coupling was
exactly the indicator of a nonzero mean (verified), and the published
off-drug coupling after exposure was 1 in every case (verified). The
mean vectors and response blocks of the published and patched files
were identical (Section 8.3). The 2024 upstream file's coupling matched
graded coupling at half-life $t_{1/2}/2$ with $g_1 = 0$, the reading of
Section 4.5, to $1.1 \times 10^{-16}$ in the same 40 cases (verified).

On the published analysis the patch restores power at the article's
half-lives to within Monte Carlo error of the no-carryover value
(Table 8): 0.638, 0.658 and 0.656 at 0.1 week against 0.668, 0.628 and
0.670 without carryover, in three independently seeded runs.

### 9.4 Why the coupling should not be computed from the means

The step rule and the 2024 ratio rule both read the coupling off the
drug-response means. That couples two parts of the model that should
be specified separately. The step rule turns a small positive mean into
full coupling. The ratio rule works only while the means compound
(Section 10), and needs a previous visit, which is why it never
couples visit 1. Computing $g_t$ from the drug state and $t_{sd}$, as in
Section 9.1, keeps the coupling correct whatever the mean model, and
leaves the decay family to be varied in one place if desired, as the
article anticipates ("the simulations could be easily adapted to
incorporate an alternative model", [1], p. 12; quoted).

## 10. How the two defects interact

The `58b32a9` file has two defects in the carryover: the compounding
recursion of the mean (`docs/53`) and the step coupling (this note).
They touch disjoint lines (83 to 92 and 134 to 137), and their effects
separate cleanly in the published process.

- **The step coupling drives the power collapse.** In the paired
  simulation of `docs/53`, Section 8.3, correcting the mean alone left
  power at 0.078 to 0.088 at every positive half-life, and changed no
  test decision at 0.1 or 0.2 week; replacing the coupling raised it to
  0.658 at 0.1 week, with or without the mean correction (quoted).
- **The recursion matters at longer half-lives,** and only for the
  means: under graded coupling, adding the mean correction changed
  power by at most 0.002 at 1 and 2 weeks (0.486 against 0.484 and
  0.344 against 0.344; quoted, same source).
- **The fixes are independent.** The step rule's trigger depends only on
  whether the means are nonzero, which is the same under either mean
  (above the underflow thresholds of Section 7.2), and the graded rule
  does not read the means at all. The two patches can be applied in
  either order and give byte-identical files (verified), and the
  `docs/53` 2x2 simulation shows the corresponding power results
  (quoted above).

The 2024 upstream code is different, because there the coupling is
computed from the means. Its ratio rule
$g_t = \mu_{BR,t}/\mu_{BR,t-1}$ reproduces a fading coupling only
because the means compound; with the anchored means of `docs/53` the
ratio becomes $2^{-s_f(t_{sd,t} - t_{sd,t-1})/t_{1/2}}$, a constant
factor per visit (derived; `docs/53`, Section 9.5, gives the derivation
and its own check). Table 11 confirms it on the upstream file itself,
with the anchored mean substituted and nothing else changed.

| Design, path | Off-drug visit ($t_{sd}$) | Upstream as is ($s_f = 2$) | Upstream, anchored means ($s_f = 2$) | Intended, stated $t_{1/2}$ |
|------|------|------|------|------|
| Hybrid, 1 | BD3 (1) | 0.250 | 0.250 | 0.500 |
| | BD4 (2) | 0.0625 | 0.250 | 0.250 |
| | COp (4) | 0.0039 | 0.0039 | 0.0625 |
| CO, 1 | COb1 (2.5) | 0.0313 | 0.0313 | 0.177 |
| | COb2 (5) | 0.00098 | 0.0313 | 0.0313 |
| | COb3 (7.5) | $3.1 \times 10^{-5}$ | 0.0313 | 0.0055 |
| | COb4 (10) | $9.5 \times 10^{-7}$ | 0.0313 | 0.00098 |

Table: Off-drug coupling $g_t$ of the 2024 upstream code at $t_{1/2} = 1$ week, as is and with only the means corrected to the anchored form, against the intended graded rule (verified)

With the means corrected and the ratio rule kept, the crossover's
coupling would sit at 0.031 for all four off-drug visits, ten weeks
after stopping as one week after, while the intended coupling falls to
0.001. The same holds at $s_f = 1$, where the plateau is 0.177
(verified; `docs/53`, Table 8). So for the current upstream code the
two fixes must go together: correcting the recursion requires
replacing the ratio rule by the direct rule of Section 9.1, which also
restores the coupling at visit 1. With the direct rule in place, the
two fixes are again independent.

## 11. How to verify independently

The snippet below needs only `data.table`, `MASS` and `corpcor`, the
packages the original code already uses, and a checkout of
`rchendrickson/pmsimstats` at `58b32a9` as the working directory. It
does not use pmsimstats-ng. It runs the original `generateData()` with
`empirical = TRUE`, for which `MASS::mvrnorm()` returns a sample whose
means and covariances equal the supplied ones exactly [3], so the
sample correlations are the code's own coupling and the regression
slope is the population slope.

```r
## Run in a checkout of github.com/rchendrickson/pmsimstats at 58b32a9.
library(data.table); library(MASS); library(corpcor)
for (f in c('utilities.R', 'buildtrialdesign.R', 'generateData.R'))
  source(file.path('R', f))
load('data/extracted_rp.rda'); load('data/extracted_bp.rda')
td <- buildtrialdesign(
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1','OL2','BD1','BD2','BD3','BD4','COd','COp'),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0)))$trialpaths[[1]]
t_half <- 1; c_bm <- 0.25
mp <- data.frame(N = 200, c.bm = c_bm, carryover_t1half = t_half,
  c.tv = .8, c.pb = .8, c.br = .8, c.cf1t = .2, c.cfct = .1)
## empirical = TRUE: sample means and covariances equal the code's own
dat <- generateData(mp, extracted_rp, extracted_bp, td,
  empirical = TRUE, makePositiveDefinite = TRUE)
br <- paste0(td$timeptnames, '.br')
code_g <- sapply(br, function(v) cor(dat$bm, dat[[v]])) / c_bm
on <- td$tod > 0
graded_g <- ifelse(on, 1, ifelse(td$tsd > 0, 2^(-td$tsd / t_half), 0))
print(data.frame(visit = td$timeptnames, on = on, tsd = td$tsd,
  code = round(code_g, 4), graded = round(graded_g, 4)), row.names = FALSE)
## E9: on-drug mean minus off-drug mean of the outcome, regressed on bm
y <- as.matrix(dat[, td$timeptnames, with = FALSE])
delta <- rowMeans(y[, on]) - rowMeans(y[, !on])
k <- c_bm * extracted_rp[cat == 'br']$sd / extracted_bp[cat == 'bm']$sd
cat(sprintf('E9 slope, code: %.4f; graded rule: %.4f\n',
  coef(lm(delta ~ dat$bm))[2], -k * (1 - mean(graded_g[!on]))))
```

Its output, run against the `58b32a9` files with the `extracted_rp.rda`
and `extracted_bp.rda` fetched from GitHub at that commit (verified):

```
 visit    on tsd code graded
   OL1  TRUE   0    1 1.0000
   OL2  TRUE   0    1 1.0000
   BD1  TRUE   0    1 1.0000
   BD2  TRUE   0    1 1.0000
   BD3 FALSE   1    1 0.5000
   BD4 FALSE   2    1 0.2500
   COd  TRUE   0    1 1.0000
   COp FALSE   4    1 0.0625
E9 slope, code: 0.0000; graded rule: -0.0949
```

The `code` column is the step rule: full coupling at the three off-drug
visits. Changing `t_half` to 0.1 leaves it unchanged; setting it to 0
makes the off-drug entries 0 and the code's slope $-0.1302$. Applying
the change of Section 9.1 before sourcing makes the `code` column equal
the `graded` column and the code's slope equal $-0.0949$.

## 12. Limitations

- **One defect, one component.** The note treats the off-drug coupling
  of the biomarker to the drug response. The carryover of the mean is
  `docs/53`; the compound-symmetry ceiling, which prevents the article's
  $c_{bm} = 0.3$ and 0.6 from being simulated as stated under either
  coupling rule, is `docs/36`.
- **The intended rule is our reading.** The article does not state an
  off-drug coupling (Section 3). Graded coupling is the reading that
  matches its description of the carryover and the pharmacology of a
  predictive biomarker; the authors may hold a different intent, which
  they are best placed to state.
- **Closed form and simulation scope.** The closed form is exact for E9
  and E9b, which reproduce the published analysis's mean estimates and
  its Hybrid power (`docs/37`, Section 6) but are not that analysis. The
  simulations cover one design (Hybrid), $N = 70$, $c_{bm} = 0.25$, no
  censoring and 500 replicates per cell; report 01 extends the
  published-code comparison to all four designs, both sample sizes and
  the censoring patterns at 250 replicates.
- **The 2024 code** was examined for its data-generating process only,
  with the published analysis. Its revised `lme_analysis()` (commit
  `325314f` and later) was not run here.
- **Section 8.4** is inferred, not tested.

## 13. Reproducibility

Every number labeled verified in Sections 1 and 4 to 11 was computed
by two scripts written for this note, in
`analysis/scripts/step-coupling/` and run from the repository root:

- `55-coupling-examples.R`: Tables 2 to 7, 9 and 11, both figures,
  the 40-case checks of the graded patch, of the package's
  `buildSigma()` and of the 2024 rule, the underflow thresholds, the
  ceilings, the marginal-distribution check and the order-of-patches
  check. It runs each version of `generateData()` with the sampler
  replaced by a recorder of the mean vector and covariance matrix, so
  that every value is the code's own. Run time about 7 minutes (421 s elapsed, one
  process).
- `55-power-check.R`: run C of Table 8 and Table 10 (the drug
  effect), the published `lme_analysis()` on data from the published,
  graded and 2024 `generateData()`, 500 replicates per cell. Run time
  11.2 minutes on 8 cores.

For example,

```bash
Rscript analysis/scripts/step-coupling/55-power-check.R --reps 500 --cores 8
Rscript analysis/scripts/step-coupling/55-coupling-examples.R out55
```

The first writes to `analysis/data/quick-sim/step-coupling/`; the
second writes the tables to `out55/`, the figures to `docs/figures/`,
and reads the stored output of the first and of
`analysis/scripts/quick-sim/hendrickson-problems/14-graded-patch-58b32a9.R`
(run A, `analysis/data/quick-sim/hendrickson-problems/graded-patch/graded-patch-reps500.csv`)
for Figure 2. Run B is the paired simulation of `docs/53`
(`analysis/scripts/carryover-recursion/53-power-check.R`); its output
is stored in `analysis/data/quick-sim/carryover-recursion/`
(`power-check-reps500.csv`, with `power-check-reps500-pairs.csv` and
`power-check-reps500-raw.csv`), and the run B columns of Table 8 match
`power-check-reps500.csv` (verified). The scripts read the
`58b32a9` sources from `analysis/data/quick-sim/hendrickson-58b32a9/R/`,
which were checked to be byte-identical to the file fetched from GitHub
at `58b32a9` (`generateData.R`, with `cmp`), and the 2024 file from
`analysis/scripts/quick-sim/hendrickson-original-comparison/vendored-hendrickson-generateData.R`,
which differs from upstream `06dac83` in one semantically identical
line (a modern `data.table` column extraction; checked with `diff`).
The parameter files `data/extracted_rp.rda` and `data/extracted_bp.rda`
in this repository are identical to those at `58b32a9` (`docs/53`,
Section 12).

## 14. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker
   validation: statistical methods and theoretical findings.
   *Frontiers in Digital Health*. 2020;2:13.
   doi:10.3389/fdgth.2020.00013.
2. Hendrickson RC. pmsimstats: R package for simulation of aggregated
   N-of-1 trials. GitHub repository rchendrickson/pmsimstats, commit
   58b32a92b99a8beb145de00c7bd8a0880f06f0f0 (2020-11-19); commit
   8609f124eef1770763ecbb4a8978a474c382f615 (2024-05-06); current
   `master` 06dac83f31ca85c994a6a2c1a0b98bc18ef233bd (2024-11-27).
   https://github.com/rchendrickson/pmsimstats.
3. Venables WN, Ripley BD. *Modern Applied Statistics with S*. 4th ed.
   New York: Springer; 2002. doi:10.1007/978-0-387-21706-2.
4. Bates D, Mächler M, Bolker B, Walker S. Fitting linear mixed-effects
   models using lme4. *Journal of Statistical Software*.
   2015;67(1):1-48. doi:10.18637/jss.v067.i01.
5. Goldberg D. What every computer scientist should know about
   floating-point arithmetic. *ACM Computing Surveys*.
   1991;23(1):5-48. doi:10.1145/103162.103163.
6. Senn S. *Cross-over Trials in Clinical Research*. 2nd ed.
   Chichester: Wiley; 2002. doi:10.1002/0470854596.

---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the numbers agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at docs/series/registry/open-items.csv -->

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
