---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# The Carryover Recursion in the Published pmsimstats Code: Compounding Decay of the Off-Drug Drug-Response Mean {.unlisted .unnumbered}
*2026-10-08 15:21 PDT*

**Author.** pmsimstats team

**Purpose.** This note describes one specific feature of the code that
accompanies Hendrickson et al. (2020) [1], the `pmsimstats` package at
commit `58b32a9` [2]. When a simulated participant stops the drug, the
code is meant to let the expected drug response decay exponentially
from its last on-drug value with a stated half-life, as the article
describes. The code instead decays the previous visit's mean, which
has already been decayed, by the cumulative time since
discontinuation. The decay therefore compounds from the second
off-drug visit of a run onward and runs faster than the stated
half-life. The note quotes the code verbatim, derives the size of the
discrepancy, works two examples with every number shown, states what
the discrepancy does and does not change downstream, and gives a
one-line correction with a check that it reproduces the intended
equation exactly. It is written to be read on its own, by readers who
have not seen any other pmsimstats-ng document, and in particular by
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
(`analysis/report/NOTATION.md`). Symbols marked *local* are defined
here for this note only. (Notation, 2026-10-09: $G(t_{od})$,
$\mu_{BR,t}$, $\mu_{\text{last}}$ and $s_f$ are now in `NOTATION.md`
and are marked canonical. The local run position $m$ is unrelated to
the Gompertz ceiling $m$ of `NOTATION.md`, which this note does not
use.)

| Symbol | Meaning | Status |
|------------|--------------------------------------------------|----------|
| $t$ | visit index within a randomization path | canonical |
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug | canonical |
| $t_{od}$ | cumulative time on drug at a visit (code `tod`); the input to the Gompertz curve | canonical |
| $t_{sd}$ | time since discontinuation (code `tsd`), cumulative within an off-drug run, zero on drug and before first exposure | canonical |
| $t_{1/2}$ | carryover half-life, in weeks (code `carryover_t1half`) | canonical |
| $\lambda$ | decay rate, $\ln 2 / t_{1/2}$ | canonical |
| $BR$ | the biological (drug) response component | canonical |
| $G(t_{od})$ | modified Gompertz mean of $BR$ on drug (code `modgompertz`) | canonical |
| $\mu_{BR,t}$ | expected $BR$ at visit $t$ (code `brmeans[p]`) | canonical |
| $\mu_{\text{last}}$ | $\mu_{BR}$ at the last on-drug visit before an off-drug run | canonical |
| $m$ | position of an off-drug visit within its run, $m = 1, \ldots, L$ | local |
| $L$ | number of off-drug visits in a run | local |
| $j$ | running index over the off-drug visits of a run, $j \le m$ | local |
| $t_{sd,m}$ | cumulative $t_{sd}$ at the $m$-th off-drug visit of a run; $t_{sd,0} = 0$ | local |
| $\mu_m$ | $\mu_{BR}$ at the $m$-th off-drug visit of a run; $\mu_0 = \mu_{\text{last}}$ | local |
| $r_m$ | ratio of the published to the intended mean at the $m$-th off-drug visit | local |
| $t_{1/2}^{\text{eff}}(m)$ | effective half-life of the published recursion at the $m$-th off-drug visit | local |
| $\Delta t$ | visit spacing within an off-drug run, when equal, in weeks | local |
| $s_f$ | the scale factor added to the exponent in 2024 (code `scalefactor`) | canonical |
| $c_{bm}$ | on-drug correlation of the biomarker with $BR$ (code `c.bm`) | canonical |
| $\hat\beta_{bm:D}$ | estimated biomarker-by-treatment interaction | canonical |
| $N$, $\alpha$ | total sample size; nominal significance level | canonical |
| $D_{bc,it}$ | continuous exposure-decayed drug indicator (code `Dbc`) | canonical |
| $\mu_p$, $\bar\mu$ | path mean and its average, in the between-path term of `docs/37` | local (`docs/37`) |

Table: Notation index for this note

**Glossary.**

- **Hybrid design.** The design that the article calls the "N-of-1
  design": eight visits at weeks 4, 8, 9, 10, 11, 12, 16 and 20, with
  open-label titration, a blinded discontinuation and a brief
  crossover, and four randomization paths (`pathA` to `pathD`, here
  paths 1 to 4).
- **OL+BDC.** Open-label treatment followed by blinded discontinuation;
  two paths.
- **CO.** Two-period crossover, eight visits 2.5 weeks apart; two
  paths.
- **Off-drug run.** A maximal sequence of consecutive off-drug visits
  after the participant has been exposed to the drug.
- **Anchored decay.** Decay of the off-drug mean from $\mu_{\text{last}}$
  by the elapsed time $t_{sd}$; the model the article describes.
- **Compounding decay.** Decay applied to an already-decayed value
  with the cumulative elapsed time, so that elapsed time is counted
  more than once.
- **The published process.** The data-generating process of the
  `58b32a9` code (`NOTATION.md`); used here only where a name is
  needed.
- **The refined process.** The pmsimstats-ng data-generating process,
  which among other changes uses anchored decay (`NOTATION.md`).

**Evidence labels.** Each claim carries one of four labels: *quoted*
(text reproduced from a source), *inspected* (read from source code),
*derived* (follows algebraically) or *verified* (computed by running
code, with the script named in Section 12).

## 1. Summary for the busy reader

**What the code does.** In `R/generateData.R` at `58b32a9`, line 88
sets the expected drug response at an off-drug visit to the previous
visit's expected response times $2^{-t_{sd}/t_{1/2}}$, where $t_{sd}$
is the *cumulative* time since the drug was stopped (inspected).
Because the previous visit's value has already been decayed, elapsed
time is counted again at every visit, and at the $m$-th off-drug visit
the mean is $\mu_{\text{last}}\,2^{-(t_{sd,1} + \cdots + t_{sd,m})/t_{1/2}}$
(derived).

**What it should do.** The article states that "the value of the BR at
the last timepoint on active drug was exponentially decayed, with the
half-life of this decay being maintained as a model parameter" [1,
p. 6] (quoted). That is $\mu_{\text{last}}\,2^{-t_{sd,m}/t_{1/2}}$.

**Where the two agree and where they differ.** They agree at every
on-drug visit and at the first off-drug visit of every run. From the
second off-drug visit on, the code's mean is smaller by the factor
$2^{-(t_{sd,1} + \cdots + t_{sd,m-1})/t_{1/2}}$ (derived). With equally
spaced visits the effective half-life at the $m$-th visit is
$2t_{1/2}/(m+1)$: two-thirds of the stated value at the second visit,
one-half at the third (derived).

**How large.** On Hybrid path 1 at $t_{1/2} = 1$ week, the mean two
weeks after stopping is 1.273 points where the intended value is 2.547
(verified, Section 6). The largest shortfall in any off-drug mean over
the three designs' paths is 1.37 points at $t_{1/2} = 1$ week and 2.50
points at 2 weeks (verified, Table 5). At the half-lives used for the
article's power results, 0.1 and 0.2 weeks, the entire carryover is
nearly gone after one week under either form, and the largest
shortfall is 0.0104 points (verified). The discrepancy is therefore
numerically negligible for those results, and grows with the half-life
and with the length of the off-drug run.

**What it does not touch.** The variances and the correlation matrix
are unchanged, so the simulated data differ only by a shift in the
mean at the affected visits (inspected, verified). A paired simulation
with the published analysis shows no detectable change in power from
correcting the mean alone (verified, Section 8.3).

**The one-line fix.** Decay by the increment in $t_{sd}$ since the
previous visit rather than by its cumulative value:

```r
brmeans[p]<-brmeans[p]+brmeans[p-1]*
  (1/2)^((d$tsd[p]-d$tsd[p-1])/modelparam$carryover_t1half)
```

The product of the increments telescopes to the intended
$2^{-t_{sd,m}/t_{1/2}}$; the corrected code reproduces the intended
means to within $1.1 \times 10^{-14}$ in all 40 path and half-life
cases checked (verified, Section 9).

## 2. Scope

This note concerns one statement: the off-drug carryover of the mean
of the drug-response component. It does not address the biomarker
coupling rule, the correlation structure, positive-definiteness repair
or the analysis model, each of which is examined elsewhere in the
project. Where those features interact with the carryover of the mean,
the interaction is stated in Section 8.

## 3. The intended model

The article describes the drug-response component as follows (quoted;
[1], Methods, p. 6):

> The BR factor was set to zero at times when participants had never
> been on active drug; however, a carryover effect was built in such
> that when a participant moved from active drug to being off active
> drug, the value of the BR at the last timepoint on active drug was
> exponentially decayed, with the half-life of this decay being
> maintained as a model parameter.

and, in the Discussion (quoted; [1], p. 12):

> In these models, we assume an exponential decay analogous to a
> pharmacokinetic half-life, although the simulations could be easily
> adapted to incorporate an alternative model.

In words: when the drug is stopped, the expected drug response starts
at its last on-drug value and halves every $t_{1/2}$ weeks. Writing
$\mu_{\text{last}}$ for that value and $t_{sd,t}$ for the time since
discontinuation at visit $t$,

$$
\mu_{BR,t} \;=\; \mu_{\text{last}}\; 2^{-t_{sd,t}/t_{1/2}}
\;=\; \mu_{\text{last}}\; e^{-\lambda\, t_{sd,t}},
\qquad \lambda = \ln 2 / t_{1/2}.
$$

This is the anchored form. It is a single exponential in elapsed time,
so its half-life is $t_{1/2}$ at every visit, and it does not depend on
how many visits fall in the off-drug run. Carryover in within-subject
designs has a long methodological literature [4]; the exponential form
is a pharmacologically motivated choice within it, and it is the one
the article's text describes.

## 4. The code as published

### 4.1 The carryover lines

`R/generateData.R` at `58b32a9`, lines 83 to 92 (quoted;
[permalink](https://github.com/rchendrickson/pmsimstats/blob/58b32a92b99a8beb145de00c7bd8a0880f06f0f0/R/generateData.R#L83-L92)):

```r
 83      brmeans<-modgompertz(d$tod,rp$max,rp$disp,rp$rate)
 84      if(nP>1){
 85        for(p in 2:nP){
 86          if(!d[p]$onDrug){
 87            if(d[p]$tsd>0){
 88              brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^(d$tsd[p]/modelparam$carryover_t1half)
 89            }
 90          }
 91        }
 92      }
```

The drug state used on line 86 is set at line 55
([permalink](https://github.com/rchendrickson/pmsimstats/blob/58b32a92b99a8beb145de00c7bd8a0880f06f0f0/R/generateData.R#L53-L55)):

```r
 53  d<-data.table(trialdesign)
 54  d$t_wk_cumulative<-cumulative(d$t_wk)
 55  d[,onDrug:=(tod>0)]
```

The file `R/generateData.R` was last changed before `58b32a9` in the
initial commit `42ac030` (2020-02-22), so these lines are the ones
that produced the article's results (inspected, GitHub commit history
for the file).

### 4.2 How the trial design supplies $t_{od}$ and $t_{sd}$

`R/buildtrialdesign.R` at `58b32a9`, lines 86 to 106 (quoted;
[permalink](https://github.com/rchendrickson/pmsimstats/blob/58b32a92b99a8beb145de00c7bd8a0880f06f0f0/R/buildtrialdesign.R#L86-L106)):

```r
 86  t_wk<-c(timepoints[1],
 87          timepoints[2:length(timepoints)]-timepoints[1:(length(timepoints)-1)])
 88  for(iTD in 1:nTD){
 89    od<-ondrug[[iTD]]
 90    od_intervals<-t_wk*od
 91    tod<-od_intervals
 92    tsd<-t_wk*(od!=1)
 93    tpb<-t_wk*(expectancies>0)
 94    for(iTP in 2:length(timepoints)){
 95      if(od_intervals[iTP]>0){
 96        tod[iTP]<-tod[iTP]+tod[iTP-1]
 97      }else{
 98        tsd[iTP]<-tsd[iTP]+tsd[iTP-1]
 99      }
100      if(tpb[iTP]>0){
101        tpb[iTP]<-tpb[iTP]+tpb[iTP-1]
102      }
103    }
104    everondrug<-(cumulative(od)>0)
105    tsd<-everondrug*tsd
106    tds[[iTD]]<-data.table(timeptnames=timeptnames,t_wk=t_wk,e=expectancies,tod=tod,tsd=tsd,tpb=tpb)
```

The modified Gompertz curve, `R/utilities.R` lines 34 to 40 (quoted;
[permalink](https://github.com/rchendrickson/pmsimstats/blob/58b32a92b99a8beb145de00c7bd8a0880f06f0f0/R/utilities.R#L34-L40)):

```r
 34  modgompertz<-function(t,maxr,disp,rate){
 35    y<-maxr*exp(-disp*exp(-rate*t))
 36    vert_offset<-maxr*exp(-disp*exp(-rate*0))
 37    y<-y-vert_offset
 38    y<-y*(maxr/(maxr-vert_offset)) # return the assymptotic max to what it was
 39    return(y)
 40  }
```

### 4.3 Line-by-line reading

All statements in this subsection are inspected.

- **`buildtrialdesign`, lines 86 to 87.** `t_wk` holds the interval
  lengths: the first visit's week, then the gaps between visits.
- **Lines 90 to 92.** Each interval is assigned to time on drug if the
  participant is on drug at its closing visit (`od = 1`, "on drug
  entire previous interval and this measurement point", in the
  package's own comment), and to time off drug otherwise.
- **Lines 94 to 99.** Moving forward through the visits, an on-drug
  interval adds to the previous `tod`; an off-drug interval adds to
  the previous `tsd`. Two consequences follow. First, `tsd` is
  *cumulative within an off-drug run*: one week, then two, then three.
  Second, `tsd` returns to 0 at every on-drug visit, and `tod` is 0 at
  every off-drug visit, so time on drug restarts at zero on
  re-exposure.
- **Lines 104 to 105.** `tsd` is set to zero before first exposure, so
  a participant who starts on placebo has no carryover.
- **`generateData`, line 83.** `brmeans` is the Gompertz mean
  $G(t_{od})$ at every visit. Because $G(0) = 0$ (line 37 subtracts the
  value at zero), it is 0 at every off-drug visit.
- **Lines 85 to 87.** The loop runs forward from the second visit and
  acts only at off-drug visits after exposure ($t_{sd} > 0$).
- **Line 88.** The off-drug mean becomes $0 + \texttt{brmeans[p-1]}
  \times 2^{-t_{sd}/t_{1/2}}$. Since the loop runs forward,
  `brmeans[p-1]` has already been overwritten at the previous
  iteration whenever visit $p-1$ was itself off drug. At the first
  off-drug visit of a run it is $\mu_{\text{last}}$; at later visits it
  is the previous, already-decayed value. The exponent, however, uses
  the cumulative `d$tsd[p]`.

The pairing on line 88 of an updated (recursive) value with a
cumulative (non-recursive) elapsed time is the source of the
discrepancy. Either element alone would be correct: a recursion with
the one-interval increment, or the cumulative time applied to
$\mu_{\text{last}}$.

### 4.4 What the 2024 revisions changed

Three later commits changed `R/generateData.R` (inspected, GitHub
commit history and per-commit patches):

- **`8609f12` (2024-05-06).** Adds an argument `scalefactor = 2` to the
  signature (line 52 of the current file) and multiplies the exponent
  on the carryover line by it (line 96 of the current file,
  [permalink](https://github.com/rchendrickson/pmsimstats/blob/06dac83f31ca85c994a6a2c1a0b98bc18ef233bd/R/generateData.R#L96)):

  ```r
  brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^(scalefactor * d$tsd[p]/modelparam$carryover_t1half)
  ```

  The argument is documented as `TODO update when understand what this
  does?` (line 43). The same commit, whose message credits edits by
  Ron Thomas, also changes the biomarker correlation rule. That rule
  is the subject of the companion note `docs/55`, but it depends on the
  carryover recursion in a way that matters for the fix (Section 9.5).
- **`f6ee86d` (2024-05-08) and `325314f` (2024-05-10).** Correct an
  index in the same-visit cross-factor correlation
  (`correlations[n2,n1]`) and adjust its spacing. Neither touches the
  carryover lines.

The current upstream code (`master` at `06dac83`, 2024-11-27) therefore
keeps the recursion with cumulative $t_{sd}$ and adds the factor
$s_f = 2$ to the exponent. Its driver `generateSimulatedResults()`
calls `generateData()` without the argument (lines 138 and 144), so
the default $s_f = 2$ applies (inspected). At the $m$-th off-drug
visit of a run the current code gives

$$
\mu_{BR} = \mu_{\text{last}}\;
  2^{-s_f\,(t_{sd,1} + \cdots + t_{sd,m})/t_{1/2}},
$$

which is the published recursion evaluated at half the stated
half-life (derived). The compounding remains, and with $s_f = 2$ the
first off-drug visit of each run is also below the intended value, by
the factor $2^{-t_{sd,1}/t_{1/2}}$ (derived). Table 2 gives the values
on Hybrid path 1. The current upstream file was run with
`empirical = TRUE` (Section 6.1) in all 40 cases of Section 9 and
matched this formula to within $1.6 \times 10^{-14}$ (verified).

| Visit | Week | $t_{sd}$ | Published, $s_f = 1$ | Current, $s_f = 2$ | Intended |
|------|------|------|------|------|------|
| BD3 | 11 | 1 | 5.094 | 2.547 | 5.094 |
| BD4 | 12 | 2 | 1.273 | 0.159 | 2.547 |
| COp | 20 | 4 | 0.268 | 0.017 | 0.268 |

Table: Off-drug drug-response means on Hybrid path 1 at $t_{1/2} = 1$ week: the 58b32a9 code, the current upstream code and the intended model (verified)

## 5. Why it compounds

Consider one off-drug run that follows a last on-drug mean
$\mu_{\text{last}}$, with off-drug visits $m = 1, \ldots, L$ at
cumulative times since discontinuation
$t_{sd,1} < t_{sd,2} < \cdots < t_{sd,L}$. Line 88 gives, with the
Gompertz term zero off drug,

$$
\mu_m = \mu_{m-1}\, 2^{-t_{sd,m}/t_{1/2}}, \quad \mu_0 = \mu_{\text{last}}
\qquad\Longrightarrow\qquad
\mu_m = \mu_{\text{last}}\;
  2^{-(t_{sd,1} + t_{sd,2} + \cdots + t_{sd,m})/t_{1/2}} ,
$$

while the intended model gives $\mu_{\text{last}}\,2^{-t_{sd,m}/t_{1/2}}$
(derived). The ratio of published to intended is

$$
r_m = 2^{-(t_{sd,1} + \cdots + t_{sd,m-1})/t_{1/2}} ,
$$

equal to 1 at $m = 1$ and strictly less than 1 for every $m \ge 2$
(derived).

**Equal spacing.** With visits $\Delta t$ weeks apart, $t_{sd,j} = j\,\Delta t$,
the sum is $\Delta t\, m(m+1)/2$, and

$$
\mu_m = \mu_{\text{last}}\; 2^{-\Delta t\, m(m+1)/(2t_{1/2})},
\qquad
r_m = 2^{-\Delta t\, m(m-1)/(2t_{1/2})}
$$

(derived). Equivalently, the code behaves at the $m$-th visit like the
intended model with an effective half-life

$$
t_{1/2}^{\text{eff}}(m) =
  \frac{t_{sd,m}}{t_{sd,1} + \cdots + t_{sd,m}}\; t_{1/2}
\;=\; \frac{2\,t_{1/2}}{m+1} \quad \text{(equal spacing)},
$$

so the effective half-life is $t_{1/2}$ at the first visit, $2t_{1/2}/3$
at the second, $t_{1/2}/2$ at the third and $2t_{1/2}/5$ at the fourth
(derived). It shrinks with every visit in the run, and the shortfall
depends on how densely the off-drug period is sampled: the same
elapsed time measured at more visits decays further. The intended
model has no such dependence.

## 6. Worked example 1: Hybrid path 1

### 6.1 Set-up

The design is the article's Hybrid ("N-of-1") design as defined in
the package's publication vignette
(`vignettes/Produce_Publication_Results_1_generate_data.Rmd`, lines 48
to 61 at `58b32a9`): visits at weeks 4, 8, 9, 10, 11, 12, 16 and 20,
path 1 (`pathA`) on drug at visits 1 to 4 and 7, off at visits 5, 6
and 8 (inspected). The drug-response parameters are the shipped
`extracted_rp` values for `br`: maximum 10.98604, displacement 5, rate
0.42 (inspected; `data/extracted_rp.rda` at `58b32a9`).

The means in this section were obtained in two independent ways that
agree to $10^{-14}$: by restating lines 83 to 92 as a function, and by
running the published `generateData()` itself with `empirical = TRUE`,
for which `MASS::mvrnorm()` returns a sample whose column means equal
the supplied mean vector exactly [3], and reading off the column means
(verified).

### 6.2 Step by step

The Gompertz means at the on-drug visits are $G(4) = 4.281307$,
$G(8) = 9.222638$, $G(9) = 9.793228$, $G(10) = 10.187065$ and, after
the restart at re-exposure, $G(4) = 4.281307$ again at week 16
(verified). At $t_{1/2} = 1$ week:

- **BD3, week 11, $t_{sd} = 1$.** Published:
  $0 + 10.1871 \times 2^{-1/1} = 5.0935$. Intended: the same, 5.0935.
- **BD4, week 12, $t_{sd} = 2$.** Published: the previous value is
  carried forward and decayed by the cumulative two weeks,
  $5.0935 \times 2^{-2/1} = 1.2734$, which is
  $10.1871 \times 2^{-(1+2)/1}$. Intended:
  $10.1871 \times 2^{-2/1} = 2.5468$. Ratio 0.500.
- **COd, week 16.** On drug; $t_{od}$ restarts at 4 weeks; 4.2813 in
  both.
- **COp, week 20, $t_{sd} = 4$.** The first visit of a new run, so both
  give $4.2813 \times 2^{-4/1} = 0.2676$.

At $t_{1/2} = 0.5$ week the same steps give 2.5468 at BD3 in both, then
$2.5468 \times 2^{-2/0.5} = 0.1592$ (published) against
$10.1871 \times 2^{-2/0.5} = 0.6367$ (intended) at BD4, a ratio of
0.250, and 0.0167 in both at COp (verified).

| Visit | Week | $D_{it}$ | $t_{od}$ | $t_{sd}$ | $G(t_{od})$ | Pub., 1 | Int., 1 | Ratio, 1 | Pub., 0.5 | Int., 0.5 | Ratio, 0.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| OL1 | 4 | 1 | 4 | 0 | 4.281 | 4.281 | 4.281 |  | 4.281 | 4.281 |  |
| OL2 | 8 | 1 | 8 | 0 | 9.223 | 9.223 | 9.223 |  | 9.223 | 9.223 |  |
| BD1 | 9 | 1 | 9 | 0 | 9.793 | 9.793 | 9.793 |  | 9.793 | 9.793 |  |
| BD2 | 10 | 1 | 10 | 0 | 10.187 | 10.187 | 10.187 |  | 10.187 | 10.187 |  |
| BD3 | 11 | 0 | 0 | 1 | 0.000 | 5.094 | 5.094 | 1.000 | 2.547 | 2.547 | 1.000 |
| BD4 | 12 | 0 | 0 | 2 | 0.000 | 1.273 | 2.547 | 0.500 | 0.159 | 0.637 | 0.250 |
| COd | 16 | 1 | 4 | 0 | 4.281 | 4.281 | 4.281 |  | 4.281 | 4.281 |  |
| COp | 20 | 0 | 0 | 4 | 0.000 | 0.268 | 0.268 | 1.000 | 0.017 | 0.017 | 1.000 |

Table: Hybrid path 1: expected drug response by visit under the published recursion (Pub.) and the intended anchored decay (Int.), at $t_{1/2} = 1$ and 0.5 week (verified)

The ratio column agrees with Section 5: BD4 is the second visit of a
run with $\Delta t = 1$, so $r_2 = 2^{-1/t_{1/2}}$, which is 0.5 at
$t_{1/2} = 1$ and 0.25 at $t_{1/2} = 0.5$ (derived). The effective
half-life at BD4 is $2t_{1/2}/3$: 0.67 week instead of 1, and 0.33
instead of 0.5.

![Hybrid path 1: expected drug response under the published recursion (dashed, triangles) and the intended anchored decay (solid, circles), at $t_{1/2} = 1$ week (left) and 0.5 week (right). Points are visit means; the lines join them and do not show the trajectory between visits. Shaded intervals are on drug. Baseline is plotted at 0. The two agree everywhere except at week 12, the second visit of the off-drug run.](figures/53-fig1-hybrid-path1.png)

## 7. Worked example 2: the crossover

The crossover design of the vignette (lines 87 to 98) has eight visits
2.5 weeks apart; path 1 (`pathA`) is on drug for the first four and
off for the last four (inspected). The off-drug run has four visits at
$t_{sd} = 2.5, 5, 7.5, 10$, so by Section 5 the ratios are
$r_m = 2^{-2.5\,m(m-1)/(2t_{1/2})}$: at $t_{1/2} = 2$ weeks,
$1, 2^{-1.25}, 2^{-3.75}, 2^{-7.5}$, that is 1, 0.420, 0.074 and 0.0055
(derived; verified in Table 4). The effective half-life at the fourth
off-drug visit is $2 \times 2/5 = 0.8$ week instead of 2.

| Visit | Week | $D_{it}$ | $t_{od}$ | $t_{sd}$ | $G(t_{od})$ | Pub., 1 | Int., 1 | Ratio, 1 | Pub., 2 | Int., 2 | Ratio, 2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| COa1 | 2.5 | 1 | 2.5 | 0 | 1.848 | 1.848 | 1.848 |  | 1.848 | 1.848 |  |
| COa2 | 5 | 1 | 5 | 0 | 5.922 | 5.922 | 5.922 |  | 5.922 | 5.922 |  |
| COa3 | 7.5 | 1 | 7.5 | 0 | 8.853 | 8.853 | 8.853 |  | 8.853 | 8.853 |  |
| COa4 | 10 | 1 | 10 | 0 | 10.187 | 10.187 | 10.187 |  | 10.187 | 10.187 |  |
| COb1 | 12.5 | 0 | 0 | 2.5 | 0.000 | 1.801 | 1.801 | 1.000 | 4.283 | 4.283 | 1.000 |
| COb2 | 15 | 0 | 0 | 5 | 0.000 | 0.056 | 0.318 | 0.177 | 0.757 | 1.801 | 0.420 |
| COb3 | 17.5 | 0 | 0 | 7.5 | 0.000 | 0.000 | 0.056 | 0.006 | 0.056 | 0.757 | 0.074 |
| COb4 | 20 | 0 | 0 | 10 | 0.000 | 0.000 | 0.010 | 0.000 | 0.002 | 0.318 | 0.006 |

Table: Crossover path 1: expected drug response by visit under the published recursion (Pub.) and the intended anchored decay (Int.), at $t_{1/2} = 1$ and 2 weeks (verified)

At $t_{1/2} = 2$ weeks the intended mean at week 15, five weeks or two
and a half half-lives after stopping, is $10.187 \times 2^{-2.5} =
1.801$ points. The published mean there is 0.757, the value the
intended model reaches only at week 17.5. At week 17.5 the published
mean, 0.056, is already below the intended value at week 20, 0.318
(verified, Table 4).

![Crossover path 1: expected drug response under the published recursion (dashed, triangles) and the intended anchored decay (solid, circles), at $t_{1/2} = 1$ week (left) and 2 weeks (right). Points are visit means joined by lines; shading marks the on-drug period. With four off-drug visits in one run, the published mean falls well below the intended one from the second off-drug visit.](figures/53-fig2-crossover-path1.png)

## 8. Consequences

### 8.1 What changes: the off-drug means

Only the expected drug response at the second and later visits of an
off-drug run changes (derived, Section 5; verified for all paths in
Table 5). Table 5 gives, for each design and half-life, the smallest
published-to-intended ratio and the largest shortfall in points over
all off-drug visits of all paths.

| Design | $t_{1/2}$ (weeks) | Smallest ratio | Largest shortfall (points) |
|------|------|------|------|
| Hybrid | 0.1 | $8.7 \times 10^{-19}$ | $9.7 \times 10^{-6}$ |
| Hybrid | 0.2 | $9.3 \times 10^{-10}$ | 0.0096 |
| Hybrid | 0.5 | $2.4 \times 10^{-4}$ | 0.478 |
| Hybrid | 1 | 0.0156 | 1.273 |
| Hybrid | 2 | 0.125 | 2.238 |
| OL+BDC | 0.1 | $9.3 \times 10^{-10}$ | $1.0 \times 10^{-5}$ |
| OL+BDC | 0.2 | $3.1 \times 10^{-5}$ | 0.0104 |
| OL+BDC | 0.5 | 0.0156 | 0.514 |
| OL+BDC | 1 | 0.125 | 1.370 |
| OL+BDC | 2 | 0.354 | 2.501 |
| CO | 0.1 | $7.0 \times 10^{-46}$ | $9.0 \times 10^{-15}$ |
| CO | 0.2 | $2.6 \times 10^{-23}$ | $3.0 \times 10^{-7}$ |
| CO | 0.5 | $9.3 \times 10^{-10}$ | 0.0096 |
| CO | 1 | $3.1 \times 10^{-5}$ | 0.262 |
| CO | 2 | 0.0055 | 1.044 |

Table: Published versus intended off-drug drug-response means: smallest ratio and largest shortfall over every off-drug visit of every path, by design and half-life (verified)

The two columns tell different stories, and both matter. The ratio
can be extreme at short half-lives, but there the intended carryover
is itself almost zero after a week: at $t_{1/2} = 0.1$ week the largest
intended off-drug mean anywhere is 0.011 points (verified). The
shortfall in points is what reaches the data, and it is negligible at
the half-lives of the article's power analyses, 0.1 and 0.2 weeks (the
`coremodelparams` grid, vignette line 318, inspected). It becomes
material, of the order of one to two and a half points on the CAPS
scale, at half-lives of one to two weeks, which the article does not
report on but which the package accepts.

### 8.2 What does not change

- **On-drug means**, and the first off-drug visit of every run, are
  identical (derived).
- **Variances.** The standard deviations are set at lines 69 to 72 and
  do not depend on the means or the half-life (inspected).
- **The correlation matrix.** The published coupling rule (lines 132
  to 138) correlates the biomarker with the drug response at every
  visit whose mean is nonzero. Both forms are strictly positive at
  exactly the same visits, so the matrix is the same (inspected; the
  pattern of nonzero means is identical in all 40 cases, verified).
  For the same seed, then, the simulated data under the two forms
  differ only by a fixed shift at the affected visits.
- **The biomarker-by-drug interaction in the generated data.** Under
  the published process the interaction is carried by the correlation
  matrix, not by the mean, so the true interaction is unchanged
  (inspected).

### 8.3 Effect on the published analysis

Because only the mean moves, any effect on inference runs through how
the analysis model absorbs a different off-drug mean profile. Two
pieces of evidence bear on this.

**A closed-form argument (quoted).** For the paired-difference
statistic E9, a two-stage proxy for the published mixed model,
`docs/37` (Section 4.4) shows that the means enter the power only
through a between-path term $(\mu_p - \bar\mu)^2$, and reports that
"In the Hybrid and OL+BDC designs the path-mean term is small: in
Hybrid at most 0.30 against residual variances of 41-49, under 0.7% of
the total, decreasing with $t_{1/2}$". The compendium's mathematics
note reaches the same conclusion and labels the effect on the mixed
model "inferred, not tested" (`analysis/report/compendium-mathematics.tex`,
section "Anchored versus compounding means in the closed form").

**A paired simulation (verified; computed for this note).** The
published `generateData()` was run in four versions: as published;
with only lines 83 to 92 replaced by the anchored form of Section 9;
with the graded biomarker coupling used elsewhere in the project
(`docs/37`, Appendix A); and with both changes. Each version fed the
published `generateSimulatedResults()` and `lme_analysis()` (the `lmer`
model with `useDE = FALSE`), with the Hybrid design, the published
parameters at $N = 70$, $c_{bm} = 0.25$, no censoring and 500
replicates per cell. All versions used the same seed at a given
half-life, so each pair of versions differing only in the carryover of
the mean shares identical noise.

| $t_{1/2}$ (weeks) | Published | Anchored only | Graded | Graded and anchored |
|------|------|------|------|------|
| 0 | 0.628 | 0.628 | 0.628 | 0.628 |
| 0.1 | 0.078 | 0.078 | 0.658 | 0.658 |
| 0.2 | 0.082 | 0.082 | 0.670 | 0.670 |
| 1 | 0.088 | 0.088 | 0.486 | 0.484 |
| 2 | 0.078 | 0.072 | 0.344 | 0.344 |

Table: Power of the published analysis for the interaction, Hybrid design, $N = 70$, $c_{bm} = 0.25$, 500 replicates per cell, by version of the 58b32a9 data generator (verified; Monte Carlo standard error at most 0.022)

| $t_{1/2}$ (weeks) | Pair | Mean change in $\hat\beta_{bm:D}$ | Largest change | Discordant decisions |
|------|------|------|------|------|
| 0.1 | anchored vs published | $3 \times 10^{-11}$ | $6.8 \times 10^{-9}$ | 0 of 500 |
| 0.2 | anchored vs published | $1.2 \times 10^{-7}$ | $6.0 \times 10^{-6}$ | 0 of 500 |
| 1 | anchored vs published | $1.2 \times 10^{-5}$ | 0.0017 | 0 of 500 |
| 2 | anchored vs published | $-7.4 \times 10^{-5}$ | 0.0049 | 3 of 500 |
| 1 | graded and anchored vs graded | $5.3 \times 10^{-6}$ | 0.0016 | 1 of 500 |
| 2 | graded and anchored vs graded | $-7.3 \times 10^{-5}$ | 0.0042 | 4 of 500 |

Table: Replicate-by-replicate effect of the anchored correction on the interaction estimate and on the test decision at $\alpha = 0.05$, same data otherwise (verified)

Three things follow (verified unless stated).

- **At $t_{1/2} = 0$ all four versions coincide**, as they must: there
  is no carryover to correct, and the shared seed gives identical data.
- **Correcting the mean alone leaves the published analysis
  essentially unchanged.** At the article's half-lives of 0.1 and 0.2
  weeks the interaction estimates move by less than $10^{-5}$ and no
  test decision changes. At 1 and 2 weeks the estimates move by at
  most 0.005, against a mean estimate near $-0.03$ to $-0.13$, and at
  most 4 of 500 decisions change; the corresponding power differences
  (0 to 0.006) are well inside one Monte Carlo standard error. In the
  published model the time term and the drug main effect absorb most
  of the change in the off-drug mean profile (inferred).
- **The power loss under carryover in the published code does not
  come from this recursion.** With the published coupling rule, power
  falls from 0.628 to about 0.08 at any positive half-life, with or
  without the correction. With graded coupling it stays near 0.66 at
  0.1 and 0.2 weeks and declines to 0.49 and 0.34 at 1 and 2 weeks,
  again with or without the correction. The coupling rule, which is a
  separate feature of the code (`docs/37`, Appendix A, where the same
  comparison at $t_{1/2} \in \{0, 0.1, 0.2, 1\}$ with different seeds
  gave published power 0.668, 0.126, 0.092 and 0.094; quoted), governs
  the power; the recursion governs only the off-drug means.

In short, for the article's analysis and parameter grid the recursion
is a discrepancy between the stated and the simulated decay that has
no material effect on the reported power. Its practical importance is
for users who set longer half-lives, who use the simulated means
themselves (for example to plot expected trajectories, as in the
article's Figure 3), or who fit analyses that model the decay.

### 8.4 Where the discrepancy could matter more

The following are inferred, not tested here.

- **Longer half-lives and longer off-drug runs.** The shortfall grows
  with both (Section 5, Table 5). Designs with several closely spaced
  off-drug visits, or half-lives of a week or more, are where the
  published process departs most from its stated parameter.
- **Analyses that model the decay.** An analysis that weights the
  off-drug visits by $2^{-t_{sd}/t_{1/2}}$ at the stated half-life (the
  continuous indicator $D_{bc,it}$ used in later work) assumes the
  anchored form. Data from
  the published process decay faster than that indicator, so the
  analysis is mildly misspecified even when its half-life equals the
  simulation's.
- **Estimating the half-life.** A procedure that estimates
  $t_{1/2}$ from simulated data would recover something closer to the
  effective half-life of Section 5 than to the stated one.

## 9. The fix

### 9.1 The minimal change (exact for exponential decay)

Because exponential decay is memoryless, the recursion is correct if
each step uses the time elapsed since the previous visit. At the first
off-drug visit of a run the previous visit is on drug, where `tsd` is
0, so the increment equals the cumulative value and the first step is
unchanged. A one-line change against `58b32a9`:

```diff
--- a/R/generateData.R
+++ b/R/generateData.R
@@ -85,7 +85,7 @@
         for(p in 2:nP){
           if(!d[p]$onDrug){
             if(d[p]$tsd>0){
-              brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^(d$tsd[p]/modelparam$carryover_t1half)
+              brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^((d$tsd[p]-d$tsd[p-1])/modelparam$carryover_t1half)
             }
           }
         }
```

The product of the increments telescopes,
$\prod_{j \le m} 2^{-(t_{sd,j} - t_{sd,j-1})/t_{1/2}} = 2^{-t_{sd,m}/t_{1/2}}$
with $t_{sd,0} = 0$, which is the intended model (derived).

### 9.2 The anchored form (general)

The one-line change relies on the exponential form. If the decay
family is to be varied, as the article anticipates ("the simulations
could be easily adapted to incorporate an alternative model"), the
decay should be applied once, to the last on-drug value:

```diff
--- a/R/generateData.R
+++ b/R/generateData.R
@@ -82,11 +82,12 @@
     if(c=="br"){
       brmeans<-modgompertz(d$tod,rp$max,rp$disp,rp$rate)
       if(nP>1){
-        for(p in 2:nP){
-          if(!d[p]$onDrug){
-            if(d[p]$tsd>0){
-              brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^(d$tsd[p]/modelparam$carryover_t1half)
-            }
+        last_on<-0
+        for(p in 1:nP){
+          if(d[p]$onDrug){
+            last_on<-brmeans[p]
+          }else if(d[p]$tsd>0){
+            brmeans[p]<-brmeans[p]+last_on*(1/2)^(d$tsd[p]/modelparam$carryover_t1half)
           }
         }
       }
```

This is the form adopted in the pmsimstats-ng package (version 0.3.0,
`R/generateData.R`, lines 276 to 291 in `buildSigma()`, added
2026-09-08; inspected). The rest of that package's data-generating
process differs from `58b32a9` in other respects, so it is not a
clean test of this change on its own; the patched `58b32a9` file is.

### 9.3 Check

Both patched versions of the `58b32a9` file were run with
`empirical = TRUE` on every path of the Hybrid, OL+BDC and CO designs
(eight paths) at $t_{1/2} \in \{0.1, 0.2, 0.5, 1, 2\}$, 40 cases in all,
and their drug-response means compared with
$\mu_{\text{last}}\,2^{-t_{sd}/t_{1/2}}$ computed directly. The largest
absolute difference was $1.2 \times 10^{-14}$ for the anchored form and
$1.1 \times 10^{-14}$ for the one-line form, floating-point rounding
(verified). On Hybrid path 1 at $t_{1/2} = 1$ week both return 5.094,
2.547 and 0.268 at BD3, BD4 and COp, the intended column of Table 3.
At $t_{1/2} = 0$ both forms, like the original, give zero carryover,
since $2^{-t_{sd}/0} = 0$ in R (inspected).

### 9.4 Re-exposure is a separate choice

Neither change alters what happens when the drug is restarted. Time on
drug restarts at zero (Section 4.3), and the on-drug mean is
$G(t_{od})$ with no residual carried in. That is a modeling choice
distinct from the decay discussed here. One interaction is worth
noting. Because the anchored form leaves a larger residual at the end
of an off-drug run, the residual can exceed the restarted Gompertz
value at the first re-exposure visit, producing a downward step in the
mean while the participant is back on drug. On Hybrid path 3 (off at
weeks 10 to 12, on at week 16) the intended mean at week 12 is 3.462 at
$t_{1/2} = 2$ weeks, below the week-16 value of 4.281, but 5.823 at
$t_{1/2} = 4$ weeks, above it (verified). The published recursion gives
1.224 and 3.462 at the same two half-lives. A restart from the residual
would remove the step; it is a separate question from the one in this
note.

### 9.5 Applying the fix to the current upstream code

The fix in Sections 9.1 and 9.2 is written against `58b32a9`. Applied
to the current upstream code it needs one companion change, because
commit `8609f12` made the biomarker correlation depend on the
recursion. In the current code the correlation between the biomarker
and the drug response at an off-drug visit is $c_{bm}$ times the ratio
of consecutive drug-response means,
$c_{bm}\,\mu_{BR,t}/\mu_{BR,t-1}$ (the block marked `RON THOMAS
VERSION`; inspected). Under the compounding recursion that ratio is
$2^{-s_f t_{sd,t}/t_{1/2}}$ with the cumulative $t_{sd}$, so the rule
reproduces a coupling that fades with time since discontinuation, at
half-life $t_{1/2}/s_f$ (derived). Under anchored means the same ratio
is $2^{-s_f (t_{sd,t} - t_{sd,t-1})/t_{1/2}}$, a constant factor per
visit: the coupling falls at the first off-drug visit and then stops
falling (derived). Table 8 shows the consequence on the crossover's
off-drug run.

| Coupling factor at off-drug visit | 1 | 2 | 3 | 4 |
|------|------|------|------|------|
| Intended, $2^{-t_{sd}/t_{1/2}}$ | 0.177 | 0.031 | 0.0055 | 0.0010 |
| Ratio rule, compounding means ($s_f = 1$) | 0.177 | 0.031 | 0.0055 | 0.0010 |
| Ratio rule, anchored means ($s_f = 1$) | 0.177 | 0.177 | 0.177 | 0.177 |

Table: Coupling factor (multiplier of $c_{bm}$) under the ratio rule of the current upstream code, crossover path 1, $t_{sd} = 2.5, 5, 7.5, 10$ weeks, $t_{1/2} = 1$ week (verified)

The two upstream changes are therefore entangled: the correlation
rule is right only while the means compound, and with the default
$s_f = 2$ it fades at half the stated half-life. Fixing the recursion
alone would leave a biomarker that stays associated with the drug
response long after the drug effect has gone. The recommended
correction computes the coupling from the time since discontinuation
directly, $c_{bm}\,2^{-t_{sd}/t_{1/2}}$ off drug after first exposure
and 0 before it, independent of the means (`docs/55`). With that rule
the two fixes are independent and can be applied in either order. It
also removes a side effect of the ratio rule: because the ratio needs
the previous visit, the current code assigns no coupling at visit 1
(the block is guarded by `if (p > 1)`), so the biomarker is
uncorrelated with the drug response at the first visit even when it
is on drug, in every design (`docs/54`, Section 3.4; inspected). The
direct rule assigns $c_{bm}$ there.

## 10. How to verify independently

The snippet below needs only `data.table`, `MASS` and `corpcor`, the
packages the original code already uses, and a checkout of
`rchendrickson/pmsimstats` at `58b32a9` as the working directory. It
does not use pmsimstats-ng. It runs the original `generateData()` with
`empirical = TRUE`, so that the sample means are the code's own mean
vector, and prints them beside the intended anchored values.

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
t_half <- 1
mp <- data.frame(N = 200, c.bm = 0, carryover_t1half = t_half,
  c.tv = .8, c.pb = .8, c.br = .8, c.cf1t = .2, c.cfct = .1)
## empirical = TRUE makes the sample means equal the code's mean vector
dat <- generateData(mp, extracted_rp, extracted_bp, td,
  empirical = TRUE, makePositiveDefinite = TRUE)
code <- colMeans(dat[, paste0(td$timeptnames, '.br'), with = FALSE])
rp <- extracted_rp[cat == 'br']
g <- modgompertz(td$tod, rp$max, rp$disp, rp$rate)
intended <- g; last <- 0
for (p in seq_along(g)) {
  if (td$tod[p] > 0) last <- g[p]
  else if (td$tsd[p] > 0) intended[p] <- last * 2^(-td$tsd[p] / t_half)
}
print(data.frame(visit = td$timeptnames, tsd = td$tsd,
  code = round(code, 3), intended = round(intended, 3)), row.names = FALSE)
```

Its output, run against the `58b32a9` files (verified):

```
 visit tsd   code intended
   OL1   0  4.281    4.281
   OL2   0  9.223    9.223
   BD1   0  9.793    9.793
   BD2   0 10.187   10.187
   BD3   1  5.094    5.094
   BD4   2  1.273    2.547
   COd   0  4.281    4.281
   COp   4  0.268    0.268
```

Changing `t_half`, the path in `ondrug`, or the design reproduces the
other entries of Tables 3 and 4. Applying either change of Section 9
before sourcing makes the two columns agree.

## 11. Limitations

- **One defect, one component.** The note treats only the off-drug
  mean of the drug response. The biomarker coupling rule, which drives
  the article's carryover power results, is a separate matter
  (Section 8.3).
- **Power check scope.** The paired simulation covers one design
  (Hybrid), one sample size ($N = 70$), one interaction strength
  ($c_{bm} = 0.25$), no censoring and 500 replicates per cell. It does
  not cover the article's full grid, its censoring patterns, or values
  of $c_{bm}$ at which the published code repairs the correlation
  matrix.
- **The article's intent** is taken from its text (Section 3). The
  article does not give the carryover as an equation, so the anchored
  equation is our reading of that text; the authors may wish to
  confirm it.
- **Section 8.4** is inferred, not tested.

## 12. Reproducibility

Every number in Sections 1 and 4 to 10 that is labeled verified was
computed by two scripts written for this note, in
`analysis/scripts/carryover-recursion/` and run from the repository
root:

- `53-carryover-examples.R`: Tables 2 to 5, both figures, the
  replica check against the published `generateData()`, the 40-case
  check of both fixes and of the current upstream file, and the
  re-exposure values.
- `53-ratio-coupling-check.R`: Table 8 (Section 9.5), the ratio
  coupling rule of the current upstream code under compounding and
  anchored means, for $s_f = 1, 2$ and $t_{1/2} = 1, 2$ weeks.
- `53-power-check.R`: the paired simulation of Section 8.3, Tables 6
  and 7 (adapted from
  `analysis/scripts/quick-sim/hendrickson-problems/14-graded-patch-58b32a9.R`).

For example,
`Rscript analysis/scripts/carryover-recursion/53-carryover-examples.R out53`
writes the tables to `out53/` and the figures to `docs/figures/`. The
paired simulation was run as

```bash
Rscript analysis/scripts/carryover-recursion/53-power-check.R --reps 500 --cores 6 --out analysis/data/quick-sim/carryover-recursion
```

which takes about 20 minutes and stores its output in
`analysis/data/quick-sim/carryover-recursion/`:
`power-check-reps500.csv` (power, mean estimate and Monte Carlo
standard error by version and half-life; Table 6),
`power-check-reps500-pairs.csv` (the paired changes; Table 7),
`power-check-reps500-raw.csv` (replicate-level estimates and
$p$-values), and the two patched generators, `gd-anchored.R` and
`gd-graded-anchored.R`. Tables 6 and 7 reproduce exactly from these
files (verified). The
scripts read the `58b32a9` sources from
`analysis/data/quick-sim/hendrickson-58b32a9/R/`, which were checked to
be byte-identical to the files fetched from GitHub at `58b32a9`, and
the 2024 file from
`analysis/scripts/quick-sim/hendrickson-original-comparison/vendored-hendrickson-generateData.R`,
which differs from upstream `06dac83` in one semantically identical
line (a modern `data.table` column extraction). The parameter files
`data/extracted_rp.rda` and `data/extracted_bp.rda` in this repository
are identical to those at `58b32a9` (checked with `identical()`).

## 13. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker
   validation: statistical methods and theoretical findings.
   *Frontiers in Digital Health*. 2020;2:13.
   doi:10.3389/fdgth.2020.00013.
2. Hendrickson RC. pmsimstats: R package for simulation of aggregated
   N-of-1 trials. GitHub repository rchendrickson/pmsimstats, commit
   58b32a92b99a8beb145de00c7bd8a0880f06f0f0 (2020-11-19); current
   `master` 06dac83f31ca85c994a6a2c1a0b98bc18ef233bd (2024-11-27).
   https://github.com/rchendrickson/pmsimstats.
3. Venables WN, Ripley BD. *Modern Applied Statistics with S*. 4th ed.
   New York: Springer; 2002. doi:10.1007/978-0-387-21706-2.
4. Senn S. *Cross-over Trials in Clinical Research*. 2nd ed.
   Chichester: Wiley; 2002. doi:10.1002/0470854596.
