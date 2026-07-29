---
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.4em}
---

# White paper: notation and terminology across the compendium

*2026-07-29 10:33 PDT*

**An audit of mathematical notation and controlled vocabulary in the
eleven manuscripts under `analysis/report/`, with the nine verified
inconsistencies and a proposed canonical convention.**

## Method and headline finding

Every `report.Rmd` was searched for the symbols carrying the
programme's substantive quantities (outcome, drug indicator, biomarker,
moderation parameter, estimand, decay rate, half-life, sample size) and
for the labels naming its designs, architectures, and simulation
studies. Counts below are literal occurrence counts in the sources.
Nothing here rests on inference about intent.

The answer to the question posed is **no**. The notation is
locally coherent within each paper and internally correct, but it is
not consistent across the compendium. Nine distinct inconsistencies
are verified below, of which three are substantive (the same symbol
carries different meanings in different papers) and six are cosmetic or
organizational (the same object carries different symbols or labels).
The substantive three are the drug indicator, the moderation
parameter, and the half-life unit. All three can be repaired mechanically
without re-running any simulation, and all three should be repaired
before the papers circulate together, because each paper cites its
neighbors and a referee reading two of them will meet the collision
directly.

## The core system, as intended

Underneath the variation there is a single coherent system, most fully
realized in papers 01 and 06.

| Object | Symbol | Definition |
|---|---|---|
| Participant, timepoint | $i$, $t$ | index over participants and weekly occasions |
| Observed outcome | $Y_{it}$ | symptom score; code column `Sx` |
| Baseline | $\mathrm{BL}_i$ | participant-specific pre-treatment level |
| Response components | $BR$, $PB$, $TV$ | pharmacological, placebo-belief, natural-history; each a modified Gompertz, each a non-negative reduction in symptoms |
| Binary drug state | $D_{it}$ | 1 on drug, 0 off drug |
| Continuous exposure | $D_{bc,it}$ | 1 on drug, $e^{-\lambda t_{sd}}$ off drug; code `Dbc` |
| Time since discontinuation | $t_{sd}$ | weeks since the last on-drug occasion; code `tsd` |
| Time on drug | $t_{od}$ | cumulative on-drug time; Gompertz input; code `tod` |
| Biomarker | $B_i$ | participant-level; standardized $b_i = (B_i - \bar{B})/s_B$ |
| Moderation, covariance channel | $c_{bm}$ | correlation between $B$ and $BR$ on drug; Architecture B |
| Moderation, mean channel | $\beta_{bm}$ | multiplier on the treatment effect; Architecture A |
| Decay rate | $\lambda$ | $\ln 2 / t_{1/2}$; code `lambda_cor` |
| Carryover half-life | $t_{1/2}$ | code `carryover_t1half` |
| AR(1) correlation | $\rho$ | within-factor serial correlation |
| Target estimand | $\beta_{bm:D}$ | the `bm:Dbc` interaction coefficient |
| Sample size | $N$ | participants |

Terminology is more stable than notation. Four design families recur
throughout, with fixed meanings: **OL** (open-label), **OL+BDC**
(open-label plus blinded discontinuation), **CO** (traditional
two-period crossover), and **Hybrid** (the Hendrickson 2020 design:
open-label titration, blinded discontinuation, then a brief crossover).
Three data-generating architectures recur with fixed meanings:
**A** (direct mean moderation), **B** (MVN differential correlation),
**C** (combined, with independent channel weights $c_{bm,A}$ and
$c_{bm,B}$). The three-component vocabulary (BR, PB, TV) is used
identically in every paper that invokes it. ADEMP terminology follows
Morris, White, and Crowther (2019) uniformly, and Monte Carlo standard
error is abbreviated MCSE throughout.

## The three substantive inconsistencies

**N1. The drug indicator carries two incompatible meanings.** Papers
01, 02, and 06 use $D_{it}$ for the *binary* on/off state and reserve
$D_{bc}$ or `Dbc` for the *continuous* exponentially decaying exposure.
Paper 03 uses $D_{it}$ for the continuous form, writing
$D_{it} = \exp(-\lambda t_{sd,it})$ off drug and 1 on drug, three
times. A reader moving from 01 to 03 encounters the same symbol with
the opposite referent, and the distinction is precisely the object of
paper 02's A1-versus-A2 comparison, so the collision sits on the
programme's central axis. Compounding this, three orthographies for
the continuous indicator coexist: the code literal `Dbc` (01, 03, 04,
05, 06, 07, 08), the math form $D_{bc}$ (02, 06, 07, 09), and $D_{it}$
(03). Papers 06 and 07 use all three within a single document.

**N2. $c_{bm}$ is used for a correlation and for a regression
coefficient.** Paper 01 is careful: $c_{bm}$ is a correlation, bounded
in $[-1,1]$, belonging to Architecture B, while $\beta_{bm}$ is the
dimensionless multiplier of Architecture A. Downstream this discipline
lapses. Paper 02's parameter table collapses the row to
"$c_{bm}$ / $\beta_{bm}$". Papers 07, 08, and 09 report Architecture-A
runs parameterized as $c_{bm} \in \{0, 0.45\}$; paper 08's simulation
is Architecture A only, yet $c_{bm}$ appears thirteen times and
$\beta_{bm}$ never. There is a defensible reason for the slippage,
namely that the Architecture-A shift is implemented as
`c.bm * bm_z * br_sd` on a standardized biomarker scaled by the BR
standard deviation, which is exactly the calibration that makes the
numeric value 0.45 comparable across architectures. But that
justification is stated only in paper 05, and it is stated in code
identifiers. Papers that use $c_{bm}$ for a mean-channel coefficient
without repeating the calibration argument are asserting an
equivalence they have not established locally.

**N3. Carryover half-life is quoted in incompatible units.** Papers 01,
02, and 07 use weeks, with the canonical grid
$t_{1/2} \in \{0, 0.5, 1.0\}$. Papers 04 and 05 use days, with 3 days
as the baseline and sweeps over $\{0.5, 1, 3, 5, 7\}$ and up to 14
days. No conversion statement appears in either group. The values are
numerically confusable: a "half-life of 1" means one week in paper 02
and one day in paper 04, a sevenfold difference in the quantity that
the compendium exists to study. Paper 05's own text mixes the two,
carrying thirteen day-denominated quantities and one week-denominated
quantity.

## The six organizational inconsistencies

**N4. Four symbols for one estimand.** The biomarker-by-treatment
interaction appears as $\beta_{bm:D}$ (03, 07, 08, 09),
$\beta_{bm}^{BR}$ and $\beta_{bm}^{\text{lumped}}$ (06, where the
superscript indexes the component, a genuine extension rather than a
conflict), $\beta_3$ (01), $\theta$ and $\hat\theta$ (02, 10), and the
bare code literal `bm:Dbc` in running prose throughout. Paper 10's use
of $\theta$ is defensible, since its argument is deliberately generic
over interaction coefficients, but paper 02 uses $\theta_{\text{true}}$
for a quantity that papers 07 to 09 call $\beta_{bm:D}$ with no bridging
sentence.

**N5. Outcome notation alternates between mathematics and code.**
$Y_{it}$ appears in 01, 02, 06, 07, and 09; the code identifier `Sx`
appears in model formulas in 01, 02, 03, 04, 05, 06, 07, and 09. Papers
08 and 10 set no outcome symbol at all. Mixing the two is acceptable in
a reproducibility section; it is not acceptable in a model definition,
where papers 02 and 06 both write the fitted model as
`Sx ~ bm + t + Dbc + bm:Dbc` immediately after defining $Y_{it}$.

**N6. Label collisions across three families.** The letters A, B, C
name architectures, but paper 02 simultaneously uses A1, A2, A3 for its
three analysis specifications, 192 occurrences of "A2" alone, in a
document that also discusses Architecture A. The prefix S names
sensitivity blocks S1 to S6 in paper 02 and sweeps S1 to S12 in paper
05, with different referents at every shared index. Simulation studies
are labeled "Study A" and "Study B" in paper 06 but "Study 1" through
"Study 5" in papers 03 and 09.

**N7. Two sample-size conventions.** Paper 01 reports $N = 35$ *per
randomization path*; papers 04 and 05 report $N = 35$ *per design*,
that is, in total. Paper 11 documents its own instance of this
explicitly, noting that its Architecture C grid "uses the earlier
driver's total-sample convention" and that matched runs are pending, so
its cells are not readable against paper 01's tables. This is the one
inconsistency the compendium has already caught in writing.

**N8. Design labels alternate between abbreviation and prose.** Papers
01, 02, and 07 use CO, Hybrid, and OL+BDC as abbreviations; papers 04,
05, and 09 predominantly spell out "traditional crossover", "hybrid",
and "open-label". Capitalization of "Hybrid" is not settled even within
papers: 01 has 21 capitalized against 3 lowercase, 05 has 12 against
128.

**N9. The sign convention is stated only in paper 06.** Paper 06 fixes
that the three components are non-negative reductions, so an increase
in any component lowers $Y$. Consequently effect sizes are negative
throughout the compendium ($-2.0$ nightmares per week in 04,
$\theta_{\text{true}} = -3.6$ in 02, $\approx -0.23$ bias in 07) while
moderation parameters are positive. No other paper states the
convention, so a reader of 07 alone cannot tell whether a negative bias
estimate means attenuation or inflation.

## Recommendation

Adopt the table in Section 2 as the compendium's canonical notation and
apply it by a single mechanical pass, in this order.

1. **Fix N1 first.** Reserve $D_{it}$ for the binary state and
   $D_{bc,it}$ for the continuous exposure everywhere, retire the bare
   `Dbc` from mathematical display (keeping it only in code listings
   and reproducibility sections), and correct paper 03's three
   occurrences. This is the only collision that can cause a reader to
   misread a result.
2. **Fix N3 next.** Convert papers 04 and 05 to weeks, or state a
   conversion in both groups' Methods. Weeks is the better target,
   since it matches the outcome's measurement interval and the
   canonical $\{0, 0.5, 1.0\}$ grid.
3. **Fix N2 by rule.** Use $c_{bm}$ only for the Architecture-B
   correlation and $\beta_{bm}$ only for the Architecture-A multiplier.
   Papers 07, 08, and 09 relabel to $\beta_{bm}$. Add one sentence
   stating the calibration $\beta_{bm} \leftrightarrow c_{bm}$ via
   $b_i$ and $\sigma_{BR}$ to any paper that reports both.
4. **Standardize N4 to $\beta_{bm:D}$**, allowing paper 10 to retain
   $\theta$ with an explicit bridging sentence, and paper 06 to retain
   its component superscripts.
5. **Relabel to remove N6.** Rename paper 02's specifications to
   **M1, M2, M3** (analysis *methods*), freeing A/B/C for
   architectures; prefix sweep labels by paper, as **02-S1** and
   **05-S1**, or rename paper 02's blocks to B1 to B6; and adopt
   "Study 1, 2, ..." uniformly.
6. **Settle N5, N7, N8, N9 editorially.** Mathematics in model
   definitions and code identifiers only in reproducibility sections;
   state the sample-size convention in every Methods section and
   resolve paper 11's pending re-run; use the abbreviations CO, Hybrid,
   OL, OL+BDC at first definition in every paper and capitalize
   Hybrid; and state the sign convention in the Methods of every paper
   that reports a signed quantity.

A single shared file, either a notation section in
`analysis/report/README.md` or a `notation.tex` include alongside
`sim-preamble.tex`, should hold the canonical table so that later
papers inherit it rather than re-deriving it. The project lexicon at
`docs/29-nof1-precision-medicine-lexicon.md` already performs this role
for prose terminology and is the natural place to cross-reference.

---
*Rendered on 2026-07-29 at 10:33 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/whitepaper-notation-audit.md*
