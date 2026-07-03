# Revision Plan: Paper 02 and 01/02 Consistency
*2026-06-16 16:40 PDT*

Consolidated, prioritised plan covering every issue raised across the
three review documents in this directory:

- `referee-report-2026-06-16.md` (top-tier referee review, M1-M7)
- `readability-review/readability-consistency-review-2026-06-16.md`
  (undergraduate-reader readability + the A2-dominance finding + 10
  consistency items)
- `sister-paper-consistency-01-02-2026-06-16.md` (cross-paper)

Each issue lists: source, severity, the concrete action, and the files
to touch. Severity: **P0** = factual error, claim is false as written;
**P1** = substantive (changes interpretation or blocks a reader); **P2**
= editorial coherence. Do P0 first; no part of the paper should be
re-rendered for circulation until P0 is clear.

---

## STATUS (updated 2026-06-16 17:14 PDT)

**Applied and verified by re-rendering all three papers (each produced a
fresh `report.pdf`).**

- **All P0 items done.** Paper 02: A1, A2, A3, A4, A5. Paper 01: E1, E2,
  E3. Paper 10: none were P0.
- **Selected P1/P2 done.** Paper 02: B1 (A3 clarification), B2 (CO
  mechanism cross-link, folded into the principal-finding rewrite), C2
  (R-version), C7 (Architecture C), C11 (S6 cross-run drift note), F1
  (CR2-scope restriction). Paper 10: F1 (limitation sentence), F3
  (reciprocal citations + bib entries), F4 (title page added, Type-I
  range harmonised). C9 effectively closed (paper 02 dropped the 40-60
  framing; paper 01 reframed the ranges as internal exploratory).

**Deferred (not yet applied) -- larger work, a re-run, or a decision:**
B3 (shared reference-cell / N-ceiling note, partially covered by E1),
B4/F2 (self-contain calibration or co-submit paper 10 -- needs a
decision), B5 (reader onboarding -- new content), B6 (wire the
corrective figure into report.Rmd), C1 (strip-claudecode -- finalisation
step), C3 (unique-cell count), C4 (Type-I multiplicity argument), C5
(rerun-design note), C6 (estimand sign), C8 (c_bm notation harmonise),
C10 (novelty demarcation), E4/E5/E6/E8 (paper 01 P1/P2). No simulation
re-run was performed; none was required for the items applied.

---

## A. P0 -- factual corrections (the paper currently claims false things)

### A1. The A2 'universal dominance' claim is false
*Source: readability review Sec. 1; refs the A2-not-universal finding.*

- **Problem.** Abstract, Section 5.1, and Conclusion 1 state A2
  'produced the highest power across all 540 cells' and the ranking held
  'across all three designs, both architectures'. Recomputed: A2 is
  strictly highest in 68/360 cells (19%); it wins in 0% of CO cells and
  0% of Architecture-A cells. Worst cell CO/B/Weibull: A2 = 0.484 vs
  A3 = 0.870.
- **Action.** Restate the claim as conditional: *A2 attains the highest
  power under Architecture B in the Hybrid and OL+BDC designs; under the
  crossover design or under Architecture A, A1/A3 are at least as
  powerful and often more so.* Relabel the N=70/Hybrid/B reference cell
  as the most favourable cell, not a representative one. Sweep every
  'across all designs / both architectures / all 540 cells' sentence.
- **Files.** `report.Rmd` Abstract, Sec. 5.1, Sec. 5.5 (architecture),
  Conclusions 1 and 5.

### A2. Conclusion 5 states A2's Architecture-A advantage is 'positive'
*Source: readability review Sec. 4 item 3.*

- **Problem.** Data show A2's mean advantage under Architecture A is
  negative in all three designs (e.g. Hybrid/A t1/2=1.0: A2 = 0.960 vs
  A1 = 0.972).
- **Action.** Correct to: under Architecture A, A2 is marginally
  *dominated*; the binary specification is preferable.
- **Files.** `report.Rmd` Conclusion 5, Sec. 5.5.

### A3. Paper 02 misdescribes the companion's carryover-loss numbers
*Source: sister-paper review Sec. 1.*

- **Problem.** Paper 02 says the companion found Architecture B loses
  40-60 pp and Architecture A 8-13 pp, and cites '0.82 at t1/2=1.0' as
  the companion's B number. Paper 01 reports B losing up to 30.6%
  relative (~24 pp, OL+BDC) and A losing 1.4-2.8%, and explicitly says
  these fall *below* the 40-60%/8-13% literature ranges. The 0.82 is
  Hendrickson's. Paper 02's own data agree with paper 01 (B 13%, A 1.6%
  at N=70 Hybrid), so the intro contradicts the results too.
- **Action.** Replace the three passages (Abstract bullet, Introduction
  orig block, Matched-decay baseline) with paper 01's actual production
  numbers. Fix percent-vs-percentage-point unit conflation. Drop the
  '0.82 = companion' attribution (keep it as Hendrickson's).
- **Files.** `report.Rmd` Abstract, Introduction, Sec. 4.2
  (matched-decay baseline). Ideally compute the quoted cell once and
  reuse in both papers.

### A4. Bias / MSE / coverage scored against the wrong estimand
*Source: referee M1; readability review Sec. 4 item 4.*

- **Problem.** Table 1 scores all three specs against one target
  (theta_true = -3.6), but A1/A3 estimate the `bm:Db` coefficient and A2
  estimates `bm:Dbc` -- different parameters (verified: `true_value`
  identical across specs in all 324 cells). A1's 'bias +1.09 / coverage
  0.83' is an artefact of scoring against A2's target.
- **Action.** Either (a) derive a per-specification population target and
  score bias/coverage/MSE against the matched target, or (b) restrict the
  cross-spec comparison to power and Type I error (both valid per spec)
  and drop or heavily caveat the bias/MSE/coverage columns. Option (b) is
  the honest minimum.
- **Files.** `report.Rmd` Sec. 3.2 (estimands), Sec. 4.1 (headline
  table), `02-summarise-grid.R` if recomputing matched targets.

### A5. Sturdevant-Lumley citation misattributed
*Source: referee M3; sister review Sec. 6.*

- **Problem.** Sec. 4.6 credits Sturdevant & Lumley (2021) with a
  'two-tier reporting convention'. That paper is a mixed-effects *test
  for carryover effects* (TROPHY), not a simulation-reporting method.
- **Action.** Remove the attribution; cite Morris et al. (2019) Sec. 7
  for the convention. Optionally engage Sturdevant-Lumley on its real
  content as a carryover-testing comparator.
- **Files.** `report.Rmd` Sec. 4.6; `references.bib`.

---

## B. P1 -- substantive (interpretation and reader comprehension)

### B1. A3 is structurally A1-plus-a-nuisance-term
*Source: referee M2; readability review Sec. 2.3.*

- **Problem.** A3 (`Sx ~ bm + t + Db + bm:Db + L`) adds `L` only as a
  main effect; the tested interaction is `bm:Db`, identical to A1. So
  'A1 ~ A3' is near-structural, not empirical, and the 'contradicts
  Jones-Kenward' framing is unsupported (the implemented A3 is not the
  parametric-vs-nonparametric contrast).
- **Action.** State plainly what A3 estimates and why it tracks A1.
  Temper all 'contradicts/does not support Jones-Kenward' language. If a
  genuine model-free carryover handling of the interaction is intended,
  implement it (`bm:L`, period-stratified contrast, or contaminated-obs
  exclusion) and re-run; otherwise present A3 honestly as A1 + nuisance.
- **Files.** `report.Rmd` Sec. 3.3.3, 4.3, 5.3, Conclusion 3;
  `simulation-core.R` only if re-specifying A3.

### B2. Add the CO mechanism cross-link to paper 01
*Source: sister review Sec. 3.*

- **Problem.** Paper 01 shows CO neutralises the Architecture-B
  correlation channel (z = 0.29, n.s.); this is exactly why A2 collapses
  under CO. Paper 02 misses the link and instead overclaims dominance.
- **Action.** Cite paper 01's CO result as the mechanism for A2's CO
  failure; pairs naturally with A1 above.
- **Files.** `report.Rmd` Sec. 5.5 / new sentence in Sec. 4 design
  discussion.

### B3. Reconcile the reference cell and N=70 ceiling effect across papers
*Source: sister review Sec. 2.*

- **Problem.** Paper 01 anchors at N=35/OL+BDC, paper 02 at N=70/Hybrid;
  paper 01 Hybrid/B at N=35 (0.944) exceeds paper 02 Hybrid/B at N=70
  (0.860), which reads as contradictory absent the ceiling-compression
  note paper 01 gives (its lines 1141-1149).
- **Action.** Add an identical short 'shared calibration / reference
  cell / N-ceiling' note to both papers so cross-cited numbers reconcile.
- **Files.** `report.Rmd` Sec. 3 (both papers).

### B4. Self-contain the calibration-conservatism result
*Source: referee M5; readability review Sec. 4 item 9.*

- **Problem.** The Type-I-conservatism interpretation and the entire S6
  CR2 recovery depend on the unpublished companion (paper 10). Reviewers
  cannot verify the 6-10% SE inflation.
- **Action.** Add a short self-contained subsection demonstrating the
  conservatism on this grid, or defer S6 and the 'mildly conservative
  lower bound' framing until paper 10 is available for linked review.
- **Files.** `report.Rmd` Sec. 4.5, 4.8 (S6); `references.bib`.

### B5. Undergraduate-reader onboarding
*Source: readability review Sec. 2.2-2.3.*

- **Problem.** Undefined jargon blocks the target reader (aggregated
  N-of-1, carryover, corCAR1/AR(1), estimand, ADEMP, MCSE, half-life,
  Architecture A/B, the OL/CO/Hybrid/OL+BDC abbreviations).
- **Action.** Add a short 'Background' box defining the setting and a
  plain-language A1/A2/A3 comparison table early; move the toy-`Dbc`
  table forward; add a one-paragraph plain explanation of the mixed
  model; add a design-abbreviation legend/schematic.
- **Files.** `report.Rmd` Sec. 1-3.

### B6. Add the corrective 'where A2 loses' exhibit
*Source: readability review Sec. 3.2, Fig. 2.*

- **Problem.** The conditional nature of A2's advantage is invisible in
  the current figures.
- **Action.** Add a figure/table of A2 advantage by design x architecture
  (draft generated at `readability-review/figures/fig2-advantage-by-
  design.png`). Regenerate inside the container for renv reproducibility
  and wire into the figure pipeline.
- **Files.** new figure script under
  `analysis/scripts/carryover-sensitivity/`; `report.Rmd` Sec. 4.

---

## C. P2 -- editorial coherence and consistency hygiene

### C1. Strip the drafting scaffold (both papers)
*Source: referee M7; readability Sec. 2.1; sister Sec. 6.* Run
`strip-claudecode` on `report.Rmd` in 01 and 02; every `rgt` block is
Lorem ipsum. Proof-read the result as continuous prose. Nothing
circulates before this.

### C2. R-version contradiction
*Source: readability Sec. 4 item 1.* Sec. 3.6 says R 4.5.3;
Reproducibility says R 4.4.x. Pick one. `report.Rmd`.

### C3. '540 cells' overstates unique coverage
*Source: referee M4; readability Sec. 4 item 6.* At t1/2=0 all five
decay forms collapse to identical DGPs (`functions.R:54`); 144/540 cells
(27%) are exact duplicates, 396 unique. Report the unique count.
`report.Rmd` Sec. 3.4-3.5.

### C4. Type-I multiplicity argument assumes independence
*Source: referee M6; readability Sec. 4 item 5.* The 540 null cells
share common random numbers within cell and duplicate across decay forms
at t1/2=0. Soften the 'expected maximum among 540 independent draws'
argument or replace with a corrected per-cell test. `report.Rmd` Sec. 4.4.

### C5. High-precision rerun is on a different design
*Source: referee M minor; readability Sec. 4 item 7.* The McNemar
p = 6.6e-17 is on OL+BDC, not the Hybrid reference cell. State this.
`report.Rmd` Sec. 4.7.

### C6. Architecture-B estimand sign and 'dimensionless' description
*Source: readability Sec. 4 item 10.* Sec. 3.2 calls it a dimensionless
correlation but scores bias in BR units with a negative sign (-3.6).
Explain the convention. `report.Rmd` Sec. 3.2.

### C7. Acknowledge Architecture C as out of scope -- DONE 2026-06-16
*Source: sister review Sec. 4.* Paper 02 said it retains 'both' DGP
architectures; the companion now has three. **Fixed** in `report.Rmd`
Sec. 3.1 (Data-generating mechanisms): bullets header and orig block now
state the two single-channel architectures (A, B) are retained and
Architecture C is out of scope. The remaining 'both architectures' usages
(L285, L1825, L2050) correctly denote the two architectures this paper
evaluates and are left as-is; L1825 is corrected separately under A1.

### C8. Harmonise c_bm notation across papers
*Source: sister review Sec. 5.* Paper 01 uses c_bm,A / c_bm,B; paper 02
uses bare c_bm. Disambiguate centrally. `report.Rmd` (both).

### C9. The uncited '40-60% / 8-13% literature' range
*Source: sister review Sec. 6.* Both papers attribute these to a
'broader N-of-1 simulation literature' with no citation; may be an
internal legacy number (archived figure5 / architecture_comparison).
Find a real source or remove from both. `report.Rmd` (both);
`references.bib`.

### C10. Novelty demarcation and missing comparators
*Source: referee M (novelty); sister Sec. 6.* Demarcate the 'no
systematic evaluation exists' claim against Hendrickson (2020) and the
2019 aggregated-N-of-1 carryover-simulation comparison (PMC6955665).
`report.Rmd` Introduction; `references.bib`.

---

## E. Paper 01 issues (line-anchored, verified against archived data)

Full §3.1 power table, Type I range, z-statistics, and the Architecture C
grid were recomputed and match the paper, except where noted. Several
items from paper 01's own 2026-06-13 referee report are already fixed (z
defined, super-additivity baseline, R version, bib entries).

### E1. P0 -- 'N = 35' is not the analysed sample size
*Verified by recompute + driver (`01-dgp-prototype.R:145-155`).* §3.1/§3.3
draw N = 35 *per randomisation path* and concatenate, so the effective
analysis sample is **70** (CO, OL+BDC; 2 paths) and **140** (Hybrid; 4
paths), disclosed only in the Reproducibility note (L2858-2862). The
Abstract (L68, L81-84), Setup (L801), and Conclusions present 'N = 35' at
face value. This is the root cause of the cross-paper power 'inversion'
(01 Hybrid B 0.944 at 140 subjects vs 02 Hybrid B 0.860 at 70).
**Fix.** State effective N wherever 'N = 35' appears; relabel tables with
effective N or add '(effective N = 70/140 after path concatenation)'.

### E2. P0 -- Architecture C 'boundary validation gate' claim is false
*Verified false by recompute.* §3.4 (L1252-1258) and Table 1 caption
(L1267-1269) state the boundary cells 'reduce the combined DGP to pure
Architecture A/B' and 'serve as internal consistency gates'. The pure-A
boundary (c_bm,a=0.45, c_bm,b=0, t1/2=1.0) gives CO 0.172 / Hybrid 0.389 /
OL+BDC 0.276 versus §3.1 Architecture A 0.725 / 0.978 / 0.730 -- a 35-59
pp gap, because §3.4 uses N=35 *total* while §3.1 uses per-path N. The
Reproducibility note (L2863-2869) now admits §3.4 'is not directly
comparable', contradicting the main-text gate language.
**Fix.** Delete the validation-gate wording; state boundary cells recover
the pure-DGP *form* but at a different sample size, so they are not
numerical checks against §3.1. Optionally re-run Architecture C at matched
per-path N to make the gate real.

### E3. P0 -- uncited '40-60% / 8-13% broader literature' ranges
*Verified: no citation anywhere.* L863, L867-868, L882-883, L935-937, and
the §5.1 table (L1866). Asserted as established literature values the
paper's results fall below; likely internal legacy numbers. This is the
upstream source of paper 02's A3 misquote.
**Fix.** Cite the specific source or delete the comparison in both papers;
if internal, attribute to the archived exploratory runs explicitly. (Same
as C9; resolve once, propagate to 01 and 02.)

### E4. P1 -- §6 'six analysis strategies' has no simulation evidence
L2027-2353; benefit table L2319-2344. The text concedes 'None has been
evaluated by simulation', but the table assigns benefit directions with no
data. **Fix.** Run the two tractable candidates (on-drug-only, weighted)
on existing infrastructure, or demote the table to prose. Note paper 02 is
effectively the simulation §6 calls for (A1/A2/A3); cross-reference it.

### E5. P1 -- sample-size convention differs from paper 02 (cross-paper)
*Verified by drivers.* §3.1/§3.3 use per-path N (effective 70/140); §3.4
and paper 02 use total N. The label 'N=35' means 35 (§3.4) vs 70-140
(§3.1) within one paper. **Fix.** Adopt one convention across both papers
or add an identical 'effective N' note to each. (Three-way item.)

### E6. P1 -- cross-paper estimand scale mismatch
*Suspected (inspected, not re-derived).* Paper 01's interaction
coefficient beta(bm:Dbc) approximately -0.234 for c_bm=0.45; paper 02
calibrates the same nominal estimand to theta_true = -3.6 (about 15x
larger). Biomarker on different scales (raw vs standardized) between
papers. **Fix.** State the biomarker scaling in each paper; note
non-comparability of coefficient magnitudes across the pair.

### E7. P1 -- Architecture C is a real production run, not a placeholder
*Verified:* C is a 27-cell production grid with data present, ahead of the
CLAUDE.md 'open work item' note. Reinforces C7 (paper 02's 'both
architectures' wording is stale).

### E8. P2 -- editorial
- Lorem ipsum in every `rgt` block (L114-116 ... L2755); run
  `strip-claudecode` (same as C1). Verify §2-3 markup is `:::` fenced.
- R version: paper 01 = R 4.6.0 (L2845) vs paper 02 = 4.5.3 / 4.4.x
  (folds into C2 as a three-way item).
- Missing methodological citations: corCAR1/AR(1) (Pinheiro & Bates
  2000), exponential PK decay, finite-mixture claim §4.2 (McLachlan &
  Peel 2000); Rizopoulos 2012 §7.2 still a stretch.
- Novelty 'Hendrickson unique in using Architecture B' (L211-215,
  L2504-2524, L2806-2814) rests on an informal survey; document the
  search or soften to 'we are not aware of'.
- §3.3 N=70 robustness table (L1106-1107) not in the verified summary
  file (separate 800-rep run); confirm archived and add MCSE columns.

---

## F. Paper 10 issues (line-anchored, verified against archived data)

**Verdict: paper 10 substantively and numerically supports what paper 02
attributes to it.** Every borrowed number reproduces from
`output/*.rds`: SE over-estimated 6-10% (kappa = 1.07/1.10/1.06 for
A1/A2/A3), realised Type I 0.029-0.039, conservatism isolated to the
corCAR1 layer (rand-int kappa 0.97-0.99 / Type I 0.049-0.055 vs +corCAR1
kappa 1.07-1.10 / Type I 0.029-0.033), CR2 restores size (kappa
0.976-1.000, Type I 0.045-0.049, df approximately 23.7), df 487 -> 23.7
(N=70) / 11.5 (N=35). Paper 10 is internally consistent and well-archived.
**No P0 issues in paper 10.** The problems are scope and sequencing.

### F1. P1 -- CR2 generalisation validated at one cell only
*Verified by recompute; paper 10 L1759-1764, L1789-1803.* Paper 10 checks
CR2 ranking-invariance and the 'mildly conservative lower bound' property
solely at Architecture B / Hybrid / N=70 / c_bm=0.45 (null + one alt
cell), and its own Limitations concede uniformity across the factorial
'was not established'. Paper 02 applies that framing to its **entire grid,
including CO and Architecture A** -- exactly the strata where A2's ranking
reverses (item A1). **Fix.** In paper 02, restrict the S6 'ranking
preserved / lower bound' claim to the cells paper 10 validated; do not
extend it to CO or Architecture A. (Couples to A1.)

### F2. P1 -- the whole pair depends on an unpublished, unfinished draft
*Inspected.* Paper 10 is `@unpublished` and itself a Lorem-ipsum draft.
Reviewers of 01/02 cannot verify the 6-10% / CR2 results. **Fix.**
Co-submit 10 with 01/02, or fold a self-contained one-paragraph
calibration result into 02 (= item B4). The S6 *table* is sound (S6
recompute matches: ref kappa_mod 1.152 -> CR2 0.992, df 23.2; rho=0.9
kappa 1.255; dropout kappa_mod 0.993 -> CR2 0.855, df 4.8, power 0.148 ->
0.126); only the generalisation sentence over-reaches.

### F3. P1 -- cross-citation is one-directional
*Verified: `references.bib` in paper 10 has hendrickson2020, morris2019,
jones2014, senn2016; no pmsimstats / carryover-sensitivity key.* Paper 02
cites 10 as 'the companion calibration study'; 10 never formally cites 01
or 02. **Fix.** Add reciprocal citations 10 -> 01 and 10 -> 02.

### F4. P2 -- editorial
- No `title-page.md` in paper 10 (01 and 02 have one); add for parity.
- Lorem ipsum in every `rgt` block (L99-145 ...); `strip-claudecode`.
- Internal: Type-I range stated as '0.029 to 0.035' (L864) and '0.029 to
  0.033' (L976) for the same cell; harmonise.

---

## G. Three-way consistency matrix (01 / 02 / 10)

The shared elements that must agree across all three papers. 'Action'
points to the item ID that fixes it.

| Shared element | Paper 01 | Paper 02 | Paper 10 | Action |
|---|---|---|---|---|
| Carryover power loss | B up to 30.6% rel, A 1.4-2.8% | own data B 13% / A 1.6%; **intro wrongly says 40-60 pp / 8-13 pp** | n/a | A3, E3 |
| '40-60% / 8-13%' literature range | uncited; results fall below | quoted as companion's findings | n/a | E3 / C9 (cite or delete in both) |
| Reference cell | OL+BDC | Hybrid | Arch B Hybrid | B3 (state relationship) |
| Sample-size convention | **per-path** N=35 (eff 70/140) in §3.1; total in §3.4 | total N=70 | total N=70 | E1, E5 (unify; state effective N) |
| Type-I conservatism / kappa / df / CR2 | reports Type I [0.010,0.074] | leans on 10; S6 table sound | **establishes** kappa 1.07-1.10, CR2->nominal, df 480->23/12 | sound; keep |
| CR2 ranking-invariance scope | n/a | claims whole grid incl CO / Arch A | **validated only Arch B/Hybrid/N70** | F1 + A1 (restrict 02) |
| A2 power, Arch B/Hybrid/N70/t1/2=1.0 | n/a (diff cell) | **0.860** (headline) vs **0.829** (S6) | **0.843** (alt cell) | C11, P10-4 (note independent MC runs + MCSE) |
| c_bm notation | c_bm,A / c_bm,B / z_B / sigma_BR | bare c_bm | uses 02's | C8 (harmonise) |
| Biomarker scale / coef magnitude | beta approx -0.234 (raw) | theta = -3.6 (standardized, sigma_BR=8) | uses 02's | E6 (state scaling; non-comparable) |
| Architecture C | real 27-cell production run | absent ('retains both') | n/a | C7 / E7 (acknowledge out of scope) |
| R version | 4.6.0 | 4.5.3 / 4.4.x | confirm | C2 (unify across three) |
| Cross-citation | cites neither 02 nor 10 | cites 10 as companion | cites neither 01 nor 02 | F3 (make reciprocal across all three) |
| Draft state | Lorem ipsum | Lorem ipsum | Lorem ipsum | C1 / E8 / F4 (strip all three) |
| Title page | present | present | **absent** | F4 |
| Venue / author | SiM / Persons / pmsimstats | SiM / Persons / pmsimstats | SiM / pmsimstats | consistent |

**Three-way coherence priorities.** (i) The N-convention split (E1/E5)
and the carryover-loss misquote (A3/E3) are the two issues that make the
01/02 pair look mutually contradictory; fix both and the headline numbers
reconcile. (ii) F1 + A1 together: paper 02 must not borrow paper 10's
lower-bound guarantee for the CO / Architecture A strata that neither
paper 10 validated nor paper 02's own data support. (iii) The 0.860 /
0.843 / 0.829 drift (C11) is within Monte Carlo noise but reads as
contradiction; quote a single shared run or annotate with MCSE.

### C11. P2 -- paper 02 internal: same cell, three power values
*Verified.* Paper 02 quotes A2 at Arch B/Hybrid/N70/c_bm=0.45/t1/2=1.0 as
0.860 (headline grid) and 0.829 (S6 table, `diag-s6-cr2`); paper 10 gives
0.843 (alt cell). The A2-A1 gap change under CR2 is +0.080 -> +0.092
(widens) in paper 02 S6 vs +0.086 -> +0.082 (narrows) in paper 10 --
opposite directions of a sub-MCSE effect. **Fix.** State these are
independent MC runs (give MCSE), or reuse one run; never describe the
within-noise gap change directionally. (Listed under section C scope.)

---

## D. Suggested execution order

Spanning all three papers. P0 false-claim corrections live only in paper
02 (A1-A5) and paper 01 (E1-E3); paper 10 has none.

1. **Paper 02 P0 batch (A1-A5)** in `report.Rmd` prose + the Table 1
   estimand fix. Re-derive every number quoted in Abstract/Conclusions
   from the archived summaries.
2. **Paper 01 P0 batch (E1-E3)**: relabel effective N, delete the
   Architecture-C validation-gate claim, resolve the uncited 40-60%/8-13%
   range. E1/E3 are prerequisites for the cross-paper batch.
3. **Decide A3's fate (B1)** -- if re-specifying, that re-run gates B6's
   figure and parts of Sec. 4. Resolve before touching the 02 results
   narrative.
4. **Three-way sync batch** (edit 01, 02, 10 together so they stay
   consistent -- section G): N convention (E1/E5/B3), carryover-loss
   numbers (A3/E3/C9), CR2-scope restriction (F1+A1), notation (C8),
   biomarker scale (E6), Architecture C (C7/E7), reciprocal citations
   (F3), R version (C2), the 0.860/0.843/0.829 drift (C11/P10-4).
5. **Reader onboarding + corrective exhibit (B5, B6)** in paper 02.
6. **Self-containment (B4 / F2)** -- fold a one-paragraph calibration
   result into 02, or co-submit 10; add 10's title page (F4).
7. **Editorial sweep** across all three (C1-C6, C10, E8, F4), finishing
   with `strip-claudecode` on 01, 02, and 10 and a clean re-render.

**Re-run scope.** No simulation re-run is required for any P0 item; all
three papers' archived results are mutually consistent once the prose is
corrected to match them (paper 02's carryover-loss intro, paper 01's N
labels, paper 02's CR2 generalisation). Re-runs are optional and only for:
A3 re-specification (B1), Architecture C at matched per-path N to make the
boundary gate real (E2), and paper 01's §6 analysis-strategy table (E4).

---
*Rendered on 2026-06-16 at 16:40 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/02-carryover-sensitivity/revision-plan-2026-06-16.md*
