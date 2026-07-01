# Architecture C extraction plan: paper 01 to paper 11
*2026-07-01 08:36 PDT*

Plan to remove Architecture C from the comprehensive paper 01
(`report.Rmd`), replace it with a short forward-pointing unifying
paragraph, and seed a dedicated companion paper (paper 11). This
sharpens paper 01 to its keystone A-versus-B contrast, removes the
manuscript's single largest referee liability (MC-1, the false
boundary-validation claim), and gives Architecture C a home where
the super-additivity result can receive the analytical grounding the
referee requested.

The same reduction applies to `report-slim.Rmd` (tracked separately
below); its `report-slim` cuts are lighter because §6 there is
already reframed.

## Rationale recap

- Architecture C is modular, not load-bearing: the A-vs-B argument
  has zero dependency on it. C is introduced only as a generalization
  after the contrast is established.
- C caused MC-1 (the one "demonstrably false" referee finding) and
  drew the "thinner ... descriptive rather than explanatory" critique
  in the Significance section. Removing it from paper 01 resolves MC-1
  at the source.
- C fits the program's spin-off pattern (papers 02-10 each drill into
  one element that surfaced while working the 01 question).

## Part 1. Cuts to `report.Rmd`

Line ranges are current as of 2026-07-01. Execute in descending line
order so earlier edits do not shift later ranges. All ranges verified
against section boundaries.

| # | What | Lines | Action |
|---|------|-------|--------|
| C1 | §2.3 "Architecture C: Combined ..." full subsection (header, 2x bullets/rgt/orig) | 476-585 | Delete; insert pointer paragraph (Part 2) in its place |
| C2 | §2.4 "Under Architecture C" block (bullets 664-679, rgt 681-689, orig 691-706) | 664-707 | Delete |
| C3 | §3.4 "Architecture C: combined-channel power results" full subsection incl. Table 1 (`tab:arch-c-power`) | 1226-1381 | Delete (retain the section-3 closing rule at 1382) |
| C4 | §4.4 "When Architecture C is appropriate" full subsection | 1790-1858 | Delete (retain the section-4 closing rule at 1859) |
| C5 | §5.1 decision table: 4-column (Criterion/A/B/C) | 1869-1892 | Rewrite to 3-column (drop the Architecture C column); see Part 2 |
| C6 | §5.1 C mixing-weight guidance (bullets 1894-1908, rgt 1910-1918, orig 1920-1931) | 1894-1931 | Replace with short A-vs-B guidance + pointer (Part 2) |
| C7 | §5.2 reporting item referencing "for Architecture C, independent evidence ..." | 1943 | Reword item to drop the C clause |
| C8 | Reproducibility (§8): Architecture C driver/data documentation | 2871, 2882-2883, 2893, 2901-2904 | Remove from paper 01; migrate to paper 11 reproducibility |

Renumbering after cuts:

- §2: removing 2.3 leaves 2.1, 2.2, and Mathematical Comparison.
  Renumber "2.4 Mathematical Comparison" -> "2.3 Mathematical
  Comparison" (reverts the earlier bump noted in CLAUDE.md).
- §3, §4: the removed subsections (3.4, 4.4) are last in their
  sections, so no renumbering.
- §5: table and guidance edited in place, no renumbering.

No change required to:

- Abstract: already A-vs-B only ("We formalize two architectures").
- Conclusions (§8 body): already A-vs-B only.

Verification after cutting: `grep -niE 'architecture[~ ]?c|c_\{bm,|c\.bm_|combined|arch-c' report.Rmd` should return only intentional pointer references.

## Part 2. Replacement text for `report.Rmd`

### 2a. Pointer paragraph (replaces C1, end of §2.2)

Place after §2.2 closes, before the "2.3 Mathematical Comparison"
header:

> A third possibility is that both channels operate at once: the
> biomarker could scale the mean drug effect (Architecture A) while
> simultaneously altering the biomarker-response correlation during
> active treatment (Architecture B). Such a combined architecture,
> carrying independent weights on the two channels, recovers A and B
> as its single-channel limits and is the natural framework when the
> operative mechanism is uncertain. Because it introduces a second
> free parameter and its power behavior warrants separate treatment,
> we develop it in a companion paper [CITE paper 11] and restrict the
> present paper to the A-versus-B contrast.

### 2b. §5.1 decision table (replaces C5), 3-column form

```
\begin{table}[htbp]
\centering
\begin{tabular}{lcc}
\toprule
Criterion & Architecture A & Architecture B \\
\midrule
Biomarker role   & Causal mediator       & Statistical predictor \\
Mechanism        & Deterministic scaling & Probabilistic association \\
Signal location  & Mean structure        & Covariance structure \\
Free parameters  & 1 ($c_{bm}$)          & 1 ($c_{bm}$) \\
Carryover impact & Modest (1--3\%)       & Design-dependent, up to 31\% \\
Appropriate when &
  \parbox[t]{3cm}{Biomarker\\determines dose\\or PK} &
  \parbox[t]{3cm}{Biomarker\\indexes a latent\\subtype} \\
\bottomrule
\end{tabular}
\end{table}
```

### 2c. §5.1 guidance (replaces C6)

Replace the C mixing-weight blocks with a short A-vs-B default rule
plus pointer:

> When the operative mechanism is unknown, Architecture B is the
> conservative default for trial-design decisions: its larger
> carryover sensitivity keeps power estimates from being optimistic
> relative to Architecture A. A combined architecture that mixes the
> two channels, and a corresponding sensitivity sweep, is developed
> in the companion paper [CITE paper 11].

## Part 3. Cuts to `report-slim.Rmd`

`report-slim` carries C in: §2.3 (431-506), the §2.4 C paragraph
(544-557), §3.4 (972-1067), the §5.1 3-column table + mixing guidance
(1217-1269), §4-biology C bullet (1181-1204 area), and the §8
reproducibility "not directly comparable" convention note
(1618-1624). Apply the same reductions; keep only the 2b/2c-style
pointer. Line ranges to be re-verified at execution time (slim is
edited less often and may drift).

## Part 4. Paper 11 skeleton

Proposed directory: `analysis/report/11-combined-dgp-architecture/`
Proposed title: "A combined data-generating architecture for
biomarker-treatment interactions: channel super-additivity and
mixing-weight sensitivity in aggregated N-of-1 trials."

Section outline:

1. Introduction. Recap A vs B from paper 01 (one paragraph, cite
   paper 01); motivate the combined architecture; state the two
   contributions (analytical super-additivity; corrected boundary
   validation).
2. The combined architecture. Migrate `report.Rmd` §2.3 (definition,
   sequential MVN-draw-then-mean-shift construction, algebraic
   boundary reductions, independence of `c.bm_a`/`c.bm_b`, the
   alpha-mixing parameterization) and the §2.4 C math paragraph.
3. Analytical results (NEW). Derive the super-additivity of the joint
   BM-BR covariance over the sum of single-channel contributions
   (referee Significance critique + Suggestion 2). This is the
   paper's reason to exist beyond a descriptive grid.
4. Simulation study. Migrate §3.4 grid, Table 1, and super-additivity
   numbers. REQUIRED FIX (MC-1): re-run the 3x3 grid under a driver
   whose Gompertz/biomarker parameters match the paper-01 §3.1
   production driver, so boundary cells reproduce pure A and pure B to
   within ~3 pp; present a boundary-validation table showing C
   boundaries alongside the paper-01 baselines (referee Suggestion 1).
   Add MCSE columns (referee Suggestion 4). Report all three
   $t_{1/2}$ panels, not only $t_{1/2}=1.0$.
5. Biological interpretation. Migrate §4.4 (two mechanistically
   distinct channels; prazosin-SBP both-channel case).
6. Guidance. Migrate §5.1 C guidance (sensitivity-analysis framework;
   C is not a default; needs independent evidence per channel).
7. Discussion. Relation to paper 01 and to the carryover-mitigation
   companion (paper 02).
8. Reproducibility. Migrate the C driver/data documentation removed
   from paper 01 §8 (C8 above): `analysis/scripts/dgp-combined/`,
   `analysis/data/dgp-combined/`, seed `set.seed(20260610)`,
   `furrr::future_map_dfr` with `plan(multicore)`.

New work beyond extraction (the gating items):

- The Section 3 analytical derivation (does not currently exist).
- The MC-1 boundary re-run against a common driver (paper 01 dodged
  this by reframing; paper 11 must actually do it).
- MCSE columns and the full three-panel grid.

## Part 5. Response-to-referee note

The referee reviewed a version containing Architecture C. The
response letter must state that C was MOVED to a dedicated companion
paper (not dropped), that the companion re-runs the boundary cells
against a common driver (resolving MC-1) and adds the analytical
super-additivity derivation (resolving the Significance-section and
Suggestion-2 critiques), and that paper 01 is thereby narrowed to the
A-versus-B contrast it is meant to make.
