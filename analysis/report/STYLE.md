# Manuscript style: voice and terminology
*2026-10-27 17:32 PDT*

**Author.** pmsimstats team

**Scope.** Every manuscript under `analysis/report/` (papers 01 to
14). Symbols are governed by `NOTATION.md`; this file governs prose.
Where the two differ on a symbol, NOTATION.md wins. Where they
differ on a word, this file wins.

## 1. Voice

- Write short declarative sentences. One idea per sentence. Aim for
  15 to 20 words. Split any sentence over 30 words.
- Use standard American scientific prose. State the fact, the
  method, the result. Do not build up to a point; state it.
- No metaphors, no rhetorical questions, no "the price is", "buys",
  "earns its place", "travels with", "the natural suspect", "holds
  its size", "the cells that have lost their signal". Say what the
  number is and what it means.
- No em dashes. No semicolons joining clauses. No parentheses holding
  a clause; use a new sentence.
- Verbs: "we simulated", "the test rejected", "power fell". Prefer the
  active voice where the actor is the authors or the method.
- Numbers: give the value and its Monte Carlo standard error or the
  replicate count. Do not characterize a number ("small", "negligible")
  without the number beside it.
- Each Results subsection opens with what was done and what was
  found, in that order. No "What fails" / "Why" / "Remedy" headings;
  say "We extended the DGP to ..." and "Under the extension ...".

## 2. Framing of paper 01

Paper 01 describes a data-generating process built on the one
published by Hendrickson et al. (2020). It is not a critique of that
paper. Each stage is an extension or enhancement of the DGP's
capabilities. Write "the original DGP" for the published process and
"the extended DGP" for the result, and describe each change as what
it adds:

- "The original DGP holds the biomarker-response correlation constant
  after first exposure. We extended it so the correlation decays at
  the carryover half-life." Not: "the published coupling rule fails".
- "The original DGP accepts any correlation and repairs a matrix that
  is not positive definite. The extended DGP checks feasibility and
  stops." Not: "silent repair is a defect".
- Report the comparison tables as what the extension changes, with
  the numbers. Do not editorialize about the original authors.
- Cite Hendrickson et al. as the source of the design machinery, the
  component decomposition, and the covariance encoding, once, in the
  Background. After that write "the original DGP".

## 3. Terminology map

Use the term on the right in prose. The term on the left may appear
only as a stored data label or code identifier, in backticks, and
must then be listed in the paper's glossary.

| Avoid | Use |
|---|---|
| coupling, biomarker coupling, $c_{bm}$ coupling | biomarker-response correlation; the correlation between the biomarker and the drug response |
| graded coupling, graded rule | decaying correlation; a correlation that decays after discontinuation at the carryover half-life |
| step coupling, step rule, published coupling rule | the original rule; constant correlation after first exposure |
| coupling gap | on/off correlation difference (the on-drug correlation minus the mean off-drug correlation) |
| ceiling, coupling ceiling, $c_{bm}^*$ | largest feasible correlation; feasibility bound |
| prognostic coupling | prognostic correlation (correlation of the biomarker with a component other than the drug response) |
| covariance binding, covariance moderation, covariance-based encoding | the correlation encoding of the interaction (define once: the interaction is a correlation between the biomarker and the drug response that is present on drug and absent off drug) |
| mean moderation | the mean-shift encoding |
| anchored carryover | carryover that decays from the last on-drug mean |
| the published recursion, compounding carryover | the original carryover rule (define) |
| the published process, RH's process | the original DGP |
| the refined process | the extended DGP (define once; keep "refined" only in code and data labels) |
| RH | Hendrickson et al.; after first use, "the original DGP" or "the original analysis" |
| separable AR(1), Kronecker | keep; define once as "the same first-order autoregressive correlation for every component, with a constant between-component correlation at every lag" |
| expectancy, expectancy response, PB | placebo response; placebo-response component; expectancy is acceptable when defined as the participant's probability of being on drug by design |
| natural course, TV | natural-history component |
| drug response, BR | drug-response component |
| gap-fill study | supplementary study |
| warp-speed bootstrap | single-resample bootstrap |
| strawman, straw-man stress test | stress test |
| calibration ratio $\kappa$ | keep with definition: mean model standard error divided by the empirical standard deviation of the estimates |
| oracle half-life | the true half-life (known to the analyst in that stage) |
| size restoration | the increase in power from correcting a conservative test |
| the four-fold interaction penalty | an interaction test needs about four times the sample of a main-effect test |
| process, machinery, construct, construction | DGP, design, correlation matrix, as appropriate |
| stage closes with a decision | the stage ends with the choice made and the reason |
| prog-PB, prog-TV-const, prog-BL, prog-TV-grow | spell out on first use (a biomarker correlated with the placebo response, with the natural-history component, with baseline severity); keep the label in backticks only where the data use it |
| Hybrid design | keep; define on first use as open-label titration, blinded discontinuation, then a brief crossover |
| OL+BDC, CO | spell out on first use |
| E9, E9b | the paired-difference statistic; define |
| G1 to G9, E0, F1, F1A0B | keep as analysis codes; the table defining them is the glossary entry |

## 4. Glossary

Every manuscript carries a "Glossary" section before the
Abbreviations. It defines, in one sentence each, every term from the
"Avoid" column that still appears in a table or data label, every
analysis code, every design name, and every project-specific quantity
($c_{bm}$, the on/off correlation difference, the feasibility bound,
$\kappa$, the half-life policy value). The glossary is the only place
a project term may be explained; the prose uses the standard term.

## 5. Checks before rendering

- `grep -n -i "coupling\|ceiling\|graded\|anchored\|gap-fill\|warp\|strawman\|RH's\|refined process\|published process"` on the Rmd. Each hit must be in a code chunk, a backticked label, or the glossary.
- `perl tools/notation-lint.pl` clean.
- No sentence over 30 words in the Abstract.
