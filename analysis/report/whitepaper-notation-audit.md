---
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.4em}
---

# White paper: notation and terminology across the compendium

*2026-07-30 20:10 PDT (second edition; supersedes the 2026-07-29
diagnosis)*

**A re-audit of mathematical notation and controlled vocabulary in the
manuscripts under `analysis/report/`, after two repair passes. Records
what was fixed, what the repairs themselves broke, and what remains.**

## Method

Sixteen manuscript sources were searched for the symbols carrying the
programme's substantive quantities and for the labels naming its
designs, architectures, specifications, and studies. The checks are now
mechanized as a linter (`notation-lint.pl`) that skips fenced R chunks
and `verbatim` blocks, strips `\texttt{}` and backtick spans, and
reports violations of the rules in `analysis/report/NOTATION.md`. Every
manuscript was re-rendered after each change; counts below are from the
sources and from `pdftotext` on the rendered output.

The corpus is smaller than in the first edition. Papers 01 and 02 have
since been reorganized: their narrative variants were promoted to
`report.Rmd` and the earlier full and slim versions moved to
per-paper `archive/` directories. This audit reflects the promoted
files.

## State after the second pass

The linter reports **clean** for all sixteen manuscripts except one
known false positive (a bare `Sx` inside a backtick span that crosses a
line break). All sixteen render. The nine inconsistencies of the first
edition now stand as follows.

| # | Finding (first edition) | State |
|---|---|---|
| N1 | Drug indicator: $D_{it}$ used for both the binary state and the continuous exposure | Fixed |
| N2 | $c_{bm}$ used for both a correlation and a regression multiplier | Fixed |
| N3 | Half-life quoted in weeks and in days with no conversion | Fixed by units paragraphs |
| N4 | Four symbols for one estimand | Standardized on $\beta_{bm:D}$ |
| N5 | Outcome alternates between $Y$ and the code identifier `Sx` | Fixed in model definitions |
| N6 | Label collisions across architectures, specifications, sweeps, studies | Fixed, after one false start |
| N7 | Two sample-size conventions (per path, total) | Documented, not resolved |
| N8 | Design labels alternate abbreviation and prose | Open (editorial) |
| N9 | Sign convention stated only in one paper | Open (editorial) |

## What the second pass found that the first pass missed

**The specification relabel was incomplete.** The first pass renamed
the three carryover specifications A1/A2/A3 to M1/M2/M3 in the
carryover manuscript and its variants, but four other papers cite those
specifications and were not touched. For a day the compendium said M2
in one paper and A2 in four others, which is worse than the original
state, because a reader cross-referencing them would conclude they were
different objects.

**The chosen replacement collided with an existing scheme.** The
calibration manuscript already uses M0 through M3 for its
working-covariance model ladder, and those labels are keys into
archived simulation output, not free text. M1/M2/M3 for specifications
therefore reintroduced exactly the class of collision the pass existed
to remove. The specifications are now **E1/E2/E3**, after E for the
exposure regressor that distinguishes them; E is unused anywhere else
in the corpus.

**The architectures were the deeper problem.** A, B, and C are opaque
(a reader must memorize which is which) and they are the reason the
specification labels collided in the first place. They are now named by
the channel that carries the interaction:

| Prose | Legends and tables | Data value |
|---|---|---|
| the mean-moderation architecture | Mean | `mean_moderation` |
| the covariance-moderation architecture | Covariance | `mvn` |
| the dual-channel architecture | Dual | `combined` |

This touched 581 label occurrences across 16 files plus the figure
drivers. Data values are unchanged, so no simulation was re-run; the
figures were regenerated from existing summaries.

**A second British-spelling family survived the first pass**, because
the first pass enumerated word stems rather than patterns: `flavour`,
`rigour`, `endeavour`, `labour`, `tumour`. Fixed in four files.

## What the repairs themselves broke

Reporting this is the point of a second audit. Three defects were
introduced by the first pass and caught here.

**Data-side labels were renamed along with display labels.** The
specification rename reached inside R chunks, changing
`method_levels <- c('A1', ...)` and three `filter(spec %in% c('A2',
'A3'))` calls. Factor levels and filters index stored data values;
renaming them silently produced empty selections and `NA` rows rather
than an error. The first repair was a display-time mapping that kept
the stored values at A1/A2/A3, which is where the third pass then
found the deeper problem (below).

**A display helper returned `NA` for unmapped inputs.** The first
version of `spec_display()` indexed a named vector, so the two
robust-inference methods (`lme+CR2`, `GEE+MD`) that are not
specifications came back `NA` and printed as `NA` in the
development-review tables. Now passes unmapped values through.

**Three manuscript files were truncated to zero bytes.** A shell loop
copied each source to a temporary file, transformed it, and moved the
result back. When the cloud file provider returned a spurious
`ENOENT` for a source that had just been moved by a concurrent
reorganization, the `cp` failed but the shell had already created an
empty output file, which `mv` then placed over the source. The files
were restored from scratchpad backups. Every subsequent script
verifies that the source is non-empty before producing output and that
the result has not shrunk before writing back. This is a sharper
version of the hazard the project guardrails already warn about: the
danger is not only in-place editors, but any pipeline whose failure
mode is an empty file.

## Third pass: retiring the mapping layer

The second pass left the specification labels stored as A1/A2/A3 and
displayed as E1/E2/E3, bridged by a `spec_display()` helper applied at
each display site. That is precisely the arrangement that had already
produced two of the three defects above, and it is fragile for a
structural reason: a mapping has to be applied everywhere, it fails
silently when a rename reaches one layer and not the other, and it
leaves a stored table that cannot be read without the manuscript
beside it.

The mapping has therefore been removed rather than documented. The
`spec` values were migrated in place: 21 `.rds` files recoded to
E1/E2/E3, with the row count and the absence of legacy values verified
on read-back for each, and pre-migration copies kept in
`analysis/.spec-migration-backup/`. The producing code was updated in
the same direction, in both the quoted literals (`'A2'`) and the
unquoted named-vector keys (`c(A2 = ...)`), which are easy to miss
because they index the same values through different syntax. The
display helpers were then deleted from the three manuscripts that
carried them.

Correctness was checked by re-running every affected figure driver and
manuscript. The development-review document verifies its own headline
ranking inline, and reports the same numbers after migration as before
(`E2=0.873, E1=0.727 -> PASS`), which is the evidence that the recode
moved labels and not values.

One mapping is deliberately retained. The `dgp_arch` values
(`mvn`, `mean_moderation`, `combined`) are package API arguments, not
free labels, so the manuscripts map them to Mean, Covariance and Dual
in the ordinary factor-level-to-label way. The distinction worth
holding onto is that a label the project owns should be renamed at
source, while a label an interface owns should be mapped at the edge.

## What remains open

**N7, the sample-size convention.** One paper reports $N$ per
randomization path and others report it in total. The combined-DGP
manuscript documents its own instance and defers to a pending
boundary re-run. This needs a simulation, not an edit.

**N8 and N9** are editorial: design labels still alternate between
abbreviation and spelled-out form, and the sign convention is stated
in only one paper. Both are best folded into the `rgt` writing pass
rather than done mechanically now.

**The `rgt` prose remains unwritten** in every paper that still carries
the three-part scaffold, so a third terminology pass will be needed
once the author's text exists. That is the argument for keeping the
linter: it converts these two one-off passes into a check that can run
on every render.

**A caution on scope.** This audit verifies internal consistency of
symbols and labels. It does not verify that any numerical result is
correct, and the label changes were deliberately confined to display
so that no stored result was touched. Two papers were promoted to
narrative variants during the pass; their archived predecessors were
not re-audited and will not match this scheme if they are ever
revived.

## The reference itself

The complete inventory follows as an appendix: mathematical symbols,
the code identifiers and stored data values they correspond to, and
the glossary of controlled vocabulary. It is not reproduced here but
bound in verbatim from `analysis/report/NOTATION.md`, which is the
canonical file the linter checks against. The two therefore cannot
drift apart: this PDF is rendered from the audit and that file
together.

\newpage
