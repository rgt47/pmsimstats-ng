# pmsimstats-ng Repository Reorganization Plan
*2026-04-22 11:30 PDT*

**Author:** pmsimstats team
**Target paradigm:** zzcollab unified research compendium
**Scope:** File-system layout only; no API or algorithmic changes.

## 1. Current-state diagnosis

The repository is a hybrid R package and research compendium that has
accumulated rapid growth in four areas without a corresponding
layout update. The principal symptoms are the following.

### 1.1 `docs/` is carrying tracked build artifacts

- 114 files total in `docs/`.
- 9 source Markdown files, 19 source LaTeX files.
- 23 rendered PDFs tracked in git.
- 59 LaTeX auxiliary files tracked in git (`.aux` x 15, `.log` x 16,
  `.out` x 15, `.toc` x 13).
- One file with an embedded space and non-standard name:
  `02-dgp-mean-moderation-vs-mvn rchcomments_rh.pdf`.

LaTeX auxiliary files and rendered PDFs should never appear in
source control in a reproducible-build workflow. They regenerate on
every render, bloat diffs, and create merge-conflict risk that
contributes no informational value.

### 1.2 `analysis/` and `manuscripts/` have overlapping namespaces

Both directories contain a `02-carryover-sensitivity/` subtree. On
inspection these are complementary rather than duplicated:
`analysis/02-carryover-sensitivity/` holds the simulation driver
scripts and a grid plan; `manuscripts/02-carryover-sensitivity/`
holds the accompanying `report.Rmd`, `report.tex`, CSL style, and
rendered PDF. The naming collision nonetheless obscures the
distinction and makes navigation inefficient.

### 1.3 `implementations/` mixes code and infrastructure at its root

The four implementation collections (`original/`,
`original-extended/`, `tidyverse/`, `simple/`) live alongside
loose parity infrastructure: `parity_baseline.rds`,
`parity_report.rds`, `parity_divergences.csv`, `parity-diff.R`,
`test-parity-extended-tidyverse.R`, and
`test-parity-original-tidyverse.R`. These six files constitute a
parity subsystem that has not been promoted to its own directory.

### 1.4 `vignettes/` contains abandoned cache artifacts

Three directories `coreresults_save1/`, `coreresults_save2/`,
`coreresults_save3/` sit under `vignettes/`. They are each roughly
4 KB (empty or near-empty) and appear to be stale vignette cache
output from earlier render attempts. They are not referenced by
any current `.Rmd`.

### 1.5 Repository root has accumulated stray files

- `Rplots.pdf`: a stray R graphics device output from an
  interactive session; not regenerable from any tracked code path.
- `.Rhistory`: 17 KB of interactive R history; should be
  gitignored.
- `GEMINI.md`: a parallel `CLAUDE.md` targeted at a different
  assistant; the two documents partially overlap and will drift
  apart without a consolidation rule.
- `scripts/` contains a single file
  (`scripts/pre-commit-parity.sh`); the directory name is generic.
- `zzcollab.yaml` contains a single setting (`default_profile:
  analysis`) and is not consulted elsewhere in the repo.

### 1.6 Package `data/` holds build outputs as tracked artifacts

`data/` contains seven `.rda` files (`CTdata.rda`,
`extracted_bp.rda`, `extracted_rp.rda`, `results_core.rda`,
`results_maxes.rda`, `results_rates.rda`,
`results_trajectories.rda`). There is no `data-raw/` directory of
scripts that regenerates them, so the binaries are effectively
primary sources rather than derived artifacts. The R-package
convention is to maintain `data-raw/*.R` producer scripts and to
treat `data/*.rda` as regenerable.

## 2. Target layout (zzcollab unified paradigm)

The zzcollab template assumes the Marwick, Boettiger, and Mullen
(2018) research compendium layout, in which a single repository
simultaneously serves as an R package (`R/`, `man/`, `tests/`,
`DESCRIPTION`, `NAMESPACE`) and as a research compendium
(`analysis/`, `manuscripts/`, `docs/`). The proposed pmsimstats-ng
layout preserves the current package bones and reorganises the
research material around them.

```
pmsimstats-ng/
|- R/                          # Package source (unchanged)
|- man/                        # Generated package documentation
|- NAMESPACE                   # Generated; do not edit by hand
|- DESCRIPTION
|- tests/                      # Package unit tests (tinytest)
|- vignettes/                  # Package vignettes, cache-free
|- inst/                       # Installed files
|- data/                       # Package datasets (.rda)
|- data-raw/                   # NEW: scripts that produce data/*.rda
|
|- analysis/                   # Research simulations
|  |- figure4/                 # Publication figure 4 pipeline
|  |- figure5/                 # Publication figure 5 pipeline
|  |- architecture-comparison/
|  |- carryover-factorial/
|  |- carryover-sensitivity/   # Renamed from 02-carryover-sensitivity
|  |- archive/                 # Historical commit snapshots
|  \- README.md                # Navigation and conventions
|
|- manuscripts/                # Manuscripts and supplementary materials
|  |- 01-biomarker-interaction-review/
|  |- 02-carryover-sensitivity/
|  |- data/                    # Manuscript-only data extracts
|  |- figures/                 # Final publication figures
|  |- tables/
|  |- references.bib           # Shared bibliography
|  \- README.md
|
|- docs/                       # Methodological documentation
|  |- source/                  # .md and .tex sources only
|  |  |- 00-documentation-index.tex
|  |  |- 01-codebase-overview.md
|  |  |- ...
|  |  \- 20-...
|  |- rendered/                # Gitignored; build output
|  |- audits/                  # Point-in-time audits
|  |  \- audit-2026-04-22.md
|  \- README.md                # Numbered-index navigation
|
|- implementations/
|  |- original/
|  |- original-extended/
|  |- tidyverse/
|  |- simple/
|  |- parity/                  # NEW: consolidated parity subsystem
|  |  |- artifacts/            # .rds, .csv baselines
|  |  |- tests/                # test-parity-*.R
|  |  |- parity-diff.R
|  |  \- README.md
|  \- README.md
|
|- tools/                      # Renamed from scripts/
|  \- pre-commit-parity.sh
|
|- .github/                    # CI workflows
|- CLAUDE.md                   # Canonical assistant-context file
|- Dockerfile
|- Makefile
|- renv.lock
|- README.md
|- LICENSE
\- .gitignore                  # Expanded: LaTeX artifacts, .Rhistory,
                               # Rplots.pdf, docs/rendered/,
                               # vignettes/*_cache/, *.rds in parity
```

## 3. Mapping from current to target

The moves are grouped by category. Items marked `(verify)` require
human judgement before executing.

### 3.1 Documentation hygiene (`docs/`)

- Delete all tracked LaTeX artifacts (`*.aux`, `*.log`, `*.out`,
  `*.toc`; 59 files) via `git rm` after adding globs to
  `.gitignore`.
- Move the 23 tracked PDFs from `docs/` to `docs/rendered/`, then
  gitignore `docs/rendered/` and `git rm --cached` the originals.
- Move source files `*.md` and `*.tex` into `docs/source/`.
- Move `docs/audit-2026-04-22.md` into `docs/audits/`.
- Merge `docs/.gitignore` into the root `.gitignore`.
- `(verify)` Relocate
  `docs/02-dgp-mean-moderation-vs-mvn rchcomments_rh.pdf` (and any
  other review-annotated PDFs) to `manuscripts/reviews/` rather
  than deleting, since these are not regenerable from source.

### 3.2 Analysis and manuscripts (`analysis/`, `manuscripts/`)

- Rename `analysis/02-carryover-sensitivity/` to
  `analysis/carryover-sensitivity/` to remove the numeric prefix
  collision with the manuscript directory.
- Retain `manuscripts/02-carryover-sensitivity/` as-is.
- Retain `analysis/archive/` as-is.
- `(verify)` Decide the disposition of `analysis/report/`: promote
  contents to `manuscripts/`, archive to `analysis/archive/`, or
  delete.

### 3.3 Parity subsystem (`implementations/`)

Consolidate the six loose parity files into
`implementations/parity/`:

- `parity_baseline.rds`, `parity_report.rds`,
  `parity_divergences.csv` to `implementations/parity/artifacts/`.
- `parity-diff.R` to `implementations/parity/`.
- `test-parity-extended-tidyverse.R` and
  `test-parity-original-tidyverse.R` to
  `implementations/parity/tests/`.
- Add `implementations/parity/README.md` describing the subsystem.

### 3.4 Root cleanup

- Rename `scripts/` to `tools/` (only file inside is
  `pre-commit-parity.sh`).
- Delete `Rplots.pdf` (`git rm`; gitignore).
- Delete `.Rhistory` (`git rm`; gitignore).
- Delete `vignettes/coreresults_save{1,2,3}/` after confirming
  they are not referenced by any vignette.
- `(verify)` Resolve `GEMINI.md`: delete, reduce to a pointer, or
  retain in parallel with `CLAUDE.md`.
- Retain `zzcollab.yaml`; consider extending with profile and
  Docker-tag pinning.

### 3.5 Package data (`data/`, new `data-raw/`)

Additive change: introduce `data-raw/` with producer scripts for
each `data/*.rda`. No `.rda` files are moved or deleted in this
reorganisation; documenting their provenance is scoped to a
follow-up sprint (see Section 7).

Two decisions require explicit input before the moves can be
executed:

1. **`GEMINI.md`**: is this still needed as a separate
   assistant-context file? If not, delete. If yes, reduce to a
   pointer that says 'see CLAUDE.md' to prevent drift.
2. **`analysis/report/`**: the contents
   (`11-gompertz-clinician-guide.Rmd`, `preamble.tex`, `report.Rmd`)
   suggest an abandoned or exploratory manuscript scaffold. Decide
   whether to promote to `manuscripts/` or archive to
   `analysis/archive/`.

## 4. `.gitignore` additions

The following globs should be added to the root `.gitignore` to
prevent future regression:

```
# LaTeX build artifacts
*.aux
*.log
*.out
*.toc
*.synctex.gz
*.fdb_latexmk
*.fls
*.bbl
*.blg

# Rendered documentation (regenerated on demand)
docs/rendered/

# R session artifacts
.Rhistory
.RData
.Ruserdata
Rplots.pdf

# Vignette and knitr caches
vignettes/*_cache/
vignettes/*_files/
vignettes/coreresults_save*/
```

Two notes on PDF handling. First, keeping source `.md` and `.tex`
tracked while excluding rendered `.pdf` makes the render a pure
build step. Second, certain PDFs with reviewer annotations (e.g.,
the `rchcomments_rh.pdf` variants) are not regenerable from the
tracked sources and should be relocated to a reviewed-artifacts
directory (proposed: `manuscripts/reviews/`) rather than simply
deleted.

## 5. Phased execution sequence

The moves decompose into five phases that can be executed
independently with checkpoint commits between each phase.

**Phase 1: LaTeX hygiene (low risk, high visual impact).**
Add gitignore entries for LaTeX artifacts, then
`git rm --cached` the 59 existing artifacts in `docs/`. Verify
that rebuilds still succeed (`make docs` or equivalent). Expected
outcome: 59-file reduction in tracked content with no behavioural
change.

**Phase 2: Rendered-output separation.**
Create `docs/rendered/`; `git mv` the 23 tracked PDFs into it;
update any render scripts (`Makefile` targets, render helpers in
`tools/`) to write into `docs/rendered/`; gitignore
`docs/rendered/`. Verify with a clean render that the rendered
PDFs appear at the new location.

**Phase 3: Docs source relocation.**
Create `docs/source/` and `docs/audits/`. `git mv` all `.md` and
`.tex` sources into `docs/source/`, and the audit file into
`docs/audits/`. Update any cross-references (between numbered
documents, from `CLAUDE.md`, from `README.md`). Update render
scripts. Verify that the rendered PDFs still reproduce.

**Phase 4: Parity subsystem consolidation.**
Create `implementations/parity/artifacts/` and
`implementations/parity/tests/`. `git mv` the six parity files
into their new locations. Update paths in `Makefile`,
`tools/pre-commit-parity.sh`, and any analysis drivers that
reference the artifacts. Add a
`implementations/parity/README.md` describing the subsystem.
Verify that the parity tests still pass.

**Phase 5: Root cleanup and minor renames.**
Delete `Rplots.pdf` and `.Rhistory` (gitignored in Phase 1).
Delete `vignettes/coreresults_save{1,2,3}/` after confirming they
are not referenced. Rename `scripts/` to `tools/`. Rename
`analysis/02-carryover-sensitivity/` to
`analysis/carryover-sensitivity/`. Resolve the `GEMINI.md` and
`analysis/report/` questions per Section 3. Update any path
references that break.

Each phase should end with: (a) a clean `R CMD check`, (b) a
successful `make test`, (c) a successful render of at least one
representative manuscript and one docs source file, and (d) a
single descriptive commit.

## 6. Risks and verification

The primary risks are broken path references. The following
inventory identifies the locations most likely to carry hard-coded
paths that require updating.

- **`Makefile`**: targets that render documents, run parity tests,
  or invoke scripts. Scan for strings matching `docs/`,
  `scripts/`, `implementations/`.
- **Analysis drivers** (`analysis/*/*.R`): `saveRDS` and
  `write.csv` calls that write into `output/` subtrees; usually
  relative so safe, but verify.
- **`tools/pre-commit-parity.sh`**: references the three parity
  artifacts directly.
- **`CLAUDE.md`**: project-level paths in the 'Documentation' and
  'Implementation' sections.
- **`README.md`**: any 'Project Structure' ASCII tree.
- **Cross-document references** in `docs/source/*.md` and
  `*.tex`: numbered cross-references (e.g., `\ref{sec:...}`) that
  point to other numbered documents, and any pandoc `[text](file)`
  links that assume a flat `docs/` layout.
- **CI workflows** in `.github/workflows/`: paths to render
  targets, test invocations.

Recommended verification after each phase:

```bash
Rscript -e 'devtools::check()'
Rscript -e 'devtools::test()'
make docs
make parity
```

A git `post-checkout` sanity script that re-runs `R CMD check`
and the parity tests on the main branch would catch regressions
during the multi-phase transition.

## 7. Out-of-scope items noted for follow-up

The following observations surfaced during the reorganisation
survey but are architectural rather than layout questions and
should be addressed in separate work.

- The four implementation collections remain a medium-term
  consolidation question (see `audit-2026-04-22.md` Phase 5). A
  clean layout does not reduce that consolidation cost but should
  precede it so the consolidation diff is not entangled with
  directory renames.
- The `data/` directory holds effectively primary-source `.rda`
  binaries without producer scripts in `data-raw/`. Documenting
  the provenance of each `.rda` is a reproducibility obligation
  that deserves its own sprint.
- The `zzcollab.yaml` at repository root is effectively vestigial
  (one setting, not consulted). Consider either expanding it to
  record profile pinning and Docker tag pinning, or removing it.

## 8. Decision request

Before executing Phase 1, confirm the following four decisions:

1. Approve the target layout in Section 2.
2. Decide the disposition of `GEMINI.md`: delete, pointer-only, or
   retain in parallel.
3. Decide the disposition of `analysis/report/`: promote to
   `manuscripts/`, move to `analysis/archive/`, or delete.
4. Confirm the reviewed-PDF handling: relocate
   `*rchcomments_rh.pdf` and similar to `manuscripts/reviews/`
   rather than deleting.

With those decisions, Phase 1 is the lowest-risk starting point
(59-file reduction with no path references to update) and can
proceed immediately.

---

*Rendered on 2026-04-22 at 11:30 PDT.*<br>
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/reorg-plan-2026-04-22.md*
