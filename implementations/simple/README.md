# simple: reduced-complexity design sandbox

## Purpose

This collection is a **reduced-complexity sandbox** for exploring
alternative trial designs, simplified data-generating processes,
and teaching/pedagogical demonstrations. It is not a parallel
full implementation of the pmsimstats-ng power-analysis pipeline.

Use this collection when you want to:

- Strip the DGP down to a minimal viable form and ask whether a
  design-level question survives the simplification.
- Prototype a new trial-design variant before investing in a full
  package-structured implementation under `original-extended/`
  or `tidyverse/`.
- Illustrate the mean-moderation (Architecture A) mechanism
  pedagogically, without the layered complexity of the full
  package.

## Contents

- `simulation.R` - Self-contained Architecture A simulation in
  base R. Generates data, fits the analysis model, and summarises
  power without external dependencies on the `R/` package or on
  the other implementation collections. Referenced in doc
  `02-dgp-mean-moderation-vs-mvn.tex` (Section 2.1) as the
  pedagogical Architecture A example.
- `docs/` - Placeholder for design-specific notes generated during
  sandbox exploration.

## Distinction from the other three collections

| Collection          | Style        | DGP architectures | Structure                          |
|---------------------|--------------|-------------------|------------------------------------|
| `original/`         | data.table   | B only            | Package (R/, tests/)               |
| `original-extended/`| data.table   | A + B             | Package (R/, tests/)               |
| `tidyverse/`        | tidyverse    | A + B             | Package (R/, tests/)               |
| `simple/`           | base R       | A only            | Single-file sandbox (no tests)     |

The other three collections are maintained in parallel and
produce comparable power estimates. `simple/` is intentionally
divergent: it prioritises transparency and modifiability over
feature parity, and should not be used for publication-grade
power calculations.

## Workflow for adding new design variants

1. Copy `simulation.R` to `simulation_<variant>.R` within this
   directory.
2. Edit the DGP and design sections to encode the new variant.
3. Keep the downstream analysis model identical to the other
   collections so that results remain comparable.
4. Record the variant's purpose, design choices, and headline
   findings in a `docs/<variant>.md` note.
5. If a variant proves useful, promote it to
   `original-extended/` or `tidyverse/` for a
   package-structured, tested implementation.
