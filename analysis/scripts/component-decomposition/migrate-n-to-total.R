## Migrate the stored `N` column of paper 06's simulation output from a
## per-path count to the total across randomization paths
## (analysis/report/NOTATION.md rule 4).
##
## Nothing about the generated data changes. These drivers passed the
## full `N` to every path, so a cell recorded as N = 150 on the
## two-path full design analyzed 300 participants. The equivalence gate
## at analysis/scripts/component-decomposition/ (see the audit trail in
## NOTATION.md) confirms that re-expressing the grid in totals and
## allocating across paths reproduces the fitted coefficients exactly,
## so the stored numbers are correct and only the label is wrong.
##
## Per-directory rules, because the studies differ in how many paths
## they run:
##   study-a-*      design is 'full' (2 paths) when the per-path count
##                  was >= 100, else 'simple' (1 path). So 100 -> 200
##                  and 150 -> 300; 35 and 70 are single-path and
##                  already totals.
##   study-contam-* same rule: 150 -> 300; 70 unchanged.
##   study-b-recovery-pilot
##                  design_hybrid_full() unconditionally, so every cell
##                  is two-path: 70 -> 140 and 150 -> 300.
##   study-b-balanced
##                  balanced-placebo design has a single path, so its
##                  N is already a total. Not migrated.
##
## Pre-migration copies go to analysis/.n-migration-backup/. Run from
## the repository root. Idempotent: a file whose N values are already
## totals is left alone.

suppressPackageStartupMessages(library(data.table))

DATA_ROOT <- file.path('analysis', 'data', 'derived_data',
                       'component-decomposition')
BACKUP_ROOT <- file.path('analysis', '.n-migration-backup',
                         'component-decomposition')

## Each rule maps an old per-path value to the corrected total.
RULES <- list(
  `study-a-alt1000-null5000`     = c(`100` = 200L, `150` = 300L),
  `study-a-N250`                 = c(`100` = 200L, `150` = 300L),
  `study-contam-alt100-null100`  = c(`150` = 300L),
  `study-contam-alt1000-null5000` = c(`150` = 300L),
  `study-b-recovery-pilot`       = c(`70` = 140L, `150` = 300L)
)
SKIP <- c('study-b-balanced')

migrate_file <- function(path, rule, backup_dir) {
  x <- readRDS(path)
  if (!('N' %in% names(x))) return(list(status = 'no-N'))
  before_n <- nrow(x)
  before_vals <- sort(unique(x$N))
  legacy <- as.integer(names(rule))
  if (!any(before_vals %in% legacy))
    return(list(status = 'already-total', vals = before_vals))

  dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(path, file.path(backup_dir, basename(path)),
            overwrite = FALSE)

  was_dt <- is.data.table(x)
  x <- as.data.table(x)
  for (old in names(rule)) x[N == as.integer(old), N := rule[[old]]]
  if (!was_dt) x <- as.data.frame(x)
  saveRDS(x, path)

  ## Read back and verify: row count preserved, no legacy value left,
  ## and every surviving value is one we expect.
  y <- readRDS(path)
  ok_rows <- nrow(y) == before_n
  ok_clean <- !any(y$N %in% legacy)
  expected <- sort(unique(c(setdiff(before_vals, legacy),
                            unname(rule[as.character(
                              intersect(before_vals, legacy))]))))
  ok_vals <- identical(sort(unique(as.integer(y$N))),
                       as.integer(expected))
  list(status = if (ok_rows && ok_clean && ok_vals) 'migrated'
                else 'FAILED',
       before = before_vals, after = sort(unique(y$N)),
       rows = before_n, ok_rows = ok_rows, ok_clean = ok_clean,
       ok_vals = ok_vals)
}

dirs <- list.dirs(DATA_ROOT, recursive = FALSE)
tally <- list(migrated = 0L, already = 0L, skipped = 0L, failed = 0L,
              no_n = 0L)

for (d in dirs) {
  nm <- basename(d)
  if (nm %in% SKIP) {
    cat(sprintf('[skip]    %-32s single-path design; N already total\n',
                nm))
    tally$skipped <- tally$skipped + 1L
    next
  }
  rule <- RULES[[nm]]
  if (is.null(rule)) {
    cat(sprintf('[skip]    %-32s no rule defined; left untouched\n', nm))
    tally$skipped <- tally$skipped + 1L
    next
  }
  files <- list.files(d, pattern = '\\.rds$', full.names = TRUE)
  cat(sprintf('\n== %s (%d files)\n', nm, length(files)))
  for (f in files) {
    r <- migrate_file(f, rule, file.path(BACKUP_ROOT, nm))
    if (r$status == 'migrated') {
      tally$migrated <- tally$migrated + 1L
    } else if (r$status == 'already-total') {
      tally$already <- tally$already + 1L
    } else if (r$status == 'no-N') {
      tally$no_n <- tally$no_n + 1L
    } else {
      tally$failed <- tally$failed + 1L
      cat(sprintf('  FAILED %s rows_ok=%s clean=%s vals_ok=%s\n',
                  basename(f), r$ok_rows, r$ok_clean, r$ok_vals))
    }
  }
  ## Report the summary file's mapping explicitly, as the visible check.
  sfile <- files[grepl('summary', basename(files))]
  for (s in sfile) {
    y <- readRDS(s)
    cat(sprintf('  %-38s N now: %s\n', basename(s),
                paste(sort(unique(y$N)), collapse = ', ')))
  }
}

cat(sprintf(paste0(
  '\nMigration complete: %d files recoded, %d already totals, ',
  '%d without an N column,\n%d directories skipped, %d FAILURES.\n'),
  tally$migrated, tally$already, tally$no_n, tally$skipped,
  tally$failed))
if (tally$failed > 0L) quit(status = 1L)
