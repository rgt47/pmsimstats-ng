## analysis/scripts/dgp-combined/02-summarise-combined.R
##
## Loads the Architecture C replicate output and produces tidy
## summary tables matching the structure of Section 3 in paper 01.
##
## Inputs:
##   analysis/data/dgp-combined/01-combined-replicates.rds
##
## Outputs (written to analysis/data/dgp-combined/):
##   02-power-grid-co.csv      -- CO design, A2 spec
##   02-power-grid-hybrid.csv  -- Hybrid design, A2 spec
##   02-power-grid-olbdc.csv   -- OL+BDC design, A2 spec
##   02-power-wide.csv         -- all designs, wide 3x3 panel

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

repo_root <- here::here()
out_dir   <- file.path(repo_root, 'analysis/data/dgp-combined')

reps <- readRDS(file.path(out_dir, '01-combined-replicates.rds'))

## -----------------------------------------------------------------
## Cell-level summaries
## -----------------------------------------------------------------

cell_summary <- reps |>
  group_by(dgp_arch, t1half, design, N, c_bm_a, c_bm_b,
           carryover_form, spec) |>
  summarise(
    n_reps     = n(),
    n_conv     = sum(converged, na.rm = TRUE),
    power      = mean(p_value[converged] < 0.05, na.rm = TRUE),
    mcse_power = sqrt(power * (1 - power) / max(n_conv, 1L)),
    mean_beta  = mean(estimate[converged], na.rm = TRUE),
    sd_beta    = sd(estimate[converged], na.rm = TRUE),
    conv_rate  = n_conv / n(),
    .groups = 'drop'
  )

## -----------------------------------------------------------------
## 3x3 power panels -- A2 spec, N = 35, t1half in {0, 0.5, 1.0}
## -----------------------------------------------------------------

make_panel <- function(design_name) {
  cell_summary |>
    dplyr::filter(spec == 'A2', N == 35L, design == design_name) |>
    dplyr::select(t1half, c_bm_a, c_bm_b, power, mcse_power) |>
    dplyr::arrange(t1half, c_bm_b, c_bm_a)
}

panel_co     <- make_panel('CO')
panel_hybrid <- make_panel('Hybrid')
panel_olbdc  <- make_panel('OLBDC')

write.csv(panel_co,    file.path(out_dir, '02-power-grid-co.csv'),
          row.names = FALSE)
write.csv(panel_hybrid, file.path(out_dir, '02-power-grid-hybrid.csv'),
          row.names = FALSE)
write.csv(panel_olbdc, file.path(out_dir, '02-power-grid-olbdc.csv'),
          row.names = FALSE)

## -----------------------------------------------------------------
## Wide format for the paper's heatmap figure (Section 3 extended)
## Power at t1half = 1.0, N = 35, A2 spec, all designs
## c_bm_a as rows, c_bm_b as columns
## -----------------------------------------------------------------

wide_panel <- cell_summary |>
  dplyr::filter(spec == 'A2', N == 35L, t1half == 1.0) |>
  dplyr::select(design, c_bm_a, c_bm_b, power) |>
  tidyr::pivot_wider(names_from = c_bm_b,
                     names_prefix = 'cbmb_',
                     values_from = power) |>
  dplyr::arrange(design, c_bm_a)

write.csv(wide_panel, file.path(out_dir, '02-power-wide.csv'),
          row.names = FALSE)

## -----------------------------------------------------------------
## Validation gate: boundary cells vs existing Architecture A/B
## -----------------------------------------------------------------

cat('\n--- Validation gate (boundary cells vs pure architectures) ---\n')

boundary <- cell_summary |>
  dplyr::filter(spec == 'A2', N == 35L,
                (c_bm_a == 0.45 & c_bm_b == 0.0) |
                (c_bm_a == 0.0  & c_bm_b == 0.45)) |>
  dplyr::mutate(
    arch_equiv = dplyr::case_when(
      c_bm_a == 0.45 & c_bm_b == 0.0 ~ 'pure_A',
      c_bm_a == 0.0  & c_bm_b == 0.45 ~ 'pure_B'
    )
  ) |>
  dplyr::select(design, t1half, arch_equiv, power, mcse_power)

print(boundary, n = Inf)

## -----------------------------------------------------------------
## Print full 3x3 panels to console
## -----------------------------------------------------------------

for (dn in c('CO', 'Hybrid', 'OLBDC')) {
  cat(sprintf('\n3x3 power panel -- %s, A2 spec, N=35:\n', dn))
  make_panel(dn) |>
    tidyr::pivot_wider(names_from = c_bm_b,
                       names_prefix = 'cbmb_',
                       values_from = c(power, mcse_power),
                       names_sep   = '_') |>
    print(n = Inf)
}

cat('\nOutputs written to', out_dir, '\n')
invisible(NULL)
