## analysis/scripts/carryover-sensitivity/spec-labels.R
##
## Display names for the analysis-model carryover specifications.
##
## The stored codes are E1/E2/E3/E7/E9 everywhere: in
## simulation-core.R, in the drivers, and in every .rds under
## analysis/data/ and output/. Only the display layer changes here,
## so archived artifacts and the superseded drafts under archive/
## stay readable without regeneration. This follows the same pattern
## as the earlier A1/A2/A3 -> E1/E2/E3 migration (NOTATION.md):
## storage is stable, display remaps.
##
## G1-G5 are a second, compact display layer added 2026-08-18 for
## table/figure headers where the full names are unwieldy (e.g.
## heatmap axis labels); prose spells out "G1 (Unadjusted)" etc. at
## first use per section, then uses the short form.
##
##   E1 -> G1 -> Unadjusted             Sx ~ bm + t + Db  + bm:Db
##   E3 -> G2 -> Lag-adjusted           Sx ~ bm + t + Db  + bm:Db + L
##   E2 -> G3 -> Exposure-weighted      Sx ~ bm + t + Dbc + bm:Dbc
##   E7 -> G4 -> AIC-selected t1/2      E2's formula, half-life chosen
##                                      by AIC over {0.25,0.5,1.0,2.0}
##   E9 -> G5 -> Paired-difference      diff ~ bm (patient-level OLS,
##                                      no repeated-measures model)
##
## Display order for the Tier 1 three is Unadjusted, Lag-adjusted,
## Exposure-weighted rather than E1/E2/E3: the first two share the
## estimand bm:Db and are directly comparable on every performance
## measure, while the third targets bm:Dbc and is comparable only on
## power and Type I error. G4/G5 (E7/E9) are Tier 2 additions (Blocks
## S7-S9) and are not part of the Tier 1 grid; spec_order/spec_labels
## below remain the Tier 1 three only. g_labels/g_order below cover
## all five, for Tier 2 S7-S9 tables and figures.

spec_order <- c('E1', 'E3', 'E2')

spec_labels <- c(E1 = 'Unadjusted',
                 E3 = 'Lag-adjusted',
                 E2 = 'Exposure-weighted')

## Five-specification G-code layer (Tier 2, Blocks S7-S9)
g_order <- c('E1', 'E3', 'E2', 'E7', 'E9')

g_code <- c(E1 = 'G1', E3 = 'G2', E2 = 'G3', E7 = 'G4', E9 = 'G5')

g_labels <- c(E1 = 'G1: Unadjusted',
             E3 = 'G2: Lag-adjusted',
             E2 = 'G3: Exposure-weighted',
             E7 = 'G4: AIC-selected t1/2',
             E9 = 'G5: Paired-difference')

g_labels_short <- c(E1 = 'G1', E3 = 'G2', E2 = 'G3', E7 = 'G4', E9 = 'G5')

g_factor <- function(x) {
  factor(x, levels = g_order, labels = unname(g_labels[g_order]))
}

spec_factor <- function(x) {
  factor(x, levels = spec_order, labels = unname(spec_labels[spec_order]))
}

spec_colours <- c('Unadjusted'        = '#1f78b4',
                  'Lag-adjusted'      = '#e31a1c',
                  'Exposure-weighted' = '#33a02c')

spec_linetypes <- c('Unadjusted'        = 'dotted',
                    'Lag-adjusted'      = 'longdash',
                    'Exposure-weighted' = 'solid')

spec_shapes <- c('Unadjusted'        = 16,
                 'Lag-adjusted'      = 4,
                 'Exposure-weighted' = 17)
