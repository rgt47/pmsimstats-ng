## analysis/scripts/carryover-sensitivity/spec-labels.R
##
## Display names for the three analysis-model carryover specifications.
##
## The stored codes are E1/E2/E3 everywhere: in simulation-core.R, in the
## drivers, and in every .rds under analysis/data/ and output/. Only the
## display layer changes here, so archived artifacts and the superseded
## drafts under archive/ stay readable without regeneration.
##
##   E1 -> Unadjusted         Sx ~ bm + t + Db  + bm:Db
##   E3 -> Lag-adjusted       Sx ~ bm + t + Db  + bm:Db + L
##   E2 -> Exposure-weighted  Sx ~ bm + t + Dbc + bm:Dbc
##
## Display order is Unadjusted, Lag-adjusted, Exposure-weighted rather
## than E1/E2/E3: the first two share the estimand bm:Db and are directly
## comparable on every performance measure, while the third targets
## bm:Dbc and is comparable only on power and Type I error.

spec_order <- c('E1', 'E3', 'E2')

spec_labels <- c(E1 = 'Unadjusted',
                 E3 = 'Lag-adjusted',
                 E2 = 'Exposure-weighted')

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
