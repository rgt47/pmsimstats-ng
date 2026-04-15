## analysis/02-carryover-sensitivity/03-render-figures.R
##
## Render manuscript-02 figures from the summarised grid. Outputs
## PDF files into manuscripts/figures/.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

repo_root <- here::here()

grid <- readRDS(file.path(repo_root,
  'manuscripts/data/02-grid-summary.rds'))$summary

fig_dir <- file.path(repo_root, 'manuscripts/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Figure 02-heatmap-matched-vs-mismatched.pdf
## Rows: DGP decay form. Columns: analysis specification. Cells:
## relative power vs the matched cell for that DGP row. One facet
## per (architecture, design, t1half).
##
## TODO: implement after simulation grid populated. Placeholder
## below writes a diagnostic empty plot so the pipeline is wired.

ggsave(
  filename = file.path(fig_dir, '02-heatmap-matched-vs-mismatched.pdf'),
  plot = ggplot() +
    annotate('text', x = 0.5, y = 0.5,
             label = 'Pending simulation run.') +
    theme_void(),
  width = 6, height = 4
)

message('Wrote placeholder figure to ', fig_dir)
