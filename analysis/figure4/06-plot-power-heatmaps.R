# 06-plot-power-heatmaps.R
#
# Generates the Figure 4 power heatmaps from Hendrickson et al.
# (2020). Can use either:
#   (a) the bundled results_core.rda (default), or
#   (b) a freshly computed .rds file from 01-generate-power-sweep.R
#
# Usage:
#   Rscript analysis/figure4/06-plot-power-heatmaps.R
#   RESULTS_FILE=analysis/output/results_core_YYYYMMDD_final.rds \
#     Rscript analysis/figure4/06-plot-power-heatmaps.R

library(devtools)
load_all(".")
library(gridExtra)

output_dir <- file.path("analysis", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------
# 1. Load results
# ---------------------------------------------------------------

results_file <- Sys.getenv("RESULTS_FILE", unset = "")

if (nzchar(results_file)) {
  cat("Loading results from:", results_file, "\n")
  simresults <- readRDS(results_file)
  source_label <- basename(results_file)
} else {
  cat("Loading bundled results_core.rda\n")
  data(results_core)
  simresults <- results_core
  source_label <- "results_core.rda (bundled)"
}

n_cells <- nrow(unique(simresults$results[,
  .(trialdesign, N, c.bm, carryover_t1half, censorparamset)]))
n_reps <- nrow(simresults$results) / n_cells
n_designs <- length(simresults$parameterselections$trialdesigns)
design_names <- paste(
  sapply(simresults$parameterselections$trialdesigns,
         function(x) x$metadata$name_shortform),
  collapse = ", ")
n_values <- paste(sort(unique(simresults$results$N)),
                  collapse = "/")
cbm_values <- paste(sort(unique(simresults$results$c.bm)),
                    collapse = ", ")
co_values <- paste(sort(unique(
  simresults$results$carryover_t1half)), collapse = ", ")

cat("Results dimensions:", dim(simresults$results), "\n")
cat("Trial designs:", design_names, "\n")
cat("Reps per cell:", n_reps, "\n\n")

provenance <- sprintf(
  paste0("Source: %s | N = %s | c.bm = {%s} | ",
         "carryover = {%s} wk | %d reps | %s"),
  source_label, n_values, cbm_values, co_values,
  n_reps, format(Sys.time(), "%Y-%m-%d")
)

# ---------------------------------------------------------------
# 2. Figure 4A: Power by design x N x c.bm x censoring
#    (carryover = 0)
# ---------------------------------------------------------------

p1 <- PlotModelingResults(
  simresults,
  param2plot = "power",
  param2vary = c("trialdesign", "N", "c.bm",
                 "censorparamset"),
  param2hold = data.table(blparamset = 1,
                          respparamset = 1,
                          carryover_t1half = 0),
  param2nothold = c("c.cfct", "modelparamset")
)

p1 <- p1 +
  ggtitle(paste("A. Effect of trial design, censoring",
                "and biomarker effect on power")) +
  labs(subtitle = provenance) +
  theme(plot.title = element_text(hjust = 0, size = 10),
        plot.subtitle = element_text(hjust = 0, size = 7,
                                     color = "grey40"))

# ---------------------------------------------------------------
# 3. Figure 4B: Power by design x N x carryover x censoring
#    (c.bm = 0.6)
# ---------------------------------------------------------------

p2 <- PlotModelingResults(
  simresults,
  param2plot = "power",
  param2vary = c("trialdesign", "N",
                 "carryover_t1half",
                 "censorparamset"),
  param2hold = data.table(blparamset = 1,
                          respparamset = 1,
                          c.bm = max(simresults$results$c.bm)),
  param2nothold = c("c.cfct", "modelparamset")
)

p2 <- p2 +
  ggtitle(paste("B. Effect of nonzero carryover term",
                "on power")) +
  labs(subtitle = provenance) +
  theme(plot.title = element_text(hjust = 0, size = 10),
        plot.subtitle = element_text(hjust = 0, size = 7,
                                     color = "grey40"))

# ---------------------------------------------------------------
# 4. Compose and save
# ---------------------------------------------------------------

ml <- arrangeGrob(grobs = list(p1, p2), nrow = 2)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

pdf_path <- file.path(output_dir,
                      paste0("figure4_", timestamp, ".pdf"))
ggsave(pdf_path, ml, width = 9, height = 9,
       units = "in", dpi = 300)
cat("Saved:", pdf_path, "\n")

png_path <- file.path(output_dir,
                      paste0("figure4_", timestamp, ".png"))
ggsave(png_path, ml, width = 9, height = 9,
       units = "in", dpi = 150)
cat("Saved:", png_path, "\n")
