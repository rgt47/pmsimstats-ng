# 02-generate-parameter-sensitivity-original.R
#
# Generates Figure 5 equivalent using the ORIGINAL code from
# the initial commit (42ac030) to match the publication.
#
# Panel A: Vary max for BR, ER, TR across {5, 8, 11}
#   Fixed: N=35, c.bm=0.6, carryover=0, no censoring
# Panel B: Vary rate for BR, ER, TR across {0.2, 0.35, 0.5}
#   Fixed: N=35, c.bm=0.3, carryover=0, no censoring
#
# Usage:
#   NREPS=200 Rscript analysis/figure5/02-generate-parameter-sensitivity-original.R

library(data.table)
library(lme4)
library(lmerTest)
library(corpcor)
library(MASS)
library(tictoc)
library(ggplot2)
library(gridExtra)

orig_dir <- file.path("analysis", "scripts", "original_R")
source(file.path(orig_dir, "utilities.R"))
source(file.path(orig_dir, "buildtrialdesign.R"))
source(file.path(orig_dir, "generateData.R"))
source(file.path(orig_dir, "censordata.R"))
source(file.path(orig_dir, "lme_analysis.R"))
source(file.path(orig_dir, "generateSimulatedResults.R"))

# Also need PlotModelingResults from current code for plotting
library(devtools)
suppressMessages(load_all("."))
# But use original generateData etc (already sourced, overrides)

output_dir <- file.path("analysis", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv("NREPS", unset = "200"))
cat("=== Figure 5 (Original Code, 42ac030) ===\n")
cat("Nreps:", Nreps, "\n\n")

# ---------------------------------------------------------------
# Trial designs
# ---------------------------------------------------------------

tdOL <- buildtrialdesign("open label", "OL",
  cumulative(rep(2.5, 8)), paste("OL", 1:8, sep = ""),
  rep(1, 8), list(pathA = rep(1, 8)))
tdOLBDC <- buildtrialdesign("OL+BDC", "OL+BDC",
  c(4, 8, 12, 16, 17, 18, 19, 20),
  c("OL1","OL2","OL3","OL4","BD1","BD2","BD3","BD4"),
  c(1, 1, 1, 1, .5, .5, .5, .5),
  list(pathA = c(1,1,1,1,1,1,0,0), pathB = c(1,1,1,1,1,0,0,0)))
tdCO <- buildtrialdesign("crossover", "CO",
  cumulative(rep(2.5, 8)),
  c(paste("COa", 1:4, sep = ""), paste("COb", 1:4, sep = "")),
  rep(.5, 8),
  list(pathA = c(1,1,1,1,0,0,0,0), pathB = c(0,0,0,0,1,1,1,1)))
tdNof1 <- buildtrialdesign("N-of-1", "N-of-1",
  c(4, 8, 9, 10, 11, 12, 16, 20),
  c("OL1","OL2","BD1","BD2","BD3","BD4","COd","COp"),
  c(1, 1, .5, .5, .5, .5, .5, .5),
  list(pathA = c(1,1,1,1,0,0,1,0), pathB = c(1,1,1,1,0,0,0,1),
       pathC = c(1,1,1,0,0,0,1,0), pathD = c(1,1,1,0,0,0,0,1)))

trialdesigns <- list(OL = tdOL, OLBDC = tdOLBDC,
                     CO = tdCO, Nof1 = tdNof1)

# ---------------------------------------------------------------
# Baseline parameters
# ---------------------------------------------------------------

extracted_bp <- data.table(
  cat = c("BL", "bm"),
  m = c(83.06897, 124.32759),
  sd = c(18.48267, 15.36159))
blparamsets <- list(list(name = "TR",
  verbaldesc = "Extracted from data", param = extracted_bp))

# ---------------------------------------------------------------
# Analysis parameters
# ---------------------------------------------------------------

analysisparams <- expand.grid(
  useDE = FALSE, t_random_slope = FALSE,
  full_model_out = FALSE)

# ---------------------------------------------------------------
# Panel A: Vary maxes (publication params: N=35, c.bm=0.6)
# ---------------------------------------------------------------

cat("--- Panel A: Varying response maxima ---\n\n")

extracted_rp <- data.table(
  cat = c("tv", "pb", "br"),
  max = c(6.50647, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.35, 0.35, 0.42),
  sd = c(10, 10, 8))

maxes <- c(5, 8, 11)
maxes_grid <- expand.grid(tv_max = maxes, pb_max = maxes,
                          br_max = maxes)

respparamsets_maxes <- vector("list", nrow(maxes_grid))
for (i in 1:nrow(maxes_grid)) {
  op <- maxes_grid[i, ]
  respparamsets_maxes[[i]] <- list(
    name = paste0("MX_tv", op$tv_max, "_pb", op$pb_max,
                  "_br", op$br_max),
    verbaldesc = "",
    param = data.table(
      cat = c("tv", "pb", "br"),
      max = as.numeric(op),
      disp = c(5, 5, 5),
      rate = c(0.35, 0.35, 0.42),
      sd = as.numeric(op)))
}

modelparams_a <- expand.grid(
  N = 35, c.bm = 0.6, carryover_t1half = 0,
  c.tv = .8, c.pb = .8, c.br = .8,
  c.cf1t = .2, c.cfct = .1)

cat("Starting Panel A at", format(Sys.time()), "\n")
simresults_maxes <- generateSimulatedResults(
  trialdesigns = list(trialdesigns$OL, trialdesigns$OLBDC,
                      trialdesigns$CO, trialdesigns$Nof1),
  respparamsets = respparamsets_maxes,
  blparamsets = blparamsets,
  censorparams = NA,
  modelparams = modelparams_a,
  simparam = list(Nreps = Nreps, progressiveSave = FALSE,
    basesavename = "fig5a_orig", nRep2save = 5,
    saveunit2start = 1, savedir = output_dir),
  analysisparams = analysisparams,
  rawdataout = FALSE)
cat("Panel A complete at", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------
# Panel B: Vary rates (publication params: N=35, c.bm=0.3)
# ---------------------------------------------------------------

cat("--- Panel B: Varying response rates ---\n\n")

rates <- c(0.2, 0.35, 0.5)
rates_grid <- expand.grid(tv_rate = rates, pb_rate = rates,
                          br_rate = rates)

respparamsets_rates <- vector("list", nrow(rates_grid))
for (i in 1:nrow(rates_grid)) {
  op <- rates_grid[i, ]
  respparamsets_rates[[i]] <- list(
    name = paste0("RT_tv", op$tv_rate, "_pb", op$pb_rate,
                  "_br", op$br_rate),
    verbaldesc = "",
    param = data.table(
      cat = c("tv", "pb", "br"),
      max = extracted_rp$max,
      disp = c(5, 5, 5),
      rate = as.numeric(op),
      sd = extracted_rp$sd))
}

modelparams_b <- expand.grid(
  N = 35, c.bm = 0.3, carryover_t1half = 0,
  c.tv = .8, c.pb = .8, c.br = .8,
  c.cf1t = .2, c.cfct = .1)

cat("Starting Panel B at", format(Sys.time()), "\n")
simresults_rates <- generateSimulatedResults(
  trialdesigns = list(trialdesigns$OL, trialdesigns$OLBDC,
                      trialdesigns$CO, trialdesigns$Nof1),
  respparamsets = respparamsets_rates,
  blparamsets = blparamsets,
  censorparams = NA,
  modelparams = modelparams_b,
  simparam = list(Nreps = Nreps, progressiveSave = FALSE,
    basesavename = "fig5b_orig", nRep2save = 5,
    saveunit2start = 1, savedir = output_dir),
  analysisparams = analysisparams,
  rawdataout = FALSE)
cat("Panel B complete at", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------
# Save
# ---------------------------------------------------------------

saveRDS(simresults_maxes, file.path(output_dir,
  paste0("results_fig5a_orig_n", Nreps, ".rds")))
saveRDS(simresults_rates, file.path(output_dir,
  paste0("results_fig5b_orig_n", Nreps, ".rds")))

# ---------------------------------------------------------------
# Plot Panel A
# ---------------------------------------------------------------

p1 <- PlotModelingResults(
  simresults_maxes,
  param2plot = "power",
  param2vary = c("trialdesign", "tv_max", "pb_max", "br_max"),
  param2hold = data.table(blparamset = 1, censorparamset = 0,
                          carryover_t1half = 0, c.bm = 0.6),
  param2nothold = c("c.cfct", "modelparamset", "respparamset"))

p1 <- p1 +
  ggtitle("A. Effect of response parameter maximum values on power") +
  labs(subtitle = sprintf(
    "Original code (42ac030) | N=35 | c.bm=0.6 | CS rho=0.8 | %d reps",
    Nreps)) +
  theme(plot.title = element_text(hjust = 0, size = 10),
        plot.subtitle = element_text(hjust = 0, size = 7,
                                     color = "grey40"))

# ---------------------------------------------------------------
# Plot Panel B
# ---------------------------------------------------------------

p2 <- PlotModelingResults(
  simresults_rates,
  param2plot = "power",
  param2vary = c("trialdesign", "tv_rate", "pb_rate", "br_rate"),
  param2hold = data.table(blparamset = 1, censorparamset = 0,
                          carryover_t1half = 0, c.bm = 0.3),
  param2nothold = c("c.cfct", "modelparamset", "respparamset"))

p2 <- p2 +
  ggtitle("B. Effect of response parameter rate values on power") +
  labs(subtitle = sprintf(
    "Original code (42ac030) | N=35 | c.bm=0.3 | CS rho=0.8 | %d reps",
    Nreps)) +
  theme(plot.title = element_text(hjust = 0, size = 10),
        plot.subtitle = element_text(hjust = 0, size = 7,
                                     color = "grey40"))

# ---------------------------------------------------------------
# Save figure
# ---------------------------------------------------------------

ml <- arrangeGrob(grobs = list(p1, p2), nrow = 2)
ggsave(file.path(output_dir, "figure5_42ac030.pdf"),
  ml, width = 9, height = 9, dpi = 300)
ggsave(file.path(output_dir, "figure5_42ac030.png"),
  ml, width = 9, height = 9, dpi = 150)

cat("Figure 5 (original) saved\n")
cat("Complete at", format(Sys.time()), "\n")
