# 01-generate-parameter-sensitivity.R
#
# Generates Figure 5 equivalent: effect of response trajectory
# parameters (max and rate) on power across trial designs.
#
# Panel A: Vary max for BR, ER, TR across {5, 8, 11}
# Panel B: Vary rate for BR, ER, TR across {0.2, 0.35, 0.5}
#
# Fixed parameters: N=35, c.bm=0.5, carryover=0, no censoring,
# AR(1) rho=0.7, lambda_cor=auto
#
# Usage:
#   Rscript analysis/figure5/01-generate-parameter-sensitivity.R
#   NREPS=200 Rscript analysis/figure5/01-generate-parameter-sensitivity.R

library(devtools)
load_all(".")
library(gridExtra)

output_dir <- file.path("analysis", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

Nreps <- as.integer(Sys.getenv("NREPS", unset = "200"))
n_cores <- as.integer(Sys.getenv("NCORES", unset = "-1"))
if (n_cores < 0) n_cores <- max(1, parallel::detectCores() - 1)

cat("=== Figure 5: Response Parameter Sensitivity ===\n")
cat("Nreps:", Nreps, "\n")
cat("Cores:", n_cores, "\n\n")

# ---------------------------------------------------------------
# Trial designs (same as Figure 4)
# ---------------------------------------------------------------

tdOL <- buildtrialdesign(
  name_longform = "open label",
  name_shortform = "OL",
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = paste("OL", 1:8, sep = ""),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)

tdOLBDC <- buildtrialdesign(
  name_longform = "open label+blinded discontinuation",
  name_shortform = "OL+BDC",
  timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptname = c("OL1", "OL2", "OL3", "OL4",
                 "BD1", "BD2", "BD3", "BD4"),
  expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

tdCO <- buildtrialdesign(
  name_longform = "traditional crossover",
  name_shortform = "CO",
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = c(paste("COa", 1:4, sep = ""),
                 paste("COb", 1:4, sep = "")),
  expectancies = rep(.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)

tdNof1 <- buildtrialdesign(
  name_longform = "primary N-of-1 design",
  name_shortform = "N-of-1",
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptname = c("OL1", "OL2", "BD1", "BD2",
                 "BD3", "BD4", "COd", "COp"),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )
)

trialdesigns <- list(OL = tdOL, OLBDC = tdOLBDC,
                     CO = tdCO, Nof1 = tdNof1)

# ---------------------------------------------------------------
# Baseline parameters
# ---------------------------------------------------------------

data(extracted_bp)
blparamsets <- list(list(name = "TR",
  verbaldesc = "Extracted from data", param = extracted_bp))

# ---------------------------------------------------------------
# Fixed model parameters
# ---------------------------------------------------------------

modelparams <- expand.grid(
  N = 35,
  c.bm = 0.5,
  carryover_t1half = 0,
  c.tv = .7, c.pb = .7, c.br = .7,
  c.cf1t = .2, c.cfct = .1
)

analysisparams <- expand.grid(
  useDE = FALSE,
  t_random_slope = FALSE,
  full_model_out = FALSE
)

# ---------------------------------------------------------------
# Panel A: Vary maxes
# ---------------------------------------------------------------

cat("--- Panel A: Varying response maxima ---\n\n")

maxes <- c(5, 8, 11)
maxes_grid <- expand.grid(tv_max = maxes, pb_max = maxes,
                          br_max = maxes)

respparamsets_maxes <- vector("list", nrow(maxes_grid))
for (i in 1:nrow(maxes_grid)) {
  op <- maxes_grid[i, ]
  respparamsets_maxes[[i]] <- list(
    name = paste0("MX_tv", op$tv_max,
                  "_pb", op$pb_max,
                  "_br", op$br_max),
    verbaldesc = paste0("max: tv=", op$tv_max,
                        " pb=", op$pb_max,
                        " br=", op$br_max),
    param = data.table(
      cat = c("tv", "pb", "br"),
      max = as.numeric(op),
      disp = c(5, 5, 5),
      rate = c(0.35, 0.35, 0.42),
      sd = as.numeric(op)
    )
  )
}

cat("Response parameter sets:", length(respparamsets_maxes), "\n")
cat("Starting Panel A simulation at", format(Sys.time()), "\n\n")

simresults_maxes <- generateSimulatedResults(
  trialdesigns = list(trialdesigns$OL, trialdesigns$OLBDC,
                      trialdesigns$CO, trialdesigns$Nof1),
  respparamsets = respparamsets_maxes,
  blparamsets = blparamsets,
  censorparams = NA,
  modelparams = modelparams,
  simparam = list(
    Nreps = Nreps,
    progressiveSave = FALSE,
    basesavename = "fig5a",
    nRep2save = 5,
    saveunit2start = 1,
    savedir = output_dir
  ),
  analysisparams = analysisparams,
  rawdataout = FALSE,
  n_cores = n_cores
)

cat("Panel A complete at", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------
# Panel B: Vary rates
# ---------------------------------------------------------------

cat("--- Panel B: Varying response rates ---\n\n")

rates <- c(0.2, 0.35, 0.5)
rates_grid <- expand.grid(tv_rate = rates, pb_rate = rates,
                          br_rate = rates)

data(extracted_rp)
respparamsets_rates <- vector("list", nrow(rates_grid))
for (i in 1:nrow(rates_grid)) {
  op <- rates_grid[i, ]
  respparamsets_rates[[i]] <- list(
    name = paste0("RT_tv", op$tv_rate,
                  "_pb", op$pb_rate,
                  "_br", op$br_rate),
    verbaldesc = paste0("rate: tv=", op$tv_rate,
                        " pb=", op$pb_rate,
                        " br=", op$br_rate),
    param = data.table(
      cat = c("tv", "pb", "br"),
      max = extracted_rp$max,
      disp = c(5, 5, 5),
      rate = as.numeric(op),
      sd = extracted_rp$sd
    )
  )
}

# Use c.bm=0.25 for Panel B (matching publication's 0.3)
modelparams_b <- expand.grid(
  N = 35,
  c.bm = 0.25,
  carryover_t1half = 0,
  c.tv = .7, c.pb = .7, c.br = .7,
  c.cf1t = .2, c.cfct = .1
)

cat("Response parameter sets:", length(respparamsets_rates), "\n")
cat("Starting Panel B simulation at", format(Sys.time()), "\n\n")

simresults_rates <- generateSimulatedResults(
  trialdesigns = list(trialdesigns$OL, trialdesigns$OLBDC,
                      trialdesigns$CO, trialdesigns$Nof1),
  respparamsets = respparamsets_rates,
  blparamsets = blparamsets,
  censorparams = NA,
  modelparams = modelparams_b,
  simparam = list(
    Nreps = Nreps,
    progressiveSave = FALSE,
    basesavename = "fig5b",
    nRep2save = 5,
    saveunit2start = 1,
    savedir = output_dir
  ),
  analysisparams = analysisparams,
  rawdataout = FALSE,
  n_cores = n_cores
)

cat("Panel B complete at", format(Sys.time()), "\n\n")

# ---------------------------------------------------------------
# Save results
# ---------------------------------------------------------------

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

saveRDS(simresults_maxes, file.path(output_dir,
  paste0("results_fig5a_n", Nreps, ".rds")))
saveRDS(simresults_rates, file.path(output_dir,
  paste0("results_fig5b_n", Nreps, ".rds")))

# ---------------------------------------------------------------
# Plot Panel A
# ---------------------------------------------------------------

p1 <- PlotModelingResults(
  simresults_maxes,
  param2plot = "power",
  param2vary = c("trialdesign", "tv_max", "pb_max", "br_max"),
  param2hold = data.table(blparamset = 1, censorparamset = 0,
                          carryover_t1half = 0, c.bm = 0.5),
  param2nothold = c("c.cfct", "modelparamset", "respparamset")
)

provenance_a <- sprintf(
  "N=35 | c.bm=0.5 | AR(1) rho=0.7 | %d reps | %s",
  Nreps, format(Sys.time(), "%Y-%m-%d"))

p1 <- p1 +
  ggtitle("A. Effect of response parameter maximum values on power") +
  labs(subtitle = provenance_a) +
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
                          carryover_t1half = 0, c.bm = 0.25),
  param2nothold = c("c.cfct", "modelparamset", "respparamset")
)

provenance_b <- sprintf(
  "N=35 | c.bm=0.25 | AR(1) rho=0.7 | %d reps | %s",
  Nreps, format(Sys.time(), "%Y-%m-%d"))

p2 <- p2 +
  ggtitle("B. Effect of response parameter rate values on power") +
  labs(subtitle = provenance_b) +
  theme(plot.title = element_text(hjust = 0, size = 10),
        plot.subtitle = element_text(hjust = 0, size = 7,
                                     color = "grey40"))

# ---------------------------------------------------------------
# Compose and save
# ---------------------------------------------------------------

ml <- arrangeGrob(grobs = list(p1, p2), nrow = 2)

ggsave(file.path(output_dir,
  paste0("figure5_", timestamp, ".pdf")),
  ml, width = 9, height = 9, units = "in", dpi = 300)
ggsave(file.path(output_dir,
  paste0("figure5_", timestamp, ".png")),
  ml, width = 9, height = 9, units = "in", dpi = 150)

cat("Figure 5 saved to:", output_dir, "\n")
cat("Complete at", format(Sys.time()), "\n")
