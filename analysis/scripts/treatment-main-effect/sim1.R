# ------------------------------------------------------------------------------
# Simulation Study: Impact of Carryover on Treatment Effect Estimates in
# Multiple Trial Designs
#
# Purpose:
#   To evaluate the bias and loss of power caused by carryover effects
#   when they are present but not adjusted for in the analysis.
#   Modular design allows comparison of different trial structures.
#
# Key Elements:
#   - Simulated data includes drug effects, carryover contamination, and noise
#   - No adjustment is made for carryover in the analysis
#   - Power is estimated by repeating simulations across a range of carryover values
#
# Dependencies: dplyr, tidyr, stringr, progress
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)

# ------------------------
# Global Parameters
# ------------------------

set.seed(2025)
n_subj <- 40              # number of subjects
base_treat_effect <-.50  # true treatment effect
sigma <- 1.0              # residual standard deviation
nsim <- 100              # simulations per delta level
alpha <- 0.05             # test significance threshold

# ------------------------
# Trial Design Definitions
# ------------------------

# Design 1: Hybrid N-of-1 Design
# - Period 1: Open-label (all on drug)
# - Periods 2-3: Blinded discontinuation (everyone reaches placebo by period 3)
# - Periods 4-5: Brief crossover
hybrid_nof1 <- list(
  name = "Hybrid N-of-1",
  sequences = list(
    S1 = c("D", "D", "P", "D", "P"),  # Late discontinuation
    S2 = c("D", "D", "P", "P", "D"),  # Late discontinuation
    S3 = c("D", "P", "P", "D", "P"),  # Early discontinuation
    S4 = c("D", "P", "P", "P", "D")   # Early discontinuation
  )
)

# Design 2: Standard Crossover Design
# - Periods 1-2: First treatment phase
# - Periods 3-5: Second treatment phase (after crossover)
standard_crossover <- list(
  name = "Standard Crossover",
  sequences = list(
    C1 = c("D", "D", "P", "P", "P"),  # Drug first, then placebo
    C2 = c("P", "P", "D", "D", "D")   # Placebo first, then drug
  )
)

# List of all designs for easy iteration
trial_designs <- list(
  hybrid_nof1 = hybrid_nof1,
  standard_crossover = standard_crossover
)

# ------------------------
# Analysis Methods
# ------------------------

# Unadjusted contrast analysis (ignores carryover)
analyze_unadjusted <- function(data) {
  contrast_df <- data %>%
    group_by(Subject) %>%
    summarise(
      mean_D = mean(Y[Treatment == "D"]),
      mean_P = mean(Y[Treatment == "P"]),
      contrast = mean_D - mean_P,
      .groups = "drop"
    )
  t.test(contrast_df$contrast, mu = 0)$p.value
}

# Template for additional analysis methods:
# analyze_new_method <- function(data) {
#   # Your analysis code here
#   # Return p-value
# }

# List of available analysis methods
analysis_methods <- list(
  "Unadjusted" = analyze_unadjusted
  # "New Method" = analyze_new_method  # Add new methods here
)

# ------------------------
# Simulation Functions
# ------------------------

# Generate simulated data for one trial
generate_trial_data <- function(delta_true, design) {
  tibble(Subject = 1:n_subj) %>%
    mutate(Sequence = sample(names(design$sequences), n(), replace = TRUE)) %>%
    rowwise() %>%
    mutate(TreatmentSeq = list(design$sequences[[Sequence]])) %>%
    unnest_longer(TreatmentSeq, values_to = "Treatment") %>%
    group_by(Subject) %>%
    mutate(
      Period = row_number(),
      Treatment = factor(Treatment, levels = c("D", "P")),
      Carryover = Treatment == "P" & lag(Treatment) == "D",
      Carryover = replace_na(Carryover, FALSE),
      mu = case_when(
        Treatment == "D" ~ base_treat_effect,
        Treatment == "P" & Carryover == TRUE ~ delta_true,
        TRUE ~ 0
      ),
      Y = mu + rnorm(n(), 0, sigma)
    ) %>%
    ungroup()
}

# Run simulation for specific delta, design, and analysis method
run_sim <- function(delta_true, design, analysis_method) {
  p_values <- numeric(nsim)

  for (i in 1:nsim) {
    trial_data <- generate_trial_data(delta_true, design)
    p_values[i] <- analysis_method(trial_data)
  }

  mean(p_values < alpha)
}

# ------------------------
# Run Power Analysis Across Carryover Values for All Designs and Methods
# ------------------------

delta_grid <- seq(0, 1, by = 0.1)

# Run simulations for all combinations of designs and analysis methods
results_list <- list()

for (design_name in names(trial_designs)) {
  design <- trial_designs[[design_name]]
  
  for (method_name in names(analysis_methods)) {
    method_func <- analysis_methods[[method_name]]
    
    cat("Running simulations for", design$name, "with", method_name, "analysis...\n")
    
    power_results <- sapply(delta_grid, function(delta) {
      cat(sprintf("  %s (%s): delta=%s\n", design$name, method_name, delta))
      run_sim(delta, design, method_func)
    })
    
    combo_name <- paste0(design_name, "_", gsub(" ", "_", method_name))
    results_list[[combo_name]] <- tibble(
      Design = design$name,
      Analysis = method_name,
      Carryover_Effect = delta_grid,
      Power = power_results
    )
    
    cat("\n")  # Add newline after progress bar
  }
}

# Combine results
all_results <- bind_rows(results_list)

# Create side-by-side comparison table
comparison_table <- all_results %>%
  select(Design, Carryover_Effect, Power) %>%
  pivot_wider(names_from = Design, values_from = Power)

# Print side-by-side results
cat("\n", "Power Comparison Across Designs (Unadjusted Analysis):\n")
cat("======================================================\n")
print(comparison_table, digits = 3)

# Print individual design results for reference
cat("\n", "Individual Design Results:\n")
cat("==========================\n")
for (design_name in names(trial_designs)) {
  cat("\n", trial_designs[[design_name]]$name, "Results:\n")
  print(filter(all_results, Design == trial_designs[[design_name]]$name) %>% 
        select(Carryover_Effect, Power), digits = 3)
}

# ------------------------
# Optional: Plot Results
# ------------------------

# Plot comparing all designs
designs_in_results <- unique(all_results$Design)
cols <- hcl.colors(length(designs_in_results), palette = "Set 2")

plot(NULL, xlim = range(delta_grid), ylim = c(0, 1),
     xlab = "Carryover Effect (delta_true)",
     ylab = "Estimated Power",
     main = "Power vs. Carryover Effect: Design Comparison (Unadjusted Analysis)")
for (j in seq_along(designs_in_results)) {
  d <- all_results[all_results$Design == designs_in_results[j], ]
  lines(d$Carryover_Effect, d$Power, col = cols[j], lwd = 2)
  points(d$Carryover_Effect, d$Power, col = cols[j], pch = 16)
}
abline(h = 0.8, lty = 2, col = "red")
legend("topright", legend = designs_in_results, col = cols, lwd = 2, pch = 16)

# Individual plots for each design
for (design_name in names(trial_designs)) {
  design_data <- filter(all_results, Design == trial_designs[[design_name]]$name)
  plot(design_data$Carryover_Effect, design_data$Power,
       type = "b", col = "steelblue", lwd = 2, pch = 16,
       xlim = range(delta_grid), ylim = c(0, 1),
       xlab = "Carryover Effect (delta_true)",
       ylab = "Estimated Power",
       main = paste("Power vs. Carryover Effect:", trial_designs[[design_name]]$name))
  abline(h = 0.8, lty = 2, col = "red")
}
