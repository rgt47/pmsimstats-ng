# ------------------------------------------------------------------------------
# Simulation Study: t-test vs GLS in 4-Sequence N-of-1 Crossover Design
# Refactored to use a `for` loop with progress bar
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(nlme)

# ---- GLOBAL PARAMETERS ----
set.seed(2025)
n_subj <- 40
base_treat_effect <- .5
sigma <- 1.0
nsim <- 1000
alpha <- 0.05

# ---- DESIGN SEQUENCES ----
sequences <- list(
  S1 = c("D", "D", "P", "D", "P", "P"),
  S2 = c("D", "D", "P", "P", "D", "D"),
  S3 = c("D", "P", "P", "D", "P", "P"),
  S4 = c("D", "P", "P", "P", "D", "D")
)

# ---- FUNCTION: Run One Simulation ----
run_one_sim <- function(delta_true) {
  subj_data <- tibble(Subject = 1:n_subj) %>%
    mutate(Sequence = sample(names(sequences), n(), replace = TRUE)) %>%
    rowwise() %>%
    mutate(TreatmentSeq = list(sequences[[Sequence]])) %>%
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

  # Paired contrast t-test
  contrast_df <- subj_data %>%
    group_by(Subject) %>%
    summarise(
      mean_D = mean(Y[Treatment == "D"]),
      mean_P = mean(Y[Treatment == "P"]),
      contrast = mean_D - mean_P,
      .groups = "drop"
    )

  p_ttest <- tryCatch(
    t.test(contrast_df$contrast, mu = 0)$p.value,
    error = function(e) NA_real_
  )

  # GLS model
  p_gls <- tryCatch({
    fit <- gls(
      Y ~ Treatment,
      data = subj_data,
      correlation = corCompSymm(form = ~1 | Subject),
      method = "REML"
    )
    summary(fit)$tTable["TreatmentP", "p-value"]
  }, error = function(e) NA_real_)

  list(p_ttest = p_ttest, p_gls = p_gls)
}

# ---- FUNCTION: Estimate Power at One Carryover Level ----
run_power_for_delta <- function(delta_true) {
  p_ttest <- numeric(nsim)
  p_gls   <- numeric(nsim)

  for (i in seq_len(nsim)) {
    sim_result <- run_one_sim(delta_true)
    p_ttest[i] <- sim_result$p_ttest
    p_gls[i]   <- sim_result$p_gls
    if (i %% 100 == 0) cat(sprintf("Delta=%s: %d/%d\n", delta_true, i, nsim))
  }

  c(
    Power_ttest = mean(p_ttest < alpha, na.rm = TRUE),
    Power_gls   = mean(p_gls < alpha, na.rm = TRUE)
  )
}

# ---- RUN POWER ANALYSIS OVER GRID ----
delta_grid <- seq(0, 1, by = 0.2)
power_results <- t(sapply(delta_grid, run_power_for_delta))

# ---- TABULATE RESULTS ----
power_table <- tibble(
  Carryover_Effect = delta_grid,
  Power_ttest = power_results[, "Power_ttest"],
  Power_gls   = power_results[, "Power_gls"]
)

print(power_table)

# ---- OPTIONAL: PLOT RESULTS ----
power_long <- tidyr::pivot_longer(power_table, cols = starts_with("Power"),
                                  names_to = "Method", values_to = "Power")

plot(NULL, xlim = range(delta_grid), ylim = c(0, 1),
     xlab = "Carryover Effect (delta_true)",
     ylab = "Estimated Power",
     main = "Power vs. Carryover Effect\nPaired Contrast t-test vs GLS")
methods <- unique(power_long$Method)
cols <- c("steelblue", "firebrick")
for (j in seq_along(methods)) {
  d <- power_long[power_long$Method == methods[j], ]
  lines(d$Carryover_Effect, d$Power, col = cols[j], lwd = 2)
  points(d$Carryover_Effect, d$Power, col = cols[j], pch = 16)
}
abline(h = 0.8, lty = 2, col = "gray40")
legend("topright", legend = methods, col = cols, lwd = 2, pch = 16)
