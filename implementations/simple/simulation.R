# simulation.R
# Simplified N-of-1 trial simulation for biomarker-moderated treatment response
# Based on Hendrickson et al. (2020) design specifications
#
# Features:
# - Four trial designs (CO, N-of-1, OL+BDC) preserved exactly from Hendrickson
# - Carryover with exponential decay (longer half-lives: 0, 0.5, 1, 2 weeks)
# - Optional carryover adjustment in analysis model
# - Dropout/censoring with constant weekly probability (MCAR)
# - ANCOVA on phase means for analysis
# - Random intercept + iid errors for data generation

library(tidyverse)

# --- Parameters -----------------------------------------------------------

n_iterations <- if (Sys.getenv("EXPLORATORY_MODE") == "TRUE") 20 else 500
base_seed <- 42

# Sample sizes (Hendrickson: 35, 70)
sample_sizes <- c(35, 70)

# Biomarker moderation effect (Hendrickson: 0-0.6)
# Full set to match screenshot layout with 6 rows per N
biomarker_mod_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

# Response parameters
br_rate <- 0.3
er_rate <- 0.15

# Variance components
sigma_between <- 1.0
sigma_within <- 0.5

# Carryover half-life in weeks (longer times to see effect)
# 0 = no carryover, 0.5/1/2 = weeks for drug effect to decay by half
carryover_halflife_levels <- c(0, 0.5, 1, 2)

# Carryover adjustment in analysis (TRUE = adjust, FALSE = no adjustment)
# Compare both to show benefit of adjustment
carryover_adjust_levels <- c(FALSE, TRUE)

# Dropout probability per week (MCAR - Missing Completely At Random)
# 0 = no dropout, 0.02/0.05/0.10 = 2%/5%/10% chance of dropping out each week
dropout_prob_levels <- c(0, 0.02, 0.05, 0.10)

# Designs to include (excluding OL which has no drug comparison)
designs <- c("CO", "N-of-1", "OL+BDC")

# --- Design Specifications ------------------------------------------------

create_design <- function(design_name) {
 switch(design_name,
    "OL" = list(
      name = "OL",
      n_weeks = 8,
      on_drug = rep(TRUE, 8),
      blinded = rep(FALSE, 8)
    ),
    "OL+BDC" = list(
      name = "OL+BDC",
      n_weeks = 8,
      on_drug_fixed = c(TRUE, TRUE, TRUE, TRUE, NA, NA, NA, NA),
      blinded = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
    ),
    "CO" = list(
      name = "CO",
      n_weeks = 8,
      on_drug_fixed = rep(NA, 8),
      blinded = rep(TRUE, 8)
    ),
    "N-of-1" = list(
      name = "N-of-1",
      n_weeks = 12,
      on_drug_fixed = c(TRUE, TRUE, TRUE, TRUE, rep(NA, 8)),
      blinded = c(FALSE, FALSE, FALSE, FALSE, rep(TRUE, 8))
    )
  )
}

randomize_assignment <- function(design) {
  if (design$name == "OL") {
    return(design$on_drug)
  }

  on_drug <- design$on_drug_fixed

  if (design$name == "OL+BDC") {
    randomized_weeks <- which(is.na(on_drug))
    if (runif(1) < 0.5) {
      on_drug[randomized_weeks] <- TRUE
    } else {
      on_drug[randomized_weeks] <- FALSE
    }
  } else if (design$name == "CO") {
    if (runif(1) < 0.5) {
      on_drug <- c(rep(TRUE, 4), rep(FALSE, 4))
    } else {
      on_drug <- c(rep(FALSE, 4), rep(TRUE, 4))
    }
  } else if (design$name == "N-of-1") {
    on_drug[1:4] <- TRUE
    cycle1_order <- if (runif(1) < 0.5) c(TRUE, TRUE, FALSE, FALSE) else
                                         c(FALSE, FALSE, TRUE, TRUE)
    cycle2_order <- if (runif(1) < 0.5) c(TRUE, TRUE, FALSE, FALSE) else
                                         c(FALSE, FALSE, TRUE, TRUE)
    on_drug[5:8] <- cycle1_order
    on_drug[9:12] <- cycle2_order
  }

  on_drug
}

# --- Carryover ------------------------------------------------------------

calc_carryover <- function(weeks_since_drug, halflife) {
  if (halflife <= 0 || weeks_since_drug <= 0) {
    0
  } else {
    0.5 ^ (weeks_since_drug / halflife)
  }
}

# --- Data Generation ------------------------------------------------------

generate_response <- function(week, on_drug, blinded, biomarker,
                              biomarker_mod, br_rate, er_rate,
                              carryover_factor = 0) {
  if (on_drug) {
    br <- br_rate * (1 + biomarker_mod * biomarker)
  } else {
    br <- carryover_factor * br_rate * (1 + biomarker_mod * biomarker)
  }

  er <- if (!blinded && on_drug) er_rate else 0

  br + er
}

generate_participant <- function(id, design, biomarker, biomarker_mod,
                                 br_rate, er_rate, sigma_between, sigma_within,
                                 carryover_halflife = 0, dropout_prob = 0) {
  on_drug <- randomize_assignment(design)
  n_weeks <- design$n_weeks
  blinded <- design$blinded

  u_i <- rnorm(1, 0, sigma_between)
  epsilon <- rnorm(n_weeks, 0, sigma_within)

  # Carryover: stateless decay from last on-drug week.
  # For weekly assessments with no drug re-initiation during
  # off-drug periods, this is equivalent to sequential
  # propagation (as used in the orig/2025 repos) because
  # each off-drug week's decay factor is (1/2)^(w - last_drug)
  # regardless of whether computed statelessly or sequentially.
  carryover_factors <- numeric(n_weeks)
  if (carryover_halflife > 0) {
    last_drug_week <- 0
    for (w in seq_len(n_weeks)) {
      if (on_drug[w]) {
        last_drug_week <- w
        carryover_factors[w] <- 0
      } else if (last_drug_week > 0) {
        weeks_since <- w - last_drug_week
        carryover_factors[w] <- calc_carryover(weeks_since, carryover_halflife)
      }
    }
  }

  dropout_week <- NA
  if (dropout_prob > 0) {
    for (w in seq_len(n_weeks)) {
      if (runif(1) < dropout_prob) {
        dropout_week <- w
        break
      }
    }
  }

  responses <- numeric(n_weeks)
  for (w in seq_len(n_weeks)) {
    signal <- generate_response(
      week = w,
      on_drug = on_drug[w],
      blinded = blinded[w],
      biomarker = biomarker,
      biomarker_mod = biomarker_mod,
      br_rate = br_rate,
      er_rate = er_rate,
      carryover_factor = carryover_factors[w]
    )
    responses[w] <- signal + u_i + epsilon[w]
  }

  result <- tibble(
    participant_id = id,
    week = seq_len(n_weeks),
    on_drug = on_drug,
    blinded = blinded,
    biomarker = biomarker,
    carryover = carryover_factors,
    response = responses
  )

  if (!is.na(dropout_week)) {
    result <- result |> filter(week < dropout_week)
  }

  result
}

generate_trial <- function(n_participants, design, biomarker_mod,
                           br_rate, er_rate, sigma_between, sigma_within,
                           carryover_halflife = 0, dropout_prob = 0) {
  biomarkers <- rnorm(n_participants, 0, 1)

  map2_dfr(seq_len(n_participants), biomarkers, function(id, bm) {
    generate_participant(
      id = id,
      design = design,
      biomarker = bm,
      biomarker_mod = biomarker_mod,
      br_rate = br_rate,
      er_rate = er_rate,
      sigma_between = sigma_between,
      sigma_within = sigma_within,
      carryover_halflife = carryover_halflife,
      dropout_prob = dropout_prob
    )
  })
}

# --- Analysis -------------------------------------------------------------

analyze_trial <- function(trial_data, design_name, adjust_carryover = FALSE) {
  if (design_name == "OL") {
    summary_data <- trial_data |>
      group_by(participant_id, biomarker) |>
      summarise(mean_response = mean(response), .groups = "drop")

    model <- lm(mean_response ~ biomarker, data = summary_data)
    p_value <- coef(summary(model))["biomarker", "Pr(>|t|)"]

  } else if (design_name == "OL+BDC") {
    # For OL+BDC: Use change from open-label baseline to blinded phase
    # Everyone gets drug in open-label (weeks 1-4)
    # Blinded phase (weeks 5-8): randomized to continue drug or switch to placebo
    # Drug effect = blinded_response - open_response
    # For drug group: should be ~0 (continuing drug)
    # For placebo group: should be negative (losing drug effect), moderated by biomarker

    summary_data <- trial_data |>
      mutate(phase = ifelse(blinded, "blinded", "open")) |>
      group_by(participant_id, biomarker, phase) |>
      summarise(mean_response = mean(response), .groups = "drop") |>
      pivot_wider(names_from = phase, values_from = mean_response) |>
      filter(!is.na(open) & !is.na(blinded))

    # Get blinded phase assignment for each participant
    blinded_assignment <- trial_data |>
      filter(blinded) |>
      group_by(participant_id) |>
      summarise(on_drug = first(on_drug), .groups = "drop")

    combined <- summary_data |>
      left_join(blinded_assignment, by = "participant_id") |>
      mutate(
        drug_effect = blinded - open,
        group = ifelse(on_drug, "drug", "placebo")
      )

    if (nrow(combined) < 10 || length(unique(combined$group)) < 2) {
      return(NA_real_)
    }

    # Test whether biomarker moderates the change from open to blinded
    # Drug group: drug_effect should be ~0 regardless of biomarker
    # Placebo group: drug_effect should be negative, MORE negative for high biomarker
    model <- lm(drug_effect ~ biomarker * group, data = combined)
    coefs <- coef(summary(model))
    if ("biomarker:groupplacebo" %in% rownames(coefs)) {
      p_value <- coefs["biomarker:groupplacebo", "Pr(>|t|)"]
    } else {
      p_value <- NA_real_
    }

  } else {
    if (adjust_carryover) {
      # NOTE: This is an oracle analysis -- br_rate is the true DGP
      # parameter, not estimable from data. Power estimates with
      # carryover adjustment are therefore upper bounds on what a
      # real analysis could achieve.
      summary_data <- trial_data |>
        mutate(adjusted_response = response - carryover * br_rate) |>
        group_by(participant_id, biomarker, on_drug) |>
        summarise(mean_response = mean(adjusted_response), .groups = "drop") |>
        pivot_wider(names_from = on_drug, values_from = mean_response,
                    names_prefix = "drug_") |>
        filter(!is.na(drug_TRUE) & !is.na(drug_FALSE)) |>
        mutate(drug_effect = drug_TRUE - drug_FALSE)
    } else {
      summary_data <- trial_data |>
        group_by(participant_id, biomarker, on_drug) |>
        summarise(mean_response = mean(response), .groups = "drop") |>
        pivot_wider(names_from = on_drug, values_from = mean_response,
                    names_prefix = "drug_") |>
        filter(!is.na(drug_TRUE) & !is.na(drug_FALSE)) |>
        mutate(drug_effect = drug_TRUE - drug_FALSE)
    }

    if (nrow(summary_data) < 5) {
      return(NA_real_)
    }

    model <- lm(drug_effect ~ biomarker, data = summary_data)
    p_value <- coef(summary(model))["biomarker", "Pr(>|t|)"]
  }

  p_value
}

# --- Main Simulation Loop -------------------------------------------------

# --- Simulation 1: Carryover Assessment ------------------------------------
param_grid_carryover <- expand_grid(
  design = designs,
  n_participants = sample_sizes,
  biomarker_mod = biomarker_mod_levels,
  carryover_halflife = carryover_halflife_levels,
  carryover_adjust = carryover_adjust_levels
)

cat("=== CARRYOVER SIMULATION ===\n")
cat("Running", n_iterations, "iterations per condition\n")
cat("Total conditions:", nrow(param_grid_carryover), "\n")
cat("Designs:", paste(designs, collapse = ", "), "\n")
cat("Carryover half-lives:", paste(carryover_halflife_levels, collapse = ", "),
    "weeks\n")
cat("Carryover adjustment:", paste(carryover_adjust_levels, collapse = ", "),
    "\n\n")

results_carryover <- param_grid_carryover |>
  mutate(
    condition_id = row_number(),
    power = pmap_dbl(
      list(design, n_participants, biomarker_mod, carryover_halflife,
           carryover_adjust, condition_id),
      function(d, n, bm_mod, co_hl, co_adj, cid) {
        set.seed(base_seed + cid)
        design_obj <- create_design(d)

        p_values <- replicate(n_iterations, {
          trial_data <- generate_trial(
            n_participants = n,
            design = design_obj,
            biomarker_mod = bm_mod,
            br_rate = br_rate,
            er_rate = er_rate,
            sigma_between = sigma_between,
            sigma_within = sigma_within,
            carryover_halflife = co_hl,
            dropout_prob = 0
          )
          analyze_trial(trial_data, d, adjust_carryover = co_adj)
        })

        mean(p_values < 0.05, na.rm = TRUE)
      }
    ),
    .progress = TRUE
  )

# --- Simulation 2: Dropout Assessment --------------------------------------
param_grid_dropout <- expand_grid(
  design = designs,
  n_participants = sample_sizes,
  biomarker_mod = biomarker_mod_levels,
  dropout_prob = dropout_prob_levels
)

cat("\n=== DROPOUT SIMULATION ===\n")
cat("Running", n_iterations, "iterations per condition\n")
cat("Total conditions:", nrow(param_grid_dropout), "\n")
cat("Dropout probabilities:", paste(dropout_prob_levels, collapse = ", "),
    "per week\n\n")

results_dropout <- param_grid_dropout |>
  mutate(
    condition_id = row_number(),
    power = pmap_dbl(
      list(design, n_participants, biomarker_mod, dropout_prob,
           condition_id),
      function(d, n, bm_mod, dp, cid) {
        set.seed(base_seed + 10000 + cid)
        design_obj <- create_design(d)

        p_values <- replicate(n_iterations, {
          trial_data <- generate_trial(
            n_participants = n,
            design = design_obj,
            biomarker_mod = bm_mod,
            br_rate = br_rate,
            er_rate = er_rate,
            sigma_between = sigma_between,
            sigma_within = sigma_within,
            carryover_halflife = 0,
            dropout_prob = dp
          )
          analyze_trial(trial_data, d, adjust_carryover = FALSE)
        })

        mean(p_values < 0.05, na.rm = TRUE)
      }
    ),
    .progress = TRUE
  )

# --- Output ---------------------------------------------------------------

cat("\n--- Carryover Power Results ---\n\n")
results_carryover |>
  mutate(power = round(power, 2)) |>
  pivot_wider(names_from = design, values_from = power) |>
  print(n = 100)

cat("\n--- Dropout Power Results ---\n\n")
results_dropout |>
  mutate(power = round(power, 2)) |>
  pivot_wider(names_from = design, values_from = power) |>
  print(n = 100)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- here::here("analysis", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

write_csv(results_carryover, file.path(output_dir,
  paste0("power_carryover_", timestamp, ".csv")))
write_csv(results_dropout, file.path(output_dir,
  paste0("power_dropout_", timestamp, ".csv")))

# --- Heatmaps -------------------------------------------------------------

create_carryover_heatmap <- function(data, title, n_iter) {
  plot_data <- data |>
    mutate(
      carryover_label = paste0("t½=", carryover_halflife, " wk"),
      carryover_label = factor(carryover_label,
        levels = paste0("t½=", sort(unique(carryover_halflife)), " wk")),
      n_label = paste0("N=", n_participants),
      power_pct = round(power * 100),
      design = factor(design, levels = c("CO", "N-of-1", "OL+BDC"))
    )

  p <- plot_data |>
    ggplot(aes(x = carryover_label, y = factor(biomarker_mod),
               fill = power)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = power_pct), size = 3, fontface = "bold") +
    facet_grid(n_label ~ design, switch = "y") +
    scale_fill_gradientn(
      colors = c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60",
                 "#1a9850"),
      values = scales::rescale(c(0, 0.25, 0.5, 0.75, 0.9, 1)),
      limits = c(0, 1),
      name = "Power (%)",
      labels = scales::percent
    ) +
    scale_y_discrete(limits = rev) +
    labs(
      title = title,
      subtitle = paste0("Hendrickson et al. 2020 approach | ", n_iter,
                        " iterations"),
      x = NULL,
      y = "Biomarker\nModeration"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      strip.placement = "outside",
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )

  p
}

create_dropout_heatmap <- function(data, title, n_iter) {
  plot_data <- data |>
    mutate(
      dropout_label = paste0(dropout_prob * 100, "%/wk"),
      dropout_label = factor(dropout_label,
        levels = paste0(sort(unique(dropout_prob)) * 100, "%/wk")),
      n_label = paste0("N=", n_participants),
      power_pct = round(power * 100),
      design = factor(design, levels = c("CO", "N-of-1", "OL+BDC"))
    )

  p <- plot_data |>
    ggplot(aes(x = dropout_label, y = factor(biomarker_mod),
               fill = power)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = power_pct), size = 3, fontface = "bold") +
    facet_grid(n_label ~ design, switch = "y") +
    scale_fill_gradientn(
      colors = c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60",
                 "#1a9850"),
      values = scales::rescale(c(0, 0.25, 0.5, 0.75, 0.9, 1)),
      limits = c(0, 1),
      name = "Power (%)",
      labels = scales::percent
    ) +
    scale_y_discrete(limits = rev) +
    labs(
      title = title,
      subtitle = paste0("Hendrickson et al. 2020 approach | ", n_iter,
                        " iterations"),
      x = "Weekly Dropout Rate",
      y = "Biomarker\nModeration"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      strip.placement = "outside",
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10)
    )

  p
}

# Generate carryover heatmaps (with/without adjustment)
results_no_adj <- results_carryover |> filter(!carryover_adjust)
results_with_adj <- results_carryover |> filter(carryover_adjust)

p_no_adj <- create_carryover_heatmap(
  results_no_adj,
  "Statistical Power: WITHOUT Carryover in Analysis Model",
  n_iterations
)

ggsave(
  file.path(output_dir, paste0("power_no_adjustment_", timestamp, ".png")),
  p_no_adj, width = 11, height = 8, dpi = 150
)

if (nrow(results_with_adj) > 0) {
  p_with_adj <- create_carryover_heatmap(
    results_with_adj,
    "Statistical Power: WITH Carryover in Analysis Model",
    n_iterations
  )

  ggsave(
    file.path(output_dir, paste0("power_with_adjustment_", timestamp, ".png")),
    p_with_adj, width = 11, height = 8, dpi = 150
  )
}

# Generate dropout heatmap
p_dropout <- create_dropout_heatmap(
  results_dropout,
  "Statistical Power: Effect of Dropout on Power",
  n_iterations
)

ggsave(
  file.path(output_dir, paste0("power_dropout_", timestamp, ".png")),
  p_dropout, width = 11, height = 8, dpi = 150
)

cat("\nResults saved to:", output_dir, "\n")
