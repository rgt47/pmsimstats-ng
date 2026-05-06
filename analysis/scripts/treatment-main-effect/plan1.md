# Detailed Plan to Expand Simulation to Cover Hendrickson Parameter Range

# 7/10/25

## Overview
Expand the existing carryover simulation to incorporate Hendrickson's three-factor response model (BR, ER, TR) and compare N-of-1 vs traditional crossover designs over a 5-period trial structure.

## 1. Effect Size Parameterization

### Current vs Target Approach
**Your current**: Single `base_treat_effect` parameter  
**Hendrickson target**: Multi-component effect sizes through BR-ER-TR decomposition

### Implementation Strategy
```r
effect_size_params <- list(
  # Biomarker-related treatment effect (BR component)
  biomarker_treatment_correlation = c(0.0, 0.3, 0.6),  # As in Hendrickson Fig 4
  br_max_response = c(5, 10, 15),  # CAPS-IV scale points
  br_rate = c(0.2, 0.35, 0.5),     # Gompertz rate parameter
  
  # Expectancy effect (ER component) 
  er_max_response = c(2, 5, 8),    # Smaller than BR typically
  er_rate = c(0.2, 0.35, 0.5),     # Similar range to BR
  
  # Time/study participation effect (TR component)
  tr_max_response = c(3, 6, 9),    # Regression to mean + study effects
  tr_rate = c(0.1, 0.2, 0.3),      # Often slower than treatment effects
  
  # Carryover characterization
  carryover_halflife = c(0, 0.1, 0.2, 0.5, 1.0)  # weeks, as in Hendrickson
)
```

### Effect Size Calculation
```r
# Total treatment effect becomes biomarker-dependent:
# For subject i with biomarker value b_i:
# Expected treatment effect = BR_max * correlation * (b_i - mean(b)) / sd(b)
# Plus expectancy and time components
```

## 2. Three-Factor Response Model (BR, ER, TR)

### Mathematical Framework
Following Hendrickson's Gompertz function approach:

```r
# Component j (BR, ER, or TR) for subject i at time t
calculate_component <- function(t, max_response, rate, displacement = 0) {
  max_response * exp(-exp(-rate * (t - displacement)))
}

generate_subject_trajectory <- function(subject_id, biomarker_value, design, params) {
  trajectory <- tibble(
    Subject = subject_id,
    Biomarker = biomarker_value,
    Period = 1:5,
    Treatment = design$sequence,
    
    # Expectancy based on known vs blinded treatment
    Expectancy = case_when(
      design$blinding[Period] == "open" ~ 1.0,      # Know on treatment
      design$blinding[Period] == "blinded" ~ 0.5,   # 50% chance on treatment  
      design$blinding[Period] == "placebo" ~ 0.0,   # Know on placebo
      TRUE ~ 0.5
    )
  ) %>%
  rowwise() %>%
  mutate(
    # BR: Biomarker-dependent biological response
    BR = ifelse(Treatment == "D", 
                calculate_component(Period, 
                                  params$br_max * biomarker_effect_size(Biomarker), 
                                  params$br_rate),
                0),
    
    # ER: Expectancy-related response  
    ER = calculate_component(Period, 
                           params$er_max * Expectancy, 
                           params$er_rate),
    
    # TR: Time-related response (study participation, regression to mean)
    TR = calculate_component(Period, params$tr_max, params$tr_rate),
    
    # Add carryover effects
    BR_with_carryover = add_carryover_effects(BR, Treatment, params$carryover_halflife),
    
    # Total observed response
    Y_true = BR_with_carryover + ER + TR,
    Y_observed = Y_true + rnorm(n(), 0, params$sigma)
  )
}

biomarker_effect_size <- function(biomarker_value) {
  # Standardized biomarker effect
  params$biomarker_treatment_correlation * (biomarker_value - mean_biomarker) / sd_biomarker
}
```

### Carryover Implementation
```r
add_carryover_effects <- function(br_values, treatments, halflife) {
  if (halflife == 0) return(br_values)
  
  br_with_carryover <- br_values
  decay_rate <- log(2) / halflife
  
  for (t in 2:length(treatments)) {
    if (treatments[t] == "P" && treatments[t-1] == "D") {
      # Add exponentially decaying carryover
      carryover <- br_values[t-1] * exp(-decay_rate * 1)  # 1 period delay
      br_with_carryover[t] <- br_with_carryover[t] + carryover
    }
  }
  return(br_with_carryover)
}
```

## 3. Trial Design Characterization

### Design Definitions for 5-Period Trials

```r
# N-of-1 Design (Hendrickson Design 4)
nof1_design <- list(
  name = "N-of-1",
  description = "Open-label → Blinded discontinuation → Brief crossover",
  sequences = list(
    # Period:        1    2    3    4    5
    N1 = list(
      treatment = c("D", "D", "P", "D", "P"),
      blinding =  c("open", "open", "blinded", "blinded", "blinded")
    ),
    N2 = list(
      treatment = c("D", "D", "P", "P", "D"), 
      blinding =  c("open", "open", "blinded", "blinded", "blinded")
    ),
    N3 = list(
      treatment = c("D", "P", "P", "D", "P"),
      blinding =  c("open", "blinded", "blinded", "blinded", "blinded") 
    ),
    N4 = list(
      treatment = c("D", "P", "P", "P", "D"),
      blinding =  c("open", "blinded", "blinded", "blinded", "blinded")
    )
  ),
  randomization_points = c("sequence", "period_2_treatment")
)

# Traditional Crossover Design (Hendrickson Design 3, adapted to 5 periods)
crossover_design <- list(
  name = "Traditional Crossover", 
  description = "Balanced crossover with equal drug/placebo exposure",
  sequences = list(
    # Period:        1    2    3    4    5
    C1 = list(
      treatment = c("D", "D", "P", "P", "P"),
      blinding =  c("blinded", "blinded", "blinded", "blinded", "blinded")
    ),
    C2 = list(
      treatment = c("P", "P", "D", "D", "D"),
      blinding =  c("blinded", "blinded", "blinded", "blinded", "blinded")
    )
  ),
  randomization_points = c("sequence")
)

trial_designs <- list(
  nof1 = nof1_design,
  crossover = crossover_design
)
```

### Design Assignment Function
```r
assign_design_sequence <- function(subject_id, design) {
  if (design$name == "N-of-1") {
    # Random assignment to one of 4 sequences
    sequence_name <- sample(names(design$sequences), 1)
    sequence <- design$sequences[[sequence_name]]
    
    # Additional randomization for period 2 in some sequences
    if (sequence_name %in% c("N3", "N4")) {
      # 50% chance of early vs late discontinuation  
      if (runif(1) < 0.5) {
        sequence$treatment[2] <- "P"  # Early discontinuation
      }
    }
  } else if (design$name == "Traditional Crossover") {
    # Simple randomization to D-first or P-first
    sequence_name <- sample(names(design$sequences), 1)
    sequence <- design$sequences[[sequence_name]]
  }
  
  return(sequence)
}
```

## 4. Integration with Existing ADEMP Structure

### Modified ADEMP Parameters
```r
sim_params <- list(
  # Existing parameters
  n_subj = 40,
  nsim = 1000, 
  alpha = 0.05,
  
  # Hendrickson-style parameters
  biomarker_mean = 140,      # Systolic BP baseline
  biomarker_sd = 20,
  baseline_caps_mean = 70,   # CAPS-IV baseline severity
  baseline_caps_sd = 15,
  
  # Three-factor model parameters
  br_max = 10,               # Max biological response (CAPS points)
  er_max = 5,                # Max expectancy response  
  tr_max = 8,                # Max time-related response
  br_rate = 0.35,            # Response rate parameters
  er_rate = 0.35,
  tr_rate = 0.2,
  
  # Biomarker-treatment interaction
  biomarker_treatment_correlation = 0.6,  # Key effect size parameter
  
  # Carryover
  carryover_halflife = 0.1,   # weeks
  
  # Noise
  sigma = 8                   # Residual SD (CAPS scale)
)
```

### Analysis Method Updates
```r
# Hendrickson-style analysis for biomarker interaction
analyze_biomarker_interaction <- function(data) {
  # For designs with both drug and placebo periods
  if (any(data$Treatment == "D") && any(data$Treatment == "P")) {
    model <- lmer(Y_observed ~ Treatment * Biomarker + Period + (1|Subject), 
                  data = data)
    interaction_coef <- fixef(model)["TreatmentP:Biomarker"]
    interaction_p <- summary(model)$coefficients["TreatmentP:Biomarker", "Pr(>|t|)"]
  } else {
    # For open-label only periods (like Hendrickson's OL design)
    model <- lmer(Y_observed ~ Biomarker * Period + (1|Subject), 
                  data = data)
    interaction_coef <- fixef(model)["Biomarker:Period"] 
    interaction_p <- summary(model)$coefficients["Biomarker:Period", "Pr(>|t|)"]
  }
  
  return(list(
    p_value = interaction_p,
    estimate = interaction_coef,
    # ... other return values
  ))
}
```

## 5. Simulation Execution Plan

### Parameter Grid
```r
# Comprehensive parameter exploration matching Hendrickson
parameter_grid <- expand_grid(
  design = c("nof1", "crossover"),
  biomarker_correlation = c(0.0, 0.3, 0.6),
  carryover_halflife = c(0, 0.1, 0.2, 0.5),
  br_max = c(5, 10, 15),
  n_subj = c(25, 35, 50)  # Sample size exploration
)
```

### Expected Outputs
- **Power curves** by design and carryover level (reproducing Hendrickson Fig 4)
- **Bias assessment** for biomarker effect estimates  
- **Component-wise analysis** showing BR vs ER vs TR contributions
- **Comparative plots** showing N-of-1 vs crossover performance

This expansion would provide direct comparison with Hendrickson's findings while maintaining your ADEMP structure and adding the missing Monte Carlo standard errors they didn't report.
