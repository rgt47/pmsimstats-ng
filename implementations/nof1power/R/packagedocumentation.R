#' nof1power: N-of-1 Clinical Trial Power Analysis
#'
#' Monte Carlo simulation framework for comparing N-of-1 and
#' parallel-group clinical trial designs for detecting
#' biomarker-treatment interactions under carryover. The package
#' mirrors the pmsimstats-ng tidyverse reference implementation for
#' all core machinery (data generation, covariance construction,
#' linear mixed-effects analysis, simulation orchestration) and
#' extends it with a carryover decay family, a sensitivity-analysis
#' workflow, and logging and parallel plumbing suited to long batch
#' runs.
#'
#' @section Trial-design and data generation:
#' \itemize{
#'   \item \code{\link{build_trial_design}}
#'   \item \code{\link{generate_data}}
#'   \item \code{\link{build_sigma_matrix}}
#' }
#'
#' @section Analysis and simulation:
#' \itemize{
#'   \item \code{\link{lme_analysis}}
#'   \item \code{\link{generate_simulated_results}}
#'   \item \code{\link{censor_data}}
#' }
#'
#' @section Carryover models:
#' \itemize{
#'   \item \code{\link{carryover_decay}}
#'   \item \code{\link{apply_carryover_to_component}}
#'   \item \code{\link{calculate_carryover}}
#'   \item \code{\link{run_carryover_sensitivity}}
#' }
#'
#' @section Utilities and infrastructure:
#' \itemize{
#'   \item \code{\link{cumulative}}, \code{\link{mod_gompertz}}
#'   \item \code{\link{modify_list}}, \code{\link{create_ar1_corr}}
#'   \item \code{\link{log_progress}},
#'     \code{\link{logging_progress_bar}}
#'   \item \code{\link{setup_parallel_processing}},
#'     \code{\link{prepare_parallel_globals}}
#' }
#'
#' @name nof1power-package
#' @aliases nof1power
#' @keywords internal
"_PACKAGE"
