#' Log progress to console and optionally to a file
#'
#' This function extends the standard cat() function to optionally write
#' the same output to a log file with appended writes.
#'
#' @param message The message to log
#' @param log_file Optional file path to write logs to. If NULL, only writes to console.
#' @param ... Additional arguments passed to cat()
#' @return Invisibly returns NULL
#' @export
log_progress <- function(message, log_file = NULL, ...) {
  # Write to stdout (console)
  cat(message, ...)
  
  # Also write to log file if specified
  if (!is.null(log_file)) {
    cat(message, ..., file = log_file, append = TRUE)
  }
  
  invisible(NULL)
}

#' Set up a text progress bar that can also log to a file
#'
#' This function creates a progress bar similar to txtProgressBar but with 
#' the option to log progress updates to a file.
#'
#' @param min Minimum value for the progress bar
#' @param max Maximum value for the progress bar
#' @param initial Initial value for the progress bar
#' @param style Style of the progress bar (see txtProgressBar)
#' @param log_file Optional file path to write progress updates to
#' @param log_frequency How often to log to file (in terms of update count)
#' @return A progress bar object with additional logging capabilities
#' @export
logging_progress_bar <- function(min = 0, max = 1, initial = 0, style = 3, 
                                log_file = NULL, log_frequency = 10) {
  # Create standard text progress bar
  pb <- txtProgressBar(min = min, max = max, initial = initial, style = style)
  
  # Add logging functionality
  update_count <- 0
  
  # Create custom update function
  update_fn <- function(value, ...) {
    # Update standard progress bar
    setTxtProgressBar(pb, value)
    
    # Log to file if needed
    if (!is.null(log_file)) {
      update_count <<- update_count + 1
      if (update_count %% log_frequency == 0 || value >= max) {
        percent <- 100 * (value - min) / (max - min)
        log_message <- sprintf("Progress: %.1f%% (%d/%d)\n", percent, value, max)
        cat(log_message, file = log_file, append = TRUE)
      }
    }
  }
  
  # Return custom progress bar with update function
  structure(
    list(
      pb = pb,
      update = update_fn,
      min = min,
      max = max,
      log_file = log_file,
      log_frequency = log_frequency
    ),
    class = "logging_progress_bar"
  )
}

#' Update a logging progress bar
#'
#' @param pb A logging progress bar object
#' @param value Current value
#' @param ... Additional arguments (not used)
#' @return The progress bar object, invisibly
#' @export
update_logging_progress_bar <- function(pb, value, ...) {
  pb$update(value)
  invisible(pb)
}

#' Close a logging progress bar
#'
#' @param pb A logging progress bar object
#' @param ... Additional arguments (not used)
#' @export
close_logging_progress_bar <- function(pb, ...) {
  close(pb$pb)
  invisible(NULL)
}