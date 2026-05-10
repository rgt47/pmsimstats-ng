#' Resolve common conflicts
#'
#' This function resolves common package conflicts in the nof1power package.
#' It explicitly declares preferences for certain functions when multiple packages
#' provide functions with the same name.
#'
#' @return NULL, invisibly
#' @export
resolve_conflicts <- function() {
  # If conflicted package is available, use it
  if (requireNamespace("conflicted", quietly = TRUE)) {
    # Prefer dplyr functions over others
    conflicted::conflict_prefer("select", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("lag", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)
    
    # Add any other conflicts that need resolution
  } else {
    # If conflicted package is not available
    message("The 'conflicted' package is not installed. Consider installing it with: install.packages('conflicted')")
  }
  
  # Return invisibly
  invisible(NULL)
}