#' Merge values from one list into another
#'
#' Overrides entries of `x` with entries of `val` and appends any
#' entries of `val` not present in `x`. A simpler analogue of
#' [utils::modifyList()] used by the analysis-options plumbing.
#'
#' @param x Original list.
#' @param val List whose entries override `x`.
#' @return A list with combined entries. If `x` is `NULL`, returns
#'   `val`; if `val` is `NULL`, returns `x`.
#' @examples
#' defaults <- list(a = 1, b = 2, c = 3)
#' override <- list(b = 20, d = 40)
#' modify_list(defaults, override)
#' @export
modify_list <- function(x, val) {
  if (is.null(x)) return(val)
  if (is.null(val)) return(x)

  vnames <- names(val)
  vnames <- vnames[nzchar(vnames)]

  for (v in vnames) {
    x[[v]] <- val[[v]]
  }

  x
}
