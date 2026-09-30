utils::globalVariables(".data")

#' Print a standardized function header banner
#'
#' Shared helper so every user-facing function opens with the same block:
#' a rule, "<name> -- <what it does>", a last-updated date, and a closing rule.
#'
#' @param name character. Function name.
#' @param what character. One-line description of what the function does.
#' @param verbose logical. If \code{TRUE} (default) the banner is printed.
#' @param updated character. Last-updated date shown in the banner.
#' @keywords internal
#' @noRd
.et_banner <- function(name, what, verbose = TRUE, updated = "09/27/2026") {
  if (!isTRUE(verbose)) return(invisible(NULL))
  bar <- strrep("-", 63)
  cat(bar, "\n",
      name, " -- ", what, "\n",
      "Last updated at ", updated, "\n",
      bar, "\n", sep = "")
  invisible(NULL)
}

#' Print a single progress/description step, gated by verbose
#'
#' @param msg character. What is being done.
#' @param verbose logical. If \code{TRUE} (default) the step is printed.
#' @keywords internal
#' @noRd
.et_step <- function(msg, verbose = TRUE) {
  if (isTRUE(verbose)) cat("  - ", msg, "\n", sep = "")
  invisible(NULL)
}
