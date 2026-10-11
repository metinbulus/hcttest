# Check before means, Winsorized SDs, or population integration are calculated.
.check_retained_sample_sizes <- function(n1, n2, trim) {
  for (n in list(n1, n2)) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n)) {
      stop("Each group sample size must be a finite numeric scalar.",
           call. = FALSE)
    }
    retained <- n - 2 * round(n * trim)
    if (n < 2 || retained < 2) {
      stop("Each group must retain at least two observations after trimming.",
           call. = FALSE)
    }
  }
  invisible(NULL)
}
