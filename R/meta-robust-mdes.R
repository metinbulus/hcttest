#' Append study-specific minimum detectable standardized effects
#'
#' @param data A nonempty data frame with unique column names.
#' @param sd1,sd2,n1,n2 Column names for positive group SDs and integer group
#'   sizes of at least two. Group 1 is the numerator of both ratios.
#' @param power Target power, strictly between zero and one; default 0.80.
#' @param alpha Significance level passed to power.hct.test; default 0.05.
#' @param alternative Two-sided or one-sided detection. Equivalence effect
#'   searches are excluded because their midpoint is not an MDES.
#' @param req.sign Direction of the effect search, '+' or '-'; default '+'.
#'   This is not inferred from observed group means.
#' @param null.d,margin Reference effect and null margin, as in power.hct.test.
#' @param method Variance estimator passed to power.hct.test; default 'hc2'.
#' @param claim.basis 'smd.ci' (default) includes standardizer uncertainty;
#'   'md.pval' uses mean-difference inference, as in power.hct.test.
#' @param test.distribution 't', 'z', 'lp', or 'lambdap', as in power.hct.test.
#' @return The input data with a numeric mdes column appended (or replaced
#'   if already present). No pooling, weights, intervals, or plot class.
#' @details Internally calls power.hct.test for each row with d = NULL,
#'   var.ratio = (sd1/sd2)^2, n.ratio = n1/n2, the row's n2, and verbose = 0.
#'   Stores the returned d without changing its sign. The selected method
#'   determines the standardizer, as in power.hct.test. Input data are not
#'   modified in place. Missing or invalid inputs produce an error identifying
#'   the row. No dependency on metafor is introduced.
#' @seealso power.hct.test, meta.robust.power
#' @examples
#' studies <- data.frame(n1 = c(107, 86), n2 = c(67, 294),
#'                       sd1 = c(2, 4), sd2 = c(11, 7))
#' meta.robust.mdes(studies, power = 0.80, method = 'hc2')
#' @export
meta.robust.mdes <- function(data, sd1 = "sd1", sd2 = "sd2",
                            n1 = "n1", n2 = "n2", power = 0.80,
                            alpha = 0.05,
                            alternative = c("two.sided", "one.sided"),
                            req.sign = "+", null.d = 0, margin = 0,
                            method = "hc2",
                            claim.basis = c("smd.ci", "md.pval"),
                            test.distribution = c("t", "z", "lp", "lambdap")) {
  alternative <- match.arg(alternative)
  claim.basis <- match.arg(claim.basis)
  test.distribution <- match.arg(test.distribution)
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= 0 || power >= 1)
    stop("'power' must be a number strictly between zero and one.", call. = FALSE)
  values <- .meta_power_columns(data, list(sd1 = sd1, sd2 = sd2, n1 = n1, n2 = n2))
  .meta_power_design(values)
  data$mdes <- .meta_power_rows(values, function(i) {
    fit <- power.hct.test(d = NULL, power = power, n2 = values$n2[i],
      n.ratio = values$n1[i] / values$n2[i],
      var.ratio = (values$sd1[i] / values$sd2[i])^2,
      alpha = alpha, alternative = alternative, req.sign = req.sign,
      null.d = null.d, margin = margin, method = method,
      claim.basis = claim.basis, test.distribution = test.distribution,
      ceil.n = FALSE, verbose = 0)
    if (length(fit$d) != 1L || !is.finite(fit$d) ||
        length(fit$power) != 1L || !is.finite(fit$power) ||
        abs(fit$power - power) > 1e-4)
      stop("The MDES solver did not attain the target power within tolerance.")
    fit$d
  })
  data
}

# Shared internal helpers for both study-level sensitivity functions.
.meta_power_columns <- function(data, fields) {
  if (!is.data.frame(data) || nrow(data) == 0L || anyDuplicated(names(data)))
    stop("'data' must be a nonempty data frame with unique column names.", call. = FALSE)
  lapply(fields, function(name) {
    if (!is.character(name) || length(name) != 1L || is.na(name) ||
        !name %in% names(data))
      stop("Each column argument must name a column in 'data'.", call. = FALSE)
    value <- data[[name]]
    if (!is.numeric(value) || !is.null(dim(value)))
      stop(sprintf("Column '%s' must be a numeric vector.", name), call. = FALSE)
    bad <- which(!is.finite(value))
    if (length(bad))
      stop(sprintf("Row %d: column '%s' must be finite and nonmissing.",
                   bad[1L], name), call. = FALSE)
    value
  })
}

.meta_power_design <- function(values) {
  bad <- which(values$sd1 <= 0 | values$sd2 <= 0 |
    values$n1 < 2 | values$n2 < 2 |
    values$n1 != floor(values$n1) | values$n2 != floor(values$n2))
  if (length(bad))
    stop(sprintf("Row %d: SDs must be positive and group sizes must be integers >= 2.",
                 bad[1L]), call. = FALSE)
}

.meta_power_rows <- function(values, calculate) {
  vapply(seq_along(values$n2), function(i) {
    tryCatch(calculate(i), error = function(e)
      stop(sprintf("Study row %d: %s", i, conditionMessage(e)), call. = FALSE))
  }, numeric(1))
}
