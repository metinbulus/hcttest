#' Append study-specific power
#'
#' @param data A nonempty data frame with unique column names.
#' @param d NULL (default), a numeric effect shared by all studies, a numeric
#'   vector with one effect per row, or the name of a numeric effect column.
#'   With NULL, calculate uncorrected standardized differences from m1 and m2.
#' @param m1,m2 Names of group-mean columns, used only when d is NULL.
#' @param sd1,sd2,n1,n2 Column names for positive group SDs and integer group
#'   sizes of at least two. Group 1 is the numerator of both ratios.
#' @param alpha Significance level; default 0.05.
#' @param alternative 'two.sided', 'one.sided', or 'two.one.sided', as in
#'   power.hct.test. One-sided direction follows the supplied effect relative
#'   to the null margin; it is not a separately fixed directional hypothesis.
#' @param null.d,margin Reference effect and null margin, as in power.hct.test.
#' @param method Variance estimator passed to power.hct.test; default 'hc2'.
#' @param claim.basis 'smd.ci' (default) or 'md.pval', as in power.hct.test.
#' @param test.distribution 't', 'z', 'lp', or 'lambdap', as in power.hct.test.
#' @return The input data with a numeric power column appended (or replaced
#'   if already present). No pooling, weights, intervals, or plot class.
#' @details Internally calls power.hct.test for each row with power = NULL,
#'   var.ratio = (sd1/sd2)^2, n.ratio = n1/n2, the row's n2, and verbose = 0.
#'   With d = NULL, d = (m1-m2)/S, where S is sd1 for method 'sd1', sd2 for
#'   method 'sd2', or the usual degrees-of-freedom-weighted pooled SD for all
#'   other methods. No Hedges correction or absolute-value transformation is
#'   applied. Supplied d values must already use the corresponding standardizer.
#'   Power based on observed means is descriptive retrospective power, not
#'   independent evidence for the observed result. For sensitivity analysis,
#'   supply a scientifically meaningful assumed effect through d.
#'   Input data are not modified in place. Invalid inputs identify the row.
#'   No dependency on metafor is introduced.
#' @seealso power.hct.test, meta.robust.mdes
#' @examples
#' studies <- data.frame(m1 = c(48, 6), m2 = c(50, 7),
#'   n1 = c(107, 86), n2 = c(67, 294), sd1 = c(2, 4), sd2 = c(11, 7))
#' meta.robust.power(studies, d = 0.30)
#' meta.robust.power(studies)
#' studies <- meta.robust.mdes(studies)
#' meta.robust.power(studies, d = 'mdes') # approximately 0.80 in each row
#' @export
meta.robust.power <- function(data, d = NULL, m1 = "m1", m2 = "m2",
                             sd1 = "sd1", sd2 = "sd2", n1 = "n1", n2 = "n2",
                             alpha = 0.05,
                             alternative = c("two.sided", "one.sided", "two.one.sided"),
                             null.d = 0, margin = 0, method = "hc2",
                             claim.basis = c("smd.ci", "md.pval"),
                             test.distribution = c("t", "z", "lp", "lambdap")) {
  alternative <- match.arg(alternative)
  claim.basis <- match.arg(claim.basis)
  test.distribution <- match.arg(test.distribution)
  values <- .meta_power_columns(data, list(sd1 = sd1, sd2 = sd2, n1 = n1, n2 = n2))
  .meta_power_design(values)
  if (is.null(d)) {
    means <- .meta_power_columns(data, list(m1 = m1, m2 = m2))
    standardizer <- if (identical(method, "sd1")) values$sd1 else
      if (identical(method, "sd2")) values$sd2 else
        sqrt(((values$n1 - 1) * values$sd1^2 +
              (values$n2 - 1) * values$sd2^2) / (values$n1 + values$n2 - 2))
    d <- (means$m1 - means$m2) / standardizer
  } else if (is.character(d)) {
    d <- .meta_power_columns(data, list(d = d))$d
  }
  if (!is.numeric(d) || !is.null(dim(d)) ||
      !length(d) %in% c(1L, nrow(data)))
    stop("'d' must be NULL, a numeric scalar, a numeric vector of nrow(data), or a column name.",
         call. = FALSE)
  d <- rep(d, length.out = nrow(data))
  bad <- which(!is.finite(d))
  if (length(bad))
    stop(sprintf("Row %d: 'd' must be finite and nonmissing.", bad[1L]), call. = FALSE)
  data$power <- .meta_power_rows(values, function(i) {
    fit <- power.hct.test(d = d[i], power = NULL, n2 = values$n2[i],
      n.ratio = values$n1[i] / values$n2[i],
      var.ratio = (values$sd1[i] / values$sd2[i])^2,
      alpha = alpha, alternative = alternative, null.d = null.d,
      margin = margin, method = method, claim.basis = claim.basis,
      test.distribution = test.distribution, ceil.n = FALSE, verbose = 0)
    if (length(fit$power) != 1L || !is.finite(fit$power) ||
        fit$power < 0 || fit$power > 1)
      stop("The power function did not return a probability between zero and one.")
    fit$power
  })
  data
}
