# Weight-mass trimming with partial boundary weights. This helper estimates
# location and descriptive scale; it does not define an analytical HC test.
.weighted_trimmed_group <- function(group, trim) {
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim < 0 || trim > 0.45) {
    stop("'trim' must be a finite number in [0, 0.45].", call. = FALSE)
  }
  ordering <- order(group$values)
  values <- group$values[ordering]
  w <- group$weights[ordering]
  # Work with proportions so the operation is invariant to weight scaling.
  mass <- w / sum(w)
  right <- cumsum(mass)
  left <- c(0, head(right, -1L))
  # Each observation occupies an interval in cumulative-weight space.
  # Retain its overlap with [trim, 1 - trim].
  retained <- pmax(0, pmin(right, 1 - trim) - pmax(left, trim))
  # Ignore negligible overlaps created by floating-point boundary roundoff.
  retained[retained < 16 * .Machine$double.eps] <- 0
  positive <- retained > 0
  if (!any(positive)) stop("Trimming leaves no retained weight.", call. = FALSE)
  trimmed_mean <- sum(retained * values) / sum(retained)
  lower <- min(values[positive])
  upper <- max(values[positive])
  winsorized <- pmin(upper, pmax(lower, values))
  winsorized_mean <- sum(mass * winsorized)
  # Match the descriptive weighted-SD definition of the untrimmed path.
  winsorized_variance <- sum(mass * (winsorized - winsorized_mean)^2) /
    (1 - sum(mass^2))
  list(mean = trimmed_mean, winsorized_sd = sqrt(winsorized_variance),
       # Preserve the existing population-style Yuen SD normalization.
       # No effective-df correction is used to claim analytical inference.
       adjusted_sd = sqrt(winsorized_variance / (1 - 2 * trim)),
       effective_n = group$effective_n,
       values = values, original_mass = mass, retained_mass = retained,
       lower = lower, upper = upper, winsorized = winsorized)
}

.weighted_trimmed_boot_statistic <- function(data, outcome, group, weights,
                                            trim, standardizer, scale,
                                            calibration = NULL) {
  levels <- levels(data[[group]])
  first <- data[[group]] == levels[1]
  second <- data[[group]] == levels[2]
  prepared <- .prepare_observation_weights(
    data[[outcome]][first], data[[outcome]][second],
    weights = list(x = data[[weights]][first], y = data[[weights]][second]),
    na.rm = FALSE)
  a <- .weighted_trimmed_group(prepared$x, trim)
  b <- .weighted_trimmed_group(prepared$y, trim)
  delta <- a$mean - b$mean
  denominator <- switch(standardizer,
    sd1 = a$adjusted_sd,
    sd2 = b$adjusted_sd,
    sdp = .pooled_sd(a$adjusted_sd, b$adjusted_sd,
                     a$effective_n, b$effective_n),
    stop("Unknown weighted effect-size standardizer.", call. = FALSE))
  if (!is.finite(denominator) || denominator <= 0) {
    stop("The weighted Winsorized standardizer is zero or undefined in a bootstrap sample.",
         call. = FALSE)
  }
  d <- delta / denominator
  if (scale) {
    if (is.null(calibration)) {
      calibration <- winsor.pop.sd(mu = 0, sigma = 1, trim = trim,
                                   distribution = "norm") / sqrt(1 - 2 * trim)
    }
    d <- d * calibration
  }
  c(delta, d)
}
