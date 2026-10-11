# Unweighted HC derivations: Bulus (2026), Research Synthesis Methods, Appendix A1.
# Observation-weight WLS support is an extension validated separately against lm/sandwich.
# Observation weights for untrimmed two-group weighted least squares.
# These are distinct from the mixture-component argument named 'weight'.
.prepare_observation_weights <- function(x, y, weights, na.rm) {
  if (is.list(weights)) {
    if (length(weights) != 2L) {
      stop("Supply weights as a concatenated vector or list(x = ..., y = ...).",
           call. = FALSE)
    }
    if (!is.null(names(weights))) {
      if (!all(c("x", "y") %in% names(weights))) {
        stop("Named weights must contain 'x' and 'y'.", call. = FALSE)
      }
      wx <- weights$x
      wy <- weights$y
    } else {
      wx <- weights[[1L]]
      wy <- weights[[2L]]
    }
  } else {
    if (!is.numeric(weights) || length(weights) != length(x) + length(y)) {
      stop("The weights vector must have length length(x) + length(y).",
           call. = FALSE)
    }
    wx <- weights[seq_along(x)]
    wy <- weights[length(x) + seq_along(y)]
  }
  prepare_group <- function(values, w) {
    if (!is.numeric(values) || !is.numeric(w) || !is.null(dim(w)) ||
        length(values) != length(w)) {
      stop("Each group needs numeric observations and a matching weights vector.",
           call. = FALSE)
    }
    missing <- is.na(values) | is.na(w)
    if (any(missing)) {
      if (!na.rm) stop("Missing observations or weights; use na.rm = TRUE.",
                        call. = FALSE)
      values <- values[!missing]
      w <- w[!missing]
    }
    if (any(!is.finite(values)) || any(!is.finite(w)) || any(w < 0)) {
      stop("Observations and weights must be finite; weights must be nonnegative.",
           call. = FALSE)
    }
    positive <- w > 0
    values <- values[positive]
    w <- w[positive]
    if (length(values) < 2L) {
      stop("Each group needs at least two positive-weight observations.",
           call. = FALSE)
    }
    total <- sum(w)
    mean <- sum(w * values) / total
    residuals <- values - mean
    effective_n <- total^2 / sum(w^2)
    variance <- sum(w * residuals^2) / (total - sum(w^2) / total)
    list(values = values, weights = w, mean = mean, sd = sqrt(variance),
         residuals = residuals, total = total, effective_n = effective_n,
         leverage = w / total, n = length(values))
  }
  list(x = prepare_group(x, wx), y = prepare_group(y, wy))
}

.weighted_hc_variance <- function(groups, method) {
  a <- groups$x
  b <- groups$y
  n <- a$n + b$n
  if (method %in% c("constant", "student")) {
    residual_variance <- (sum(a$weights * a$residuals^2) +
                            sum(b$weights * b$residuals^2)) / (n - 2)
    return(c(residual_variance / a$total, residual_variance / b$total))
  }
  if (method %in% c("sd1", "sd2")) {
    variance <- if (method == "sd1") a$sd^2 else b$sd^2
    return(c(variance / a$effective_n, variance / b$effective_n))
  }
  h <- c(a$leverage, b$leverage)
  score <- c(a$weights * a$residuals, b$weights * b$residuals)
  exponent <- switch(method,
    hc0 = 0, hc1 = 0, hc2 = 1, welch = 1, hc3 = 2,
    hc4 = pmin(4, n * h / 2),
    hc4m = pmin(1, n * h / 2) + pmin(1.5, n * h / 2),
    hc5 = pmin(n * h / 2, pmax(4, n * 0.7 * max(h) / 2)) / 2,
    stop("Unknown weighted HC method.", call. = FALSE)
  )
  omega <- score^2 / (1 - h)^exponent
  if (method == "hc1") omega <- omega * n / (n - 2)
  c(sum(omega[seq_len(a$n)]) / a$total^2,
    sum(omega[a$n + seq_len(b$n)]) / b$total^2)
}
