# Internal helper for untrimmed mixture-normal summary statistics.
.mixture_mean_sd <- function(mu, sigma, weight) {
  valid_numeric <- function(x) {
    is.numeric(x) && is.null(dim(x)) && length(x) > 0L &&
      all(is.finite(x))
  }
  if (!valid_numeric(mu) || !valid_numeric(sigma) || !valid_numeric(weight)) {
    stop("Mixture means, SDs, and weights must be finite numeric vectors.",
         call. = FALSE)
  }
  if (length(weight) == 1L) weight <- c(weight, 1 - weight)
  if (any(weight < 0) || abs(sum(weight) - 1) > 1e-10) {
    stop("Mixture weights must be nonnegative and sum to one.", call. = FALSE)
  }
  # Match the scalar convenience behavior of the existing trimmed helpers.
  if (length(mu) == 1L) mu <- rep(mu, length(weight))
  if (length(sigma) == 1L) sigma <- rep(sigma, length(weight))
  if (length(mu) != length(weight) || length(sigma) != length(weight)) {
    stop("Mixture means, SDs, and weights must have matching component counts.",
         call. = FALSE)
  }
  if (any(sigma < 0)) {
    stop("Mixture component SDs must be nonnegative.", call. = FALSE)
  }
  mean_mix <- sum(weight * mu)
  # Within + between is algebraically identical to
  # sum(weight * (sigma^2 + mu^2)) - mean_mix^2, and avoids cancellation.
  variance_mix <- sum(weight * sigma^2) +
    sum(weight * (mu - mean_mix)^2)
  list(mean = mean_mix, sd = sqrt(variance_mix))
}
