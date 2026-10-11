# Internal Borenstein et al. (2009), Chapter 24, equations 24.4-24.5, p.230.
# Equal-weight composite: mean(yi), with sampling variance sum(V)/k^2.
# rho is a WORKING COMMON SAMPLING-ERROR correlation, not raw outcome cor.
# No composite t degrees-of-freedom rule is assumed; inference is normal.
.validate_pool_rho <- function(rho, k) {
  if (is.null(rho)) return(invisible(NULL))
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) ||
      rho > 1 || rho < if (k > 1L) -1 / (k - 1L) else -1)
    stop("'pool.rho' must be one finite common sampling-error correlation in ",
         if (k > 1L) sprintf("[%.6g, 1].", -1 / (k - 1L)) else "[-1, 1].",
         call. = FALSE)
  invisible(rho)
}

# Validate a supplied sampling-error correlation matrix in outcome order.
.validate_pool_cor <- function(R, outcomes) {
  if (is.null(R)) return(NULL)
  k <- length(outcomes)
  if (!is.matrix(R) || !is.numeric(R) ||
      !identical(dim(R), c(as.integer(k), as.integer(k))) || any(!is.finite(R)))
    stop("'pool.cor' must be a finite k by k numeric matrix.", call. = FALSE)
  named <- !is.null(rownames(R)) || !is.null(colnames(R))
  if (named) {
    if (is.null(rownames(R)) || is.null(colnames(R)) ||
        anyDuplicated(rownames(R)) || anyDuplicated(colnames(R)) ||
        !setequal(rownames(R), outcomes) || !setequal(colnames(R), outcomes))
      stop("'pool.cor' row and column names must match all outcomes.", call. = FALSE)
    R <- R[outcomes, outcomes, drop = FALSE]
  }
  if (max(abs(R - t(R))) > 1e-8 || any(abs(diag(R) - 1) > 1e-8) ||
      any(abs(R) > 1 + 1e-8) ||
      min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) < -1e-8)
    stop("'pool.cor' must be symmetric, positive semidefinite, with unit diagonal.",
         call. = FALSE)
  R <- (R + t(R)) / 2
  dimnames(R) <- list(outcomes, outcomes)
  R
}

.borenstein_smd_pool <- function(x, pool.rho, alpha, pool.cor = NULL) {
  if (is.null(pool.rho) && is.null(pool.cor)) return(NULL)
  if (!is.null(pool.rho) && !is.null(pool.cor))
    stop("Supply only one of pool.rho and pool.cor.", call. = FALSE)
  k <- nrow(x)
  .validate_pool_rho(pool.rho, k)
  if (k < 1L || any(!is.finite(x$yi)) || any(!is.finite(x$vi) | x$vi <= 0))
    stop("SMD pooling requires finite effects and positive finite sampling variances.",
         call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) stop("Pooling alpha must lie in (0, 1).", call. = FALSE)
  if (is.null(pool.cor)) {
    R <- matrix(pool.rho, k, k)
    diag(R) <- 1
  } else {
    # Stored names allow safe subsetting/reordering of already fitted tables.
    if (!is.null(rownames(pool.cor)) && !is.null(colnames(pool.cor)) &&
        all(x$outcome %in% rownames(pool.cor)) && all(x$outcome %in% colnames(pool.cor)))
      pool.cor <- pool.cor[x$outcome, x$outcome, drop = FALSE]
    R <- .validate_pool_cor(pool.cor, x$outcome)
  }
  V <- R * outer(sqrt(x$vi), sqrt(x$vi))
  diag(V) <- x$vi
  variance <- sum(V) / k^2
  if (!is.finite(variance) || variance <= 0)
    stop("The assumed correlation gives a nonpositive pooled variance.", call. = FALSE)
  estimate <- mean(x$yi)
  se <- sqrt(variance)
  critical <- stats::qnorm(1 - alpha / 2)
  data.frame(effect = "smd", estimate = estimate, variance = variance, se = se,
             lcl = estimate - critical * se, ucl = estimate + critical * se,
             p = 2 * stats::pnorm(-abs(estimate / se)), k = k,
             alpha = alpha, pool.rho = if (is.null(pool.rho)) NA_real_ else pool.rho, pool.test = "z",
             pool.method = "Borenstein equal-weight composite", row.names = NULL)
}
