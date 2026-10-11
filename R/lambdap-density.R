# CRAN release: import these Lambda-prime functions from the pwrss package
# instead of maintaining local definitions. They are internal in this draft.
# Application reference: Bulus (2026), Research Synthesis Methods,
# confidence-interval discussion. Distribution reference: Lecoutre (2007).
# Lambda-prime density: replacement for the supplied density definitions.
# plambdap() is defined in lambdap.R in this package.
# Preserve the original numerical density calculation; apply log afterwards.
.dlambdap <- function(x, df, ncp, log = FALSE) {
    gradf <- function(q) plambdap(q, df, ncp, lower.tail = TRUE, log.p = FALSE)

    x <- as.numeric(x)
    d <- stats::numericDeriv(quote(gradf(x)), "x")
    gradient <- as.numeric(attr(d, "gradient"))

    if (gradient > 0) {
        density <- gradient
    } else {
        tol <- sqrt(.Machine$double.eps)
        if (d < tol || d > (1 - tol)) {
            density <- 0
        } else {
            epsf <- function(e) {
                as.numeric(attr(
                    stats::numericDeriv(quote(gradf(x)), "x", eps = e),
                    "gradient"
                ))
            }
            epsv <- vapply(seq(0, 3), function(e) {
                vapply(seq(tol, 10 ^ -e, length.out = 1e3), epsf, numeric(1))
            }, numeric(1e3))
            epsv <- epsv[, colSums(epsv == 0) < 1e2, drop = FALSE]
            epsv <- epsv[, which.min(apply(epsv, 2, stats::sd))]
            density <- stats::median(epsv[epsv != 0 & !is.na(epsv)])
        }
    }

    if (log) base::log(density) else density
}

# CRAN release: import dlambdap from pwrss; keep internal in this draft.
dlambdap <- Vectorize(.dlambdap,
                     vectorize.args = c("x", "df", "ncp"),
                     SIMPLIFY = TRUE)
