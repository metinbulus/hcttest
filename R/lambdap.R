# CRAN release: import these Lambda-prime functions from the pwrss package
# instead of maintaining local definitions. They are internal in this draft.
# Adapted from user-supplied metinbulus/pwrss R/lambdap.R:
# https://github.com/metinbulus/pwrss/blob/main/R/lambdap.R
# The supplied density functions are in lambdap-density.R, with the approved
# correction to log-density handling. These functions use stats only.

# CRAN release: import plambdap from pwrss; keep internal in this draft.
plambdap <- function(q, df, ncp, lower.tail = TRUE, log.p = FALSE) {
    suppressWarnings(stats::pt(q = ncp, df = df, ncp = q,
                              lower.tail = !lower.tail, log.p = log.p))
}

.qlambdap <- function(p, df, ncp, lower.tail = TRUE, log.p = FALSE) {
    if (log.p) {
        if (p == -Inf) return(ifelse(lower.tail, -Inf, Inf))
        if (p == 0)    return(ifelse(lower.tail, Inf, -Inf))
        if (p > 0)    return(NA)
    } else {
        if (p == 1) return(ifelse(lower.tail, Inf, -Inf))
        if (p == 0) return(ifelse(lower.tail, -Inf, Inf))
        if (p < 0 || p > 1) return(NA)
    }

    if (lower.tail) {
        zerof <- function(q) plambdap(q, df = df, ncp = ncp,
                                     lower.tail = lower.tail, log.p = log.p) - p
    } else {
        zerof <- function(q) p - plambdap(q, df = df, ncp = ncp,
                                         lower.tail = lower.tail, log.p = log.p)
    }
    zmax <- 2 * max(stats::qnorm(p, lower.tail = TRUE, log.p = log.p),
                    stats::qnorm(p, lower.tail = FALSE, log.p = log.p))
    flim <- c(min(-1, ncp - zmax), max(1, ncp + zmax))
    for (i in seq(2)) {
        while (sign(zerof(flim[i])) == ifelse(i == 1, 1, -1) &&
               abs(flim[i]) < 1e8) {
            flim[i] <- 2 * flim[i]
        }
    }
    stats::uniroot(zerof, flim, tol = 1e-12)$root
}

# CRAN release: import qlambdap from pwrss; keep internal in this draft.
qlambdap <- Vectorize(.qlambdap, vectorize.args = c("p", "df", "ncp"),
                     SIMPLIFY = TRUE)

# CRAN release: import rlambdap from pwrss; keep internal in this draft.
rlambdap <- function(n, df, ncp) {
    stats::rnorm(n) + ncp * sqrt(stats::rchisq(n, df = df) / df)
}
