# Primary reference: Bulus (2026), Research Synthesis Methods,
# Implications for Meta-Analysis and Table 3. Batch interfaces are convenience extensions.
# Study effects use the existing core; pooling delegates to metafor.
.batch_smd_column <- function(data, column, argument, numeric = FALSE) {
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !column %in% names(data)) {
    stop(sprintf("'%s' must name a column in 'data'.", argument), call. = FALSE)
  }
  value <- data[[column]]
  if (numeric && (!is.numeric(value) || !is.null(dim(value)))) {
    stop(sprintf("Column '%s' must be a numeric vector.", column), call. = FALSE)
  }
  value
}

.batch_smd_settings <- function(dots, reserved) {
  if (length(dots) && (is.null(names(dots)) || any(is.na(names(dots))) ||
                      any(names(dots) == "") || anyDuplicated(names(dots)))) {
    stop("Additional hct.test arguments must have unique names.", call. = FALSE)
  }
  conflict <- intersect(names(dots), reserved)
  if (length(conflict)) {
    stop(paste("Use the explicit batch arguments instead of supplying:",
               paste(conflict, collapse = ", ")), call. = FALSE)
  }
  dots
}

.batch_smd_result <- function(fit, ss.adjust, ci, scale) {
  effect <- robust.smd(fit, ss.adjust = ss.adjust, ci = ci,
                       scale = scale, verbose = FALSE)
  data.frame(m1 = fit$parms$m1, m2 = fit$parms$m2,
             sd1 = fit$parms$sd1, sd2 = fit$parms$sd2,
             n1 = fit$parms$n1, n2 = fit$parms$n2,
             delta = fit$test$delta, se.delta = fit$test$se,
             t = fit$test$t, df = fit$test$df, p = fit$test$p,
             lcl.delta = fit$test$lcl, ucl.delta = fit$test$ucl,
             yi = effect$es$d, vi = effect$es$se^2,
             se = effect$es$se, lcl = effect$es$lcl, ucl = effect$es$ucl,
             ci = effect$parms$ci, alpha = fit$parms$alpha,
             two.tailed = fit$parms$two.tailed,
             method = fit$parms$method, trim = fit$parms$trim,
             ss.adjust = ss.adjust, scale = scale,
             weighted = isTRUE(fit$parms$weighted),
             n_eff1 = if (isTRUE(fit$parms$weighted)) fit$parms$n_eff1 else NA_real_,
             n_eff2 = if (isTRUE(fit$parms$weighted)) fit$parms$n_eff2 else NA_real_)
}


# Pool using metafor rather than introducing a new heterogeneity estimator.
.meta_smd_pool <- function(x, pool.method, pool.test, alpha) {
  if (nrow(x) < 2L) return(NULL)
  if (any(!is.finite(x$vi) | x$vi <= 0) ||
      any(!is.finite(x$se.delta) | x$se.delta <= 0))
    stop("Pooling requires positive finite sampling variances.", call. = FALSE)
  models <- list(
    smd = metafor::rma.uni(yi = x$yi, vi = x$vi, method = pool.method,
                           test = pool.test, level = 100 * (1 - alpha)),
    raw = metafor::rma.uni(yi = x$delta, vi = x$se.delta^2, method = pool.method,
                           test = pool.test, level = 100 * (1 - alpha)))
  summary <- do.call(rbind, lapply(names(models), function(effect) {
    fit <- models[[effect]]
    data.frame(effect = effect, estimate = as.numeric(fit$b), se = fit$se,
               lcl = fit$ci.lb, ucl = fit$ci.ub, p = fit$pval,
               tau2 = fit$tau2, I2 = fit$I2, Q = fit$QE, Q.p = fit$QEp,
               k = fit$k, pool.method = pool.method, pool.test = pool.test,
               row.names = NULL)
  }))
  list(summary = summary, models = models)
}

# Equivalence intervals use the fitted SMD CI method; raw inference stays central t.
.multi_smd_equivalence <- function(x, margins, equivalence.effect, ci, fits) {
  if (any(!is.finite(x$alpha) | x$alpha <= 0 | x$alpha >= 0.5))
    stop("Equivalence alpha must be between zero and 0.5.", call. = FALSE)
  if (!is.null(margins) && (!is.numeric(margins) || length(margins) != 2L ||
      any(!is.finite(margins)) || margins[1] >= margins[2]))
    stop("'margins' must contain two increasing finite equivalence bounds.", call. = FALSE)
  x$lcl.equivalence.raw <- x$delta - stats::qt(1 - x$alpha, x$df) * x$se.delta
  x$ucl.equivalence.raw <- x$delta + stats::qt(1 - x$alpha, x$df) * x$se.delta
  x$lcl.outer.raw <- x$delta - stats::qt(1 - x$alpha / 2, x$df) * x$se.delta
  x$ucl.outer.raw <- x$delta + stats::qt(1 - x$alpha / 2, x$df) * x$se.delta
  # Reuse robust.smd, including its standardizer, correction and trim calibration.
  intervals <- lapply(seq_len(nrow(x)), function(i) {
    fit <- fits[[i]]
    fit$parms$two.tailed <- TRUE
    fit$parms$alpha <- 2 * x$alpha[i]
    inner <- robust.smd(fit, ss.adjust = x$ss.adjust[i], ci = ci,
                        scale = x$scale[i], verbose = FALSE)$es
    fit$parms$alpha <- x$alpha[i]
    outer <- robust.smd(fit, ss.adjust = x$ss.adjust[i], ci = ci,
                        scale = x$scale[i], verbose = FALSE)$es
    # The estimate for an observed t of one is the positive multiplier c
    # in d = c * lambda. This avoids dividing by an observed zero estimate.
    fit$test$delta <- fit$test$se
    multiplier <- robust.smd(fit, ss.adjust = x$ss.adjust[i], ci = FALSE,
                             scale = x$scale[i], verbose = FALSE)$es$d
    c(inner$lcl, inner$ucl, outer$lcl, outer$ucl, multiplier)
  })
  intervals <- do.call(rbind, intervals)
  x$lcl.equivalence.smd <- intervals[, 1]
  x$ucl.equivalence.smd <- intervals[, 2]
  x$lcl.outer.smd <- intervals[, 3]
  x$ucl.outer.smd <- intervals[, 4]
  estimate <- if (equivalence.effect == "raw") x$delta else x$yi
  se <- if (equivalence.effect == "raw") x$se.delta else x$se
  if (any(!is.finite(se) | se <= 0))
    stop("Equivalence inference requires positive finite standard errors.", call. = FALSE)
  x$equivalence.effect <- equivalence.effect
  x$equivalence.level <- 1 - 2 * x$alpha
  x$lcl.equivalence <- x[[paste0("lcl.equivalence.", equivalence.effect)]]
  x$ucl.equivalence <- x[[paste0("ucl.equivalence.", equivalence.effect)]]
  x$margin.lower <- if (is.null(margins)) NA_real_ else margins[1]
  x$margin.upper <- if (is.null(margins)) NA_real_ else margins[2]
  x$p.lower <- x$p.upper <- x$p.max <- NA_real_
  x$equivalent <- NA
  if (!is.null(margins)) {
    if (equivalence.effect == "raw" || ci == "ct") {
      x$p.lower <- stats::pt((estimate - margins[1]) / se, x$df, lower.tail = FALSE)
      x$p.upper <- stats::pt((estimate - margins[2]) / se, x$df)
    } else if (ci == "z") {
      x$p.lower <- stats::pnorm((estimate - margins[1]) / se, lower.tail = FALSE)
      x$p.upper <- stats::pnorm((estimate - margins[2]) / se)
    } else if (ci == "nct") {
      # Same noncentral-t inversion as the SMD interval, at each bound.
      x$p.lower <- suppressWarnings(stats::pt(x$t, x$df,
                         ncp = margins[1] / intervals[, 5], lower.tail = FALSE))
      x$p.upper <- suppressWarnings(stats::pt(x$t, x$df,
                         ncp = margins[2] / intervals[, 5], lower.tail = TRUE))
    }
    x$p.max <- pmax(x$p.lower, x$p.upper)
    x$equivalent <- x$p.max < x$alpha
  }
  x
}

#' Compute standardized differences for study tables or multiple outcomes
#'
#' @description Batch interfaces to hct.test and robust.smd. Each input study or outcome produces one result row.
#' @param data Nonempty data frame with unique column names.
#' @param m1,m2 Column names containing first/second group means. Each row represents a study. Values must be numeric and finite.
#' @param sd1,sd2 Column names containing first/second group SDs. Values must be numeric, finite, and nonnegative.
#' @param n1,n2 Column names containing first/second group sample sizes. Values must be integers of at least two.
#' @param study Optional column name for study identifiers; otherwise row numbers are used. Duplicate identifiers are retained.
#' @param outcomes Character vector of unique numeric outcome column names.
#' @param group Column name identifying the two groups.
#' @param weights NULL, the name of a numeric observation-weight column, or a numeric vector aligned with data rows. Available only in multi.robust.smd.
#' @param method Testing variance method passed to hct.test; the default matches hct.test.
#' @param df Degrees-of-freedom convention passed to hct.test.
#' @param alpha Significance level.
#' @param two.tailed Logical; use two-sided intervals and tests.
#' @param trim Trimming proportion for raw data in multi.robust.smd; analytical weighted trimming remains disabled.
#' @param na.rm Logical; remove missing raw observations separately for each outcome. With weights, remove incomplete observation-weight pairs.
#' @param ss.adjust Logical; apply the existing small-sample correction in robust.smd.
#' @param ci Effect-size interval approximation passed to robust.smd: nct, ct, z, or FALSE. Raw-difference inference always uses central t intervals. p is the raw-difference p value.
#' @param scale Logical; apply the existing normal calibration of trimmed effect sizes.
#' @param pool Logical; pool independent study effects using metafor. FALSE returns the effect table without fitting a pooled model. At least two rows are required.
#' @param pool.method Meta-analysis estimator, default REML; alternatives include DL, PM, ML, and EE (equal-effect model). Distinct tau-squared estimates can change the pooled effect.
#' @param pool.test Pooled-inference method passed to metafor, default knha (Knapp-Hartung); z, t, and adhoc are also supported.
#' @param margins NULL or two increasing finite equivalence bounds for multi.robust.smd. Bounds are in equivalence.effect units; NULL leaves equivalence p values undefined.
#' @param equivalence.effect Effect tested against the bounds: smd (default) or raw. Use raw for outcomes already divided by fixed reference SDs.
#' @param ... Additional uniquely named hct.test arguments. Arguments controlled by the batch interface cannot be overridden here.
#'
#' @details The difference is first group minus second group. In multi.robust.smd, group ordering follows the factor levels (or the factor conversion of the group column). Study effects delegate to hct.test and robust.smd; no new within-study variance estimator is introduced. Pooling uses metafor.
#'
#' meta.robust.smd computes study effects and, by default, pools them through metafor::rma.uni. The pooled attribute has raw and smd summary rows; meta.models retains the corresponding metafor fits. The weight column contains unnormalized SMD inverse-variance weights 1/(vi + tau2). Pooling assumes independent rows; correlated outcomes require a separate model. pool = FALSE skips pooling and omits the weight column. It supports untrimmed scalar summaries only. Observation weights, raw data, mixture-component list columns, and analytical trimming are not accepted by this interface. Raw group means and SDs alone do not determine empirical trimmed summaries.
#'
#' multi.robust.smd performs a separate analysis for each outcome, retaining the requested order. With na.rm = TRUE, sample sizes can differ across outcomes. It returns individual outcome estimates without a pooled effect. No multiplicity correction is applied.
#'
#' SMD equivalence inherits ci: ct and z use estimate/SE TOST, while nct uses the same noncentral-t inversion, standardizer, small-sample correction and trim calibration as robust.smd. Raw-difference equivalence retains central-t inference. p.lower tests the lower bound, p.upper tests the upper bound, and p.max is their maximum. Equivalence requires p.max < alpha. lcl.equivalence and ucl.equivalence form the corresponding 100*(1-2*alpha) percent CI (90 percent at alpha = 0.05), irrespective of two.tailed. A 100*(1-alpha) percent outer interval supplies ribbon fading. Both raw and SMD intervals are retained for plotting, but p values test only equivalence.effect. With ci = FALSE, SMD equivalence limits and p values are NA. No multiplicity correction is applied. Weighted inference conditions on supplied weights.
#'
#' The yi column is the d estimate returned by robust.smd (small-sample corrected when ss.adjust = TRUE). vi is its squared reported SE. These inherit the existing effect-size approximations, including the limitations for weighted effect-size SEs. They are not newly validated meta-analysis variance formulas. Undefined standardized estimates or variances raise an error. Other errors identify the study row or outcome; no unsuccessful row is silently removed. Printing is suppressed during each analysis.
#'
#' @return Meta-analysis pooling results are available as attr(x, "pooled"), with pooling options retained in attr(x, "pool.settings"). Multi-outcome results have no pooled effect. Meta-analysis fitted models are available as attr(x, "meta.models"). Multi-outcome tables include bounds, p.lower, p.upper, p.max, equivalent, and equivalence CI columns. A data frame of class meta.robust.smd or multi.robust.smd, with a plot method, containing identifiers, analyzed group moments and sizes, raw difference, raw SE, t statistic, df, p value, raw interval, yi, vi, effect-size SE and interval, method, trim, correction/scaling indicators, and weight information. Weighted raw-data rows include Kish effective sample sizes. Raw trimmed sample sizes are the retained counts used by the core. The analysis attribute retains the original input data frame and evaluated arguments, including whether the source function was meta.robust.smd or multi.robust.smd.
#'
#' @examples studies <- data.frame(study = c("A", "B"),
#'                       m1 = c(1, 2), m2 = c(0.5, 1.5),
#'                       sd1 = c(1, 1.2), sd2 = c(1.1, 1.5),
#'                       n1 = c(30, 40), n2 = c(35, 45))
#' meta.robust.smd(studies, study = "study", method = "hc3", ci = "ct")
#'
#' set.seed(1)
#' dat <- data.frame(group = factor(rep(c("A", "B"), each = 30)),
#'                   anxiety = rnorm(60), wellbeing = rnorm(60),
#'                   w = runif(60, 0.5, 2))
#' multi.robust.smd(dat, outcomes = c("anxiety", "wellbeing"),
#'                  group = "group", weights = "w", method = "hc3", ci = "ct")
#'
#' @references Bulus, M. (2026). Robust standard errors and confidence intervals for standardized mean difference. \emph{Research Synthesis Methods}. Advance online publication.
#'
#' Lakens, D. (2017). Equivalence tests: A practical primer for t tests, correlations, and meta-analyses. Social Psychological and Personality Science, 8(4), 355-362. doi:10.1177/1948550617697177. Pooling uses metafor; batch helpers do not establish agreement with all values in the earlier proof study table.
#'
#'
#' @rdname meta.robust.smd
#' @export meta.robust.smd
meta.robust.smd <- function(data, m1 = "m1", m2 = "m2",
                            sd1 = "sd1", sd2 = "sd2", n1 = "n1", n2 = "n2",
                            study = NULL, method = "constant", df = "regression",
                            alpha = 0.05, two.tailed = TRUE,
                            ss.adjust = FALSE, ci = "nct", scale = TRUE,
                            pool = TRUE, pool.method = "REML", pool.test = "knha", ...) {
  analysis <- .record_analysis("meta.robust.smd", match.call(), environment())
  if (!is.logical(pool) || length(pool) != 1L || is.na(pool))
    stop("'pool' must be TRUE or FALSE.", call. = FALSE)
  pool.method <- match.arg(pool.method, c("REML", "DL", "PM", "ML", "HE", "HS", "HSk",
                                        "SJ", "EB", "PMM", "EE", "FE"))
  pool.test <- match.arg(pool.test, c("knha", "z", "t", "adhoc"))
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be between zero and one.", call. = FALSE)
  if (!is.data.frame(data) || nrow(data) < 1L || anyDuplicated(names(data))) {
    stop("'data' must be a nonempty data frame with unique column names.", call. = FALSE)
  }
  columns <- list(m1 = m1, m2 = m2, sd1 = sd1, sd2 = sd2, n1 = n1, n2 = n2)
  values <- lapply(names(columns), function(name) {
    .batch_smd_column(data, columns[[name]], name, numeric = TRUE)
  })
  names(values) <- names(columns)
  id <- if (is.null(study)) seq_len(nrow(data)) else {
    .batch_smd_column(data, study, "study")
  }
  dots <- .batch_smd_settings(list(...),
    c("x", "y", names(columns), "weights", "method", "df", "alpha",
      "two.tailed", "ci", "verbose", "trim", "na.rm"))
  results <- lapply(seq_len(nrow(data)), function(i) {
    tryCatch({
      row <- lapply(values, function(value) value[i])
      if (any(!is.finite(unlist(row))) || row$sd1 < 0 || row$sd2 < 0 ||
          row$n1 < 2 || row$n2 < 2 || row$n1 != floor(row$n1) ||
          row$n2 != floor(row$n2)) {
        stop("Means/SDs must be finite, SDs nonnegative, and sample sizes integers >= 2.")
      }
      fit <- do.call(hct.test, c(row, list(method = method, df = df,
                     alpha = alpha, two.tailed = two.tailed,
                     ci = "ct", trim = 0, verbose = FALSE), dots))
      result <- .batch_smd_result(fit, ss.adjust, ci, scale)
      if (!is.finite(result$yi) || !is.finite(result$vi)) {
        stop("The standardized estimate or its variance is undefined.")
      }
      data.frame(row = i, study = id[i], result, check.names = FALSE)
    }, error = function(e) {
      stop(sprintf("Study row %d: %s", i, conditionMessage(e)), call. = FALSE)
    })
  })
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  class(out) <- c("meta.robust.smd", "data.frame")
  attr(out, "analysis") <- analysis
  attr(out, "pool.settings") <- list(pool = pool, pool.method = pool.method,
                                     pool.test = pool.test, alpha = alpha)
  if (pool) {
    pooled <- .meta_smd_pool(out, pool.method, pool.test, alpha)
    if (is.null(pooled)) warning("At least two rows are required for pooling; returning the study table only.", call. = FALSE)
    attr(out, "pooled") <- pooled$summary
    attr(out, "meta.models") <- pooled$models
    out$weight <- if (is.null(pooled)) 1 / out$vi else 1 / (out$vi + pooled$models$smd$tau2)
  }
  out
}

#' @rdname meta.robust.smd
#' @export multi.robust.smd
multi.robust.smd <- function(data, outcomes, group, weights = NULL,
                             method = "constant", df = "regression",
                             alpha = 0.05, two.tailed = TRUE, trim = 0,
                             na.rm = FALSE, ss.adjust = FALSE,
                             ci = "nct", scale = TRUE, margins = NULL,
                             equivalence.effect = c("smd", "raw"),
                             ...) {
  analysis <- .record_analysis("multi.robust.smd", match.call(), environment())
  if (any(c("pool.rho", "pool.cor") %in% names(list(...))))
    stop("Pooling has been removed from multi.robust.smd; omit pool.rho and pool.cor.", call. = FALSE)
  equivalence.effect <- match.arg(equivalence.effect)
  ci <- match.arg(as.character(ci), c("nct", "ct", "z", "FALSE"))
  if ("equivalence.ci" %in% names(list(...)))
    stop("The equivalence.ci argument has been removed; specify ci instead.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 0.5)
    stop("Equivalence alpha must be between zero and 0.5.", call. = FALSE)
  if (!is.null(margins) && (!is.numeric(margins) || length(margins) != 2L ||
      any(!is.finite(margins)) || margins[1] >= margins[2]))
    stop("'margins' must contain two increasing finite equivalence bounds.", call. = FALSE)
  if (!is.data.frame(data) || nrow(data) < 1L || anyDuplicated(names(data))) {
    stop("'data' must be a nonempty data frame with unique column names.", call. = FALSE)
  }
  if (!is.character(outcomes) || !length(outcomes) || anyNA(outcomes) ||
      anyDuplicated(outcomes)) {
    stop("'outcomes' must contain unique outcome column names.", call. = FALSE)
  }
  labels <- .batch_smd_column(data, group, "group")
  for (outcome in outcomes) .batch_smd_column(data, outcome, "outcomes", TRUE)
  w <- if (is.character(weights)) {
    .batch_smd_column(data, weights, "weights", TRUE)
  } else weights
  if (!is.null(w) && (!is.numeric(w) || !is.null(dim(w)) || length(w) != nrow(data))) {
    stop("'weights' must name a numeric column or be a vector aligned with data.", call. = FALSE)
  }
  dots <- .batch_smd_settings(list(...),
    c("x", "y", "m1", "m2", "sd1", "sd2", "n1", "n2", "formula", "data",
      "weights", "method", "df", "alpha", "two.tailed", "ci", "trim", "na.rm", "verbose"))
  results <- lapply(outcomes, function(outcome) {
    tryCatch({
      # A fixed temporary formula supports outcome names containing spaces.
      current <- data.frame(.outcome = data[[outcome]], .group = labels)
      keep <- !is.na(labels) & !is.na(data[[outcome]])
      if (!is.null(w)) {
        keep <- keep & !is.na(w) & w > 0
        keep[is.na(keep)] <- FALSE
      }
      group_levels <- levels(droplevels(as.factor(labels[keep])))
      fit <- do.call(hct.test, c(list(formula = .outcome ~ .group, data = current,
                       weights = w, method = method, df = df, alpha = alpha,
                       two.tailed = two.tailed, ci = "ct", trim = trim,
                       na.rm = na.rm, verbose = FALSE), dots))
      result <- .batch_smd_result(fit, ss.adjust, ci, scale)
      if (!is.finite(result$yi) || !is.finite(result$vi)) {
        stop("The standardized estimate or its variance is undefined.")
      }
      row <- data.frame(outcome = outcome, group1 = group_levels[1], group2 = group_levels[2],
                        result, check.names = FALSE)
      attr(row, "equivalence.fit") <- fit
      row
    }, error = function(e) {
      stop(sprintf("Outcome '%s': %s", outcome, conditionMessage(e)), call. = FALSE)
    })
  })
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out <- .multi_smd_equivalence(out, margins, equivalence.effect, ci,
                                lapply(results, attr, "equivalence.fit"))
  attr(out, "equivalence.fit") <- NULL
  class(out) <- c("multi.robust.smd", "data.frame")
  attr(out, "analysis") <- analysis
  out
}
