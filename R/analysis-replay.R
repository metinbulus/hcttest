.capture_analysis_rng <- function() {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    invisible(stats::runif(1))
  list(kind = RNGkind(), seed = get(".Random.seed", envir = .GlobalEnv))
}

.record_analysis <- function(fun, call, frame, overrides = list(), rng = NULL) {
  form <- formals(get(fun, envir = environment(.record_analysis)))
  args <- list()
  for (name in setdiff(names(form), "...")) {
    # Summary inputs without defaults are absent for raw-data calls.
    absent <- eval(substitute(missing(ARG), list(ARG = as.name(name))), frame)
    if (absent && identical(form[[name]], quote(expr = ))) next
    args[name] <- if (name %in% names(overrides)) overrides[name] else
      list(get(name, envir = frame, inherits = FALSE))
  }
  if ("..." %in% names(form)) args <- c(args, eval(quote(list(...)), frame))
  if (is.null(rng)) rng <- .capture_analysis_rng()
  list(fun = fun, call = call, args = args, rng = rng)
}

.replay_analysis <- function(x, alpha) {
  analysis <- attr(x, "analysis")
  if (is.null(analysis))
    stop("Recreate this object to retain the original data and analysis arguments.")
  args <- analysis$args
  if (analysis$fun == "robust.smd") args$x$parms$alpha <- alpha else args$alpha <- alpha
  if ("verbose" %in% names(args)) args$verbose <- FALSE
  old.kind <- RNGkind()
  had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old.kind))
    if (had.seed) assign(".Random.seed", old.seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  do.call(RNGkind, as.list(analysis$rng$kind))
  assign(".Random.seed", analysis$rng$seed, envir = .GlobalEnv)
  do.call(get(analysis$fun, envir = environment(.replay_analysis)), args)
}

.analysis_alpha <- function(x) {
  analysis <- attr(x, "analysis")
  if (is.null(analysis)) stop("Recreate this object to retain its original analysis inputs.")
  if (analysis$fun == "robust.smd") analysis$args$x$parms$alpha else analysis$args$alpha
}

.plot_replayed_intervals <- function(x, effect, boot.ci) {
  if (inherits(x, "multi.robust.smd")) {
    required <- paste0(c("lcl.equivalence.", "ucl.equivalence.", "lcl.outer.", "ucl.outer."), effect)
    if (!all(required %in% names(x)))
      stop("Recreate the multi.robust.smd object to retain equivalence intervals.")
    return(data.frame(lower.inner = x[[required[1]]], upper.inner = x[[required[2]]],
                      lower.outer = x[[required[3]]], upper.outer = x[[required[4]]],
                      inner.level = x$equivalence.level, outer.level = 1 - x$alpha))
  }
  extract <- function(result) {
    if (inherits(result, "hct.test")) return(cbind(result$test$lcl, result$test$ucl))
    if (inherits(result, "robust.smd")) return(cbind(result$es$lcl, result$es$ucl))
    if (inherits(result, "boot.hct.test")) {
      index <- match(if (effect == "raw") "Raw Scores" else "Standardized", result$est.type)
      suffix <- if (boot.ci == "percentile") "perc" else "bca"
      return(cbind(result[[paste0("lcl.boot.", suffix)]][index],
                   result[[paste0("ucl.boot.", suffix)]][index]))
    }
    index <- if (inherits(x, "meta.robust.smd")) match(x$row, result$row) else
      match(x$outcome, result$outcome)
    if (anyNA(index)) stop("Replayed analysis rows do not match the original object.")
    if (effect == "raw") cbind(result$lcl.delta[index], result$ucl.delta[index]) else
      cbind(result$lcl[index], result$ucl[index])
  }
  original.alpha <- .analysis_alpha(x)
  result95 <- if (isTRUE(all.equal(original.alpha, 0.05))) x else .replay_analysis(x, 0.05)
  result99 <- if (isTRUE(all.equal(original.alpha, 0.01))) x else .replay_analysis(x, 0.01)
  ci95 <- extract(result95)
  ci99 <- extract(result99)
  data.frame(lower95 = ci95[, 1], upper95 = ci95[, 2],
             lower99 = ci99[, 1], upper99 = ci99[, 2],
             lower.inner = ci95[, 1], upper.inner = ci95[, 2],
             lower.outer = ci99[, 1], upper.outer = ci99[, 2],
             inner.level = 0.95, outer.level = 0.99)
}
