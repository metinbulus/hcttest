.draw_hct_ribbon <- function(estimate, lower.inner, upper.inner, lower.outer, upper.outer, y, pretty.height, color, pretty.lwd) {
      ribbon.color <- color
      half.height <- pretty.height / 2
      if (upper.outer > lower.outer) {
        # One raster avoids seams between narrow rectangle slices.
        # Solid inner band; fade only between inner and outer limits.
        edges <- seq(lower.outer, upper.outer, length.out = 1201)
        mid <- (utils::head(edges, -1) + utils::tail(edges, -1)) / 2
        alpha <- rep(0.35, length(mid))
        left.tail <- mid < lower.inner
        right.tail <- mid > upper.inner
        if (lower.inner > lower.outer)
          alpha[left.tail] <- 0.35 * (mid[left.tail] - lower.outer) /
            (lower.inner - lower.outer)
        if (upper.outer > upper.inner)
          alpha[right.tail] <- 0.35 * (upper.outer - mid[right.tail]) /
            (upper.outer - upper.inner)
        alpha <- pmin(0.35, pmax(0, alpha))
        rgba <- grDevices::col2rgb(ribbon.color) / 255
        colors <- grDevices::rgb(rgba[1], rgba[2], rgba[3], alpha = alpha)
        ribbon <- grDevices::as.raster(matrix(colors, nrow = 1))
        graphics::rasterImage(ribbon, edges[1], y - half.height,
                              utils::tail(edges, 1), y + half.height,
                              interpolate = FALSE)
        # Full-height vertical boundaries, thinner/lighter than the estimate.
        graphics::segments(c(lower.inner, upper.inner), y - half.height,
                           c(lower.inner, upper.inner), y + half.height,
                           col = grDevices::adjustcolor(ribbon.color, alpha.f = 0.65),
                           lwd = pretty.lwd / 3, lend = "butt")
      }
      graphics::segments(estimate, y - half.height,
                         estimate, y + half.height,
                         col = ribbon.color, lwd = pretty.lwd, lend = "butt")
}

.as_hctplot <- function(x, effect = NULL,
                       boot.ci = c("percentile", "bca"), labels = NULL) {
  boot.ci <- match.arg(boot.ci)
  if (!is.null(effect)) effect <- match.arg(effect, c("raw", "smd"))

  if (inherits(x, "hct.test")) {
    if (is.null(effect)) effect <- "raw"
    if (effect != "raw")
      stop("Use robust.smd(x) first to choose the SMD interval/settings.")
    two.tailed <- x$parms$two.tailed
    d <- data.frame(label = "", estimate = x$test$delta,
                    lower = x$test$lcl, upper = x$test$ucl, se = x$test$se)
  } else if (inherits(x, "robust.smd")) {
    if (is.null(effect)) effect <- "smd"
    if (effect != "smd") stop("A robust.smd object contains standardized results.")
    two.tailed <- x$parms$two.tailed
    d <- data.frame(label = "", estimate = x$es$d,
                    lower = x$es$lcl, upper = x$es$ucl, se = x$es$se)
  } else if (is.data.frame(x) &&
             all(c("est.type", "est.data", "lcl.boot.perc",
                   "ucl.boot.perc", "lcl.boot.bca", "ucl.boot.bca") %in% names(x))) {
    if (is.null(effect)) effect <- "smd"
    selected <- if (effect == "raw") "Raw Scores" else "Standardized"
    i <- which(x$est.type == selected)
    if (length(i) != 1L) stop("Expected one bootstrap row for the selected effect.")
    suffix <- if (boot.ci == "percentile") "perc" else "bca"
    two.tailed <- TRUE
    d <- data.frame(label = "",
                    estimate = x$est.data[i],
                    lower = x[[paste0("lcl.boot.", suffix)]][i],
                    upper = x[[paste0("ucl.boot.", suffix)]][i], se = x$se.boot[i])
    # The point is the observed estimate, rather than the bootstrap mean.
  } else if (is.data.frame(x) &&
             all(c("yi", "lcl", "ucl", "delta", "lcl.delta", "ucl.delta") %in% names(x)) &&
             any(c("study", "outcome") %in% names(x))) {
    if (is.null(effect)) effect <- "smd"
    row.labels <- if ("outcome" %in% names(x)) x$outcome else x$study
    if (!all(c("ci", "two.tailed") %in% names(x)))
      stop("Recreate the batch object to retain its fitted confidence-interval settings.")
    two.tailed <- x$two.tailed
    if (effect == "smd") {
      d <- data.frame(label = as.character(row.labels), estimate = x$yi,
                      lower = x$lcl, upper = x$ucl, se = x$se)
    } else {
      d <- data.frame(label = as.character(row.labels), estimate = x$delta,
                      lower = x$lcl.delta, upper = x$ucl.delta, se = x$se.delta)
    }
  } else {
    stop("Supply output from hct.test, robust.smd, boot.hct.test, meta.robust.smd or multi.robust.smd.")
  }

  if (!nrow(d)) stop("There are no estimates to plot.")
  if (!is.null(labels)) {
    if (length(labels) != nrow(d) || anyNA(labels))
      stop("labels must contain one nonmissing label per plotted row.")
    d$label <- as.character(labels)
  }
  if (anyNA(d$label)) stop("Supply labels for rows with missing identifiers.")
  if (!all(vapply(d[c("estimate", "lower", "upper")], is.numeric, logical(1))))
    stop("Estimates and confidence limits must be numeric.")
  if (any(!is.finite(d$estimate))) stop("Point estimates must be finite.")
  known <- !is.na(d$lower) & !is.na(d$upper)
  if (any(d$lower[known] > d$upper[known])) stop("Lower limits exceed upper limits.")
  if (any(!known) && !inherits(x, "multi.robust.smd")) warning("Some confidence limits are missing; those rows show points only.")
  bands <- .plot_replayed_intervals(x, effect = effect, boot.ci = boot.ci)
  bands$ci <- if (inherits(x, "boot.hct.test")) boot.ci else
    if (inherits(x, "hct.test") || inherits(x, "robust.smd")) x$parms$ci else
      if (effect == "raw") "ct" else x$ci
  if (inherits(x, "multi.robust.smd")) {
    d$lower <- bands$lower.inner
    d$upper <- bands$upper.inner
    bands$ci <- if (effect == "raw") "ct" else x$ci
    two.tailed <- TRUE
  }
  d <- cbind(d, bands)
  d$two.tailed <- rep(two.tailed, length.out = nrow(d))
  reference <- 0
  if (inherits(x, "meta.robust.smd")) {
    settings <- attr(x, "pool.settings")
    models <- attr(x, "meta.models")
    reference <- NULL
    if (!is.null(settings) && isTRUE(settings$pool) && nrow(x) >= 2L) {
      model <- models[[effect]]
      estimates <- if (effect == "smd") x$yi else x$delta
      variances <- if (effect == "smd") x$vi else x$se.delta^2
      if (is.null(model) || !identical(as.numeric(model$yi), as.numeric(estimates)) ||
          !identical(as.numeric(model$vi), as.numeric(variances))) {
        # Subsetting a data frame retains attributes: refit the displayed subset.
        pooled <- .meta_smd_pool(x, settings$pool.method, settings$pool.test, settings$alpha)
        model <- pooled$models[[effect]]
      }
      reference <- as.numeric(model$b)
    }
  }
  attr(d, "reference") <- reference
  attr(d, "reference.series") <- if (length(reference)) 1L else integer(0)
  attr(d, "reference.colored") <- inherits(x, "meta.robust.smd")
  attr(d, "pooled") <- .pooled_hctplot(x, effect)
  attr(d, "effect") <- effect
  class(d) <- c("hctplot", "data.frame")
  d
}

# The shared S3 method. Row 1 appears at the top; equal-sized points are used.
# Infinite or out-of-range intervals have arrows at the plotting boundary.
#' @noRd
#' @export
plot.hctplot <- function(x, margins = NULL, reference = "auto", labels = NULL,
                         xlim = NULL, xlab = NULL, main = NULL,
                         legend.labels = NULL, legend.position = "topright",
                         col = c("#001F3F", "#B35400"), pch = 19, cex = 1,
                         ci.lwd = 1.5, cap = 0.08, label.padding = 1,
                         pretty = FALSE, pretty.height = 0.18, pretty.lwd = 3, ...) {
  series <- if (".series" %in% names(x)) x$.series else rep(1L, nrow(x))
  row <- if (".row" %in% names(x)) x$.row else seq_len(nrow(x))
  n.series <- max(series)
  if (!is.null(legend.labels)) {
    if (!is.character(legend.labels) ||
        length(legend.labels) != n.series ||
        anyNA(legend.labels)) {
      stop("Supply one nonmissing legend label per analysis object.")
    }
  }
  pooled <- attr(x, "pooled")
  n.pooled <- if (is.null(pooled)) 0L else nrow(pooled)
  pooled.labels <- if (n.pooled) if (n.series == 1L) "Pooled" else
    paste("Pooled", pooled$.series) else character(0)
  if (!length(col) || anyNA(col) || length(col) > 2L ||
      (n.series == 2L && length(col) != 2L))
    stop("Supply one color for a single object or two colors for two objects.")
  tryCatch(grDevices::col2rgb(col), error = function(e) stop("Invalid col value."))
  row.colors <- col[series]
  axis.labels <- x$label[match(seq_len(max(row)), row)]
  if (!is.null(margins) && (!is.numeric(margins) || any(!is.finite(margins))))
    stop("margins must be NULL or finite numeric x-axis positions.")
  automatic.reference <- identical(reference, "auto")
  if (automatic.reference) reference <- attr(x, "reference")
  if (!is.null(reference) && (!is.numeric(reference) || !length(reference) ||
      length(reference) > n.series || any(!is.finite(reference))))
    stop("reference must be auto, NULL, or one finite number per object.")
  if (is.null(xlab)) xlab <- if (attr(x, "effect") == "raw")
    "Mean Difference" else "Standardized Mean Difference"
  if (!is.logical(pretty) || length(pretty) != 1L || is.na(pretty))
    stop("pretty must be TRUE or FALSE.")
  if (!is.numeric(pretty.height) || length(pretty.height) != 1L ||
      !is.finite(pretty.height) || pretty.height <= 0 || pretty.height >= 1)
    stop("pretty.height must be between 0 and 1 (row-spacing units).")
  if (!is.numeric(pretty.lwd) || length(pretty.lwd) != 1L ||
      !is.finite(pretty.lwd) || pretty.lwd <= 0)
    stop("pretty.lwd must be positive.")
  # Both styles use the same outer bounds for identical axes.
  range.lower <- c(x$lower, x$lower.outer, pooled$lower, pooled$lower.outer)
  range.upper <- c(x$upper, x$upper.outer, pooled$upper, pooled$upper.outer)
  if (pretty) {
    if (any(!x$two.tailed))
      stop("Shaded ribbons require two-sided inference; use pretty = FALSE for one-sided fits.")
    if (any(!is.finite(as.matrix(x[c("lower.inner", "upper.inner", "lower.outer", "upper.outer")]))))
      stop("The fitted ribbon intervals are unavailable or infinite. Recreate older bootstrap objects, or use pretty = FALSE.")
    if (any(x$lower.outer > x$lower.inner | x$upper.outer < x$upper.inner |
            x$lower.inner > x$upper.inner)) stop("Fitted confidence intervals are not nested.")
    attr(x, "pretty.interval") <- unique(x$ci)
  }
  if (is.null(xlim)) {
    values <- c(x$estimate, pooled$estimate, range.lower, range.upper, margins, reference)
    values <- values[is.finite(values)]
    finite.ci <- is.finite(range.lower) & is.finite(range.upper)
    widths <- range.upper[finite.ci] - range.lower[finite.ci]
    # One full CI width of empty space on BOTH sides, then center at zero.
    # With several rows, use the widest finite CI for a shared plotting range.
    # Missing/infinite/zero-width CIs use the finite span as a fallback.
    padding <- if (length(widths)) max(widths) else 0
    if (padding <= 0) padding <- max(diff(range(values)), 0.1)
    half.range <- max(abs(values)) + padding
    xlim <- c(-half.range, half.range)
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || any(!is.finite(xlim)) ||
      xlim[1] >= xlim[2]) stop("xlim must contain two increasing finite numbers.")
  if (any(c(x$estimate, pooled$estimate) < xlim[1] | c(x$estimate, pooled$estimate) > xlim[2]))
    stop("xlim must include every point estimate.")

  if (!is.numeric(label.padding) || length(label.padding) != 1L ||
      !is.finite(label.padding) || label.padding < 0)
    stop("label.padding must be one finite nonnegative number (margin lines).")

  # Make room for labels, and restore the user's graphics settings afterwards.
  old <- graphics::par(c("mar", "las", "xpd"))
  on.exit(graphics::par(old), add = TRUE)
  has.labels <- any(nzchar(axis.labels))
  if (has.labels || n.pooled) {
    label.width <- max(graphics::strwidth(c(axis.labels, pooled.labels), units = "inches", cex = cex))
    left <- label.width / (graphics::par("csi") * graphics::par("mex")) + 1.25 + label.padding
  } else {
    left <- 0.75
  }
  pooled.lines <- if (n.pooled) 2 + (seq_len(n.pooled) - 1) * 1.15 else numeric(0)
  xlab.line <- if (n.pooled) max(pooled.lines) + 1.4 else 2.8
  graphics::par(mar = c(max(4, xlab.line + 1.2), left, 2.5, 0.75), las = 1, xpd = FALSE)
  y <- max(row) + 1L - row
  row.height <- if (pretty) pretty.height else 2 * cap
  half.offset <- 0
  if (n.series == 2L) {
    half.offset <- (row.height + 0.08) / 2
    y <- y + ifelse(series == 1L, half.offset, -half.offset)
  }
  y.padding <- max(0.5, half.offset + row.height / 2 + 0.1)
  graphics::plot(NA_real_, NA_real_, xlim = xlim, ylim = c(1 - y.padding, max(row) + y.padding),
                 xaxs = "i", yaxs = "i", yaxt = "n", ylab = "", bty = "n",
                 xlab = "", main = main, ...)
  graphics::mtext(xlab, side = 1, line = xlab.line)
  if (has.labels) {
    separate.labels <- n.series == 2L &&
      any(x$label[series == 1L] != x$label[series == 2L])
    
    if (separate.labels) {
      graphics::axis(
        2, at = y, labels = x$label,
        tick = FALSE, cex.axis = cex, gap.axis = -1
      )
    } else {
      graphics::axis(
        2, at = rev(seq_len(max(row))), labels = axis.labels,
        tick = FALSE, cex.axis = cex
      )
    }
  }
  if (pretty) {
    for (i in seq_len(nrow(x)))
      .draw_hct_ribbon(x$estimate[i], x$lower.inner[i], x$upper.inner[i],
                       x$lower.outer[i], x$upper.outer[i], y[i], pretty.height,
                       row.colors[i], pretty.lwd)
  } else {
  for (i in seq_len(nrow(x))) {
    if (!is.na(x$lower[i]) && !is.na(x$upper[i])) {
      lower <- max(xlim[1], x$lower[i])
      upper <- min(xlim[2], x$upper[i])
      if (lower <= upper) {
        graphics::segments(lower, y[i], upper, y[i], col = row.colors[i], lwd = ci.lwd)
        if (is.finite(x$lower[i]) && x$lower[i] >= xlim[1] && x$lower[i] <= xlim[2])
          graphics::segments(lower, y[i] - cap, lower, y[i] + cap, col = row.colors[i], lwd = ci.lwd)
        if (is.finite(x$upper[i]) && x$upper[i] >= xlim[1] && x$upper[i] <= xlim[2])
          graphics::segments(upper, y[i] - cap, upper, y[i] + cap, col = row.colors[i], lwd = ci.lwd)
        if (x$lower[i] < xlim[1])
          graphics::arrows(xlim[1] + 0.04 * diff(xlim), y[i], xlim[1], y[i],
                           length = 0.07, col = row.colors[i], lwd = ci.lwd)
        if (x$upper[i] > xlim[2])
          graphics::arrows(xlim[2] - 0.04 * diff(xlim), y[i], xlim[2], y[i],
                           length = 0.07, col = row.colors[i], lwd = ci.lwd)
      }
    }
  }
  graphics::points(x$estimate, y, pch = pch, cex = cex, col = row.colors)
  }
  # Reference/margin lines remain visible above both styles.
  if (!is.null(reference)) {
    reference.colors <- if (automatic.reference && isTRUE(attr(x, "reference.colored")))
      col[attr(x, "reference.series")] else rep("grey65", length(reference))
    for (j in seq_along(reference))
      graphics::abline(v = reference[j], col = reference.colors[j], lty = 3, lwd = 1.5)
  }
  if (length(margins)) graphics::abline(v = margins, col = "red", lty = 2, lwd = 1.5)
  if (n.pooled) attr(x, "pooled") <- .draw_hct_pooled(pooled, pooled.labels, pooled.lines,
      xlim, col, pretty, pretty.height, pretty.lwd, ci.lwd, pch, cex)
  if (!is.null(legend.labels)) {
    graphics::legend(
      x = legend.position,
      legend = legend.labels,
      col = col[seq_len(n.series)],
      pch = rep_len(pch, n.series),
      bty = "n",
      cex = cex
    )
  }
  x$.y <- y
  attr(x, "reference.drawn") <- reference
  attr(x, "xlim") <- xlim
  invisible(x)
}

# Convenient entry point for every existing output type.
.hctplot <- function(
    x, y = NULL, effect = NULL,
    boot.ci = c("percentile", "bca"),
    labels = NULL, ...
) {
  first <- .as_hctplot(
    x, effect = effect, boot.ci = boot.ci
  )
  
  if (is.null(y)) {
    if (is.list(labels))
      stop("For one object, labels must be a character vector.")
    
    if (!is.null(labels)) {
      if (length(labels) != nrow(first) || anyNA(labels))
        stop("Supply one nonmissing label per row.")
      
      first$label <- as.character(labels)
    }
    
    first$.series <- 1L
    first$.row <- seq_len(nrow(first))
    
    return(plot(first, ...))
  }
  
  second <- .as_hctplot(
    y, effect = attr(first, "effect"), boot.ci = boot.ci
  )
  
  if (nrow(first) != nrow(second))
    stop("The two objects must have the same number of rows.")
  
  if (is.list(labels)) {
    if (length(labels) != 2L ||
        any(lengths(labels) != nrow(first)) ||
        anyNA(labels[[1]]) ||
        anyNA(labels[[2]])) {
      stop("A labels list must contain one complete label vector per object.")
    }
    
    # Explicit label lists pair rows by position.
    first$label <- as.character(labels[[1]])
    second$label <- as.character(labels[[2]])
    
  } else {
    kind <- function(z) {
      if (is.data.frame(z) && "outcome" %in% names(z)) {
        "outcome"
      } else if (is.data.frame(z) && "study" %in% names(z)) {
        "study"
      } else {
        "single"
      }
    }
    
    if (kind(x) != kind(y))
      stop("Match two study tables or two outcome tables, or supply labels as a list.")
    
    keys1 <- first$label
    keys2 <- second$label
    
    if (anyDuplicated(keys1) || anyDuplicated(keys2))
      stop("Study/outcome identifiers must be unique.")
    
    if (!setequal(keys1, keys2))
      stop("The two objects must have matching study/outcome names.")
    
    second <- second[match(keys1, keys2), , drop = FALSE]
    
    if (!is.null(labels)) {
      if (length(labels) != nrow(first) || anyNA(labels))
        stop("Supply one nonmissing label per matched row.")
      
      first$label <- second$label <- as.character(labels)
    }
  }
  
  first$.series <- 1L
  second$.series <- 2L
  first$.row <- second$.row <- seq_len(nrow(first))
  
  combined <- rbind(
    as.data.frame(first),
    as.data.frame(second)
  )
  
  pooled.first <- attr(first, "pooled")
  pooled.second <- attr(second, "pooled")
  
  if (!is.null(pooled.second))
    pooled.second$.series <- 2L
  
  attr(combined, "pooled") <-
    if (is.null(pooled.first)) {
      pooled.second
    } else if (is.null(pooled.second)) {
      pooled.first
    } else {
      rbind(pooled.first, pooled.second)
    }
  
  attr(combined, "reference") <- c(
    attr(first, "reference"),
    attr(second, "reference")
  )
  
  attr(combined, "reference.series") <- c(
    rep(1L, length(attr(first, "reference"))),
    rep(2L, length(attr(second, "reference")))
  )
  
  attr(combined, "reference.colored") <-
    isTRUE(attr(first, "reference.colored")) &&
    isTRUE(attr(second, "reference.colored"))
  
  attr(combined, "effect") <- attr(first, "effect")
  class(combined) <- c("hctplot", "data.frame")
  
  # Forward pretty, margins, legend settings, colors, etc.
  plot(combined, ...)
}

#' Plot mean differences and confidence intervals
#'
#' @description Forest-style plots for individual estimates, studies and multiple outcomes.
#' @param y Optional second analysis object to overlay. Rows are matched by unique study/outcome names or by custom labels supplied as a list.
#' @param x An object returned by hct.test, robust.smd, boot.hct.test, meta.robust.smd or multi.robust.smd.
#' @param effect NULL selects raw differences for hct.test and standardized differences for other outputs. Use "raw" or "smd" for bootstrap and batch outputs. For an hct.test object, first call robust.smd to plot an SMD.
#' @param boot.ci Stored bootstrap interval used by the ordinary style: "percentile" or "bca".
#' @param margins NULL or finite numeric x-axis positions for red dashed vertical lines. Values are in the plotted effect units.
#' @param reference Default auto uses the pooled effect for meta-analysis objects, with one line in each object color. Other objects use zero. NULL omits reference lines; numeric values override the defaults.
#' @param labels Optional display labels in first-object order. For two objects, a list of two unique label vectors can explicitly define row matching. A single common label vector renames rows after their original identifiers have been matched. Defaults to study/outcome names, or blank labels for single estimates.
#' @param xlim Two increasing finite x-axis limits. NULL gives a common zero-centered range for both styles. Explicit limits must contain all point estimates; confidence intervals can be clipped.
#' @param xlab Bottom axis title. Defaults to "Mean Difference" or "Standardized Mean Difference".
#' @param main Optional plot title.
#' @param col Colors for estimates and intervals in both styles. Defaults to dark navy and dark orange. With two objects, supply two colors, one per object. Ribbons use lighter transparent versions of these colors.
#' @param pch Point symbol for the ordinary style.
#' @param cex Point size and row-label text scaling.
#' @param ci.lwd Confidence-interval line width for the ordinary style.
#' @param cap Half-height of ordinary interval endpoint ticks, in row-spacing units.
#' @param label.padding Extra outer space before row labels, in margin lines; default one.
#' @param pretty Logical; FALSE draws intervals and points. TRUE draws shaded ribbons. Multi-outcome plots use equivalence intervals (90 percent at alpha = 0.05) and fade to 95 percent; other plots use 95/99 percent fitted intervals.
#' @param pretty.height Ribbon height in row-spacing units; default 0.18.
#' @param pretty.lwd Estimate-line width in the shaded style. The inner interval boundary lines have one-third this width.
#' @param ... Additional arguments passed to the initial base graphics plot.
#'
#' @details The first row is shown at the top, with signed estimates and equal-sized points. Bootstrap points use observed estimates, rather than bootstrap means. Meta-analysis reference lines use pooled estimates from the selected metafor model.
#'
#' Multi-outcome SMD plots inherit the specified ci method for equivalence intervals; raw-difference plots retain central-t inference. Both styles show 100*(1-2*alpha) percent intervals (90 percent by default), with shaded fading to 100*(1-alpha) percent (95 percent by default). Other ordinary plots use their stored intervals. Missing limits show a point only; infinite or clipped ordinary intervals have boundary arrows.
#'
#' For other objects, pretty = TRUE uses both levels from the fitted interval method: z uses normal quantiles; ct uses central-t quantiles and the fitted df; nct reruns the original analysis function using retained inputs and arguments, changing only alpha to 0.01 for 99 percent limits or 0.05 for 95 percent limits. SMD plots call robust.smd with the retained fitted hct.test object and original settings; batch plots rerun the original meta.robust.smd or multi.robust.smd call. Bootstrap plots rerun boot.hct.test with the original data and settings, using its percentile or BCa limits. Unavailable BCa bands, infinite limits, disabled CIs and one-sided inference cannot be displayed as shaded ribbons; ordinary plots remain available.
#'
#' The ribbon uses a lighter version of its object color inside the inner interval, with thin vertical lines at its limits. It fades to the outer interval, allowing asymmetric intervals. The estimate line stays within the ribbon height. No explanatory annotation is added.
#'
#' Both styles share a zero-centered default x-axis range based on finite stored and outer ribbon limits, estimates, margins and reference. Padding on either side is at least the widest finite stored or outer ribbon interval width. Missing or infinite widths use the finite span as a fallback. xlim overrides this rule. Computing the common range also replays the analyses for ordinary plots. The original interval is reused when its alpha already matches the requested level.
#'
#' Two objects are matched by unique identifiers and reordered to the first object. Unequal counts, duplicate identifiers or unmatched names raise errors. A labels list can supply unique matching identifiers. The first object is placed slightly above each row label and the second slightly below, with a small gap based on the ribbon or endpoint-tick height. Axes and effect scales are shared; each object retains its own fitted CI method. No inference about the difference between objects is calculated.
#'
#' Meta-analysis pooled summaries appear below the tick labels and above the axis title, using the fitted metafor model and its CI method. Two meta-analysis objects have separate vertically spaced summaries in matching colors. Displayed study subsets are repooled. Multi-outcome plots show only individual outcomes and use zero as the default reference. Both styles use the same axis limits.
#'
#' Replay retains evaluated data and arguments, rather than looking up variable names in the user workspace. The original starting RNG state is reused during replay, and the user RNG state is restored afterwards. Sequential bootstrap replay is reproducible; parallel reproducibility depends on the selected backend. Plotting bootstrap objects can take as long as the original bootstrap. Recreate objects from older versions to retain these inputs.
#'
#' Plots use base graphics and restore the changed graphics parameters afterwards. Long labels or many rows may require a larger graphics device.
#'
#' @return For meta-analysis plots, the pooled attribute contains displayed pooled summaries, CI levels, object-series identifiers and margin positions. Invisibly, a data frame of the plotted estimates, stored limits and SEs, with .series, .row and .y columns identifying objects, matched rows and drawn positions, and the x-axis limits in its xlim attribute. The lower.inner, upper.inner, lower.outer and upper.outer columns and inner.level/outer.level describe the displayed bands. Meta-analysis objects retain lower95/upper95/lower99/upper99 aliases.
#'
#' @examples fit <- hct.test(m1 = 0.4, m2 = 0, sd1 = 1, sd2 = 1.5,
#'                 n1 = 60, n2 = 60, method = "hc2", ci = "z", verbose = FALSE)
#' plot(fit, margins = c(-0.5, 0.5))
#' effect <- robust.smd(fit, ci = "z", verbose = FALSE)
#' plot(effect, pretty = TRUE, margins = c(-0.2, 0.2))
#' 
#' studies <- data.frame(study = c("Study A", "Study B"),
#'                       m1 = c(0.4, -0.2), m2 = c(0, 0),
#'                       sd1 = c(1, 1.2), sd2 = c(1.5, 1),
#'                       n1 = c(60, 50), n2 = c(60, 55))
#' tab <- meta.robust.smd(studies, study = "study", ci = "z")
#' ordinary <- plot(tab, margins = c(-0.2, 0.2))
#' shaded <- plot(tab, pretty = TRUE, margins = c(-0.2, 0.2))
#' # Compare the identical default x-axis limits numerically:
#' rbind(ordinary = attr(ordinary, "xlim"), shaded = attr(shaded, "xlim"))
#' plot(tab, pretty = TRUE, xlim = c(-1.5, 1.5), label.padding = 2)
#'
#' studies2 <- studies
#' studies2$m1 <- studies2$m1 + 0.2
#' tab2 <- meta.robust.smd(studies2, study = "study", ci = "z")
#' plot(tab, tab2, pretty = TRUE, col = c("navy", "darkorange3"))
#'
#' @rdname hcttest-plot
#' @method plot hct.test
#' @export
plot.hct.test <- function(x, y = NULL, effect = NULL, boot.ci = c("percentile", "bca"),
                          margins = NULL, reference = "auto", labels = NULL,
                          xlim = NULL, xlab = NULL, main = NULL,
                          legend.labels = NULL, legend.position = "topright",
                          col = c("#001F3F", "#B35400"), pch = 19, cex = 1,
                          ci.lwd = 1.5, cap = 0.08, label.padding = 1,
                          pretty = FALSE, pretty.height = 0.18, pretty.lwd = 3, ...) {
  .hctplot(x, y = y, effect = effect, boot.ci = boot.ci, labels = labels,
           legend.labels = legend.labels, legend.position = legend.position,
            margins = margins, reference = reference, xlim = xlim,
            xlab = xlab, main = main, col = col, pch = pch, cex = cex,
            ci.lwd = ci.lwd, cap = cap, label.padding = label.padding,
            pretty = pretty, pretty.height = pretty.height,
            pretty.lwd = pretty.lwd, ...)
}

#' @rdname hcttest-plot
#' @method plot robust.smd
#' @export
plot.robust.smd <- plot.hct.test

#' @rdname hcttest-plot
#' @method plot boot.hct.test
#' @export
plot.boot.hct.test <- plot.hct.test

#' @rdname hcttest-plot
#' @method plot meta.robust.smd
#' @export
plot.meta.robust.smd <- plot.hct.test

#' @rdname hcttest-plot
#' @method plot multi.robust.smd
#' @export
plot.multi.robust.smd <- plot.hct.test
