#' Overlay Funnel Plots for Study-Level Effect Estimates
#' 
#' @description Draws a funnel plot for one or two meta.robust.smd objects.
#' Each object has its own pooled reference estimate and colored contours.
#' 
#' @param x An object returned by \code{meta.robust.smd()}.
#' @param y An optional second object. Both objects must contain the same unique, nonmissing study names. The second object is aligned to the first by name.
#' @param effect \code{"smd"} plots \code{yi} against \code{sqrt(vi)}; \code{"raw"} plots \code{delta} against \code{se.delta}.
#' @param col One color, or one color per object, used for points, reference lines, and contour shading. Defaults to dark navy and dark orange.
#' @param pch One plotting symbol, or one per object.
#' @param cex A positive point-size multiplier.
#' @param level Pseudo confidence levels as probabilities strictly between zero and one. Levels are sorted and duplicates removed.
#' @param shade Logical; draw translucent, separately shaded contour bands.
#' @param connect Logical; connect corresponding points when two objects are supplied.
#' @param contour.legend Logical; show the contour legend when \code{legend = TRUE} and shading or contour boundaries are displayed.
#' @param jitter Logical; apply small, deterministic horizontal display offsets. Defaults to \code{FALSE}.
#' @param jitter.width Maximum absolute horizontal displacement as a fraction of the x-axis range, from zero to 0.05. Default is 0.008.
#' @param refline One finite reference value shared by both objects, or one per object. \code{NULL} extracts the selected pooled estimate from each object's \code{pooled} attribute. Use \code{refline = 0} for zero-centered contours.
#' @param xlim Two increasing finite x-axis limits. The default includes estimates and the widest contours, with padding.
#' @param ylim Standard-error limits as \code{c(0, maximum)}. The axis is reversed so zero appears at the top. Defaults to 1.08 times the largest standard error.
#' @param legend Logical master switch. \code{FALSE} suppresses both legends.
#' @param legend.labels One model label per object. Defaults to model numbers and, when available, estimator names.
#' @param xlab Optional x-axis label.
#' @param main Optional plot title.
#' @param model.legend Logical; show model labels when \code{legend = TRUE}.
#' @param legend.position Position keyword passed to \code{graphics::legend()} for the model legend.
#' @param contour.position Position keyword for the contour legend.
#' @param legend.cex Positive size multiplier for both legends.
#' @param contours Logical; draw contour boundary lines.
#' @param show.refline Logical; draw the vertical reference lines. Turning these lines off does not change the contour centers.
#' @param grid Logical; draw horizontal grid lines.
#' @param shade.alpha Opacity in [0, 1], either one value or one per sorted level. \code{NULL} uses decreasing opacity from the inner to outer bands.
#' @param point.alpha Point opacity in [0, 1].
#' @param contour.lty Line types recycled across sorted contour levels.
#' @param contour.lwd Positive contour line width.
#' @param label \code{FALSE} suppresses study labels; \code{TRUE} labels every point; \code{"out"} labels points outside the largest pseudo confidence region, using original estimates before jitter.
#' @param label.cex Positive study-label size multiplier.
#' @param ... Additional graphical arguments passed to the initial base plot, such as \code{cex.axis}, \code{cex.lab}, \code{axes}, and \code{font.lab}.
#' 
#' @details At standard error \eqn{s}, contours for level \eqn{C} are
#' \deqn{r \pm \Phi^{-1}((1+C)/2)s,}
#' where \eqn{r} is the selected reference estimate. These are normal-based pseudo
#' confidence regions, not confidence intervals for the pooled estimate. They
#' ignore uncertainty in the reference estimate and do not include between-study
#' heterogeneity in their widths. They do not inherit the study-level confidence
#' interval distribution or alpha setting.
#' 
#' For levels 0.90, 0.95, and 0.99, critical values are approximately 1.645, 1.960,
#' and 2.576. With \code{refline = 0}, these contours correspond to two-sided
#' normal-test p-value thresholds of 0.10, 0.05, and 0.01. With pooled centers,
#' they indicate deviations from the pooled estimates instead.
#' 
#' Shaded bands within each funnel do not overlap. Shading from two funnels
#' may overlap and blend colors. All study estimates must be finite and sampling
#' standard errors positive and finite.
#' 
#' Jitter changes only plotted x-coordinates. It leaves original estimates,
#' standard errors, contour geometry, and pooled centers unchanged. Corresponding
#' studies receive the same horizontal offset in both models. Offsets are stable
#' under row reordering and do not change the random-number generator state.
#' 
#' Funnel asymmetry is not sufficient evidence of publication bias. Dependence
#' among multiple effect estimates from the same study should be considered when
#' interpreting the plot or conducting associated statistical tests.
#' 
#' @return Invisibly returns a list containing \code{data} (one data frame per object,
#' with \code{study}, \code{estimate}, \code{se}, and \code{estimate.plot}),
#' \code{refline}, \code{level}, \code{xlim}, and \code{ylim}.
#' 
#' @examples studies <- readRDS(system.file(
#'   "extdata", "table3-studies.rds", package = "hcttest"
#' ))
#' a <- meta.robust.smd(studies, study = "study", method = "constant",
#'   ss.adjust = TRUE, ci = "ct", pool = TRUE, pool.method = "REML")
#' b <- meta.robust.smd(studies, study = "study", method = "hc2",
#'   ss.adjust = TRUE, ci = "ct", pool = TRUE, pool.method = "REML")
#' funnel(a, b, legend.labels = c("Constant", "HC2"))
#' funnel(a, b, refline = 0, jitter = TRUE, legend = FALSE)
#' funnel(b, level = 0.95, shade = FALSE, legend = FALSE)
#' 
#' @references Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., and Rushton, L. (2008).
#' Contour-enhanced meta-analysis funnel plots help distinguish publication bias
#' from other causes of asymmetry. \emph{Journal of Clinical Epidemiology},
#' \bold{61}(10), 991--996. doi:10.1016/j.jclinepi.2007.11.010.
#' 
#' Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package.
#' \emph{Journal of Statistical Software}, \bold{36}(3), 1--48.
#' doi:10.18637/jss.v036.i03.
#' 
#' @seealso \code{\link{meta.robust.smd}}, \code{metafor::funnel()}
#' @export
funnel <- function(x, y = NULL, effect = c("smd", "raw"),
                   col = c("#001F3F", "#B35400"),
                   pch = c(16, 17), cex = 1.1,
                   level = c(0.90, 0.95, 0.99), shade = TRUE, connect = FALSE,
                   contour.legend = TRUE, jitter = FALSE, jitter.width = 0.008,
                   refline = NULL, xlim = NULL, ylim = NULL,
                   legend = TRUE, legend.labels = NULL,
                   xlab = NULL, main = NULL,
                   model.legend = TRUE, legend.position = "topright",
                   contour.position = "topleft", legend.cex = 0.8,
                   contours = TRUE, show.refline = TRUE, grid = TRUE,
                   shade.alpha = NULL, point.alpha = 0.8,
                   contour.lty = c(1, 2, 3), contour.lwd = 0.9,
                   label = FALSE, label.cex = 0.65, ...) {
  effect <- match.arg(effect)
  objects <- if (is.null(y)) list(x) else list(x, y)
  k <- length(objects)
  fail <- function(message) stop(message, call. = FALSE)
  flags <- list(model.legend = model.legend, contours = contours,
                show.refline = show.refline, grid = grid)
  for (nm in names(flags))
    if (!is.logical(flags[[nm]]) || length(flags[[nm]]) != 1L || is.na(flags[[nm]]))
      fail(paste0("'", nm, "' must be TRUE or FALSE."))
  if (!(identical(label, FALSE) || identical(label, TRUE) || identical(label, "out")))
    fail("'label' must be FALSE, TRUE, or 'out'.")
  for (v in list(legend.cex, label.cex, contour.lwd))
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0)
      fail("Legend size, label size, and contour width must be positive.")
  if (!is.numeric(point.alpha) || length(point.alpha) != 1L ||
      !is.finite(point.alpha) || point.alpha < 0 || point.alpha > 1)
    fail("'point.alpha' must be between zero and one.")
  if (!length(contour.lty) || anyNA(contour.lty)) fail("Supply contour line types.")
  if (!is.numeric(level) || !length(level) ||
      any(!is.finite(level) | level <= 0 | level >= 1))
    fail("'level' must contain probabilities between zero and one.")
  level <- sort(unique(level))
  if (!is.logical(contour.legend) || length(contour.legend) != 1L || is.na(contour.legend))
    fail("'contour.legend' must be TRUE or FALSE.")
  if (!is.logical(shade) || length(shade) != 1L || is.na(shade) ||
      !is.logical(connect) || length(connect) != 1L || is.na(connect) ||
      !is.logical(legend) || length(legend) != 1L || is.na(legend))
    fail("'shade', 'connect', and 'legend' must be TRUE or FALSE.")
  if (!is.logical(jitter) || length(jitter) != 1L || is.na(jitter))
    fail("'jitter' must be TRUE or FALSE.")
  if (!is.numeric(jitter.width) || length(jitter.width) != 1L ||
      !is.finite(jitter.width) || jitter.width < 0 || jitter.width > 0.05)
    fail("'jitter.width' must be between zero and 0.05 (fraction of x-axis range).")
  if (!is.numeric(cex) || length(cex) != 1L || !is.finite(cex) || cex <= 0)
    fail("'cex' must be a positive finite number.")
  if (!length(col) || !length(pch)) fail("Supply colors and plotting symbols.")
  col <- rep_len(col, k)
  pch <- rep_len(pch, k)
  grDevices::col2rgb(col)
  extract <- function(object) {
    if (!inherits(object, "meta.robust.smd") || !is.data.frame(object))
      fail("Supply objects returned by meta.robust.smd().")
    fields <- if (effect == "smd") c("study", "yi", "vi") else
      c("study", "delta", "se.delta")
    if (!all(fields %in% names(object))) fail("Required study fields are missing.")
    estimate <- object[[fields[2L]]]
    if (effect == "smd" && (!is.numeric(object$vi) ||
        any(!is.finite(object$vi) | object$vi <= 0)))
      fail("Sampling variances must be positive and finite.")
    se <- if (effect == "smd") sqrt(object$vi) else object$se.delta
    labels <- as.character(object$study)
    if (!length(labels) || anyNA(labels) || any(!nzchar(labels)) ||
        anyDuplicated(labels))
      fail("Each object must have unique, nonmissing study names.")
    if (!is.numeric(estimate) || !is.numeric(se) ||
        any(!is.finite(estimate)) || any(!is.finite(se) | se <= 0))
      fail("Study estimates must be finite and standard errors positive and finite.")
    data.frame(study = labels, estimate = estimate, se = se)
  }
  dat <- lapply(objects, extract)
  if (k == 2L) {
    if (!setequal(dat[[1]]$study, dat[[2]]$study))
      fail("The two objects must contain the same study names.")
    dat[[2]] <- dat[[2]][match(dat[[1]]$study, dat[[2]]$study), , drop = FALSE]
    rownames(dat[[2]]) <- NULL
  }
  if (is.null(refline)) {
    refline <- vapply(objects, function(object) {
      pooled <- attr(object, "pooled")
      if (!is.data.frame(pooled) ||
          !all(c("effect", "estimate") %in% names(pooled)))
        fail("No pooled estimate found; use pool = TRUE or supply 'refline'.")
      value <- pooled$estimate[pooled$effect == effect]
      if (length(value) != 1L || !is.finite(value))
        fail("A finite pooled estimate is required for the selected effect.")
      value
    }, numeric(1))
  }
  if (!is.numeric(refline) || !length(refline) ||
      !length(refline) %in% c(1L, k) || any(!is.finite(refline)))
    fail("'refline' must contain one finite value or one per object.")
  refline <- rep_len(refline, k)
  if (is.null(legend.labels)) {
    legend.labels <- vapply(seq_len(k), function(i) {
      methods <- unique(as.character(objects[[i]]$method))
      paste0("Model ", i, if (length(methods) == 1L) paste0(" (", toupper(methods), ")") else "")
    }, character(1))
  }
  if (length(legend.labels) != k || anyNA(legend.labels))
    fail("Supply one legend label per object.")
  all.se <- unlist(lapply(dat, function(z) z$se))
  if (is.null(ylim)) ylim <- c(0, 1.08 * max(all.se))
  if (!is.numeric(ylim) || length(ylim) != 2L || any(!is.finite(ylim)) ||
      ylim[1L] != 0 || ylim[2L] <= 0)
    fail("'ylim' must be c(0, positive_maximum_standard_error).")
  critical <- stats::qnorm((1 + level) / 2)
  if (is.null(xlim)) {
    limits <- range(c(unlist(lapply(dat, function(z) z$estimate)),
                      refline - max(critical) * ylim[2L],
                      refline + max(critical) * ylim[2L]))
    xlim <- limits + c(-1, 1) * 0.06 * diff(limits)
  }
  if (!is.numeric(xlim) || length(xlim) != 2L || any(!is.finite(xlim)) ||
      xlim[1L] >= xlim[2L]) fail("'xlim' must contain two increasing finite values.")
  offsets <- rep(0, nrow(dat[[1]]))
  if (jitter && length(offsets) > 1L) {
    # Reproducible offsets; shared by models, independent of RNG state.
    order.id <- match(dat[[1]]$study, sort(dat[[1]]$study))
    offsets <- (order.id * ((sqrt(5) - 1) / 2)) %% 1
    offsets <- offsets - mean(offsets)
    offsets <- offsets / max(abs(offsets)) * jitter.width * diff(xlim)
  }
  for (i in seq_len(k)) dat[[i]]$estimate.plot <- dat[[i]]$estimate + offsets
  if (is.null(xlab)) xlab <- if (effect == "smd")
    "Standardized Mean Difference" else "Mean Difference"
  graphics::plot(NA_real_, NA_real_, type = "n", xlim = xlim,
                 ylim = rev(ylim), xaxs = "i", yaxs = "i",
                 xlab = xlab, ylab = "Standard Error", main = main,
                 bty = "l", las = 1, ...)
  ticks <- pretty(ylim)
  ticks <- ticks[ticks >= 0 & ticks <= ylim[2L]]
  if (grid) graphics::abline(h = ticks, col = "grey92", lwd = 0.7)
  alpha <- if (length(level) == 1L) 0.14 else
    seq(0.24, 0.06, length.out = length(level))
  if (!is.null(shade.alpha)) {
    if (!is.numeric(shade.alpha) || !length(shade.alpha) ||
        !length(shade.alpha) %in% c(1L, length(level)) ||
        any(!is.finite(shade.alpha) | shade.alpha < 0 | shade.alpha > 1))
      fail("'shade.alpha' must be in [0, 1], with one value or one per sorted level.")
    alpha <- rep_len(shade.alpha, length(level))
  }
  line.types <- rep_len(contour.lty, length(level))
  if (shade) for (i in seq_len(k)) {
    for (j in seq_along(level)) {
      if (j == 1L) {
        graphics::polygon(
          c(refline[i], refline[i] - critical[j] * ylim[2L],
            refline[i] + critical[j] * ylim[2L]),
          c(0, ylim[2L], ylim[2L]),
          col = grDevices::adjustcolor(col[i], alpha.f = alpha[j]), border = NA)
      } else {
        for (side in c(-1, 1)) graphics::polygon(
          c(refline[i], refline[i] + side * critical[j - 1L] * ylim[2L],
            refline[i] + side * critical[j] * ylim[2L]),
          c(0, ylim[2L], ylim[2L]),
          col = grDevices::adjustcolor(col[i], alpha.f = alpha[j]), border = NA)
      }
    }
  }
  if (connect && k == 2L) graphics::segments(
    dat[[1]]$estimate.plot, dat[[1]]$se, dat[[2]]$estimate.plot, dat[[2]]$se,
    col = "grey75", lwd = 0.7)
  for (i in seq_len(k)) {
    if (contours) for (j in seq_along(level)) graphics::segments(
      rep(refline[i], 2), rep(0, 2),
      refline[i] + c(-1, 1) * critical[j] * ylim[2L],
      rep(ylim[2L], 2),
      col = grDevices::adjustcolor(col[i], alpha.f = 0.65),
      lty = line.types[j], lwd = contour.lwd)
    if (show.refline) graphics::abline(v = refline[i], col = col[i], lty = 3, lwd = 1.3)
    graphics::points(dat[[i]]$estimate.plot, dat[[i]]$se,
                     col = grDevices::adjustcolor(col[i], alpha.f = point.alpha),
                     pch = pch[i], bg = grDevices::adjustcolor(col[i], alpha.f = point.alpha),
                     cex = cex)
    if (!identical(label, FALSE)) {
      keep <- if (identical(label, "out"))
        abs(dat[[i]]$estimate - refline[i]) > max(critical) * dat[[i]]$se else
        rep(TRUE, nrow(dat[[i]]))
      graphics::text(dat[[i]]$estimate.plot[keep], dat[[i]]$se[keep],
                     labels = dat[[i]]$study[keep], pos = 4,
                     cex = label.cex, col = col[i])
    }
  }
  if (legend && model.legend) graphics::legend(legend.position, legend = legend.labels,
                               col = col, pt.bg = col, pch = pch, bty = "n", cex = legend.cex)
  if (legend && contour.legend && (contours || shade)) graphics::legend(contour.position,
    legend = paste0(format(100 * level, trim = TRUE), "%"),
    title = "Pseudo CI levels", col = "grey45", lty = line.types,
    lwd = contour.lwd, bty = "n", cex = legend.cex)
  invisible(list(data = dat, refline = refline, level = level,
                 xlim = xlim, ylim = ylim))
}
