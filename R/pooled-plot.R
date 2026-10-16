# Internal preparation and rendering of batch pooled summaries.
.pooled_hctplot <- function(x, effect) {
  if (!inherits(x, "meta.robust.smd")) return(NULL)
  settings <- attr(x, "pool.settings")
  if (is.null(settings)) return(NULL)
  alpha <- settings$alpha
  outer.alpha <- min(alpha, 0.01)
  if (inherits(x, "meta.robust.smd") && isTRUE(settings$pool) && nrow(x) >= 2L) {
    inner <- .meta_smd_pool(x, settings$pool.method, settings$pool.test, alpha)$summary
    outer <- .meta_smd_pool(x, settings$pool.method, settings$pool.test, outer.alpha)$summary
    inner <- inner[inner$effect == effect, , drop = FALSE]
    outer <- outer[outer$effect == effect, , drop = FALSE]
  } else return(NULL)
  data.frame(estimate = inner$estimate, se = inner$se,
             lower = inner$lcl, upper = inner$ucl,
             lower.inner = inner$lcl, upper.inner = inner$ucl,
             lower.outer = outer$lcl, upper.outer = outer$ucl,
             inner.level = 1 - alpha, outer.level = 1 - outer.alpha, .series = 1L)
}

.draw_hct_pooled <- function(pooled, labels, lines, xlim, col, pretty,
                              pretty.height, pretty.lwd, ci.lwd, pch, cex) {
  line.height <- graphics::par("csi") * graphics::par("mex")
  bottom <- graphics::grconvertY(graphics::par("usr")[3], "user", "inches")
  y <- graphics::grconvertY(bottom - lines * line.height, "inches", "user")
  half <- abs(graphics::grconvertY(bottom + pretty.height * line.height, "inches", "user") -
                graphics::grconvertY(bottom, "inches", "user"))
  old.xpd <- graphics::par("xpd")
  old.usr <- graphics::par("usr")
  on.exit({
    graphics::par(xpd = old.xpd)
    graphics::clip(old.usr[1], old.usr[2], old.usr[3], old.usr[4])
  }, add = TRUE)
  graphics::par(xpd = NA)
  graphics::text(xlim[1], y, labels = labels, pos = 2, cex = cex, col = col[pooled$.series])
  for (i in seq_len(nrow(pooled))) {
    color <- col[pooled$.series[i]]
    graphics::clip(xlim[1], xlim[2], y[i] - half * 1.5, y[i] + half * 1.5)
    if (pretty) {
      .draw_hct_ribbon(pooled$estimate[i], pooled$lower.inner[i], pooled$upper.inner[i],
                       pooled$lower.outer[i], pooled$upper.outer[i], y[i], 2 * half,
                       color, pretty.lwd)
    } else {
      lower <- max(xlim[1], pooled$lower[i]); upper <- min(xlim[2], pooled$upper[i])
      graphics::segments(lower, y[i], upper, y[i], col = color, lwd = ci.lwd)
      if (pooled$lower[i] >= xlim[1])
        graphics::segments(lower, y[i] - half, lower, y[i] + half, col = color, lwd = ci.lwd)
      if (pooled$upper[i] <= xlim[2])
        graphics::segments(upper, y[i] - half, upper, y[i] + half, col = color, lwd = ci.lwd)
      if (pooled$lower[i] < xlim[1])
        graphics::arrows(xlim[1] + 0.04 * diff(xlim), y[i], xlim[1], y[i],
                         length = 0.07, col = color, lwd = ci.lwd)
      if (pooled$upper[i] > xlim[2])
        graphics::arrows(xlim[2] - 0.04 * diff(xlim), y[i], xlim[2], y[i],
                         length = 0.07, col = color, lwd = ci.lwd)
      graphics::points(pooled$estimate[i], y[i], pch = pch, cex = cex, col = color)
    }
  }
  pooled$.y <- y
  pooled
}
