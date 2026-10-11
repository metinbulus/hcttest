library(hcttest)

# Multi-outcome calculations and both plot styles omit pooled summaries.
dat <- subset(iris, Species != "virginica")
outcomes <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
a <- multi.robust.smd(dat, outcomes, "Species", ci = "ct")
b <- multi.robust.smd(dat, rev(outcomes), "Species", method = "hc2", ci = "ct")
stopifnot(is.null(attr(a, "pooled")), is.null(attr(a, "pool.settings")),
          !any(c("pool.rho", "pool.cor") %in% names(formals(multi.robust.smd))))
file <- tempfile(fileext = ".pdf")
grDevices::pdf(file, width = 9, height = 6)
ordinary <- plot(a, b, pretty = FALSE)
shaded <- plot(a, b, pretty = TRUE)
# Old saved objects must not reintroduce a pooled footer.
old <- a
attr(old, "pool.settings") <- list(pool.rho = 0.6, alpha = 0.05)
attr(old, "pooled") <- data.frame(estimate = 1, lcl = 0, ucl = 2)
old.plot <- plot(old, pretty = TRUE)
stopifnot(is.null(attr(ordinary, "pooled")), is.null(attr(shaded, "pooled")),
          is.null(attr(old.plot, "pooled")),
          identical(attr(ordinary, "xlim"), attr(shaded, "xlim")),
          all(attr(shaded, "reference.drawn") == 0))

# Independent-study pooled plots must match the actual metafor fits.
studies <- data.frame(study = LETTERS[1:4], m1 = c(.1, .7, -.2, .4), m2 = 0,
                       sd1 = c(1, 1.4, 1.2, .8), sd2 = c(1.3, 1, .9, 1.2),
                       n1 = c(40, 60, 80, 35), n2 = c(45, 55, 90, 40))
m1 <- meta.robust.smd(studies, study = "study", method = "hc2", ci = "ct", alpha = .1)
m2 <- meta.robust.smd(studies[4:1, ], study = "study", method = "hc3", ci = "ct")
meta.plot <- plot(m1, m2, pretty = TRUE)
pooled <- attr(meta.plot, "pooled")
for (i in 1:2) {
  model <- attr(list(m1, m2)[[i]], "meta.models")$smd
  stopifnot(isTRUE(all.equal(pooled$estimate[i], as.numeric(model$b))),
            isTRUE(all.equal(pooled$lower[i], model$ci.lb)),
            isTRUE(all.equal(pooled$upper[i], model$ci.ub)))
}
grDevices::dev.off()
unlink(file)
