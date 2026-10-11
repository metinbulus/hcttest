library(hcttest)

# Pooling is checked against the established metafor implementation.
studies <- data.frame(study = c("A", "B", "C", "D"),
                      m1 = c(0.1, 1.2, -0.5, 0.8), m2 = 0,
                      sd1 = c(1, 1.8, 1.2, 0.8), sd2 = c(1.3, 1, 0.9, 1.2),
                      n1 = c(25, 60, 150, 35), n2 = c(30, 55, 140, 40))
for (estimator in c("REML", "DL", "PM")) {
  object <- meta.robust.smd(studies, study = "study", method = "hc2",
                            ci = "ct", pool.method = estimator)
  reference <- metafor::rma.uni(object$yi, object$vi, method = estimator, test = "knha")
  pooled <- attr(object, "pooled")
  pooled <- pooled[pooled$effect == "smd", ]
  stopifnot(isTRUE(all.equal(pooled$estimate, as.numeric(reference$b))),
            isTRUE(all.equal(pooled$tau2, reference$tau2)),
            isTRUE(all.equal(object$weight, 1 / (object$vi + reference$tau2))))
}

# Raw HC2/Welch equivalence is checked against base R's independent t tests.
set.seed(2026)
a <- rnorm(90, 0.1, 1)
b <- rnorm(65, 0, 1.5)
dat <- data.frame(group = factor(c(rep("A", length(a)), rep("B", length(b)))),
                  outcome = c(a, b))
multi <- multi.robust.smd(dat, "outcome", "group", method = "hc2",
                          df = "welch-satterthwaite", ci = "nct",
                          margins = c(-0.5, 0.5), equivalence.effect = "raw")
lower.test <- t.test(a, b, mu = -0.5, alternative = "greater")
upper.test <- t.test(a, b, mu = 0.5, alternative = "less")
interval <- t.test(a, b, conf.level = 0.9)$conf.int
stopifnot(isTRUE(all.equal(multi$p.lower, lower.test$p.value)),
          isTRUE(all.equal(multi$p.upper, upper.test$p.value)),
          isTRUE(all.equal(unname(c(multi$lcl.equivalence, multi$ucl.equivalence)),
                           as.numeric(interval))),
          identical(multi$equivalent, multi$lcl.equivalence > -0.5 &
                    multi$ucl.equivalence < 0.5))

# Plot return values verify references, matching, and consistent axis ranges.
second <- meta.robust.smd(studies[4:1, ], study = "study", method = "hc3", ci = "ct")
plot.file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot.file, width = 9, height = 6)
ordinary <- plot(object, second, pretty = FALSE)
shaded <- plot(object, second, pretty = TRUE)
balance <- plot(multi, effect = "raw", pretty = TRUE, margins = c(-0.5, 0.5))
subset.plot <- plot(object[1:3, ], pretty = TRUE)
grDevices::dev.off()
unlink(plot.file)
expected <- c(attr(object, "pooled")$estimate[1], attr(second, "pooled")$estimate[1])
subset.model <- metafor::rma.uni(object$yi[1:3], object$vi[1:3], method = "PM", test = "knha")
stopifnot(isTRUE(all.equal(attr(shaded, "reference.drawn"), expected)),
          identical(attr(ordinary, "xlim"), attr(shaded, "xlim")),
          isTRUE(all.equal(balance$lower, multi$lcl.equivalence.raw)),
          identical(attr(balance, "reference.drawn"), 0),
          isTRUE(all.equal(attr(subset.plot, "reference.drawn"), as.numeric(subset.model$b))))

# Equivalence inherits the SMD method, including trim calibration and correction.
for (choice in c("z", "ct", "nct")) {
  fit <- hct.test(a, b, method = "hc2", trim = 0.2, alpha = 0.1,
                  verbose = FALSE)
  reference <- robust.smd(fit, ci = choice, ss.adjust = TRUE, verbose = FALSE)$es
  object <- multi.robust.smd(dat, "outcome", "group", method = "hc2",
                             trim = 0.2, ci = choice, ss.adjust = TRUE,
                             margins = c(-0.5, 0.5))
  stopifnot(isTRUE(all.equal(object$lcl.equivalence, reference$lcl)),
            isTRUE(all.equal(object$ucl.equivalence, reference$ucl)),
            identical(object$equivalent, object$lcl.equivalence > -0.5 &
                        object$ucl.equivalence < 0.5))
  if (choice == "nct") {
    multiplier <- object$yi / object$t
    stopifnot(isTRUE(all.equal(object$p.lower,
                pt(object$t, object$df, ncp = -0.5 / multiplier, lower.tail = FALSE))),
              isTRUE(all.equal(object$p.upper,
                pt(object$t, object$df, ncp = 0.5 / multiplier))),
              abs(pt(object$t, object$df, ncp = object$lcl.equivalence / multiplier,
                       lower.tail = FALSE) - 0.05) < 1e-8,
              abs(pt(object$t, object$df, ncp = object$ucl.equivalence / multiplier,
                       lower.tail = TRUE) - 0.05) < 1e-8)
  }
}
stopifnot(!"equivalence.ci" %in% names(formals(multi.robust.smd)))
# Zero observed contrast must not cause division by zero in nct bound tests.
zero <- data.frame(group = factor(rep(c("A", "B"), each = 10)),
                   outcome = rep(seq_len(10), 2))
zero.fit <- multi.robust.smd(zero, "outcome", "group", ci = "nct",
                             margins = c(-0.5, 0.5))
stopifnot(is.finite(zero.fit$p.lower), is.finite(zero.fit$p.upper),
          isTRUE(all.equal(zero.fit$p.lower, zero.fit$p.upper)))
