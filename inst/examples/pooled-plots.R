library(hcttest)

# Dependent outcomes from the same people: individual outcome estimates.
dat <- subset(iris, Species != "virginica")
outcomes <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
constant <- multi.robust.smd(dat, outcomes, "Species", method = "constant",
                             ci = "ct")
hc2 <- multi.robust.smd(dat, outcomes, "Species", method = "hc2",
                        ci = "ct")
print(constant[, c("outcome", "yi", "lcl.equivalence", "ucl.equivalence")], row.names = FALSE)
print(hc2[, c("outcome", "yi", "lcl.equivalence", "ucl.equivalence")], row.names = FALSE)
plot(constant, hc2, pretty = TRUE, col = c("#001F3F", "#B35400"))

# Independent studies: pooled summaries come from metafor, retaining model CIs.
studies <- data.frame(study = c("A", "B", "C", "D"),
                       m1 = c(0.1, 0.7, -0.2, 0.4), m2 = 0,
                       sd1 = c(1, 1.4, 1.2, 0.8), sd2 = c(1.3, 1, 0.9, 1.2),
                       n1 = c(40, 60, 80, 35), n2 = c(45, 55, 90, 40))
meta1 <- meta.robust.smd(studies, study = "study", method = "constant", ci = "ct")
meta2 <- meta.robust.smd(studies, study = "study", method = "hc2", ci = "ct")
plot(meta1, meta2, pretty = TRUE, col = c("#001F3F", "#B35400"))
print(attr(meta1, "pooled"), row.names = FALSE)
print(attr(meta2, "pooled"), row.names = FALSE)
