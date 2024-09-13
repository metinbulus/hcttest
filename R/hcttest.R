# install.packages("sadists")
library(sadists)

############################
# Welch degrees of freedom #
# Validation: TRUE         #
############################
.welch_df <- function(gamma1, gamma2, f1, f2, sd1, sd2) {
  (gamma1 * sd1^2 + gamma2 * sd2^2)^2 /
    (gamma1^2 * sd1^4 / f1 + gamma2^2 * sd2^4 / f2 )
}


# Lecoutre, B. (2007). Another look at the confidence intervals for the noncentral T distribution. 
# Journal of Modern Applied Statistical Methods, 6, 107-116.
.confint_nct <- function(ncp, df, alpha) {
  
  lcl.t <- sadists::qlambdap(p = alpha, df = df, t = ncp)
  ucl.t <- sadists::qlambdap(p = 1 - alpha, df = df, t = ncp)
  c(lcl.t, ucl.t)
  
} # .confint_nct

# pooled standard deviation
.pooled_sd <- function(sd1, sd2, n1, n2) {
  
  sqrt(
    ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  )
  
}


##################################
# Welch's t independent samples  # 
# Validation: TRUE               #
##################################
# constructs confidence intervals 
# from non-central t distribution
hc.t.test.default <- function(x = NULL, y = NULL, 
                              m1, m2, sd1, sd2, n1, n2,
                              alpha = 0.05, two.tailed = TRUE, 
                              ci = c("ct", "nct", "FALSE"),
                              method = c("constant", "student", "sd1", "sd2",
                                         "hc0", "hc1", "hc2", "welch",
                                         "hc3", "hc4", "hc4m", "hc5"),
                              df = c("regression", "welch-satterthwaite"),
                              verbose = TRUE, ...) { 
  
  method <- match.arg(method)
  df <- match.arg(df)
  ci <- as.character(ci)
  ci <- match.arg(ci)
  
  if(!is.null(x) & !is.null(y)) {
    
    m1 = mean(x)
    m2 = mean(y)
    sd1 = sd(x)
    sd2 = sd(y)
    n1 = length(x)
    n2 = length(y)
    
  }
  
  # global parameters
  delta <- (m1 - m2) # raw mean difference
  ifelse(two.tailed,
         alpha.for.test <- alpha / 2,
         alpha.for.test <- alpha)
  alpha.for.power <- alpha
  
  # sampling variance 
  if(method == "constant" | method == "student") {
    
    p <- 2
    var.pooled <- ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - p)
    gamma1 <- 1 / n1
    gamma2 <- 1 / n2
    var.delta <- var.pooled * (gamma1 + gamma2)
    
  } else if(method == "hc0") {
    
    gamma1 <- (n1 - 1) / n1^2
    gamma2 <- (n2 - 1) / n2^2
    var.delta <- gamma1 * sd1^2 + gamma2 * sd2^2
    
  } else if(method == "hc1") { 
    
    p <- 2
    gamma1 <- ((n1 + n2) / (n1 + n2 - p)) * ((n1 - 1) / n1^2)
    gamma2 <- ((n1 + n2) / (n1 + n2 - p)) * ((n2 - 1) / n2^2)
    var.delta <- gamma1 * sd1^2 + gamma2 * sd2^2
    
  } else if(method == "hc2" | method == "welch") { 
    
    gamma1 <- 1 / n1
    gamma2 <- 1 / n2
    var.delta <- gamma1 * sd1^2 + gamma2 * sd2^2
    
  } else if(method == "hc3") {
    
    gamma1 <- 1 / (n1 - 1)
    gamma2 <- 1 / (n2 - 1)
    var.delta <- gamma1 * sd1^2 + gamma2 * sd2^2
    
  } else if(method == "hc4") {
    
    p <- 2
    delta1 <- min(4, (n1 + n2) / (p*n1))
    delta2 <- min(4, (n1 + n2) / (p*n2))
    gamma1 <- ((n1 - 1) / n1)^(1 - delta1) * (1 / n1)
    gamma2 <- ((n2 - 1) / n2)^(1 - delta2) * (1 / n2)
    var.delta <- gamma1 * sd1^2 + gamma2 * sd2^2
    
  } else if(method == "hc4m") {
    
    p <- 2
    xi1 <- 1
    xi2 <- 1.5
    delta1 <- min(xi1, (n1 + n2) / (n1 * p)) + min(xi2, (n1 + n2) / (n1 * p))
    delta2 <- min(xi1, (n1 + n2) / (n2 * p)) + min(xi2, (n1 + n2) / (n2 * p))
    gamma1 <- ((n1 - 1) / n1)^(1 - delta1) * (1 / n1)
    gamma2 <- ((n2 - 1) / n2)^(1 - delta2) * (1 / n2)
    var.delta <- gamma1 * sd1^2 + gamma2 * sd2^2
    
  } else if(method == "hc5") {
    
    p <- 2
    k <- 0.7
    delta1 <- min((n1 + n2) / (p * n1), 
                  max(4, k * max(1 / n1, 1 / n2) / (p / (n1 + n2)))) 
    delta2 <- min((n1 + n2) / (p * n2), 
                  max(4, k * max(1 / n1, 1 / n2) / (p / (n1 + n2)))) 
    gamma1 <- ((n1 - 1) / n1)^(1 - delta1 / 2) * (1 / n1)
    gamma2 <- ((n2 - 1) / n2)^(1 - delta2 / 2) * (1 / n2)
    var.delta <- gamma1 * sd1^2 + gamma2 * sd2^2
    
  } else if (method == "sd1") {
    
    gamma1 <- 1 / n1
    gamma2 <- 1 / n2
    var.delta <- sd1^2 * (gamma1 + gamma2) 
    
  } else if(method == "sd2") {
    
    gamma1 <- 1 / n1
    gamma2 <- 1 / n2
    var.delta <- sd2^2 * (gamma1 + gamma2) 
    
  } else {
    
    stop("Unknown method request", call. = FALSE)
    
  }
  
  # degrees of freedom 
  p <- 2
  if(df == "welch-satterthwaite") {
    def <- .welch_df(gamma1 = gamma1, gamma2 = gamma2,
                     f1 = n1 - 1, f2 = n2 - 1, 
                     sd1 = sd1, sd2 = sd2)
    if(method == "sd1" | method == "sd2") {
      stop("Welch-Satterthwaite degrees of freedom is not applicable to 'sd1' or 'sd2' testing methods", call. = FALSE)
    }
  } else {
    if(method == "welch") {
      stop("Welch-Satterthwaite test requires Welch-Satterthwaite degrees of freedom", call. = FALSE)
    }
    if (method == "sd1") {
      def <- n1 - 1
    } else if(method == "sd2") {
      def <- n2 - 1
    } else {
      def <- n1 + n2 - p
    }
  }
  
  
  # confidence intervals 
  if(ci == "nct") {
    
    V <- delta / sqrt(var.delta)
    confint.t <- .confint_nct(ncp = V, df = def, alpha.for.test)
    confint.delta <- confint.t * sqrt(var.delta)
    
  } else if(ci == "ct") {
    
    ct <- qt(p = alpha.for.test, ncp = 0, df = def, log.p = FALSE, lower.tail = FALSE)
    confint.delta <- c(delta - ct * sqrt(var.delta), delta + ct * sqrt(var.delta))
    
  } else {
    
    confint.delta <- c(NA, NA)
    
  } 
  
  # get p values
  t.observed <- delta / sqrt(var.delta)
  if(ci == "ct") {
    
    p.value <- pt(q = abs(t.observed), ncp = 0, df = def, lower.tail = FALSE)
    if(two.tailed) p.value <- 2 * p.value
    
  } else {
    
    p.value <- NA
    
  }
  
  test <- data.frame(delta = delta,
                     se = sqrt(var.delta),
                     t = t.observed,
                     df = def,
                     p = p.value,
                     lcl = confint.delta[[1]],
                     ucl = confint.delta[[2]])
  
  parms <- data.frame(m1 = m1, m2 = m2,
                      sd1 = sd1, sd2 = sd2,
                      n1 = n1, n2 = n2, 
                      alpha = alpha, df = df,
                      two.tailed = two.tailed, 
                      ci = ci, method = method)
  
  if(verbose) {
    
    print(
      
      round(test, 3)
      
    )
    
    if(is.na(p.value)) message("P-value is not calculated, use confidence interval instead")
    
    if(any(is.na(confint.delta))) message("Confidence interval is not calculated")
    
  }
  
  invisible(
    
    structure(
      list(test = test, parms = parms),
      class = "hc.t.test")
    
  )
  
} # hc.t.test.default()


hc.t.test.formula <- function(formula, data, ...) {
  
  formula <- as.formula(formula)
  vars <- all.vars(formula)
  y <- vars[1]
  x <- vars[2]
  
  x_levels <- unique(data[[x]])
  if (length(x_levels) > 2) {
    stop(paste("The predictor variable '", x, "' has more than two levels. It should be binary."))
  }
  
  # Ensure x is a factor
  if (!is.factor(data[[x]])) {
    data[[x]] <- as.factor(data[[x]])
  }
  
  # Split y based on levels of x
  y1 <- data[[y]][data[[x]] == levels(data[[x]])[1]]
  y2 <- data[[y]][data[[x]] == levels(data[[x]])[2]]
  
  parms <- list(m1 = mean(y1), m2 = mean(y2),
                sd1 = sd(y1), sd2 = sd(y2),
                n1 = length(y1), n2 = length(y2),
                ...)
  
  do.call("hc.t.test", parms)
  
}


hc.t.test <- function(x, ...) {
  UseMethod("hc.t.test")
}

smd.hc.t.test <- function(x = NULL,
                          ss.adjust = FALSE, 
                          ci = c("nct", "z", "ct", "FALSE"),
                          verbose = TRUE) {
  
  if(!inherits(x, "hc.t.test")) 
    stop("Please provide an object retruned from 'hc.t.test()' function",
         call. = FALSE)
  
  ci <- as.character(ci)
  ci <- match.arg(ci)
  
  n1 <- x$parms$n1
  n2 <- x$parms$n2
  m1 <- x$parms$m1
  m2 <- x$parms$m2
  sd1 <- x$parms$sd1
  sd2 <- x$parms$sd2
  alpha <- x$parms$alpha
  two.tailed <- x$parms$two.tailed
  method <- x$parms$method
  
  df <- x$test$df
  delta <- x$test$delta
  se.delta <- x$test$se
  var.delta <- se.delta^2
  
  
  ifelse(two.tailed,
         alpha.for.test <- alpha / 2,
         alpha.for.test <- alpha)
  alpha.for.power <- alpha
  
  
  if(tolower(method) == "sd1") {
    std <- sd1 
    var.std <- sd1^2 / (2 * (n1 - 1))
  } else if(tolower(method) == "sd2") {
    std <- sd2
    var.std <- sd2^2 / (2 * (n2 - 1))
  } else {
    std <- .pooled_sd(sd1 = sd1, sd2 = sd2, n1 = n1, n2 = n2)
    var.std <- ((n1 - 1) * sd1^4 + (n2 - 1) * sd2^4) / ((2 * (n1 - 1) * sd1^2 + 2 * (n2 - 1) * sd2^2) * (n1 + n2 - 2))
  }
  
  ifelse(ss.adjust, 
         ss.adj <- exp(lgamma(df / 2) - log(sqrt(df / 2)) - lgamma((df - 1) / 2)),
         ss.adj <- 1)
  
  d <- (delta / std) * ss.adj
  var.d <- var.delta / std^2 + d^2 * var.std / std^2
  se.d <- sqrt(var.d)
  
  
  # confidence intervals 
  if(tolower(ci) == "nct") {
    
    V <- delta / se.delta
    confint.t <- .confint_nct(ncp = V, df = df, alpha.for.test)
    confint.d <- confint.t * (se.delta / std) * ss.adj 
    
  } else if(tolower(ci) == "ct") {
    
    ct <- qt(p = alpha.for.test, ncp = 0, df = df, log.p = FALSE, lower.tail = FALSE)
    confint.d <- c(d - ct * se.d, d + ct * se.d)
    
  } else if(tolower(ci) == "z") {
    
    cz <- qnorm(p = alpha.for.test, mean = 0, sd = 1, lower.tail = FALSE)
    confint.d <- c(d - cz * se.d, d + cz * se.d)
    
  } else {
    
    confint.d <- c(NA, NA)
    
  } 
  
  
  es <- data.frame(d = d,
                   se = se.d,
                   lcl = confint.d[[1]],
                   ucl = confint.d[[2]])
  
  parms <- x$parms
  
  if(verbose) {
    print(
      round(es, 3)
    )
  }
  
  invisible(list(es = es, parms = parms))
  
} # smd.hc.t.test()
