# Primary methodological reference:
# Bulus, M. (2026). Robust standard errors and confidence intervals for standardized mean difference. Research Synthesis Methods. Advance online publication.
# Full package citation is available through citation("hcttest").

# rmixnorm <- function(n, means, sds, probs) {
#  components <- sample(length(means), size = n, replace = TRUE, prob = probs)
#  rnorm(n, mean = means[components], sd = sds[components])
# }

# winsorize trim % from each side
winsor <- function(x, trim = 0.10, type = 7, na.rm = FALSE) {
  n <- length(x)
  q.trim.left <- quantile(x = x, probs = trim, type = type, na.rm = na.rm)
  q.trim.right <- quantile(x = x, probs = 1 - trim, type = type, na.rm = na.rm)
  x[x < q.trim.left] <- q.trim.left
  x[x > q.trim.right] <- q.trim.right
  x
} # winsor 

# remove trim % from each side
trim <- function(x, trim = 0.10, na.rm = FALSE) {
  if (!is.numeric(trim) || length(trim) != 1L ||
      !is.finite(trim) || trim < 0 || trim >= 0.5) {
    stop("'trim' must be a finite number in [0, 0.5).", call. = FALSE)
  }
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  n.trim <- round(n * trim)
  if (n.trim == 0L) return(x)
  if (2 * n.trim >= n) {
    stop("Trimming would leave no observations.", call. = FALSE)
  }
  sort.idx <- order(x, na.last = TRUE)
  keep <- seq.int(n.trim + 1L, n - n.trim)
  x[sort.idx[keep]]
} # trim


dmixnormal <- function(x, mu, sigma, weight) {
  # sum(weight * dnorm(x, mu, sigma))
  KScorrect::dmixnorm(x = x, mean = mu, sd = sigma, pro = weight)
}
# dmixnormal_vec <- Vectorize(dmixnormal, vectorize.args = "x")
pmixnormal <- function(q, mu, sigma, weight) {
  # sum(weight * pnorm(q, mu, sigma))
  KScorrect::pmixnorm(q = q, mean = mu, sd = sigma, pro = weight)
}
# pmixnormal_vec <- Vectorize(pmixnormal, vectorize.args = "q")
qmixnormal <- function(p, mu, sigma, weight) {
  ## Use numerical search to find quantile
  # objective <- function(x) pmixnormal(x, mu, sigma, weight) - p
  ## Find reasonable bounds
  # lower <- min(qnorm(max(p - 0.01, 0.0001), mu, sigma))
  # upper <- max(qnorm(min(p + 0.01, 0.9999), mu, sigma))
  # tryCatch({
  #   uniroot(objective, c(lower, upper), tol = 1e-8)$root
  # }, error = function(e) {
  ## If uniroot fails, use optimization
  #   optimize(function(x) abs(objective(x)), c(lower, upper))$minimum
  # })
  # 
  KScorrect::qmixnorm(p = p, mean = mu, sd = sigma, pro = weight)
}
# qmixnormal_vec <- Vectorize(qmixnormal, vectorize.args = "p")
rmixnormal <- function(n, mu, sigma, weight) {
  # components <- sample(length(mu), size = n, replace = TRUE, prob = weight)
  # rnorm(n, mean = mu[components], sd = sigma[components])
  as.vector(KScorrect::rmixnorm(n = n, mean = mu, sd = sigma, pro = weight))
}



qbbinom <- function(p, size, alpha, beta) {
  # Handle cases where p is a vector
  sapply(p, function(prob) {
    if (is.na(prob)) return(NA)
    if (prob < 0 || prob > 1) stop("p must be between 0 and 1")
    
    # Generate all possible scores for this test size
    x_range <- 0:size
    
    # Calculate the cumulative probabilities
    probs <- extraDistr::pbbinom(x_range, size = size, alpha = alpha, beta = beta)
    
    # Find the smallest x where the cumulative probability is >= p
    return(min(x_range[probs >= prob]))
  })
}

# var cannot be larger than n_items * p_bar * (1 - p_bar)

rbetabinom_mean_rel <- function(n_students = 30, n_items = 10, target_mean, target_rel, draw_data = FALSE) {
  
  # 1. Target proportion (mu)
  mu <- target_mean / n_items
  
  # 2. Solve for S (alpha + beta)
  # Based on: rel = [n / (n-1)] * [1 - (mu*(1-mu) / var_p)]
  # For Beta-Binomial, this simplifies to:
  S <- (n_items * (1 - target_rel)) / target_rel - 1
  
  # 3. Calculate Alpha and Beta
  alpha_param <- mu * S
  beta_param <- (1 - mu) * S
  
  # 4. Calculate Rho (Intra-item correlation)
  # rho = 1 / (S + 1)
  rho_calc <- 1 / (S + 1)
  
  # 5. Calculate the Theoretical SD
  # Formula for Beta-Binomial Variance: 
  # Var = n * mu * (1-mu) * [1 + (n-1) * rho]
  theo_var <- n_items * mu * (1 - mu) * (1 + (n_items - 1) * rho_calc)
  theo_sd <- sqrt(theo_var)
  
  if (alpha_param <= 0 || beta_param <= 0 || rho_calc <= 0 || rho_calc >= 1) {
    stop("Infeasable beta-binomial parameters for this mean, SD, and test length.")
  }
  
  # 6. Simulate
  if(draw_data) {
    true_abilities <- rbeta(n_students, alpha_param, beta_param)
    scores <- rbinom(n_students, n_items, true_abilities)
  } else {
    scores <- NA
  }
  
  
  return(list(
    x = scores, 
    mean = target_mean,
    sd = theo_sd,
    alpha = alpha_param, 
    beta = beta_param,
    rho = rho_calc,
    reliability = target_rel
  ))
} # rbetabinom_mean_rel


rbetabinom_mean_sd <- function(n_students = 30, n_items = 10, target_mean, target_sd, draw_data = FALSE) {
  
  # 1. Convert targets to proportion scale (0 to 1)
  mu <- target_mean / n_items
  target_var_p <- (target_sd^2) / (n_items^2)
  
  # 2. Calculate Binomial noise (minimum possible variance)
  # If target_sd is lower than this, the simulation is mathematically impossible
  min_var_p <- (mu * (1 - mu)) / n_items
  
  if (target_var_p <= min_var_p) {
    stop(paste("Target SD is too low. Minimum SD for this mean is:", 
               round(sqrt(min_var_p) * n_items, 3)))
  }
  
  # 3. Solve for Rho (Overdispersion)
  # Rho is the correlation between items, also the "intra-class correlation"
  rho_calc <- (target_var_p - min_var_p) / ((mu * (1 - mu) / n_items) * (n_items - 1))
  
  # 4. Calculate Alpha and Beta
  alpha_param <- mu * (1 / rho_calc - 1)
  beta_param <- (1 - mu) * (1 / rho_calc - 1)
  
  # 6. Calculate Reliability (KR-21/Alpha equivalent)
  reliability <- (n_items / (n_items - 1)) * (1 - (target_mean * (n_items - target_mean)) / (n_items * target_sd^2))
  
  if (alpha_param <= 0 || beta_param <= 0 || rho_calc <= 0 || reliability <= 0 || rho_calc >= 1 || reliability >= 1) {
    stop("Infeasable beta-binomial parameters for this mean, SD, and test length.")
  }
  
  if(draw_data) { 
    # 5. Simulate
    # Generate true abilities from Beta
    true_abilities <- rbeta(n_students, shape1 = alpha_param, shape2 = beta_param)
    # Generate observed scores from Binomial
    scores <- rbinom(n_students, size = n_items, prob = true_abilities)
  } else {
    scores <- NA
  }
  
  
  return(list(
    x = scores, 
    mean = target_mean,
    sd = target_sd,
    alpha = alpha_param, 
    beta = beta_param,
    rho = rho_calc,
    reliability = reliability
  ))
} # rbetabinom_mean_sd

simparms_betabinom <- function(delta = 1, mean = 7, n_items = 10, ss = 30, ss.ratio = 1, rel = 0.60) {
  
  parms1 <- rbetabinom_mean_rel(n_items = n_items , target_mean = mean + delta, target_rel = rel)
  parms2 <- rbetabinom_mean_rel(n_items = n_items , target_mean = mean, target_rel = rel)
  
  sd1 <- parms1$sd
  sd2 <- parms2$sd
  var.ratio <- sd1^2 / sd2^2
  
  n2 <- ss / (ss.ratio + 1)
  n1 <- ss - n2
  
  d  <- delta / .pooled_sd(sd1 = sd1, sd2 = sd2, n1 = n1, n2 = n2) # rel 0.50
  return(data.frame(es = d, ss = ss, ss.ratio = ss.ratio, var.ratio = var.ratio, 
                    mu1 = mean + delta, mu2 = mean, sd1 = sd1, sd2 = sd2, n1 = n1, n2 = n2, rel = rel))
  
  
}

# Wrapper function to solve for delta
solve_for_delta <- function(target_d, mean = 7, n_items = 10, ss = 30, ss.ratio = 1, rel = 0.60) {
  
  if(target_d == 0) {
    
    return(0)
    
  } else {
    
    # We search for delta in a reasonable range (e.g., between 0 and 10)
    # 10 is usually safe for a 10-item test since delta can't exceed n_items - mean
    result <- uniroot(
      f = function(x) simparms_betabinom(delta = x,
                                         mean = mean,
                                         n_items =  n_items,
                                         ss = ss,
                                         ss.ratio = ss.ratio,
                                         rel = rel)$es - target_d, 
      interval = c(0.0001, n_items - mean - 0.0001) # range for delta; upper bound depends on test length
    )
    
    return(result$root)
    
  }
  
} # solve_for_delta



trimmed_mean_bbinom <- function(size, alpha, beta, trim = 0.10) {
  
  if (trim == 0) return(size * (alpha / (alpha + beta)))
  
  x <- 0:size
  weights <- dbbinom(x, size, alpha, beta)
  
  # Find indices to trim
  cum_probs <- cumsum(weights)
  low_idx <- min(which(cum_probs > trim))
  high_idx <- max(which(cum_probs <= (1 - trim))) + 1
  
  # Adjust the 'partial' bars at the cutoffs
  # Lower tail adjustment
  weights[low_idx] <- cum_probs[low_idx] - trim
  if (low_idx > 1) weights[1:(low_idx-1)] <- 0
  
  # Upper tail adjustment
  weights[high_idx] <- (1 - trim) - (cum_probs[high_idx] - weights[high_idx])
  if (high_idx < length(x)) weights[(high_idx+1):length(x)] <- 0
  
  # Calculate weighted mean of the remaining mass
  # We divide by (1 - 2*trim) to re-normalize the remaining 80% (or whatever remains)
  t_mean <- sum(x * weights) / (1 - 2 * trim)
  
  return(t_mean)
} # trimmed_mean_bbinom


winsorized_sd_bbinom <- function(size, alpha, beta, trim = 0.10) {
  
  x <- 0:size
  probs <- dbbinom(x, size, alpha, beta)
  k <- trim
  
  if (k == 0) {
    # Standard SD if no trimming
    mu <- sum(x * probs)
    v <- sum((x - mu)^2 * probs)
    return(list(win_mean = mu, win_sd = sqrt(v)))
  }
  
  cum_probs <- cumsum(probs)
  
  # Identify indices of the quantiles
  low_idx  <- min(which(cum_probs > k))
  high_idx <- min(which(cum_probs >= (1 - k)))
  
  # Initialize Winsorized Weights
  win_weights <- probs
  
  # 1. Collapse the Lower Tail into low_idx
  # Everything below low_idx is set to 0, and its mass moves to low_idx
  if (low_idx > 1) {
    mass_below <- cum_probs[low_idx - 1]
    win_weights[1:(low_idx - 1)] <- 0
    win_weights[low_idx] <- win_weights[low_idx] + mass_below
  }
  
  # 2. Collapse the Upper Tail into high_idx
  # Everything above high_idx is set to 0, and its mass moves to high_idx
  if (high_idx < length(x)) {
    mass_above <- 1 - cum_probs[high_idx]
    win_weights[(high_idx + 1):length(x)] <- 0
    win_weights[high_idx] <- win_weights[high_idx] + mass_above
  }
  
  # 3. Calculate Winsorized Mean (Theoretical)
  win_mean <- sum(x * win_weights)
  
  # 4. Calculate Winsorized Variance
  # Var = E[X^2] - (E[X])^2 using Winsorized weights
  win_var <- sum((x - win_mean)^2 * win_weights)
  
  # 5. The Winsorized SD is the square root of the variance
  # Note: In sample statistics, a correction factor is sometimes applied
  # but for the theoretical distribution, we use the root variance.
  return(sqrt(win_var))
  
} # winsorized_sd_bbinom

# 1. Trimmed Mean for Negative Binomial
trimmed_mean_nbinom <- function(mu, size, trim = 0.10) {
  
  if (trim == 0) return(mu)
  
  max_x <- qnbinom(0.99999, mu = mu, size = size)
  x <- 0:max_x
  probs <- dnbinom(x, mu = mu, size = size)
  cum_probs <- cumsum(probs)
  
  k <- trim
  
  # Find indices
  low_idx  <- min(which(cum_probs > k))
  high_idx <- min(which(cum_probs >= (1 - k)))
  
  weights <- probs
  
  # Adjust Lower Tail
  weights[low_idx] <- cum_probs[low_idx] - k
  if (low_idx > 1) weights[1:(low_idx - 1)] <- 0
  
  # Adjust Upper Tail
  if (high_idx > 1) {
    weights[high_idx] <- (1 - k) - cum_probs[high_idx - 1]
  } else {
    weights[high_idx] <- 1 - k
  }
  if (high_idx < length(x)) weights[(high_idx + 1):length(x)] <- 0
  
  # Normalize by the remaining mass
  return(sum(x * weights) / (1 - 2 * k))
} # trimmed_mean_nbinom

# 2. Winsorized SD for Negative Binomial
winsorized_sd_nbinom <- function(mu, size, trim = 0.10) {
  
  max_x <- qnbinom(0.99999, mu = mu, size = size)
  x <- 0:max_x
  probs <- dnbinom(x, mu = mu, size = size)
  cum_probs <- cumsum(probs)
  
  if (trim == 0) {
    std_sd <- sqrt(mu + mu^2/size) # Theoretical NB variance
    return(list(win_mean = mu, win_sd = std_sd))
  }
  
  k <- trim
  low_idx  <- min(which(cum_probs > k))
  high_idx <- min(which(cum_probs >= (1 - k)))
  
  win_weights <- probs
  
  # Collapse Lower Tail
  if (low_idx > 1) {
    mass_below <- cum_probs[low_idx - 1]
    win_weights[1:(low_idx - 1)] <- 0
    win_weights[low_idx] <- win_weights[low_idx] + mass_below
  }
  
  # Collapse Upper Tail
  if (high_idx < length(x)) {
    mass_above <- 1 - cum_probs[high_idx]
    win_weights[(high_idx + 1):length(x)] <- 0
    win_weights[high_idx] <- win_weights[high_idx] + mass_above
  }
  
  win_mean <- sum(x * win_weights)
  win_var <- sum((x - win_mean)^2 * win_weights)
  
  return(sqrt(win_var))
} # winsorized_sd_nbinom 

# 1. Trimmed Mean for Poisson
# trimmed_mean_pois(lambda = sd1^2, offset = sd1^2 - mu1, trim = trim.prop)
# winsorized_sd_pois(lambda = sd1^2, offset = sd1^2 - mu1, trim = trim.prop)

trimmed_mean_pois <- function(lambda, offset, trim = 0.10) {
  
  max_x <- max(10, qpois(0.99999, lambda = lambda) ) 
  x <- 0:max_x
  weights <- dpois(x, lambda = lambda)
  cum_probs <- cumsum(weights)
  
  # Find indices where the thresholds fall
  low_idx  <- min(which(cum_probs > trim))
  high_idx <- min(which(cum_probs >= (1 - trim)))
  
  # Proportional Adjustment: Lower Tail
  weights[low_idx] <- cum_probs[low_idx] - trim
  if (low_idx > 1) weights[1:(low_idx - 1)] <- 0
  
  # Proportional Adjustment: Upper Tail
  if (high_idx > 1) {
    weights[high_idx] <- (1 - trim) - cum_probs[high_idx - 1]
  } else {
    weights[high_idx] <- 1 - trim
  }
  if (high_idx < length(x)) weights[(high_idx + 1):length(x)] <- 0
  
  # Normalize by the remaining mass (1 - 2k)
  return(sum((x - offset) * weights) / (1 - 2 * trim))
} # trimmed_mean_pois

# 2. Winsorized Stats for Poisson
winsorized_sd_pois <- function(lambda, offset, trim = 0.10) {
  
  max_x <- max(10, qpois(0.99999, lambda = lambda)) 
  x <- 0:max_x
  weights <- dpois(x, lambda = lambda)
  cum_probs <- cumsum(weights)
  
  low_idx  <- min(which(cum_probs > trim))
  high_idx <- min(which(cum_probs >= (1 - trim)))
  
  # Collapse Lower Tail into the low cutoff
  if (low_idx > 1) {
    mass_below <- cum_probs[low_idx - 1]
    weights[1:(low_idx - 1)] <- 0
    weights[low_idx] <- weights[low_idx] + mass_below
  }
  
  # Collapse Upper Tail into the high cutoff
  if (high_idx < length(x)) {
    mass_above <- 1 - cum_probs[high_idx]
    weights[(high_idx + 1):length(x)] <- 0
    weights[high_idx] <- weights[high_idx] + mass_above
  }
  
  win_mean <- sum((x - offset) * weights)
  win_var <- sum((x - offset - win_mean)^2 * weights)
  
  return(sqrt(win_var))
} # winsorized_sd_pois



trim.pop.mean <- function(mu,
                          sigma, 
                          trim, 
                          size, # for beta-binomial
                          weight = 0.50,
                          distribution = c("norm", "mixnorm", "betabin", "negbin", "pois", "exp", "uni", "lnorm")) {
  
  distribution <- match.arg(distribution)
  
  # The trimmed region goes from Q(trim) to Q(1-trim)
  
  # generate variables depending on the distribution type
  # given means and standard deviations
  if(distribution == "normal" || distribution == "norm") {
    
    lower <- qnorm(p = trim, mean = mu, sd = sigma)
    upper <- qnorm(p = 1 - trim, mean = mu, sd = sigma)
    
    integrand <- function(x) {
      x * dnorm(x, mean = mu, sd = sigma)
    }
    
    trimmed_mean <- (1 / (1 - 2*trim)) * integrate(integrand, lower = lower, upper = upper)$value
    
  } else if(distribution == "mixnorm" || distribution == "mnorm") {
    
    if(length(mu) == 1) mu <- c(mu, mu)
    if(length(sigma) == 1) sigma <- c(sigma, sigma)
    if(length(weight) == 1) weight <- c(weight, 1 - weight)
    
    lower <- qmixnormal(p = trim, mu = mu, sigma = sigma, weight = weight)
    upper <- qmixnormal(p = 1 - trim, mu = mu, sigma = sigma, weight = weight)
    
    integrand <- function(x) {
      x * dmixnormal(x = x, mu = mu, sigma = sigma, weight = weight)
    }
    
    trimmed_mean <- (1 / (1 - 2*trim)) * integrate(integrand, lower = lower, upper = upper)$value
    
  } else if(distribution == "negbin" || distribution == "nb") {
    
    trimmed_mean <- trimmed_mean_nbinom(mu = mu, size = mu^2 / (sigma^2 - mu), trim = trim)
    
  } else if(distribution == "poisson" || distribution == "pois") {
    
    trimmed_mean <- trimmed_mean_pois(lambda = sigma^2, offset = sigma^2 - mu, trim = trim)
    
  } else if(distribution == "betabinom" || distribution == "betabin" || distribution == "bb") {
    
    bb_parms <- rbetabinom_mean_sd(n_items = size, target_mean = mu, target_sd = sigma, draw_data = FALSE)
    alpha <- bb_parms$alpha
    beta <- bb_parms$beta 
    
    trimmed_mean <- trimmed_mean_bbinom(size = size, alpha = alpha, beta = beta, trim = trim)
    
  } else if(distribution == "exponential" || distribution == "exp") {
    
    lower <- qexp(p = trim, rate  = 1 / sigma) - sigma + mu
    upper <- qexp(p = 1 - trim, rate  = 1 / sigma) - sigma + mu
    
    integrand <- function(x) {
      x * dexp(x + sigma - mu, rate  = 1 / sigma) 
    }
    
    trimmed_mean <- (1 / (1 - 2*trim)) * integrate(integrand, lower = lower, upper = upper)$value
    
  } else if(distribution == "uniform" || distribution == "uni") {
    
    a <- mu - sigma * sqrt(3)
    b <- mu + sigma * sqrt(3)
    lower <- qunif(p = trim, min = a, max = b)
    upper <- qunif(p = 1 - trim, min = a, max = b)
    
    integrand <- function(x) {
      x * dunif(x, min = a, max = b)
    }
    
    trimmed_mean <- (1 / (1 - 2*trim)) * integrate(integrand, lower = lower, upper = upper)$value
    
  } else if(distribution == "lognormal" || distribution == "lognorm" || distribution == "lnorm") {
    
    meanlog <- 2 * log(mu) - log(mu^2 + sigma^2) / 2
    sdlog <- sqrt(-2 * log(mu) + log(mu^2 + sigma^2))
    lower <- qlnorm(p = trim, meanlog = meanlog, sdlog = sdlog)
    upper <- qlnorm(p = 1 - trim, meanlog = meanlog, sdlog = sdlog)
    
    integrand <- function(x) {
      x * dlnorm(x, meanlog = meanlog, sdlog = sdlog)
    }
    
    trimmed_mean <- (1 / (1 - 2*trim)) * integrate(integrand, lower = lower, upper = upper)$value
    
  } else {
    
    stop("Distribution type not supported", call. = FALSE)
    
  } # end of distribution 
  
  # Integrate over the trimmed region
  # (1 / (1 - 2*trim)) is scaling factor
  # if(distribution == "betabinom" || distribution == "betabin" || distribution == "bb" || distribution == "poisson" || distribution == "pois" || distribution == "negbin" || distribution == "nb") {
  #  trimmed_mean <- (1 / (1 - 2*trim)) * sum(sapply(round(lower):round(upper), integrand))
  # } else {
  #   trimmed_mean <- (1 / (1 - 2*trim)) * integrate(integrand, lower = lower, upper = upper)$value
  # }
  
  return(trimmed_mean)
  
} # trim.pop.mean

winsor.pop.sd <- function(mu,
                          sigma,
                          size, # for beta-binomial
                          trim,
                          weight = 0.50,
                          distribution = c("norm", "mixnorm", "betabin", "negbin", "pois", "exp", "uni", "lnorm")) {
  
  distribution <- match.arg(distribution)
  
  # Find Winsorization points
  # generate variables depending on the distribution type
  # given means and standard deviations
  if(distribution == "normal" || distribution == "norm") {
    
    lower <- qnorm(p = trim, mean = mu, sd = sigma)
    upper <- qnorm(p = 1 - trim, mean = mu, sd= sigma)
    
    mean_integrand <- function(x) {
      x * dnorm(x, mean = mu, sd = sigma)
    }
    
    var_integrand <- function(x) {
      wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
      (x - wins_mean)^2 * dnorm(x, mean = mu, sd = sigma)
    }
    
    wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
    wins_variance <- trim * (lower - wins_mean)^2 + integrate(var_integrand, lower, upper)$value + trim * (upper - wins_mean)^2
    wins_sigma <- sqrt(wins_variance)
    
  } else if(distribution == "mixnorm" || distribution == "mnorm") {
    
    if(length(mu) == 1) mu <- c(mu, mu)
    if(length(sigma) == 1) sigma <- c(sigma, sigma)
    if(length(weight) == 1) weight <- c(weight, 1 - weight)
    
    lower <- qmixnormal(p = trim, mu = mu, sigma = sigma, weight = weight)
    upper <- qmixnormal(p = 1 - trim, mu = mu, sigma = sigma, weight = weight)
    
    mean_integrand <- function(x) {
      x * dmixnormal(x, mu = mu, sigma = sigma, weight = weight)
    }
    
    var_integrand <- function(x) {
      wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
      (x - wins_mean)^2 * dmixnormal(x = x, mu = mu, sigma = sigma, weight = weight)
    }
    
    wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
    wins_variance <- trim * (lower - wins_mean)^2 + integrate(var_integrand, lower, upper)$value + trim * (upper - wins_mean)^2
    wins_sigma <- sqrt(wins_variance)
    
  } else if(distribution == "negbin" || distribution == "nb") {
    
    wins_sigma <- winsorized_sd_nbinom(mu = mu, size = mu^2 / (sigma^2 - mu), trim = trim)
    
  } else if(distribution == "poisson" || distribution == "pois") {
    
    wins_sigma <- winsorized_sd_pois(lambda = sigma^2, offset = sigma^2 - mu, trim = trim)
    
  } else if(distribution == "betabinom" || distribution == "betabin" || distribution == "bb") {
    
    bb_parms <- rbetabinom_mean_sd(n_items = size, target_mean = mu, target_sd = sigma, draw_data = FALSE)
    alpha <- bb_parms$alpha
    beta <- bb_parms$beta 
    
    wins_sigma <- winsorized_sd_bbinom(size = size, alpha = alpha, beta = beta, trim = trim)
    
  } else if(distribution == "exponential" || distribution == "exp") {
    
    lower <- qexp(p = trim, rate  = 1 / sigma) - sigma + mu
    upper <- qexp(p = 1 - trim, rate  = 1 / sigma) - sigma + mu
    
    mean_integrand <- function(x) {
      x * dexp(x + sigma - mu, rate  = 1 / sigma) 
    }
    
    var_integrand <- function(x) {
      wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
      (x - wins_mean)^2 * dexp(x + sigma - mu, rate  = 1 / sigma) 
    }
    
    wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
    wins_variance <- trim * (lower - wins_mean)^2 + integrate(var_integrand, lower, upper)$value + trim * (upper - wins_mean)^2
    wins_sigma <- sqrt(wins_variance)
    
  } else if(distribution == "uniform" || distribution == "uni") {
    
    a <- mu - sigma * sqrt(3)
    b <- mu + sigma * sqrt(3)
    lower <- qunif(p = trim, min = a, max = b)
    upper <- qunif(p = 1 - trim, min = a, max = b)
    
    mean_integrand <- function(x) {
      x * dunif(x, min = a, max = b)
    }
    
    var_integrand <- function(x) {
      wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
      (x - wins_mean)^2 * dunif(x, min = a, max = b)
    }
    
    wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
    wins_variance <- trim * (lower - wins_mean)^2 + integrate(var_integrand, lower, upper)$value + trim * (upper - wins_mean)^2
    wins_sigma <- sqrt(wins_variance)
    
  } else if(distribution == "lognormal" || distribution == "lognorm" || distribution == "lnorm") {
    
    meanlog <- 2 * log(mu) - log(mu^2 + sigma^2) / 2
    sdlog <- sqrt(-2 * log(mu) + log(mu^2 + sigma^2))
    lower <- qlnorm(p = trim, meanlog = meanlog, sdlog = sdlog)
    upper <- qlnorm(p = 1 - trim, meanlog = meanlog, sdlog = sdlog)
    
    mean_integrand <- function(x) {
      x * dlnorm(x, meanlog = meanlog, sdlog = sdlog)
    }
    
    var_integrand <- function(x) {
      wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
      (x - wins_mean)^2 * dlnorm(x, meanlog = meanlog, sdlog = sdlog)
    }
    
    wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
    wins_variance <- trim * (lower - wins_mean)^2 + integrate(var_integrand, lower, upper)$value + trim * (upper - wins_mean)^2
    wins_sigma <- sqrt(wins_variance)
    
  } else {
    
    stop("Distribution type not supported", call. = FALSE)
    
  } # end of distribution 
  
  # winsorized mean, variance, and standardized deviation
  # if(distribution == "betabinom" || distribution == "betabin" || distribution == "bb" || distribution == "poisson" || distribution == "pois" || distribution == "negbin" || distribution == "nb") {
  #  wins_mean <- lower * trim + sum(sapply(round(lower):round(upper), mean_integrand)) + upper * trim
  #  wins_variance <- trim * (lower - wins_mean)^2 + sum(sapply(round(lower):round(upper), var_integrand)) + trim * (upper - wins_mean)^2
  #  wins_sigma <- sqrt(wins_variance)
  # } else {
  #  wins_mean <- lower * trim + integrate(mean_integrand, lower, upper)$value + upper * trim
  #  wins_variance <- trim * (lower - wins_mean)^2 + integrate(var_integrand, lower, upper)$value + trim * (upper - wins_mean)^2
  #  wins_sigma <- sqrt(wins_variance)
  #}
  
  return(wins_sigma)
  
} # winsor.pop.sd


############################
# Welch degrees of freedom #
# Validation: TRUE         #
############################
.welch_df <- function(gamma1, gamma2, f1, f2, sd1, sd2) {
  (gamma1 * sd1^2 + gamma2 * sd2^2)^2 /
    (gamma1^2 * sd1^4 / f1 + gamma2^2 * sd2^4 / f2)
}



# Lecoutre, B. (2007). Another look at the confidence intervals for the noncentral T distribution. 
# Journal of Modern Applied Statistical Methods, 6, 107-116.
# for 95% CI use alpha = 0.025
# Confidence limits are Lambda-prime quantiles; alpha is a tail probability.
.confint_nct <- function(ncp, df, alpha) {
  qlambdap(p = c(alpha, 1 - alpha), df = df, ncp = ncp)
}


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
# Primary reference: Bulus (2026), Research Synthesis Methods;
# Table 1 and Appendix A1 (mean-difference HC variances), Method (TW methods).
#' @rdname hct.test
#' @method hct.test default
#' @export
hct.test.default <- function(x = NULL, y = NULL, 
                             m1, m2, # define length of two for each for mixture 
                             sd1, sd2, # define length of two for each for mixture 
                             n1, n2, 
                             alpha = 0.05, two.tailed = TRUE, 
                             ci = c("ct", "nct", "z", "FALSE"), 
                             method = c("constant", "student", "sd1", "sd2",
                                        "hc0", "hc1", "hc2", "welch",
                                        "hc3", "hc4", "hc4m", "hc5"),
                             df = c("regression", "welch-satterthwaite"),
                             # distribution is needed when trim > 0 and m1, m2, sd1, and sd2 are provided
                             trim = 0, na.rm = FALSE, weights = NULL,
                             distribution = c("norm", "norm"), # or "mixnorm", "pois", "exp", "uni", "lnorm" # to get pop diff and sd
                             mixnorm.weights = 0.50, # mixture proportions
                             betabin.size = 20, # beta-binomial size
                             norm.transform.g = NULL, 
                             norm.transform.h = NULL,
                             verbose = TRUE, ...) { 
  
  analysis <- .record_analysis("hct.test.default", match.call(), environment())

  # Local names used by internal distribution helpers.
  weight <- mixnorm.weights
  size <- betabin.size
  g <- norm.transform.g
  h <- norm.transform.h

  method <- match.arg(method)
  # distribution <- match.arg(distribution)
  df <- tolower(match.arg(df))
  ci <- as.character(ci)
  ci <- tolower(match.arg(ci))

  if (!is.numeric(trim) || length(trim) != 1L ||
      !is.finite(trim) || trim < 0 || trim > 0.45) {
    stop("'trim' must be a finite number in [0, 0.45].", call. = FALSE)
  }
  
  weighted.groups <- NULL
  if (!is.null(weights)) {
    if (is.null(x) || is.null(y)) {
      stop("Observation weights require raw data in both groups.", call. = FALSE)
    }
    if (trim > 0) {
      stop("Analytical weighted trimmed inference is pending review; use trim = 0 or boot.hct.test() for weighted trimming.",
           call. = FALSE)
    }
    if (df != "regression" || method %in% c("welch", "sd1", "sd2")) {
      stop("This weighting draft supports Student and HC0-HC5/HC4m with regression degrees of freedom; other weighted inference rules are pending review.",
           call. = FALSE)
    }
    weighted.groups <- .prepare_observation_weights(x, y, weights, na.rm)
    x <- weighted.groups$x$values
    y <- weighted.groups$y$values
  }
  if(na.rm && is.null(weighted.groups)) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
  }
  
  if(!is.null(x) & !is.null(y)) {
    
    n1 <- length(x)
    n2 <- length(y)
    .check_retained_sample_sizes(n1, n2, trim)
    
    if (!is.null(weighted.groups)) {
      m1 <- weighted.groups$x$mean
      m2 <- weighted.groups$y$mean
      sd1 <- weighted.groups$x$sd
      sd2 <- weighted.groups$y$sd
    } else if(trim == 0) {
      m1 <- mean(x)
      m2 <- mean(y)
      sd1 <- sd(x)
      sd2 <- sd(y)
    } else {
      if(trim > 0.45) stop("Trim Proportion > 0.45 is not allowed", call. = FALSE)
      m1 <- mean(trim(x, trim = trim, na.rm = na.rm))
      m2 <- mean(trim(y, trim = trim, na.rm = na.rm))
      h1 <- n1 - 2*round(trim * n1)
      h2 <- n2 - 2*round(trim * n2)
      # should we multiply sd1 and sd2 with sqrt((n1 - 1) / (h1 - 1)) 
      # to get Yuen (1973, 1974) definition of Winsorized standard deviation?! YES
      sd1 <- sd(winsor(x, trim = trim, na.rm = na.rm)) * sqrt((n1 - 1) / (h1 - 1))
      sd2 <- sd(winsor(y, trim = trim, na.rm = na.rm)) * sqrt((n2 - 1) / (h2 - 1))
      n1 <- n1 - 2*round(trim * n1)
      n2 <- n2 - 2*round(trim * n2)
      
    }
    
  } else {
    .check_retained_sample_sizes(n1, n2, trim)
    user.parms.names <- names(as.list(match.call()))
    use.dist <- "distribution" %in% user.parms.names & !is.null(distribution)
    use.gh <- !is.null(g) & !is.null(h) 
    
    if(length(distribution) == 1) distribution <- c(distribution, distribution)

    if(trim > 0) {
      
      if(trim > 0.45) stop("Trim Proportion > 0.45 is not allowed", call. = FALSE)
      
      if(use.dist || !use.gh) {
        
        if(!all(distribution %in% c("norm", "mixnorm", "betabin", "negbin", "pois", "exp", "uni", "lnorm" ))) 
          stop("The distribution argument requires a character string (or a pair of strings) identifying the distribution(s) to use.\n Please select from the following valid options: 'norm' (Normal), 'pois' (Poisson), 'exp' (Exponential), 'uni' (Uniform), or 'lnorm' (Lognormal).", call. = FALSE)
        
        if(distribution[1] != "mixnorm") {
          if(length(m1) != 1) stop("'m1' should be a scalar value", call. = FALSE)
          if(length(sd1) != 1) stop("'sd1' should be a scalar value", call. = FALSE)
        }
        
        if(distribution[2] != "mixnorm") {
          if(length(m2) != 1) stop("'m2' should be a scalar value", call. = FALSE)
          if(length(sd1) != 1) stop("'sd1' should be a scalar value", call. = FALSE)
        }
        
        m1.trimmed <- trim.pop.mean(mu = m1, sigma = sd1, trim = trim, size = size, weight = weight, distribution = distribution[1])
        m2.trimmed <- trim.pop.mean(mu = m2, sigma = sd2, trim = trim, size = size, weight = weight, distribution = distribution[2])
        # should we multiply sd1 and sd2 with sqrt(1 / (1 - 2*trim))  
        # to get Yuen (1973, 1974) definition of Winsorized standard deviation?! 
        sd1.winsored <- winsor.pop.sd(mu = m1, sigma = sd1, trim = trim, size = size, weight = weight, distribution = distribution[1]) * sqrt(1 / (1 - 2*trim)) 
        sd2.winsored <- winsor.pop.sd(mu = m2, sigma = sd2, trim = trim, size = size, weight = weight, distribution = distribution[2]) * sqrt(1 / (1 - 2*trim)) 
        n1 <- n1 - 2*round(trim * n1)
        n2 <- n2 - 2*round(trim * n2)
        
        m1 <- m1.trimmed
        m2 <- m2.trimmed
        
        sd1 <- sd1.winsored
        sd2 <- sd2.winsored
        
      } else if(!use.dist & use.gh) {
        
        if(length(g) == 1) g <- c(g, g)
        if(length(h) == 1) h <- c(h, h)
       
        x1 <- rnorm(2e6, mean = m1, sd = sd1)
        x2 <- rnorm(2e6, mean = m2, sd = sd2)
        v1 <- norm.gh(x1, g = g[1], h = h[1])
        v2 <- norm.gh(x2, g = g[2], h = h[2])
        m1 <- mean(trim(v1, trim = trim))
        m2 <- mean(trim(v2, trim = trim))
        sd1 <- sd(winsor(v1, trim = trim)) * sqrt(1 / (1 - 2*trim)) 
        sd2 <- sd(winsor(v2, trim = trim)) * sqrt(1 / (1 - 2*trim)) 
        n1 <- n1 - 2*round(trim * n1)
        n2 <- n2 - 2*round(trim * n2)
        
      } else {
        
        stop("Invalid argument combination. You must specify only ONE of the following two options:\n",
             "  1. The 'distribution' argument (e.g., 'normal'), OR\n",
             "  2. BOTH 'norm.transform.g' and 'norm.transform.h' arguments.")
        
      } # use.method
        
    } # trim 
    
    if(trim == 0 & distribution[1] == "mixnorm") {
      mixture <- .mixture_mean_sd(mu = m1, sigma = sd1, weight = weight)
      m1 <- mixture$mean
      sd1 <- mixture$sd
    } # trim = 0 for mixture normal
    
    if(trim == 0 & distribution[2] == "mixnorm") {
      mixture <- .mixture_mean_sd(mu = m2, sigma = sd2, weight = weight)
      m2 <- mixture$mean
      sd2 <- mixture$sd
    } # trim = 0 for mixture normal
    
  } # data
  
  # global parameters
  delta <- (m1 - m2) # raw mean difference
  ifelse(two.tailed,
         alpha.for.test <- alpha / 2,
         alpha.for.test <- alpha)
  alpha.for.power <- alpha
  

  if(!is.na(sd1) & sd1 == 0) {
    sd1 <- 0.0000001 # to avoid division by zero for method sd1
    warning("'sd1' is 0, using 0.0000001", call. = FALSE)
  }
  if(!is.na(sd2) & sd2 == 0) {
    sd2 <- 0.0000001 # to avoid division by zero for method sd2
    warning("'sd1' is 0, using 0.0000001", call. = FALSE)
  }
  
  # sampling variance 
  if (!is.null(weighted.groups)) {
    weighted.variances <- .weighted_hc_variance(weighted.groups, method)
    var.delta <- sum(weighted.variances)
  } else if(method == "constant" | method == "student") {
    
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
    if (!is.null(weighted.groups)) {
      def <- sum(weighted.variances)^2 /
        sum(weighted.variances^2 /
              (c(weighted.groups$x$effective_n,
                 weighted.groups$y$effective_n) - 1))
    } else {
      def <- .welch_df(gamma1 = gamma1, gamma2 = gamma2,
                     f1 = n1 - 1, f2 = n2 - 1, 
                     sd1 = sd1, sd2 = sd2)
    }
    if(method == "sd1" | method == "sd2") {
      stop("Welch-Satterthwaite degrees of freedom is not applicable to 'sd1' or 'sd2' testing methods", call. = FALSE)
    }
  } else {
    if(method == "welch") {
      stop("Welch-Satterthwaite test requires Welch-Satterthwaite degrees of freedom", call. = FALSE)
    }
    if (method == "sd1") {
      def <- if (is.null(weighted.groups)) n1 - 1 else weighted.groups$x$effective_n - 1
    } else if(method == "sd2") {
      def <- if (is.null(weighted.groups)) n2 - 1 else weighted.groups$y$effective_n - 1
    } else {
      def <- n1 + n2 - p
    }
  }
  
  # confidence intervals 
  if(ci == "nct") {
    
    V <- delta / sqrt(var.delta)
    
    confint.t <- .confint_nct(ncp = V, df = def, alpha.for.test)
    confint.delta <- confint.t * sqrt(var.delta)
    
    # if(V <= 37) { # numerical issues with non-central t-distribution starts roughly around 37, check MBESS package 
    #   confint.t <- .confint_nct(ncp = V, df = def, alpha.for.test)
    #   confint.delta <- confint.t * sqrt(var.delta)
    # } else {
    #   warning("Switching to central t-distribution for CI because NCP is large", call. = FALSE)
    #   ct <- qt(p = alpha.for.test, ncp = 0, df = def, log.p = FALSE, lower.tail = FALSE)
    #   confint.delta <- c(delta - ct * sqrt(var.delta), delta + ct * sqrt(var.delta))
    # }
   
  } else if(ci == "ct") {
    
    ct <- qt(p = alpha.for.test, ncp = 0, df = def, log.p = FALSE, lower.tail = FALSE)
    confint.delta <- c(delta - ct * sqrt(var.delta), delta + ct * sqrt(var.delta))
    
  } else if(ci == "z"){
    
    cz <- qnorm(p = alpha.for.test, mean = 0, sd = 1, lower.tail = FALSE)
    confint.delta <- c(delta - cz * sqrt(var.delta), delta + cz * sqrt(var.delta))
    
  } else {
    
    confint.delta <- c(NA, NA)
    
  }
  
  # get p values
  t.observed <- delta / sqrt(var.delta)
  if(ci == "ct") {
    
    p.value <- pt(q = abs(t.observed), ncp = 0, df = def, lower.tail = FALSE)
    if(two.tailed) p.value <- 2 * p.value
  } else if(ci == "z") {
    p.value <- pnorm(abs(t.observed), lower.tail = FALSE)
    if(two.tailed) p.value <- 2 * p.value
  } else {
    # if(two.tailed) p.value <- 2 * p.value
    
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
                      trim = trim, 
                      alpha = alpha, df = df,
                      two.tailed = two.tailed, 
                      ci = ci, method = method)

  if (!is.null(weighted.groups)) {
    parms$weighted <- TRUE
    parms$n_eff1 <- weighted.groups$x$effective_n
    parms$n_eff2 <- weighted.groups$y$effective_n
  }
  
  if(verbose) {
    
    print(
      
      round(test, 4)
      
    )
    
    if(is.na(p.value)) message("P-value is not calculated, use confidence interval instead")
    
    if(any(is.na(confint.delta))) message("Confidence interval is not calculated")
    
  }
  
  invisible(
    
    structure(
      list(test = test, parms = parms),
      class = "hct.test", analysis = analysis)
    
  )
  
} # hct.test.default()


#' @rdname hct.test
#' @method hct.test formula
#' @export
hct.test.formula <- function(formula, data, weights = NULL, ...) {
  
  formula <- as.formula(formula)
  vars <- all.vars(formula)
  y <- vars[1]
  x <- vars[2]

  observation.weights <- eval(substitute(weights), envir = data,
                              enclos = parent.frame())
  if (!is.null(observation.weights)) {
    if (!is.numeric(observation.weights) ||
        length(observation.weights) != nrow(data)) {
      stop("Formula weights must be a numeric vector with one value per data row.",
           call. = FALSE)
    }
    dots <- list(...)
    missing.group <- is.na(data[[x]])
    if (any(missing.group)) {
      if (!isTRUE(dots$na.rm)) {
        stop("Missing group labels; use na.rm = TRUE.", call. = FALSE)
      }
      data <- data[!missing.group, , drop = FALSE]
      observation.weights <- observation.weights[!missing.group]
    }
    data[[x]] <- droplevels(as.factor(data[[x]]))
    if (nlevels(data[[x]]) != 2L) {
      stop("Weighted analyses require exactly two observed groups.", call. = FALSE)
    }
  }
  
  x_levels <- unique(data[[x]])
  if (length(x_levels) > 2) {
    stop(paste("The predictor variable '", x, "' has more than two levels. It should be binary."))
  }
  
  # ensure x is a factor
  if (!is.factor(data[[x]])) {
    data[[x]] <- as.factor(data[[x]])
  }
  
  # split y based on levels of x
  y1 <- data[[y]][data[[x]] == levels(data[[x]])[1]]
  y2 <- data[[y]][data[[x]] == levels(data[[x]])[2]]
  
  # n1 <- length(y1)
  # n2 <- length(y2)
  # m1 <- mean(y1)
  # m2 <- mean(y2)
  # sd1 <- sd(y1)
  # sd2 <- sd(y2)
  
  # if(trim == 0) {
  #   m1 <- mean(y1)
  #   m2 <- mean(y2)
  #   sd1 <- sd(y1)
  #   sd2 <- sd(y2)
  # } else {
  #   m1 <- mean(trim(x, trim = trim))
  #   m2 <- mean(trim(y, trim = trim))
  #   h1 <- n1 - 2*round(trim * n1)
  #   h2 <- n2 - 2*round(trim * n2)
  #   # should we multiply sd1 and sd2 with sqrt((n1 - 1) / (h1 - 1)) 
  #   # to get Yuen (1973, 1974) definition of Winsorized standard deviation?! YES
  #   sd1 <- sd(winsor(x, trim = trim)) * sqrt((n1 - 1) / (h1 - 1))
  #   sd2 <- sd(winsor(y, trim = trim)) * sqrt((n2 - 1) / (h2 - 1))
  #   n1 <- n1 - 2*round(trim * n1)
  #   n2 <- n2 - 2*round(trim * n2)
  # }
  
  parms <- list(x = y1, y = y2, 
                # m1 = m1, m2 = m2,
                # sd1 = sd1, sd2 = sd2,
                # n1 = n1, n2 = n2,
                ...)
  if (!is.null(observation.weights)) {
    parms$weights <- list(
      x = observation.weights[data[[x]] == levels(data[[x]])[1]],
      y = observation.weights[data[[x]] == levels(data[[x]])[2]])
  }
  
  do.call("hct.test", parms)
  
}


#' Two-group heteroscedasticity-consistent tests
#'
#' @description Two-group heteroscedasticity-consistent tests.
#' @param x,y Numeric observations for the first and second groups. For the formula interface, the first factor level defines the first group. The difference is first minus second.
#' @param m1,m2 Group means for summary-data analysis; vectors of component means for mixture summaries.
#' @param sd1,sd2 Group SDs, or component SDs for mixture summaries.
#' @param n1,n2 Group sample sizes for summary-data analysis.
#' @param mixnorm.weights Mixture-component proportions for summary-data distributions; distinct from observation weights.
#' @param alpha Significance level.
#' @param verbose Logical; print results.
#' @param two.tailed Logical; use a two-sided test. One-sided tests use the magnitude of the observed difference.
#' @param ci Confidence-interval approximation: central t (ct), noncentral t inversion (nct), normal (z), or FALSE to omit intervals.
#' @param method Student/constant variance, a named HC estimator, Welch, or the sd1/sd2 variance convention.
#' @param df Regression or Welch-Satterthwaite degrees of freedom.
#' @param betabin.size Size parameter for beta-binomial summary distributions.
#' @param trim Fraction trimmed from each tail, between 0 and 0.45. Raw unweighted data trim round(n * trim) observations per tail.
#' @param na.rm Logical; remove incomplete observations. Weighted observations and their weights are removed together.
#' @param weights Observation weights for raw data only. Supply a concatenated vector for x then y, or list(x = ..., y = ...). The formula method accepts a data column expression or a numeric vector aligned with data.
#' @param distribution Two group distributions used for trimmed summary-data calculations.
#' @param norm.transform.g,norm.transform.h Optional g-and-h parameters for summary-data moment calculations.
#' @param formula A formula of the form outcome ~ group.
#' @param data Data frame containing outcome and group; exactly two observed groups are required.
#' @param ... Additional arguments passed from the formula method to the default method.
#'
#' @details Observation weights must be finite and nonnegative; zero-weight rows are excluded. At least two positive-weight observations per group are required. Weighted analytical inference currently supports trim = 0, regression df, and Student/constant or HC0, HC1, HC2, HC3, HC4, HC4m, HC5. Weighted trimmed analytical tests, weighted Welch-Satterthwaite df, and weighted sd1/sd2 tests raise errors. Use \code{boot.hct.test} for supported weighted trimmed bootstrap estimates. Weighted Student inference treats weights as precision weights; weighted HC inference uses WLS residual scores and leverages. The validation vignette compares untrimmed mean-difference inference with lm and sandwich; it does not externally validate weighted effect-size SEs. For ci = "z", p uses the standard-normal upper tail of the absolute estimate/SE ratio, doubled when two.tailed = TRUE. The one-sided convention follows the observed direction, as for the existing central-t test.
#'
#' @return Invisibly, an object of class hct.test containing a test data frame and analysis parameters. The analysis attribute retains evaluated group data or summaries and arguments for replay at a different alpha.
#'
#' @examples set.seed(1)
#' a <- rnorm(30)
#' b <- rnorm(45, 0.5, 2)
#' fit <- hct.test(a, b, method = "hc3", verbose = FALSE)
#' fit$test
#' dat <- data.frame(y = c(a, b), group = factor(rep(c("A", "B"), c(30, 45))),
#'                   w = runif(75, 0.5, 2))
#' hct.test(y ~ group, dat, weights = w, method = "hc3", verbose = FALSE)
#'
#' @references Bulus, M. (2026). Robust standard errors and confidence intervals for standardized mean difference. \emph{Research Synthesis Methods}. Advance online publication.
#'
#' See Table 1 and Appendix A1 for mean-difference variance estimators, and the Method section for trimmed/Winsorized extensions. Observation-weight support is a later implementation extension; this reference does not establish validation of that extension.
#'
#' @rdname hct.test
#' @export hct.test
hct.test <- function(x, ...) {
  UseMethod("hct.test")
}

# Primary reference: Bulus (2026), Research Synthesis Methods;
# Method, Table 1 and Appendix A2 (SMD and standardizer sampling variances).
#' Standardized mean differences from an HC test
#'
#' @description Standardized mean differences from an HC test.
#' @param x An hct.test object.
#' @param ss.adjust Logical; apply a small-sample Hedges correction.
#' @param ci Interval approximation: nct, ct, z, or FALSE.
#' @param scale Logical; apply normal-population calibration to trimmed effect sizes.
#' @param verbose Logical; print results.
#'
#' @details Uses group moments and the mean-difference variance stored in x. It does not consume raw data and therefore has no weights argument. Weighted descriptive SDs use denominator sum(w)-sum(w^2)/sum(w); Kish effective sample sizes enter descriptive pooling and effect-size corrections. Weighted effect-size SE corrections are approximations and are not validated by the lm/sandwich mean-difference benchmark. The fitted testing method determines the standardizer. When ci = "z", the es data frame includes p from the standard-normal upper tail of abs(d/se), doubled for two-sided inference; p is NA for other interval choices. This SMD test uses the SMD SE and can differ from the raw-difference test.
#'
#' @return Invisibly, an object of class robust.smd with an es data frame and analysis parameters. The analysis attribute retains the original fitted test and evaluated SMD arguments for replay at other confidence levels.
#'
#' @examples set.seed(1)
#' fit <- hct.test(rnorm(30), rnorm(45, 0.5), method = "hc2", verbose = FALSE)
#' robust.smd(fit, ci = FALSE, verbose = FALSE)$es
#'
#' @references Bulus, M. (2026). Robust standard errors and confidence intervals for standardized mean difference. \emph{Research Synthesis Methods}. Advance online publication.
#'
#' See the Method section, Table 1, and Appendix A2 for SMD sampling variance and standardizer variance. The paper supplies the primary theoretical reference; implementation agreement and weighted extensions require separate checks.
#'
#' @rdname robust.smd
#' @export robust.smd
robust.smd <- function(x = NULL,
                         ss.adjust = FALSE, 
                         ci = c("nct", "ct", "z", "FALSE"),
                         scale = TRUE, # only works when trim != 0
                         verbose = TRUE) {
  
  analysis <- .record_analysis("robust.smd", match.call(), environment())

  if(!inherits(x, "hct.test")) 
    stop("Please provide an object retruned from 'hct.test()' function",
         call. = FALSE)
  
  ci <- as.character(ci)
  ci <- match.arg(ci)
  
  n1 <- x$parms$n1
  n2 <- x$parms$n2
  # The input object already contains weighted moments and HC variance.
  # Effective sample sizes supply the moment-based standardizer corrections.
  if (isTRUE(x$parms$weighted)) {
    n1 <- x$parms$n_eff1
    n2 <- x$parms$n_eff2
  }
  m1 <- x$parms$m1
  m2 <- x$parms$m2
  sd1 <- x$parms$sd1
  sd2 <- x$parms$sd2
  trim <- x$parms$trim
  alpha <- x$parms$alpha
  two.tailed <- x$parms$two.tailed
  method <- x$parms$method
  
  df <- x$parms$df <- x$test$df
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
    
    # if(V <= 37) {
    #   confint.t <- .confint_nct(ncp = V, df = df, alpha.for.test)
    #   confint.d <- confint.t * (se.delta / std) * ss.adj 
    # } else {
    #   warning("Switching to normal approximation for CI because NCP is large", call. = FALSE)
    #   cz <- qnorm(p = alpha.for.test, mean = 0, sd = 1, lower.tail = FALSE)
    #   confint.d <- c(d - cz * se.d, d + cz * se.d)
    # }
    
  } else if(tolower(ci) == "ct") {
    
    ct <- qt(p = alpha.for.test, ncp = 0, df = df, log.p = FALSE, lower.tail = FALSE)
    confint.d <- c(d - ct * se.d, d + ct * se.d)
    
  } else if(tolower(ci) == "z") {
    
    cz <- qnorm(p = alpha.for.test, mean = 0, sd = 1, lower.tail = FALSE)
    confint.d <- c(d - cz * se.d, d + cz * se.d)
    
  } else {
    
    confint.d <- c(NA, NA)
    
  } 
  
  # normal-theory calibration of robust effect size 
  # does not change the test statistics 
  if(scale) {
    
    if(trim == 0) {
      sd.scaling.factor <- 1
    } else {
      # sqrt(1 / (1 - 2*trim))  factor preserves the definition of Yuen (1974) SD winsorization, which is slightly difference from Algina et al. (2005)
      sd.scaling.factor <- winsor.pop.sd(mu = 0, sigma = 1, trim = trim, distribution = "norm") * sqrt(1 / (1 - 2*trim)) 
    }
    
    es <- data.frame(d = d * sd.scaling.factor,
                     se = se.d * sd.scaling.factor,
                     lcl = confint.d[[1]] * sd.scaling.factor,
                     ucl = confint.d[[2]] * sd.scaling.factor)
      
  } else {
    
    es <- data.frame(d = d,
                     se = se.d,
                     lcl = confint.d[[1]],
                     ucl = confint.d[[2]])
    
  } # scale 
  
  
  # The normal test uses the same estimate and SE as the reported SMD interval.
  es$p <- if (ci == "z") {
    p <- pnorm(abs(es$d / es$se), lower.tail = FALSE)
    if (two.tailed) 2 * p else p
  } else NA_real_

  parms <- x$parms
  parms$ci <- ci
  parms$ss.adjust <- ss.adjust
  parms$scale <- scale

  if(verbose) {
    print(
      round(es, 4)
    )
  }
  
  invisible(
    
    structure(
      list(es = es, parms = parms),
      class = "robust.smd", analysis = analysis)
    
  )
  
} # robust.smd()





# Simulation context: Bulus (2026), Research Synthesis Methods, Method.
#' Generate two-group data
#'
#' @description Generate two-group data.
#' @param mu1,mu2 Target means; component means for mixture normals.
#' @param sd1,sd2 Target SDs; component SDs for mixture normals.
#' @param n1,n2 Group sample sizes.
#' @param mixnorm.weights Mixture-normal component proportions.
#' @param betabin.size Number of items for beta-binomial generation.
#' @param format Return two outcome vectors, a data frame, or outcomes and group indicators.
#' @param distribution Generating distribution; norm or its alias normal is required when norm.transform.g or norm.transform.h is supplied.
#' @param norm.transform.g,norm.transform.h Optional g-and-h shape parameters: a scalar for both groups or a pair of group values. Omitted parameters are zero. norm.transform.h must be nonnegative.
#'
#' @details Distribution-specific parameter constraints apply. Exponential and Poisson draws are shifted to obtain the requested mean; uniform endpoints are obtained from the mean and SD. Normal mixtures use component parameters rather than overall moments. When norm.transform.g or norm.transform.h is supplied, normal draws are standardized by their sample moments, transformed as expm1(g*z)/g times exp(h*z^2/2) (using z when g is zero), and rescaled to preserve the original normal draws' realized sample mean and SD. This sample calibration creates dependent observations within each transformed sample; it is a simulation convention rather than independent draws with analytically calibrated population moments. Zero g and h reproduce the original normal draws exactly. At least two observations per group are required; extreme parameters that produce nonfinite scales raise an error. This function has no observation-weight argument.
#'
#' @return A list of group vectors or a data frame, according to format.
#'
#' @examples set.seed(2)
#' gendat.hct.test(0, 0.5, 1, 2, 30, 40, distribution = "norm",
#'                 norm.transform.g = c(0, 0.3), norm.transform.h = 0.1)
#'
#' set.seed(1)
#' gendat.hct.test(mu1 = 0, mu2 = 0.5, sd1 = 1, sd2 = 2, n1 = 20, n2 = 30)
#'
#' @references Bulus, M. (2026). Robust standard errors and confidence intervals for standardized mean difference. \emph{Research Synthesis Methods}. Advance online publication.
#'
#' See the Monte Carlo design in the Method section for the simulation context of data generation.
#'
#' @rdname gendat.hct.test
#' @export gendat.hct.test
gendat.hct.test <- function(mu1, mu2, sd1, sd2, n1, n2, 
                            mixnorm.weights = 0.50, # for mixture normal
                            group.levels = c("Group 1", "Group 2"),
                            format = c("list.y", "data", "list.yx"),
                            distribution = c("norm", "mixnorm", "lnorm", "uni", "exp", "pois", "negbin", "betabin"),
                            betabin.size = ceiling(max(mu1, mu2) + 3 * max(sd1, sd2)),
                            norm.transform.g = NULL, norm.transform.h = NULL) {
  
  size <- betabin.size
  g <- norm.transform.g
  h <- norm.transform.h
  
  if (length(group.levels) != 2L ||
      anyNA(group.levels) ||
      anyDuplicated(as.character(group.levels))) {
    stop("'group.levels' must contain two distinct, nonmissing labels.")
  }
  group.levels <- as.character(group.levels)

  if (identical(distribution, "normal")) distribution <- "norm"
  distribution <- tolower(match.arg(distribution))
  fmt <- tolower(match.arg(format))
  use.gh <- !is.null(g) || !is.null(h)
  if (use.gh) {
    if (distribution != "norm") {
      stop("'norm.transform.g' and 'norm.transform.h' are available only with distribution = 'norm' (or 'normal').",
           call. = FALSE)
    }
    shape <- function(value, name) {
      if (is.null(value)) return(c(0, 0))
      if (!is.numeric(value) || !is.null(dim(value)) ||
          !length(value) %in% c(1L, 2L) || any(!is.finite(value))) {
        stop(sprintf("'%s' must be a finite number or a pair of group values.", name),
             call. = FALSE)
      }
      rep(value, length.out = 2L)
    }
    g <- shape(g, "norm.transform.g")
    h <- shape(h, "norm.transform.h")
    if (any(h < 0)) stop("'norm.transform.h' must be nonnegative.", call. = FALSE)
    if (any(!is.finite(c(n1, n2))) || any(c(n1, n2) < 2) ||
        any(c(n1, n2) != floor(c(n1, n2)))) {
      stop("g-and-h generation requires integer group sample sizes of at least two.",
           call. = FALSE)
    }
  }
  # generate variables depending on the distribution type
  # given means and standard deviations
  if(distribution == "normal" || distribution == "norm") {
    
    y1 <- rnorm(n = n1, mean = mu1, sd = sd1)
    y2 <- rnorm(n = n2, mean = mu2, sd = sd2) 
    
  } else if(distribution == "mixture-normal" || distribution == "mixnorm" || distribution == "mnorm") {
    
    y1 <- rmixnormal(n = n1, mu = mu1, sigma = sd1, weight = mixnorm.weights)
    y2 <- rmixnormal(n = n2, mu = mu2, sigma = sd2, weight = mixnorm.weights)
    
  } else if(distribution == "negbin" || distribution == "nb") {
    
    y1 <- rnbinom(n = n1, mu = mu1, size = mu1^2 / (sd1^2 - mu1)) # prob = size/(size+mu) by default
    y2 <- rnbinom(n = n2, mu = mu2, size = mu2^2 / (sd2^2 - mu2)) # prob = size/(size+mu) by default
    
  } else if(distribution == "poisson" || distribution == "pois") {
    
    y1 <- rpois(n = n1, lambda = sd1^2) - (sd1^2 - mu1)
    y2 <- rpois(n = n2, lambda = sd2^2) - (sd2^2 - mu2)
    
  } else if(distribution == "betabinom" || distribution == "betabin" || distribution == "bb") {
    
    y1 <- rbetabinom_mean_sd(n_students = n1, n_items = size, target_mean = mu1, target_sd = sd1, draw_data = TRUE)$x
    y2 <- rbetabinom_mean_sd(n_students = n2, n_items = size, target_mean = mu2, target_sd = sd2, draw_data = TRUE)$x
    
  } else if(distribution == "exponential" || distribution == "exp") {
    
    y1 <- rexp(n = n1, rate  = 1 / sd1) - sd1 + mu1
    y2 <- rexp(n = n2, rate  = 1 / sd2) - sd2 + mu2
    
  } else if(distribution == "uniform" || distribution == "uni") {
    
    a1 <- mu1 - sd1 * sqrt(3)
    b1 <- mu1 + sd1 * sqrt(3)
    y1 <- runif(n = n1, min = a1, max = b1)
    
    a2 <- mu2 - sd2 * sqrt(3)
    b2 <- mu2 + sd2 * sqrt(3)
    y2 <- runif(n = n2, min = a2, max = b2)
    
  } else if(distribution == "lognormal" || distribution == "lognorm" || distribution == "lnorm") {
    
    meanlog1 <- 2 * log(mu1) - log(mu1^2 + sd1^2) / 2
    sdlog1 <- sqrt(-2 * log(mu1) + log(mu1^2 + sd1^2))
    y1 <- rlnorm(n = n1, meanlog = meanlog1, sdlog = sdlog1)
    
    meanlog2 <- 2 * log(mu2) - log(mu2^2 + sd2^2) / 2
    sdlog2 <- sqrt(-2 * log(mu2) + log(mu2^2 + sd2^2))
    y2 <- rlnorm(n = n2, meanlog = meanlog2, sdlog = sdlog2)
    
  } else {
    
    stop("Distribution type not supported", call. = FALSE)
    
  } # end of distribution

  if (use.gh) {
    # Reuse the existing transformation and its sample-moment rescaling.
    y1 <- norm.gh(y1, g = g[1], h = h[1], rescale = TRUE)
    y2 <- norm.gh(y2, g = g[2], h = h[2], rescale = TRUE)
  }
  if (fmt == "list.y") {
    
    return(list(y1 = y1, y2 = y2))
    
  } else {
    
    if (length(group.levels) != 2L ||
        anyNA(group.levels) ||
        anyDuplicated(as.character(group.levels))) {
      stop("'group.levels' must contain two distinct, nonmissing labels.")
    }
    
    group.levels <- as.character(group.levels)
    
    x1 <- factor(rep(group.levels[1], n1), levels = group.levels)
    x2 <- factor(rep(group.levels[2], n2), levels = group.levels)
    
    if (fmt == "list.yx") {
      
      return(list(y1 = y1, y2 = y2, x1 = x1, x2 = x2))
      
    } else {
      
      return(data.frame(
        y = c(y1, y2),
        x = factor(
          c(rep(group.levels[1], n1), rep(group.levels[2], n2)),
          levels = group.levels
        )
      ))
    }
    
  }
  
} # gendat.hct.test


# Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages, 
# American Statistician 50, 361-365. doi: 10.2307/2684934.
# combine raw and standardized boots
# Primary reference: Bulus (2026), Research Synthesis Methods;
# Method and Results (bootstrap comparisons). Observation-weight paths are extensions.
#' Bootstrap raw and standardized group differences
#'
#' @description Bootstrap raw and standardized group differences.
#' @param formula A formula of the form outcome ~ group.
#' @param data Data frame with exactly two observed groups.
#' @param alpha Significance level.
#' @param verbose Logical; print results.
#' @param trim Fraction of each tail to remove, between 0 and 0.45.
#' @param scale Logical; apply normal calibration to trimmed standardized differences.
#' @param weights Observation weights as a data column expression or vector aligned with rows; finite, nonnegative values are required.
#' @param standardizer Pooled group SD (sdp), first group SD (sd1), or second group SD (sd2).
#' @param distribution Distribution used for parametric generation; weighted parametric generation supports norm only.
#' @param boot.type Ordinary stratified row bootstrap or parametric bootstrap.
#' @param parallel Parallel mode passed to boot::boot.
#' @param ncpus Number of worker processes.
#' @param ci.quantile.type Legacy argument retained in the interface; currently unused by boot::boot.ci.
#' @param boot.rep Number of bootstrap replicates.
#' @param na.rm Logical; jointly remove missing outcomes, groups, and weights.
#' @param ... Additional hct.test arguments for untrimmed statistics or parametric generating fits.
#'
#' @details Ordinary resampling draws rows uniformly within each group and carries weights with their rows. Observation weights are not importance-resampling probabilities. Weighted normal parametric generation holds weights fixed and uses untrimmed weighted means with row SD sqrt(sum(w * residual^2)/(n_group-1))/sqrt(w_i).
#'
#' With weights and positive trim, trimming removes the requested weight mass from each tail, retaining partial boundary weights and recomputing cutoffs in each replicate. The weighted Winsorized SD uses denominator sum(w)-sum(w^2)/sum(w), then divides by sqrt(1-2*trim). Pooled standardization uses original Kish effective sample sizes. scale = TRUE additionally multiplies the standardized difference by the corresponding adjusted standard-normal Winsorized SD. These are descriptive scale conventions; inference uses bootstrap variability. Count trimming without weights may differ from weight-mass trimming with unit weights at fractional boundaries. A zero or undefined standardizer raises an error, including in a replicate. Analytical weighted trimmed HC inference remains unavailable. Extra hct.test arguments do not change the dedicated weighted trimmed statistic.
#'
#' @return Invisibly, a two-row data frame of class boot.hct.test, with a plot method, containing observed estimates, bootstrap means and SEs, percentile intervals, and BCa intervals when available. Parametric BCa fields are NA. The analysis attribute retains the original data, evaluated arguments and starting RNG state. Plotting at other confidence levels reruns boot.hct.test with the changed alpha.
#'
#' @examples \dontrun{
#' set.seed(1)
#' dat <- data.frame(y = rnorm(100), group = factor(rep(c("A", "B"), each = 50)),
#'                   w = runif(100, 0.5, 2))
#' boot.hct.test(y ~ group, dat, weights = w, trim = 0.1, boot.rep = 999,
#'               ncpus = 1, parallel = "no", verbose = FALSE)
#' }
#'
#' @references Bulus, M. (2026). Robust standard errors and confidence intervals for standardized mean difference. \emph{Research Synthesis Methods}. Advance online publication.
#'
#' See the bootstrap comparisons in the Method and Results sections. The observation-weight bootstrap paths and partial-weight trimming in this implementation are later extensions; their regression checks do not reproduce the paper's coverage study.
#'
#' @rdname boot.hct.test
#' @export boot.hct.test
boot.hct.test <- function(formula, data, alpha = 0.05, 
                          trim = 0, scale = TRUE, weights = NULL,
                          standardizer = c("sdp", "sd1", "sd2"),
                          distribution = c("norm", "mixnorm", "lnorm", "uni", "exp", "pois", "negbin", "betabin"), # for parametric bootstrap
                          boot.type = c("ordinary", "parametric"), # 'permutation' (disabled)
                          parallel = c("no", "multicore", "snow"),
                          ncpus = parallel::detectCores() - 1,
                          ci.quantile.type = 8, # Hyndman and Fan (1996) recommended type 8, which is median-unbiased and does not depend on the distribution of x
                          boot.rep = 2000,
                          na.rm = FALSE,
                          verbose = TRUE,
                          ...) {

  distribution <- tolower(match.arg(distribution))
  boot.type <- tolower(match.arg(boot.type))
  standardizer <- tolower(match.arg(standardizer))
  boot.type <- tolower(match.arg(boot.type))
  parallel <- tolower(match.arg(parallel))
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) ||
      trim < 0 || trim > 0.45) {
    stop("'trim' must be a finite number in [0, 0.45].", call. = FALSE)
  }

  vars <- all.vars(formula)
  y.name <- vars[1]
  group.name <- vars[2]
  if (length(vars) != 2L || !all(vars %in% names(data))) {
    stop("Use a formula of the form outcome ~ group with columns in 'data'.",
         call. = FALSE)
  }
  weight.column <- NULL
  observation.weights <- eval(substitute(weights), envir = data,
                              enclos = parent.frame())
  analysis <- .record_analysis("boot.hct.test", match.call(), environment(),
                               overrides = list(weights = observation.weights))
  if (!is.null(observation.weights)) {
    if (boot.type == "parametric" && distribution != "norm") {
      stop("Weighted parametric generation currently supports distribution = 'norm' only.",
           call. = FALSE)
    }
    if (!is.numeric(observation.weights) || !is.null(dim(observation.weights)) ||
        length(observation.weights) != nrow(data)) {
      stop("Bootstrap weights must be a numeric vector with one value per data row.",
           call. = FALSE)
    }
    observed <- observation.weights[!is.na(observation.weights)]
    if (any(!is.finite(observed)) || any(observed < 0)) {
      stop("Observation weights must be finite and nonnegative.", call. = FALSE)
    }
    # Store weights in the resampled data so rows and weights stay paired.
    weight.column <- ".hcttest_observation_weights"
    while (weight.column %in% names(data)) weight.column <- paste0(weight.column, "_")
    data[[weight.column]] <- observation.weights
  }
  incomplete <- is.na(data[[y.name]]) | is.na(data[[group.name]])
  if (!is.null(weight.column)) {
    incomplete <- incomplete | is.na(data[[weight.column]])
  }
  if (any(incomplete)) {
    if (!na.rm) {
      stop("Outcome, group, or observation weights contain missing values; use na.rm = TRUE to remove incomplete rows.",
           call. = FALSE)
    }
    data <- data[!incomplete, , drop = FALSE]
  }
  if (!is.null(weight.column)) {
    data <- data[data[[weight.column]] > 0, , drop = FALSE]
  }
  data[[group.name]] <- droplevels(as.factor(data[[group.name]]))
  levels <- levels(data[[group.name]])
  if (length(levels) != 2L) {
    stop("The group variable must have exactly two observed levels.",
         call. = FALSE)
  }
  .check_retained_sample_sizes(sum(data[[group.name]] == levels[1]),
                               sum(data[[group.name]] == levels[2]),
                               if (is.null(weight.column)) trim else 0)
  n.total <- nrow(data)
  
  if(verbose) { 
    # starting time
    start.time <- Sys.time()
  }
  
  if(tolower(standardizer) == "sdp") {
    method <- "constant"
  } else {
    method <- tolower(standardizer)
  }
  trim.calibration <- NULL
  if (!is.null(weight.column) && trim > 0 && scale) {
    trim.calibration <- winsor.pop.sd(mu = 0, sigma = 1, trim = trim,
                                      distribution = "norm") / sqrt(1 - 2 * trim)
  }
  
  
  if (boot.type == "parametric") {
    
    # Estimate generating parameters from untrimmed analysis data.
    # The statistic below still applies the user's trim and scale settings.
    generating.weights <- if (is.null(weight.column)) NULL else data[[weight.column]]
    mle.parms <- hct.test(formula = formula, data = data,
                              weights = generating.weights,
                              verbose = FALSE, trim = 0, na.rm = na.rm, ...)
    
    mle <- list(
      dist = distribution,
      m1 = mle.parms$parms$m1,
      m2 = mle.parms$parms$m2,
      sd1 = mle.parms$parms$sd1,
      sd2 = mle.parms$parms$sd2,
      n1 = sum(data[[group.name]] == levels[1]),
      n2 = sum(data[[group.name]] == levels[2])
    )
    mle$precision.weighted <- !is.null(weight.column)
    if (mle$precision.weighted) {
      first <- data[[group.name]] == levels[1]
      second <- data[[group.name]] == levels[2]
      w1 <- data[[weight.column]][first]
      w2 <- data[[weight.column]][second]
      # Group-specific residual variance scales for precision-weighted normals.
      # Their ratio to observation weights is invariant to weight rescaling.
      mle$sd1 <- sqrt(sum(w1 * (data[[y.name]][first] - mle$m1)^2) / (mle$n1 - 1))
      mle$sd2 <- sqrt(sum(w2 * (data[[y.name]][second] - mle$m2)^2) / (mle$n2 - 1))
      mle$observation.sd1 <- mle$sd1 / sqrt(w1)
      mle$observation.sd2 <- mle$sd2 / sqrt(w2)
    }
    
    # Generate each group from the selected distribution and untrimmed moments.
    parboot.gendat <- function(data, mle) {
      if (mle$precision.weighted) {
        data[[y.name]][data[[group.name]] == levels[1]] <-
          rnorm(n = mle$n1, mean = mle$m1, sd = mle$observation.sd1)
        data[[y.name]][data[[group.name]] == levels[2]] <-
          rnorm(n = mle$n2, mean = mle$m2, sd = mle$observation.sd2)
        return(data)
      }
      gendat <- gendat.hct.test(n1 = mle$n1, 
                                n2 = mle$n2,
                                mu1 =  mle$m1,
                                mu2 =  mle$m2,
                                sd1 =  mle$sd1,
                                sd2 =  mle$sd2,
                                distribution = mle$dist)
      data[[y.name]][data[[group.name]] == levels[1]] <- gendat$y1
      data[[y.name]][data[[group.name]] == levels[2]] <- gendat$y2
      # data[[y.name]][data[[group.name]] == levels[1]] <- rnorm(n = mle$n1, mean = mle$m1, sd = mle$sd1)
      # data[[y.name]][data[[group.name]] == levels[2]] <- rnorm(n = mle$n2, mean = mle$m2, sd = mle$sd2)
      return(data)
    }
    
  } else {
    
    parboot.gendat <-  function(data, mle) data
    mle <- NULL
    
  } # if parametric
  
  stat.fun <- function(data, idx = NULL, ...) {
    # permutation ci (disabled) assume two groups roughly has same variance and only location shifts
    # not useful as much as non-parametric and parametric 
    if (boot.type == "permutation") {
      data[[y.name]] <- data[[y.name]][idx] 
    } else if (boot.type == "ordinary") {
      data <- data[idx, , drop = FALSE]
    }
    
    if (!is.null(weight.column) && trim > 0) {
      return(.weighted_trimmed_boot_statistic(
        data, outcome = y.name, group = group.name, weights = weight.column,
        trim = trim, standardizer = standardizer, scale = scale,
        calibration = trim.calibration))
    }
    stat.weights <- if (is.null(weight.column)) NULL else data[[weight.column]]
    # Bootstrap SD1/SD2 selects the effect-size denominator. The weighted
    # raw statistic can use the supported Student variance in either case.
    test.method <- if (!is.null(stat.weights) && method %in% c("sd1", "sd2")) {
      "constant"
    } else method
    obj.rmd <- hct.test(formula = formula, data = data, method = test.method,
                        weights = stat.weights, verbose = FALSE,
                        trim = trim, na.rm = na.rm, ...)
    if (!is.null(stat.weights)) obj.rmd$parms$method <- method
    obj.smd <- robust.smd(obj.rmd, scale = scale, ci = FALSE, verbose = FALSE)
    return(c(obj.rmd$test$delta, obj.smd$es$d))
    
  } # stat.fun
  
  # standardized mean difference
  boot.obj <- boot::boot(data = data,
                         statistic = stat.fun,
                         R = boot.rep, 
                         sim = boot.type,
                         strata = data[[group.name]],
                         ran.gen = parboot.gendat, mle = mle,
                         parallel = parallel,
                         ncpus = ncpus)
  data.est <- boot.obj$t0
  # sweep(x = boot.obj$t, MARGIN = 1, STATS = boot.obj$t0, FUN = "+")
  if(boot.type == "permutation") boot.obj$t <- sweep(x = boot.obj$t, MARGIN = 2, STATS = boot.obj$t0, FUN = "+") # for 'permutation' (disabled)
  boot.est <- colMeans(boot.obj$t, na.rm = na.rm)
  boot.se <- apply(X = boot.obj$t, FUN = sd, MARGIN = 2, na.rm = na.rm)

  # boot.lcl <- apply(X = boot.obj$t, FUN = function(x) quantile(x = x, probs = alpha / 2, type = ci.quantile.type, na.rm = na.rm), MARGIN = 2)
  # boot.ucl <-  apply(X = boot.obj$t, FUN = function(x) quantile(x = x, probs = 1 - alpha / 2, type = ci.quantile.type, na.rm = na.rm), MARGIN = 2)
  ci.obj.idx1.perc <- boot::boot.ci(boot.obj, conf = 1 - alpha, type = "perc", index = 1)$percent
  ci.obj.idx2.perc  <- boot::boot.ci(boot.obj, conf = 1 - alpha, type = "perc", index = 2)$percent
  boot.lcl.perc <- c(ci.obj.idx1.perc [4], ci.obj.idx2.perc [4])
  boot.ucl.perc <- c(ci.obj.idx1.perc [5], ci.obj.idx2.perc [5])
  if(boot.type == "ordinary") {
    ci.obj.idx1.bca  <- boot::boot.ci(boot.obj, conf = 1 - alpha, type = "bca", index = 1)$bca
    ci.obj.idx2.bca <- boot::boot.ci(boot.obj, conf = 1 - alpha, type = "bca", index = 2)$bca
    boot.lcl.bca <- c(ci.obj.idx1.bca[4], ci.obj.idx2.bca[4])
    boot.ucl.bca <- c(ci.obj.idx1.bca[5], ci.obj.idx2.bca[5])
  } else {
    boot.lcl.bca <- c(NA, NA)
    boot.ucl.bca <- c(NA, NA)
  }
 
  boot.out <- data.frame(est.type = c("Raw Scores", "Standardized"),
                         est.data = c(data.est[1], data.est[2]),
                         est.boot = c(boot.est[1], boot.est[2]),
                         se.boot = c(boot.se[1], boot.se[2]),
                         lcl.boot.perc = c(boot.lcl.perc[1], boot.lcl.perc[2]),
                         ucl.boot.perc = c(boot.ucl.perc[1], boot.ucl.perc[2]),
                         lcl.boot.bca = c(boot.lcl.bca[1], boot.lcl.bca[2]),
                         ucl.boot.bca = c(boot.ucl.bca[1], boot.ucl.bca[2]))
  
  if(verbose) {
    # tiend time
    end.time <- Sys.time()
    run.time <- as.numeric(difftime(end.time, start.time, units = "secs"))
    hours   <- floor(run.time / 3600)
    minutes <- floor((run.time %% 3600) / 60)
    seconds <- round(run.time %% 60)
    runtime_string <- sprintf("%02d hours : %02d minutes : %02d seconds\n\n",
                              hours, minutes, seconds)
    
    cat(runtime_string)
    
    print(boot.out)
    
  }
  
  attr(boot.out, "analysis") <- analysis
  class(boot.out) <- c("boot.hct.test", "data.frame")
  return(invisible(boot.out))
  
} # boot.hct.test
# when permutation is used, the confidence interval is derived from null distribution + the estimate
# this is same as using p-values, and also similar to central t approach to ci construction


#####################################
# Variance Ratio Constraints        # 
# Validation: TRUE                  #
#####################################
vrc <- function(var.ratio, n1 = 30, n2 = 30,
                standardizer = c("sdp", "sd1", "sd2", "avg")) {
  
  kappa <- n1 / n2
  standardizer <- match.arg(standardizer)
  
  if(standardizer == "sdp") {
    
    sd2  <- sqrt(
      (n2 * (kappa + 1) - 2) / (n2 - 1 + var.ratio * (n2 * kappa - 1))
    )
    
    sd1 <- sqrt(var.ratio * sd2^2)
    
  } else if(standardizer == "sd2") {
    
    sd2 <- 1
    sd1 <- sqrt(var.ratio * sd2^2)
    
  } else if(standardizer == "sd1") {
    
    sd1 <- 1
    sd2 <- sqrt(sd1^2 / var.ratio)
    
  } else if(standardizer == "avg") {
    
    sd2 = sqrt(2 / (var.ratio + 1))
    sd1 <- sqrt(var.ratio * sd2^2)
    
  } else {
    
    stop("Not a valid standardizer for this function", call. = FALSE)
    
  } # else if
  
  data.frame(sd1 = sd1, sd2 = sd2)
  
} # variance ratio constraint 



##############################################
#  Power based on lambda-prime distribution  #
##############################################
#' Prospective power for independent groups
#'
#' @description Calculate power, sample size, or a detectable standardized effect
#'   using constant-variance or heteroskedasticity-consistent methods.
#' @param d Planned standardized mean difference; NULL to solve for the effect. Uses the pooled within-group SD, except sd1 and sd2 use the corresponding group SD. No Hedges correction is applied.
#' @param null.d Reference standardized effect. Margins are relative to this value.
#' @param margin Scalar null margin for one.sided or two.sided. For two.one.sided, a nonzero scalar expands to symmetric bounds; alternatively supply two distinct bounds.
#' @param req.sign Direction for an effect-size search: + or -. For two.one.sided, 0 searches inside the equivalence bounds and reports the attainable range in a warning.
#' @param var.ratio Positive ratio of group 1 variance to group 2 variance.
#' @param n.ratio Positive allocation ratio n1/n2.
#' @param n2 Second-group sample size; NULL to solve for sample size. A supplied value must be an integer of at least two.
#' @param power Target power between 0.01 and 0.99; NULL to calculate power. Exactly one of d, n2, and power must be NULL.
#' @param alpha Significance level, strictly between zero and one.
#' @param alternative two.sided, one.sided, or two.one.sided. One-sided direction follows the effect relative to the null margin. two.one.sided uses equivalence testing inside the bounds and minimum-effect testing outside them.
#' @param claim.basis md.pval uses the raw mean-difference SE with a unit standardizer. smd.ci also includes uncertainty in the estimated standardizer using the delta-method SMD variance.
#' @param ceil.n Round calculated n2 and reported n1 upward. As in pwrss, power uses n1 = n.ratio*n2, which may differ from the reported rounded n1.
#' @param verbose 0 suppresses printing; 1 prints the calculation header, hypotheses, and results; 2 also prints key parameters and definitions, following pwrss. Logical values are accepted as 0 or 1.
#' @param method hc2 (default), constant, student (alias for constant), sd1, sd2, hc0, hc1, welch (alias for hc2), hc3, hc4, hc4m, or hc5.
#' @param test.distribution t (default), z, lp, or lambdap (alias for lp). Lambda-prime calculations reuse the internal package functions.
#'
#' @details
#' Planning is for independent, unweighted, untrimmed groups. The interface,
#' root-search intervals, equivalence search, and rounding conventions follow
#' pwrss::power.t.welch(). Raw-data, object-based, and cluster-design arguments
#' from the earlier power function are not supported.
#' 
#' The planned standardizer is constrained to one while preserving var.ratio.
#' For md.pval the calculation uses Var(MD). For smd.ci it uses
#' Var(d) = Var(MD) + d^2 Var(S), where S is the selected unit standardizer.
#' The standardizer-variance expressions follow Bulus (2026), Equations (5)-(8).
#' 
#' hc2 and welch use Welch-Satterthwaite degrees of freedom, enabling comparison
#' with power.t.welch(). Other HC methods and constant use n1+n2-2. sd1 and sd2
#' use n1-1 and n2-1, respectively. The z distribution ignores degrees of freedom.
#' The t and Lambda-prime calculations require at least three degrees of freedom.
#' These plug-in calculations are distribution-based planning approximations;
#' matching software calculations does not establish their sampling accuracy for
#' every HC method or nonnormal population.
#' 
#' For an equivalence effect-size search (req.sign = 0), d is the midpoint of the
#' attainable effect range, with endpoints reported in a warning. Power at that
#' midpoint can exceed the requested power, as in power.t.welch(). Effects exactly
#' on a two.one.sided boundary raise an error because the requested testing branch
#' is ambiguous. A continuous sample-size search starts at n2 = 3 and ends at 1e10;
#' if the lower endpoint implies invalid group sizes or degrees of freedom, the
#' search reports an infeasible design, following the reference search convention.
#'
#' @return
#' Invisibly, a list of class hct.power. It contains parms (original arguments),
#' test, df, ncp, null.ncp, t.alpha (critical values), d, power, n (n1 and n2), and
#' n.total. Additional fields are requested, method, test.distribution, claim.basis,
#' se (SE used), se.md, se.smd, variance.md, variance.smd, variance.standardizer,
#' and sd (planned group SDs). The new interface returns a list rather than the
#' numeric power value returned by the earlier object/noncentrality-based function.
#'
#' @references Bulus, M. (2026). Robust standard errors and confidence intervals for standardized mean difference. Research Synthesis Methods. Advance online publication.
#'
#' @seealso hct.test, robust.smd
#'
#' @examples
#' # Calculate power.
#' power.hct.test(d = 0.5, n2 = 60, var.ratio = 2,
#'                method = "hc2", verbose = 0)
#' 
#' # Calculate sample size.
#' power.hct.test(d = 0.5, power = 0.80, var.ratio = 2,
#'                method = "hc2", claim.basis = "smd.ci", verbose = 0)
#' 
#' # Calculate a detectable effect.
#' power.hct.test(n2 = 60, power = 0.80, method = "constant", verbose = 0)
#' 
#' # Equivalence planning using the normal approximation.
#' power.hct.test(d = 0, margin = c(-0.3, 0.3), n2 = 150,
#'                alternative = "two.one.sided", test.distribution = "z",
#'                verbose = 0)
#' 
#' # Lambda-prime power using existing internal distribution functions.
#' power.hct.test(d = 0.5, n2 = 60, test.distribution = "lp", verbose = 0)
#'
#' @export power.hct.test
power.hct.test <- function(d = NULL, null.d = 0, margin = 0, req.sign = "+",
                           var.ratio = 1, n.ratio = 1, n2 = NULL,
                           power = NULL, alpha = 0.05,
                           alternative = c("two.sided", "one.sided", "two.one.sided"),
                           claim.basis = c("md.pval", "smd.ci"),
                           ceil.n = TRUE, verbose = 1,
                           method = c("hc2", "constant", "student", "sd1", "sd2",
                                      "hc0", "hc1", "welch", "hc3", "hc4", "hc4m", "hc5"),
                           test.distribution = c("t", "z", "lp", "lambdap")) {
  alternative <- match.arg(alternative)
  claim.basis <- match.arg(claim.basis)
  method <- match.arg(method)
  test.distribution <- match.arg(test.distribution)
  func.parms <- as.list(environment())
  if (method == "student") method <- "constant"
  if (test.distribution == "lambdap") test.distribution <- "lp"
  
  scalar <- function(x, name, lower = -Inf, upper = Inf) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x <= lower || x >= upper)
      stop("'", name, "' must be a finite scalar between ", lower,
           " and ", upper, " (exclusive).", call. = FALSE)
  }
  if (!is.null(d)) scalar(d, "d")
  scalar(null.d, "null.d")
  scalar(var.ratio, "var.ratio", 0)
  scalar(n.ratio, "n.ratio", 0)
  scalar(alpha, "alpha", 0, 1)
  if (!is.null(power)) {
    scalar(power, "power", 0, 1)
    if (power < 0.01 || power > 0.99)
      stop("'power' must be between 0.01 and 0.99, as in pwrss.", call. = FALSE)
  }
  if (!is.null(n2)) {
    scalar(n2, "n2", 1)
    if (abs(n2 - round(n2)) > sqrt(.Machine$double.eps))
      stop("A supplied 'n2' must be an integer, as in pwrss.", call. = FALSE)
  }
  if (!is.logical(ceil.n) || length(ceil.n) != 1L || is.na(ceil.n))
    stop("'ceil.n' must be TRUE or FALSE.", call. = FALSE)
  if (is.logical(verbose) && length(verbose) == 1L && !is.na(verbose)) {
    verbose <- as.integer(verbose)
  } else if (!is.numeric(verbose) || length(verbose) != 1L ||
             !is.finite(verbose) || !(verbose %in% c(-1, 0, 1, 2))) {
    verbose <- 1
  }
  if (!is.numeric(margin) || !length(margin) || any(!is.finite(margin)))
    stop("'margin' must contain finite numbers.", call. = FALSE)
  if (alternative == "two.one.sided") {
    if (length(margin) == 1L && margin != 0) margin <- abs(margin) * c(-1, 1)
    if (length(margin) != 2L || min(margin) == max(margin))
      stop("Supply a nonzero scalar margin or two distinct bounds.", call. = FALSE)
  } else if (length(margin) != 1L) {
    stop("Supply one margin for one.sided or two.sided.", call. = FALSE)
  }
  missing.inputs <- c(is.null(d), is.null(n2), is.null(power))
  if (sum(missing.inputs) != 1L)
    stop("Exactly one of 'd', 'n2', and 'power' must be NULL.", call. = FALSE)
  requested <- c("es", "n", "power")[missing.inputs]
  
  # Direction conventions and root-search limits follow pwrss.
  check.pos_sign <- function(req.sign) {
    if (length(req.sign) != 1L || is.na(req.sign))
      stop("Invalid 'req.sign'.", call. = FALSE)
    if (req.sign %in% c("+", 1, "1", "+1", "positive", "pozitive")) return(TRUE)
    if (req.sign %in% c("-", -1, "-1", "negative")) return(FALSE)
    stop("'req.sign' must be '+' or '-'.", call. = FALSE)
  }
  check.null_sign <- function(req.sign, alternative) {
    if (length(req.sign) != 1L || is.na(req.sign))
      stop("Invalid 'req.sign'.", call. = FALSE)
    if (req.sign %in% c(" ", 0, "0", "", "null")) {
      if (alternative != "two.one.sided")
        stop("'req.sign = 0' requires two.one.sided.", call. = FALSE)
      return(TRUE)
    }
    check.pos_sign(req.sign)
    FALSE
  }
  
  # Retrieve the package functions once, and only when lp is requested.
  if (test.distribution == "lp") {
    lp.cdf <- plambdap
    lp.quantile <- qlambdap
  }
  cdf <- function(q, parameter, df, lower.tail = TRUE) {
    switch(test.distribution,
           z = stats::pnorm(q, mean = parameter, lower.tail = lower.tail),
           t = suppressWarnings(stats::pt(q, df = df, ncp = parameter, lower.tail = lower.tail)),
           lp = lp.cdf(q, df = df, ncp = parameter, lower.tail = lower.tail))
  }
  quantile <- function(p, parameter, df, lower.tail = TRUE) {
    switch(test.distribution,
           z = stats::qnorm(p, mean = parameter, lower.tail = lower.tail),
           t = suppressWarnings(stats::qt(p, df = df, ncp = parameter, lower.tail = lower.tail)),
           lp = lp.quantile(p, df = df, ncp = parameter, lower.tail = lower.tail))
  }
  
  planned.moments <- function(n1, n2, var.ratio, method, df.method) {
    total <- n1 + n2
    if (method == "sd1") {
      v1 <- 1
      v2 <- 1 / var.ratio
      var.standardizer <- 1 / (2 * (n1 - 1))
    } else if (method == "sd2") {
      v1 <- var.ratio
      v2 <- 1
      var.standardizer <- 1 / (2 * (n2 - 1))
    } else {
      # Scale the planned pooled standardizer to one, as in pwrss Welch power.
      denominator <- (n1 - 1) / (total - 2) * var.ratio +
        (n2 - 1) / (total - 2)
      v1 <- var.ratio / denominator
      v2 <- 1 / denominator
      var.standardizer <- ((n1 - 1) * v1^2 + (n2 - 1) * v2^2) /
        (2 * (total - 2)^2)
    }
    leverage <- c(1 / n1, 1 / n2)
    sizes <- c(n1, n2)
    if (method == "constant") {
      variance <- 1 / n1 + 1 / n2
      gamma <- c(1 / n1, 1 / n2)
    } else if (method %in% c("hc2", "welch")) {
      gamma <- 1 / sizes
      variance <- sum(gamma * c(v1, v2))
    } else if (method == "hc0") {
      gamma <- (sizes - 1) / sizes^2
      variance <- sum(gamma * c(v1, v2))
    } else if (method == "hc1") {
      gamma <- total / (total - 2) * (sizes - 1) / sizes^2
      variance <- sum(gamma * c(v1, v2))
    } else if (method == "hc3") {
      gamma <- 1 / (sizes - 1)
      variance <- sum(gamma * c(v1, v2))
    } else if (method %in% c("hc4", "hc4m", "hc5")) {
      relative <- total / (2 * sizes)
      exponent <- if (method == "hc4") {
        pmin(4, relative)
      } else if (method == "hc4m") {
        pmin(1, relative) + pmin(1.5, relative)
      } else {
        pmin(relative, max(4, 0.7 * max(leverage) / (2 / total))) / 2
      }
      gamma <- (1 - leverage)^(1 - exponent) / sizes
      variance <- sum(gamma * c(v1, v2))
    } else {
      # sd1 and sd2 use one group's variance for the mean-difference SE.
      variance <- 1 / n1 + 1 / n2
      gamma <- c(1 / n1, 1 / n2)
    }
    df <- if (df.method == "welch-satterthwaite") {
      if (method %in% c("sd1", "sd2"))
        stop("Welch-Satterthwaite df are unavailable for sd1/sd2.", call. = FALSE)
      sum(gamma * c(v1, v2))^2 / sum((gamma * c(v1, v2))^2 / (sizes - 1))
    } else if (method == "sd1") n1 - 1 else if (method == "sd2") n2 - 1 else total - 2
    if (any(!is.finite(c(v1, v2, variance, var.standardizer, df))) ||
        variance <= 0 || df <= 0)
      stop("The planned moments or degrees of freedom are invalid.", call. = FALSE)
    list(variance = variance, var.standardizer = var.standardizer,
         df = df, sd1 = sqrt(v1), sd2 = sqrt(v2))
  }
  
  
  planned.power <- function(d, null.d, margin, var.ratio, n2, n.ratio,
                            alpha, alternative, claim.basis) {
    n1 <- n.ratio * n2
    if (!is.finite(n1) || n1 <= 1 || n2 <= 1)
      stop("Each planned group size must exceed one.", call. = FALSE)
    df.method <- if (method %in% c("hc2", "welch"))
      "welch-satterthwaite" else "regression"
    moments <- planned.moments(n1, n2, var.ratio, method, df.method)
    df <- moments$df
    if (test.distribution != "z" && df < 3)
      stop("Degrees of freedom cannot be smaller than 3, as in pwrss.", call. = FALSE)
    # Standardizer is constrained to one. Equations (5)-(8), pp. 11-12.
    variance.smd <- moments$variance + d^2 * moments$var.standardizer
    se <- sqrt(if (claim.basis == "md.pval") moments$variance else variance.smd)
    ncp <- (d - null.d) / se
    null.ncp <- margin / se
    if (alternative == "two.sided") {
      critical <- c(quantile(alpha / 2, null.ncp, df),
                    quantile(alpha / 2, null.ncp, df, lower.tail = FALSE))
      probability <- cdf(critical[1], ncp, df) +
        cdf(critical[2], ncp, df, lower.tail = FALSE)
    } else if (alternative == "one.sided") {
      lower.tail <- ncp < null.ncp
      critical <- quantile(alpha, null.ncp, df, lower.tail)
      probability <- cdf(critical, ncp, df, lower.tail)
    } else if (ncp > min(null.ncp) && ncp < max(null.ncp)) {
      critical <- c(quantile(alpha, min(null.ncp), df, lower.tail = FALSE),
                    quantile(alpha, max(null.ncp), df))
      probability <- cdf(critical[2], ncp, df) +
        cdf(critical[1], ncp, df, lower.tail = FALSE) - 1
    } else if (ncp < min(null.ncp) || ncp > max(null.ncp)) {
      critical <- c(quantile(alpha / 2, min(null.ncp), df),
                    quantile(alpha / 2, max(null.ncp), df, lower.tail = FALSE))
      probability <- cdf(critical[1], ncp, df) +
        cdf(critical[2], ncp, df, lower.tail = FALSE)
    } else {
      stop("An effect exactly on a two.one.sided boundary is ambiguous; ",
           "specify an effect inside or outside the bounds.", call. = FALSE)
    }
    if (!is.finite(probability)) stop("Nonfinite power calculation.", call. = FALSE)
    list(power = min(1, max(0, probability)), t.alpha = critical,
         ncp = ncp, null.ncp = null.ncp, df = df, se = se,
         se.md = sqrt(moments$variance), se.smd = sqrt(variance.smd),
         variance.md = moments$variance, variance.smd = variance.smd,
         variance.standardizer = moments$var.standardizer,
         sd = c(sd1 = moments$sd1, sd2 = moments$sd2))
  }
  
  min.planned.power <- function(d, n2, power) {
    power - planned.power(d = d, null.d = null.d, margin = margin, var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio,
                          alpha = alpha, alternative = alternative, claim.basis = claim.basis)$power
  } # min.planned.power (for uniroot)
  
  if (requested == "n") {
    
    n2 <- try(stats::uniroot(function(n2) min.planned.power(d, n2, power), interval = c(3, 1e10))$root, silent = TRUE)
    if (inherits(n2, "try-error") || n2 == 1e10) stop("Design is not feasible. Possibly very large effect size.", call. = FALSE)
    
    n2 <- ifelse(ceil.n, ceiling(n2), n2)
    
  } else if (requested == "es") {
    
    if (alternative == "two.one.sided" && check.null_sign(req.sign, alternative)) {
      
      lower.int <- c(min(margin) + null.d, mean(margin) + null.d) + c(+1e-7, 0)
      upper.int <- c(mean(margin) + null.d, max(margin) + null.d) + c(0, -1e-7)
      d.lower <- stats::optimize(f = function(d) min.planned.power(d, n2, power) ^ 2, interval = lower.int, tol = 1e-12)$minimum
      d.upper <- stats::optimize(f = function(d) min.planned.power(d, n2, power) ^ 2, interval = upper.int, tol = 1e-12)$minimum
      
      d <- mean(c(d.lower, d.upper))
      
      pwr.lower <- planned.power(d = d.lower, null.d = null.d, margin = margin, var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio,
                                 alpha = alpha, alternative = alternative, claim.basis = claim.basis)$power
      pwr.upper <- planned.power(d = d.upper, null.d = null.d, margin = margin, var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio,
                                 alpha = alpha, alternative = alternative, claim.basis = claim.basis)$power
      
      if (round(pwr.lower, 3) >= power && round(pwr.upper, 3) >= power) {
        
        warning("Target effect ranges from ", round(d.lower, 4), " to ", round(d.upper, 4), " within the null bounds.",
                call. = FALSE)
        
      } else {
        
        warning("The target power rate cannot be achieved within the null bounds.", call. = FALSE)
        
      }
      
    } else {
      
      if (check.pos_sign(req.sign)) {
        d.int <- c(max(margin) + null.d, 10) + c(+1e-7, -1e-7)
      } else {
        d.int <- c(-10, min(margin) + null.d) + c(+1e-7, -1e-7)
      }
      
      d <- try(stats::uniroot(f = function(d) min.planned.power(d, n2, power), interval = d.int, tol = 1e-12)$root, silent = TRUE)
      if (inherits(d, "try-error")) stop("Design is not feasible.", call. = FALSE)
      
    } # two.one.sided?
    
  } # ss or es?
  
  n1 <- ifelse(ceil.n, ceiling(n.ratio * n2), n.ratio * n2)
  n <- c(n1 = n1, n2 = n2)
  
  # calculate power (if requested == "power") or update it (if requested == "n")
  pwr.obj <- planned.power(d = d, null.d = null.d, margin = margin, var.ratio = var.ratio, n2 = n2, n.ratio = n.ratio,
                           alpha = alpha, alternative = alternative, claim.basis = claim.basis)
  
  result <- list(parms = func.parms, test = test.distribution,
                 df = pwr.obj$df, ncp = pwr.obj$ncp, null.ncp = pwr.obj$null.ncp,
                 t.alpha = pwr.obj$t.alpha, d = d, power = pwr.obj$power,
                 n = n, n.total = sum(n), requested = requested,
                 method = method, test.distribution = test.distribution,
                 claim.basis = claim.basis, se = pwr.obj$se,
                 se.md = pwr.obj$se.md, se.smd = pwr.obj$se.smd,
                 variance.md = pwr.obj$variance.md,
                 variance.smd = pwr.obj$variance.smd,
                 variance.standardizer = pwr.obj$variance.standardizer,
                 sd = pwr.obj$sd)
  if (verbose > 0) {
    # All printing helpers remain local. ASCII formatting matches pwrss 1.3.1.
    topic <- function(label) paste0(strrep("-", 52), "\n", label, "\n",
                                    strrep("-", 52), "\n")
    format.value <- function(value, signed = FALSE) {
      integer <- is.numeric(value) &&
        !any(abs(value - round(value)) > .Machine$double.eps^(2/3))
      text <- sprintf(paste0("%", if (signed) "+" else "",
                             ".", if (integer) 0 else 3, "f"), value)
      text <- paste(text, collapse = if (length(value) < 3) " and " else ", ")
      if (signed) text <- gsub("^\\+| \\+", " ", text)
      text
    }
    heading <- switch(requested,
                      n = "             SAMPLE SIZE CALCULATION              ",
                      power = "                POWER CALCULATION                 ",
                      es = "      MINIMUM DETECTABLE EFFECT CALCULATION       ")
    cat("+", strrep("-", 50), "+\n|", heading, "|\n+",
        strrep("-", 50), "+\n\n", sep = "")
    test.name <- if (method == "constant") "Student's" else
      if (method %in% c("hc2", "welch")) "Welch's" else toupper(method)
    distribution.name <- switch(test.distribution, t = "T-Test", z = "Z-Test",
                                lp = "Lambda-Prime Test")
    cat(test.name, " ", distribution.name, " (Independent Samples)\n\n", sep = "")
    
    if (alternative == "two.sided") {
      rhs <- if (margin == 0) "0" else "margin"
      h0 <- paste("d - null.d", " =", rhs)
      h1 <- paste("d - null.d", "!=", rhs)
    } else if (alternative == "one.sided") {
      rhs <- if (margin == 0) "0" else "margin"
      less <- result$ncp < result$null.ncp
      h0 <- paste("d - null.d", if (less) ">=" else "<=", rhs)
      h1 <- paste("d - null.d", if (less) " <" else " >", rhs)
    } else if (result$ncp > min(result$null.ncp) &&
               result$ncp < max(result$null.ncp)) {
      h0 <- c("d - null.d <= min(margin) or", "d - null.d >= max(margin)")
      h1 <- c("d - null.d  > min(margin) and", "d - null.d  < max(margin)")
    } else {
      h0 <- c("d - null.d >= min(margin) and", "d - null.d <= max(margin)")
      h1 <- c("d - null.d  < min(margin) or", "d - null.d  > max(margin)")
    }
    continuation <- paste0("\n", strrep(" ", 21))
    cat(topic("Hypotheses"))
    cat(sprintf("  H0 (Null)        : %s\n", paste(h0, collapse = continuation)))
    cat(sprintf("  H1 (Alternative) : %s\n\n", paste(h1, collapse = continuation)))
    
    if (verbose == 2) {
      keys <- list(Margin = margin, `Degrees of Freedom` = result$df,
                   `Non-centrality of Alt.` = result$ncp,
                   `Non-centrality of Null` = result$null.ncp,
                   `Critical Value` = result$t.alpha)
      signed <- any(unlist(keys) < 0)
      width <- max(nchar(names(keys)))
      cat(topic("Key Parameters"))
      for (name in names(keys))
        cat("  ", name, strrep(" ", width - nchar(name)), " = ",
            format.value(keys[[name]], signed), "\n", sep = "")
      cat("\n")
    }
    cat(topic("Results"))
    marker <- function(quantity) if (requested == quantity) "  <<" else ""
    cat("  Effect Size (d)      = ", format.value(d), " (vs. null.d = ",
        format.value(null.d), ")", marker("es"), "\n", sep = "")
    cat("  Sample Size          = ", paste(round(n, 3), collapse = " and "),
        marker("n"), "\n", sep = "")
    cat(sprintf("  Type 1 Error (alpha) = %.3f\n", alpha))
    cat(sprintf("  Type 2 Error (beta)  = %.3f\n", 1 - result$power))
    cat(sprintf("  Statistical Power    = %.3f%s\n\n", result$power, marker("power")))
    if (verbose == 2) {
      cat(topic("Definitions"))
      cat("  Margin : Smallest d - null.d difference that matters\n\n")
    }
  }
  invisible(structure(result, class = "hct.power"))
} # power.hct.test()
