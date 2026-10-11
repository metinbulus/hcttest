# skew and kurtosis for Tukey g-and-h transformation 
# were calculated following Headrick et al. (2008)

# Headrick, T. C., Kowalchuk, R. K., & Sheng, Y. (2008). 
# Parametric Probability Densities and Distribution Functions 
# for Tukey g-and-h Transformations and their Use for Fitting Data
# Applied Mathematical Sciences, 2(9), 449-462. 
# https://api.semanticscholar.org/CorpusID:59493469

skew.gh <- function(g, h) {
  
  if(g != 0 & h != 0) {
    
    m1.gh <- (exp(g^2 / (2 - 2*h)) - 1) / (g * sqrt(1 - h))
    m2.gh <- (1 - 2 * exp(g^2 / (2 - 4*h)) + exp(2*g^2 / (1 - 2*h))) / (g^2 * sqrt(1 - 2*h))
    m3.gh <- (3 * exp(g^2 / (2 - 6*h)) + exp(9*g^2 / (2 - 6*h)) - 3 * exp(2*g^2 / (1 - 3*h)) - 1) / (g^3 * sqrt(1 - 3*h))
    m4.gh <- (exp(8*g^2 / (1 - 4*h)) * (1 + 6 * exp(6*g^2 / (4*h - 1)) + exp(8*g^2 / (4*h - 1)) - 4 * exp(7*g^2 / (8*h - 2)) - 4 * exp(15*g^2 / (8*h - 2)))) / (g^4 * sqrt(1 - 4*h))
    
    (m3.gh - 3 * m2.gh * m1.gh + 2 * m1.gh^3) / (m2.gh - m1.gh^2)^(3/2)
    
  } else if(g != 0 & h == 0) {
    
    (3 * exp(2*g^2) + exp(3*g^2) - 4)^(1/2)
    
  } else if(g == 0 & h != 0) {
    
    0
    
  } else if(g == 0 & h == 0) {
    
    0
    
  }
  
} # skew.gh


kurt.gh <- function(g, h) {
  
  if(g != 0 & h != 0) {
    
    m1.gh <- (exp(g^2 / (2 - 2*h)) - 1) / (g * sqrt(1 - h))
    m2.gh <- (1 - 2 * exp(g^2 / (2 - 4*h)) + exp(2*g^2 / (1 - 2*h))) / (g^2 * sqrt(1 - 2*h))
    m3.gh <- (3 * exp(g^2 / (2 - 6*h)) + exp(9*g^2 / (2 - 6*h)) - 3 * exp(2*g^2 / (1 - 3*h)) - 1) / (g^3 * sqrt(1 - 3*h))
    m4.gh <- (exp(8*g^2 / (1 - 4*h)) * (1 + 6 * exp(6*g^2 / (4*h - 1)) + exp(8*g^2 / (4*h - 1)) - 4 * exp(7*g^2 / (8*h - 2)) - 4 * exp(15*g^2 / (8*h - 2)))) / (g^4 * sqrt(1 - 4*h))
    
    (m4.gh - 4 * m3.gh * m1.gh - 3 * m2.gh^2 + 12 * m2.gh * m1.gh^2 - 6 * m1.gh^4) / (m2.gh - m1.gh^2)^2
    
  } else if(g != 0 & h == 0) {
    
    3 * exp(2*g^2) + 2 * exp(3*g^2) + exp(4*g^2) - 6
    
  } else if(g == 0 & h != 0) {
    
    3 * (1 - 2*h)^3 * (1 / (1 - 4*h)^(5/2) + 1 / (2*h - 1)^3)
    
  } else if(g == 0 & h == 0) {
    
    0
    
  }
  
} # kurt.gh

# get g-and-h parameters given skew and kurt
skew_kurtosis_to_gh <- function(skew, kurt) {
  
  optim.out <- optim(par = c(g = 0, h = 0), 
                     fn = function(x) {
                       # if(x[2] > 0.25) x[2] <- 0.25 - 1e-6
                       skew0 <- skew.gh(g = x[1], h = x[2])
                       kurt0 <- kurt.gh(g = x[1], h = x[2])
                       abs(skew - skew0) + abs(kurt - kurt0)
                     }, 
                     control = list(abstol = 1e-9)
  )
  
  # cat("Convergence:", optim.out$convergence,
  #    "\nMessage:", optim.out$message, "\n")
  
  optim.out$par
  
} # skew_kurtosis_to_gh



# transform from z to g-and-h distribution
# Hoaglin, 1985?
tukey_gh_scaler <- function(x, g = 0, h = 0,
                            rescale = FALSE) {
  if (g == 0 && h == 0) return(x)
  a <- mean(x)
  b <- sd(x)
  if (is.finite(b) && b == 0) return(x)
  if (!is.finite(b) || b <= 0) {
    stop("g-and-h transformation requires a finite, positive sample SD.", call. = FALSE)
  }
  z <- (x - mean(x)) / sd(x)
  transformed <- if (g == 0) z else expm1(g * z) / g
  transformed <- transformed * exp(h * z^2 / 2)
  if (any(!is.finite(transformed)) || !is.finite(sd(transformed)) ||
      sd(transformed) <= 0) {
    stop("g-and-h parameters produced an undefined sample scale; use less extreme parameters.",
         call. = FALSE)
  }
  v <- a + b * transformed
  
  if(rescale) {
    
    v <- a + b * (v - mean(v)) / sd(v)
    
  }
  
  v
  
} # tukey_gh_scaler()

norm.gh <- function(x, g, h, skewness = NULL, kurtosis = NULL, rescale = TRUE) {
  
  if(!is.null(skewness) & !is.null(kurtosis)) {
    gh <- skew_kurtosis_to_gh(skew = skewness, kurt = kurtosis)
    v <- tukey_gh_scaler(x, g = gh[[1]], h = gh[[2]], rescale = rescale)
  } else {
    v <- tukey_gh_scaler(x, g = g, h = h, rescale = rescale)
  }
  
  return(invisible(v))
  
}


mixnorm_moments <- function(mu, sigma, weight) {
  # Validate inputs
  if (length(mu) != 2 || length(sigma) != 2 || length(weight) != 2) {
    stop("All inputs must be length 2")
  }
  if (abs(sum(weight) - 1) > 1e-10) stop("weights must sum to 1")
  
  # Extract
  mu1 <- mu[1]; mu2 <- mu[2]
  sigma1 <- sigma[1]; sigma2 <- sigma[2]
  w1 <- weight[1]; w2 <- weight[2]
  
  # MEAN
  mean_mix <- w1*mu1 + w2*mu2
  
  # VARIANCE
  E_Y2 <- w1*(sigma1^2 + mu1^2) + w2*(sigma2^2 + mu2^2)
  var_mix <- E_Y2 - mean_mix^2
  sigma_mix <- sqrt(var_mix)
  
  # THIRD CENTRAL MOMENT
  delta1 <- mu1 - mean_mix
  delta2 <- mu2 - mean_mix
  M3 <- w1*(delta1^3 + 3*delta1*sigma1^2) + w2*(delta2^3 + 3*delta2*sigma2^2)
  
  # FOURTH CENTRAL MOMENT
  M4 <- w1*(delta1^4 + 6*delta1^2*sigma1^2 + 3*sigma1^4) + 
    w2*(delta2^4 + 6*delta2^2*sigma2^2 + 3*sigma2^4)
  
  # STANDARDIZED MOMENTS
  skew <- M3 / var_mix^(3/2)
  kurt <- M4 / var_mix^2
  
  return(c(mean = mean_mix, sd = sigma_mix, skew = skew, kurtosis = kurt))
} # mixnorm_moments 
