\name{hc.t.test}
\alias{hc.t.test}
\alias{hc.t.test.default}
\alias{hc.t.test.formula}


\title{Heteroscedasticity-Consistent (HC) T-test}

\description{
Heteroscedasticity-consistent (HC) t-test procedures.  
}

\usage{

hc.t.test(x, \dots)

\method{hc.t.test}{default}(x = NULL, y = NULL, 
                  m1, m2, sd1, sd2, n1, n2,
                  alpha = 0.05, two.tailed = TRUE, 
                  ci = c("ct", "nct", "FALSE"),
                  method = c("constant", "student", "sd1", "sd2",
                             "hc0", "hc1", "hc2", "welch",
                             "hc3", "hc4", "hc4m", "hc5"),
                  df = c("regression", "welch-satterthwaite"),
                  verbose = TRUE, \dots)
                  
\method{hc.t.test}{formula}(formula, data, \dots)

}

\arguments{
  \item{x}{vector of observations for the first group}
  \item{y}{vector of observations for the second group}
  \item{m1}{mean of the first group}
  \item{m2}{mean of the second group}
  \item{sd1}{standard deviation of the first group}
  \item{sd2}{standard deviation of the second group}
  \item{n1}{sample size for the first group}
  \item{n2}{sample size for the second group}
  \item{alpha}{type 1 error rate}
  \item{two.tailed}{logical; a two sided test? Default is TRUE}
  \item{ci}{character; distribution to construct confidence intervals. "ct" for central t-distribution, "nct" for non-central t-distribution}
  \item{method}{method to estimate standard errors and margins. Any of the "constant", "student", "sd1", "sd2", "hc0", "hc1", "hc2", "welch", "hc3", "hc4", "hc4m", "hc5"}
  \item{df}{character; whether degrees of freedom should be based on linear regression ("regression") or Welch-Satterthwaite ("welch-satterthwaite")}
  \item{verbose}{logical; print results? Default is TRUE}
  \item{formula}{formula; in the form of outcome ~ group}
  \item{data}{data; data that includes 'outcome and group}
  \item{...}{for generic consistency}
  
}

\examples{

y1 <- rnorm(20, 0, 5)
y2 <- rnorm(100, 0.50, 1)
x <- rep(c(0,1), c(length(y1), length(y2)))
dat <- data.frame(y = c(y1,y2), x = x)

hc.t.test(y ~ x, dat, method = "hc2", df = "welch")
hc.t.test(y1, y2, method = "hc2", df = "welch")
hc.t.test(m1 = mean(y1), m2 = mean(y2),
          sd1 = sd(y1), sd2 = sd(y2),
          n1 = length(y1), n2 = length(y2),
          method = "hc2", df = "welch")

}
