\name{smd.hc.t.test}
\alias{smd.hc.t.test}


\title{Standardized Mean Difference for Heteroscedasticity-Consistent (HC) T-test}

\description{
Standardized mean difference for heteroscedasticity-consistent (HC) t-test procedures. Parameters are extracted from hc.t.test() function.  
}

\usage{
smd.hc.t.test(x = NULL,
              ss.adjust = FALSE, 
              ci = c("nct", "z", "ct", "FALSE"),
              verbose = TRUE)
}

\arguments{
  \item{x}{an object returned from hc.t.test() function}
  \item{ss.adjust}{logical; apply Hedges small sample correction?  Default is FALSE}
  \item{ci}{character; distribution to construct confidence intervals. "ct" for central t-distribution, "nct" for non-central t-distribution, "z" for normal distribution}
  \item{verbose}{logical; print results? Default is TRUE}

}

\examples{

y1 <- rnorm(20, 0, 5)
y2 <- rnorm(100, 0.50, 1)

results <- hc.t.test(y1, y2, method = "hc2", df = "welch")

smd.hc.t.test(results)

}
