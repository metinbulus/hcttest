\name{meta.robust.power}
\alias{meta.robust.power}
\title{Append study-specific power}
\description{Run \code{power.hct.test()} internally for each study and append its calculated power to the input data.}
\usage{
meta.robust.power(data, d = NULL, m1 = "m1", m2 = "m2",
                             sd1 = "sd1", sd2 = "sd2", n1 = "n1", n2 = "n2",
                             alpha = 0.05,
                             alternative = c("two.sided", "one.sided", "two.one.sided"),
                             null.d = 0, margin = 0, method = "hc2",
                             claim.basis = c("smd.ci", "md.pval"),
                             test.distribution = c("t", "z", "lp", "lambdap"))
}
\arguments{
\item{data}{A nonempty data frame with unique column names.}
\item{sd1, sd2, n1, n2}{Names of numeric columns containing positive finite group standard deviations and integer sample sizes of at least two. The variance ratio is group 1 divided by group 2; the allocation ratio is n1/n2.}
\item{alpha}{Significance level, strictly between zero and one; default 0.05.}
\item{null.d, margin}{Reference standardized effect and null margin passed unchanged to \code{power.hct.test()}. Margins are relative to null.d. For equivalence testing, a scalar margin expands to symmetric bounds, or supply two distinct bounds.}
\item{method}{Variance estimator passed to \code{power.hct.test()}: hc2 (default), constant, student (alias for constant), sd1, sd2, hc0, hc1, welch, hc3, hc4, hc4m, or hc5.}
\item{claim.basis}{\code{"smd.ci"} (default) uses the delta-method standardized-effect variance, including standardizer uncertainty. \code{"md.pval"} uses mean-difference inference with a unit standardizer. Passed unchanged to \code{power.hct.test()}.}
\item{test.distribution}{Reference distribution: t (default), z, lp, or lambdap (alias for lp), as supported by \code{power.hct.test()}.}
\item{d}{NULL (default), a finite numeric effect shared by all rows, a finite numeric vector with one effect per row, or the name of a numeric effect-size column. NULL calculates effects from the group means and SDs. Supplied effects must use the standardizer appropriate to method.}
\item{m1, m2}{Names of finite numeric group-mean columns. Used only when d is NULL.}
\item{alternative}{Two-sided (default), one-sided, or two.one.sided, as in \code{power.hct.test()}. One-sided direction follows the supplied effect relative to the null margin; this is not a separately fixed directional hypothesis.}
}
\details{
Each row calls \code{power.hct.test(power = NULL, ...)} with its effect, n2, variance ratio \code{(sd1/sd2)^2}, and allocation ratio \code{n1/n2}. With d = NULL, the effect is (m1-m2)/S. S is sd1 for method sd1, sd2 for method sd2, and the degrees-of-freedom-weighted pooled SD for all other methods. No Hedges correction or absolute-value transformation is used.

Power calculated from observed means is descriptive retrospective power, not independent evidence for a finding. For sensitivity analysis, supply a scientifically meaningful assumed effect through d.

All calls use \code{ceil.n = FALSE} and \code{verbose = 0}. Rows and existing columns are retained. An existing power column is replaced. Missing or invalid row inputs produce an error identifying the row; no rows are silently dropped. The input object is not modified in place. There is no pooling, metafor call, meta-analytic weighting, interval calculation, or dedicated plotting method.
}
\value{The input data frame with one numeric column, \code{power}, appended or replaced.}
\seealso{\code{\link{power.hct.test}}, \code{\link{meta.robust.mdes}}}
\examples{
studies <- data.frame(m1 = c(48, 6), m2 = c(50, 7),
  n1 = c(107, 86), n2 = c(67, 294), sd1 = c(2, 4), sd2 = c(11, 7))
meta.robust.power(studies, d = 0.30)
meta.robust.power(studies)
studies <- meta.robust.mdes(studies)
meta.robust.power(studies, d = 'mdes') # approximately 0.80 in each row
}
