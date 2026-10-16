\name{meta.robust.mdes}
\alias{meta.robust.mdes}
\title{Append study-specific minimum detectable standardized effects}
\description{Run \code{power.hct.test()} internally for each study and append its solved effect to the input data.}
\usage{
meta.robust.mdes(data, sd1 = "sd1", sd2 = "sd2",
                            n1 = "n1", n2 = "n2", power = 0.80,
                            alpha = 0.05,
                            alternative = c("two.sided", "one.sided"),
                            req.sign = "+", null.d = 0, margin = 0,
                            method = "hc2",
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
\item{power}{Target power, strictly between zero and one; default 0.80.}
\item{alternative}{Two-sided (default) or one-sided detection. Equivalence effect-size searches are excluded because the reported midpoint is not an MDES.}
\item{req.sign}{Search direction, \code{"+"} (default) or \code{"-"}. Direction is not inferred from group means.}
}
\details{
Each row calls \code{power.hct.test(d = NULL, ...)} with the requested target power, the row's n2, variance ratio \code{(sd1/sd2)^2}, and allocation ratio \code{n1/n2}. The returned d is stored unchanged, including its sign. Means are not required. The selected method determines the standardizer, as in \code{power.hct.test()}. Solver results must attain target power within an absolute tolerance of 0.0001.

All calls use \code{ceil.n = FALSE} and \code{verbose = 0}. Rows and existing columns are retained. An existing mdes column is replaced. Missing or invalid row inputs produce an error identifying the row; no rows are silently dropped. The input object is not modified in place. There is no pooling, metafor call, meta-analytic weighting, interval calculation, or dedicated plotting method.
}
\value{The input data frame with one numeric column, \code{mdes}, appended or replaced.}
\seealso{\code{\link{power.hct.test}}, \code{\link{meta.robust.power}}}
\examples{
studies <- data.frame(n1 = c(107, 86), n2 = c(67, 294),
                      sd1 = c(2, 4), sd2 = c(11, 7))
meta.robust.mdes(studies, power = 0.80, method = 'hc2')
}
