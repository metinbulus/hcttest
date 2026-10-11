% Maintained help file; source documentation is in R/power-hct-test.R.
\name{power.hct.test}
\alias{power.hct.test}
\title{Prospective power for independent groups}
\description{Calculate power, sample size, or a detectable standardized effect using constant-variance or heteroskedasticity-consistent methods.}
\usage{
power.hct.test(d = NULL, null.d = 0, margin = 0, req.sign = "+",
                           var.ratio = 1, n.ratio = 1, n2 = NULL,
                           power = NULL, alpha = 0.05,
                           alternative = c("two.sided", "one.sided", "two.one.sided"),
                           claim.basis = c("md.pval", "smd.ci"),
                           ceil.n = TRUE, verbose = 1,
                           method = c("hc2", "constant", "student", "sd1", "sd2",
                                      "hc0", "hc1", "welch", "hc3", "hc4", "hc4m", "hc5"),
                           test.distribution = c("t", "z", "lp", "lambdap"))
}
\arguments{
  \item{d}{Planned standardized mean difference; NULL to solve for the effect. Uses the pooled within-group SD, except sd1 and sd2 use the corresponding group SD. No Hedges correction is applied.}
  \item{null.d}{Reference standardized effect. Margins are relative to this value.}
  \item{margin}{Scalar null margin for one.sided or two.sided. For two.one.sided, a nonzero scalar expands to symmetric bounds; alternatively supply two distinct bounds.}
  \item{req.sign}{Direction for an effect-size search: + or -. For two.one.sided, 0 searches inside the equivalence bounds and reports the attainable range in a warning.}
  \item{var.ratio}{Positive ratio of group 1 variance to group 2 variance.}
  \item{n.ratio}{Positive allocation ratio n1/n2.}
  \item{n2}{Second-group sample size; NULL to solve for sample size. A supplied value must be an integer of at least two.}
  \item{power}{Target power between 0.01 and 0.99; NULL to calculate power. Exactly one of d, n2, and power must be NULL.}
  \item{alpha}{Significance level, strictly between zero and one.}
  \item{alternative}{two.sided, one.sided, or two.one.sided. One-sided direction follows the effect relative to the null margin. two.one.sided uses equivalence testing inside the bounds and minimum-effect testing outside them.}
  \item{claim.basis}{md.pval uses the raw mean-difference SE with a unit standardizer. smd.ci also includes uncertainty in the estimated standardizer using the delta-method SMD variance.}
  \item{ceil.n}{Round calculated n2 and reported n1 upward. As in pwrss, power uses n1 = n.ratio*n2, which may differ from the reported rounded n1.}
  \item{verbose}{0 suppresses printing; 1 prints the calculation header, hypotheses, and results; 2 also prints key parameters and definitions, following pwrss. Logical values are accepted as 0 or 1.}
  \item{method}{hc2 (default), constant, student (alias for constant), sd1, sd2, hc0, hc1, welch (alias for hc2), hc3, hc4, hc4m, or hc5.}
  \item{test.distribution}{t (default), z, lp, or lambdap (alias for lp). Lambda-prime calculations reuse the internal package functions.}
}
\details{
Planning is for independent, unweighted, untrimmed groups. The interface,
root-search intervals, equivalence search, and rounding conventions follow
pwrss::power.t.welch(). Raw-data, object-based, and cluster-design arguments
from the earlier power function are not supported.

The planned standardizer is constrained to one while preserving var.ratio.
For md.pval the calculation uses Var(MD). For smd.ci it uses
Var(d) = Var(MD) + d^2 Var(S), where S is the selected unit standardizer.
The standardizer-variance expressions follow Bulus (2026), Equations (5)-(8).

hc2 and welch use Welch-Satterthwaite degrees of freedom, enabling comparison
with power.t.welch(). Other HC methods and constant use n1+n2-2. sd1 and sd2
use n1-1 and n2-1, respectively. The z distribution ignores degrees of freedom.
The t and Lambda-prime calculations require at least three degrees of freedom.
These plug-in calculations are distribution-based planning approximations;
matching software calculations does not establish their sampling accuracy for
every HC method or nonnormal population.

For an equivalence effect-size search (req.sign = 0), d is the midpoint of the
attainable effect range, with endpoints reported in a warning. Power at that
midpoint can exceed the requested power, as in power.t.welch(). Effects exactly
on a two.one.sided boundary raise an error because the requested testing branch
is ambiguous. A continuous sample-size search starts at n2 = 3 and ends at 1e10;
if the lower endpoint implies invalid group sizes or degrees of freedom, the
search reports an infeasible design, following the reference search convention.
}
\value{
Invisibly, a list of class hct.power. It contains parms (original arguments),
test, df, ncp, null.ncp, t.alpha (critical values), d, power, n (n1 and n2), and
n.total. Additional fields are requested, method, test.distribution, claim.basis,
se (SE used), se.md, se.smd, variance.md, variance.smd, variance.standardizer,
and sd (planned group SDs). The new interface returns a list rather than the
numeric power value returned by the earlier object/noncentrality-based function.
}
\references{
Bulus, M. (2026). Robust standard errors and confidence intervals for standardized mean difference. Research Synthesis Methods. Advance online publication.
}
\seealso{\code{\link{hct.test}}, \code{\link{robust.smd}}}
\examples{
# Calculate power.
power.hct.test(d = 0.5, n2 = 60, var.ratio = 2,
               method = "hc2", verbose = 0)

# Calculate sample size.
power.hct.test(d = 0.5, power = 0.80, var.ratio = 2,
               method = "hc2", claim.basis = "smd.ci", verbose = 0)

# Calculate a detectable effect.
power.hct.test(n2 = 60, power = 0.80, method = "constant", verbose = 0)

# Equivalence planning using the normal approximation.
power.hct.test(d = 0, margin = c(-0.3, 0.3), n2 = 150,
               alternative = "two.one.sided", test.distribution = "z",
               verbose = 0)

# Lambda-prime power using existing internal distribution functions.
power.hct.test(d = 0.5, n2 = 60, test.distribution = "lp", verbose = 0)
}
