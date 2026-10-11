\name{funnel}
\alias{funnel}
\title{Overlay Funnel Plots for Study-Level Effect Estimates}
\description{
Draws a funnel plot for one or two objects returned by
\code{meta.robust.smd()}. Each object has its own pooled reference estimate,
colored contours, and plotting symbols. Both objects share the same axes.
}
\usage{
funnel(x, y = NULL, effect = c("smd", "raw"),
  col = c("#001F3F", "#B35400"), pch = c(16, 17), cex = 1.1,
  level = c(0.90, 0.95, 0.99), shade = TRUE, connect = FALSE,
  contour.legend = TRUE, jitter = FALSE, jitter.width = 0.008,
  refline = NULL, xlim = NULL, ylim = NULL,
  legend = TRUE, legend.labels = NULL, xlab = NULL, main = NULL,
  model.legend = TRUE, legend.position = "topright",
  contour.position = "topleft", legend.cex = 0.8,
  contours = TRUE, show.refline = TRUE, grid = TRUE,
  shade.alpha = NULL, point.alpha = 0.8,
  contour.lty = c(1, 2, 3), contour.lwd = 0.9,
  label = FALSE, label.cex = 0.65, ...)
}
\arguments{
\item{x}{An object returned by \code{meta.robust.smd()}.}
\item{y}{An optional second object. Both objects must contain the same unique,
nonmissing study names. The second object is aligned to the first by name.}
\item{effect}{\code{"smd"} plots \code{yi} against \code{sqrt(vi)};
\code{"raw"} plots \code{delta} against \code{se.delta}.}
\item{col}{One color, or one color per object, used for points, reference lines,
and contour shading. Defaults to dark navy and dark orange.}
\item{pch}{One plotting symbol, or one per object.}
\item{cex}{A positive point-size multiplier.}
\item{level}{Pseudo confidence levels as probabilities strictly between zero
and one. Levels are sorted and duplicates removed.}
\item{shade}{Logical; draw translucent, separately shaded contour bands.}
\item{connect}{Logical; connect corresponding points when two objects are supplied.}
\item{contour.legend}{Logical; show the contour legend when \code{legend = TRUE}
and shading or contour boundaries are displayed.}
\item{jitter}{Logical; apply small, deterministic horizontal display offsets.
Defaults to \code{FALSE}.}
\item{jitter.width}{Maximum absolute horizontal displacement as a fraction of
the x-axis range, from zero to 0.05. Default is 0.008.}
\item{refline}{One finite reference value shared by both objects, or one per
object. \code{NULL} extracts the selected pooled estimate from each object's
\code{pooled} attribute. Use \code{refline = 0} for zero-centered contours.}
\item{xlim}{Two increasing finite x-axis limits. The default includes estimates
and the widest contours, with padding.}
\item{ylim}{Standard-error limits as \code{c(0, maximum)}. The axis is reversed
so zero appears at the top. Defaults to 1.08 times the largest standard error.}
\item{legend}{Logical master switch. \code{FALSE} suppresses both legends.}
\item{legend.labels}{One model label per object. Defaults to model numbers and,
when available, estimator names.}
\item{xlab}{Optional x-axis label.}
\item{main}{Optional plot title.}
\item{model.legend}{Logical; show model labels when \code{legend = TRUE}.}
\item{legend.position}{Position keyword passed to \code{graphics::legend()}
for the model legend.}
\item{contour.position}{Position keyword for the contour legend.}
\item{legend.cex}{Positive size multiplier for both legends.}
\item{contours}{Logical; draw contour boundary lines.}
\item{show.refline}{Logical; draw the vertical reference lines. Turning these
lines off does not change the contour centers.}
\item{grid}{Logical; draw horizontal grid lines.}
\item{shade.alpha}{Opacity in [0, 1], either one value or one per sorted level.
\code{NULL} uses decreasing opacity from the inner to outer bands.}
\item{point.alpha}{Point opacity in [0, 1].}
\item{contour.lty}{Line types recycled across sorted contour levels.}
\item{contour.lwd}{Positive contour line width.}
\item{label}{\code{FALSE} suppresses study labels; \code{TRUE} labels every
point; \code{"out"} labels points outside the largest pseudo confidence region,
using original estimates before jitter.}
\item{label.cex}{Positive study-label size multiplier.}
\item{...}{Additional graphical arguments passed to the initial base plot,
such as \code{cex.axis}, \code{cex.lab}, \code{axes}, and \code{font.lab}.}
}
\details{
At standard error \eqn{s}, contours for level \eqn{C} are
\deqn{r \pm \Phi^{-1}((1+C)/2)s,}
where \eqn{r} is the selected reference estimate. These are normal-based pseudo
confidence regions, not confidence intervals for the pooled estimate. They
ignore uncertainty in the reference estimate and do not include between-study
heterogeneity in their widths. They do not inherit the study-level confidence
interval distribution or alpha setting.

For levels 0.90, 0.95, and 0.99, critical values are approximately 1.645, 1.960,
and 2.576. With \code{refline = 0}, these contours correspond to two-sided
normal-test p-value thresholds of 0.10, 0.05, and 0.01. With pooled centers,
they indicate deviations from the pooled estimates instead.

Shaded bands within each funnel do not overlap. Shading from two funnels
may overlap and blend colors. All study estimates must be finite and sampling
standard errors positive and finite.

Jitter changes only plotted x-coordinates. It leaves original estimates,
standard errors, contour geometry, and pooled centers unchanged. Corresponding
studies receive the same horizontal offset in both models. Offsets are stable
under row reordering and do not change the random-number generator state.

Funnel asymmetry is not sufficient evidence of publication bias. Dependence
among multiple effect estimates from the same study should be considered when
interpreting the plot or conducting associated statistical tests.
}
\value{
Invisibly returns a list containing \code{data} (one data frame per object,
with \code{study}, \code{estimate}, \code{se}, and \code{estimate.plot}),
\code{refline}, \code{level}, \code{xlim}, and \code{ylim}.
}
\seealso{\code{\link{meta.robust.smd}}, \code{metafor::funnel()}}
\examples{
studies <- readRDS(system.file(
  "extdata", "table3-studies.rds", package = "hcttest"
))
a <- meta.robust.smd(studies, study = "study", method = "constant",
  ss.adjust = TRUE, ci = "ct", pool = TRUE, pool.method = "REML")
b <- meta.robust.smd(studies, study = "study", method = "hc2",
  ss.adjust = TRUE, ci = "ct", pool = TRUE, pool.method = "REML")
funnel(a, b, legend.labels = c("Constant", "HC2"))
funnel(a, b, refline = 0, jitter = TRUE, legend = FALSE)
funnel(b, level = 0.95, shade = FALSE, legend = FALSE)
}
\references{
Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., and Rushton, L. (2008).
Contour-enhanced meta-analysis funnel plots help distinguish publication bias
from other causes of asymmetry. \emph{Journal of Clinical Epidemiology},
\bold{61}(10), 991--996. doi:10.1016/j.jclinepi.2007.11.010.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package.
\emph{Journal of Statistical Software}, \bold{36}(3), 1--48.
doi:10.18637/jss.v036.i03.
}
\keyword{hplot}
