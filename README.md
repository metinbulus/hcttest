# hcttest

Two-group mean differences, robust standardized differences, bootstrap inference,
data generation, and power calculations.

## Public functions

- `hct.test()` analyzes raw data or group summaries.
- `robust.smd()` computes a standardized difference from a fitted test.
- `boot.hct.test()` bootstraps raw and standardized differences.
- `gendat.hct.test()` generates two-group data.
- `meta.robust.smd()` computes effect-size rows from a study summary table.
- `multi.robust.smd()` analyzes multiple outcome columns.
- `power.hct.test()` computes power from an object or noncentrality parameter.

All other functions, including Lambda-prime helpers, are internal.

## Observation and mixture weights

`weights` (plural) supplies per-observation weights to hct.test(), boot.hct.test(),
and multi.robust.smd(). The raw x/y interface accepts a concatenated vector or
list(x = ..., y = ...); the formula and bootstrap interfaces accept a data column
expression or vector. multi.robust.smd() accepts a column name or vector.
Missing values are removed jointly when na.rm = TRUE; zero-weight rows are
excluded. Weights must be finite and nonnegative.

`mixnorm.weight` in summary-distribution calculations and `wt` in data generation are
mixture-component proportions, not participant weights. meta.robust.smd() uses
study summaries and returns effects and their variances; it pools independent rows by default using metafor (REML and Knapp-Hartung). Set pool = FALSE for an effect table only. robust.smd() inherits weighting from its fitted object.

Analytical observation weighting supports untrimmed Student and HC0-HC5/HC4m
with regression df. Weighted analytical trimming, Welch-Satterthwaite df,
sd1/sd2 tests, and nonnormal parametric bootstrap remain unavailable.
Weighted trimming is supported in boot.hct.test() by total weight mass with
partial boundary weights. Existing weighted SMD SE corrections are approximations.

## g-and-h generation

```r
gendat.hct.test(0, 0.5, 1, 2, 30, 40, distribution = "norm",
                norm.transform.g = c(0, 0.3), norm.tranform.h = 0.1)
```

Shape parameters are optional scalars or group-specific pairs and are accepted
only for a normal base distribution. The existing sample-moment rescaling
preserves the normal draws' realized mean and SD, rather than forcing exact
population moments. Zero shapes preserve the original draws.

## Normal p values

For ci = "z", hct.test() reports a normal-approximation p value for the raw
difference, and robust.smd() reports one using d/se. Batch output distinguishes
raw p from p.smd. The existing observed-direction convention is used for one-sided
inference.

## Examples and validation

The study-tables vignette reproduces the supplied 29-row Table 3 dataset and
separately shows recalculated results; its reported weights differ from the
current formulas. The validation vignette contains short, self-contained checks,
one at a time, that execute when knitted. It calls no validation/test scripts.
Only the optional HC reference comparisons require sandwich; other examples use
base R, hcttest, and its existing boot dependency. knitr/rmarkdown build the
vignettes. R execution and rendering have not yet been verified in this environment.

## Reference

Bulus, M. (2026). Robust standard errors and confidence intervals for standardized
mean difference. Research Synthesis Methods. Advance online publication.

Use `citation("hcttest")`. Observation-weight extensions require validation
separate from the paper's original derivations and simulations.

## Updating documentation with roxygen2

Function documentation now lives in roxygen comments in R/. Run from this
package's source directory:

```r
devtools::document(roclets = c("rd", "collate", "namespace"))
devtools::check()
```

Edit those source comments rather than man/ or NAMESPACE. The S3 methods are
registered through their roxygen tags and remain internal. Encoding is UTF-8.
The shipped help/namespace files are marked for roxygen ownership and can be
regenerated. Roxygen generation has not been executed in the preparation
environment because R is unavailable.

Extract this archive into a **fresh folder**, rather than merging it into an
older source tree. Older drafts contained tests/, inst/validation/, and
man/lambdap.Rd; they are deliberately absent now. Retaining those old files
causes checks to run obsolete APIs and internal-function examples. The dataset
uses RDS version 2 to avoid the automatic R >= 3.5.0 dependency warning from a
version-3 serialized object. The unrelated quarto/TMPDIR message in the supplied
Windows log is outside the package source changes.

If updating an existing source folder, remove these two obsolete help files
before regenerating documentation (run from the package source directory):

```r
unlink(c("man/lambdap.Rd", "man/batch.robust.smd.Rd"))
devtools::document(roclets = c("rd", "collate", "namespace"))
devtools::check()
```

The shared help page for meta.robust.smd() and multi.robust.smd() is now
man/meta.robust.smd.Rd. Lambda-prime helpers have no public help page.

## Plotting

All five analysis outputs support plot():

```r
fit <- hct.test(m1 = 0.4, m2 = 0, sd1 = 1, sd2 = 1.5,
                n1 = 60, n2 = 60, ci = "z", verbose = FALSE)
effect <- robust.smd(fit, ci = "z", verbose = FALSE)
plot(effect, margins = c(-0.2, 0.2))
plot(effect, pretty = TRUE, margins = c(-0.2, 0.2))
```

The default range is identical for the ordinary and shaded styles, symmetric
about zero, with a full interval width of padding on either side. Override it
with xlim. Study/outcome names appear on the left; label.padding controls outer
space before those labels.

For individual and meta-analysis objects, the ordinary style uses stored confidence limits. Shaded 95%/99% bands follow
the fitted z, ct or nct method. Analyses are rerun at alpha = 0.01 using their original retained inputs and
arguments, so the original functions perform the noncentral-t transformation.
Bootstrap plotting also reruns boot.hct.test, reusing its starting RNG state
and restoring the current RNG afterwards. It can take as long as fitting. Shading is solid inside 95%, with thin
vertical boundaries, then fades to transparent at 99%.
See help("hcttest-plot", package = "hcttest") for all settings.

Two objects can be overlaid directly, without par(new = TRUE):

```r
plot(first, second, pretty = TRUE, col = c("navy", "darkorange3"))
```

The second object's rows are reordered to match the first by unique study or
outcome names. Both styles honor col; default colors are dark navy and dark
orange. A labels list can provide explicit matching identifiers. Both objects
share axes, with the first estimate slightly above each label and the second
slightly below, separated by a small gap. Recreate analysis objects after
installing this update so the fitted CI settings are retained.

Retained inputs can be inspected directly:

```r
attr(meta_result, "analysis")$args$data
attr(multi_result, "analysis")$args$data
attr(bootstrap_result, "analysis")$args$data
attr(smd_result, "analysis")$args$x  # original fitted hct.test object
```

### Pooling and equivalence

meta.robust.smd() accepts pool.method (default "REML") and pool.test
(default "knha"). attr(x, "pooled") contains raw and SMD pooled summaries,
including tau2, I2 and Q; attr(x, "meta.models") contains the metafor fits.
The weight column contains unnormalized SMD random-effects weights.
A different tau2 estimator can change the pooled point estimate. Pooling
assumes independent rows. Meta-analysis plot reference lines mark each
pooled estimate in its specified color; reference = NULL suppresses them.

multi.robust.smd() accepts margins = c(lower, upper), equivalence.effect
("smd" or "raw"). SMD equivalence uses the specified ci (ct, z or nct),
with 100*(1-2*alpha) percent intervals and matching one-sided bound tests.
The nct choice reuses robust.smd's noncentral-t inversion and scaling.
Raw-difference equivalence retains central-t inference. With ci = FALSE,
SMD equivalence inference is undefined. p.max is max(p.lower, p.upper);
equivalence requires p.max < alpha.

Install metafor with the other package dependencies. The equivalence and pooling checks run in R; vignette rebuilding requires the
suggested vignette packages.

## Batch plots

multi.robust.smd returns individual outcome estimates and equivalence results.
Its plots retain a zero reference and have no pooled summary.
Meta-analysis plots display metafor pooled summaries below the tick labels,
using each fitted model's confidence interval method and matching colors.
Both plot styles use the same axis limits.
