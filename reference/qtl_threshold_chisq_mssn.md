# Minimum Sample Size for a Threshold-Selected Genotype Chi-Square Test

Calculates the minimum selected case and control sample sizes after
deriving their genotype distributions from a single-trait Falconer
threshold model.

## Usage

``` r
qtl_threshold_chisq_mssn(
  power,
  alpha,
  qtl_var,
  tau,
  pd,
  x_upper,
  x_lower,
  k = 1,
  verbose = TRUE
)
```

## Arguments

- power:

  Target power in \\(0,1)\\.

- alpha:

  Significance level in \\(0,1)\\.

- qtl_var:

  Variance in the standardized trait attributable to the QTL. Must be in
  \\(0,1)\\.

- tau:

  Finite dominance-to-additivity ratio \\\delta/a\\.

- pd:

  Frequency of the increaser allele. Must be in \\(0,1)\\.

- x_upper:

  Percentage selected from the upper population tail in percentile mode.

- x_lower:

  Percentage selected from the lower population tail in percentile mode.

- k:

  Positive control-to-case ratio \\N\_{control}/N\_{case}\\.

- verbose:

  Logical. If `TRUE`, prints a concise summary.

## Value

An object of class `"qtl_threshold_chisq_mssn"` containing selected MSSN
values, target NCP, internal `S`, thresholds, penetrances, prevalences,
conditional frequencies, Falconer parameters, and separately labelled
expected population screening counts.

## Details

The primary MSSN is the selected case-control sample size used by the
association test. The separately reported screening quantities are
expected population counts needed to obtain those selected samples under
simple population sampling; they are not statistical MSSN values.
Genotype-specific distributions and threshold selection follow Eqs. 6.1
and 6.3–6.5 (pp. 324–326). The downstream genotype-test NCP is Eq. 1.22
(p. 26). The function solves for selected cases and controls; expected
screening counts are reported separately and are not rounded MSSNs.

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer, Eqs. 1.22, 6.1, and 6.3–6.5, pp. 26
and 324–326.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`qtl_threshold_chisq_power`](https://akilanthony.github.io/pawh/reference/qtl_threshold_chisq_power.md)
and
[`qtl_falconer_threshold_parameters`](https://akilanthony.github.io/pawh/reference/qtl_falconer_threshold_parameters.md).

## Examples

``` r
qtl_threshold_chisq_mssn(
  power = 0.8, alpha = 0.0001,
  qtl_var = 0.025, tau = 0.5, pd = 0.15,
  x_upper = 5, x_lower = 5, verbose = FALSE
)
```
