# Power for a Threshold-Selected Genotype Chi-Square Test

Calculates power for the 2-df genotype chi-square test after a
quantitative trait is converted into upper-tail cases and lower-tail
controls.

## Usage

``` r
qtl_threshold_chisq_power(
  N_case,
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

- N_case:

  Positive number of selected cases.

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

An object of class `"qtl_threshold_chisq_power"` containing selected
sample sizes, power, non-centrality parameter, internal effect component
`S`, thresholds, penetrances, prevalences, conditional genotype
frequencies, and Falconer parameters.

## Details

Genotype-specific normal means and the common residual variance follow
the mixture in Eq. 6.1. Upper-tail cases and lower-tail controls are
constructed with Eqs. 6.3–6.5 (pp. 324–326); the middle of the trait
distribution is excluded. The resulting conditional genotype
probabilities are compared by the two-degree-of-freedom genotype
chi-square test using Eq. 1.22 (p. 26). `N_case` and `k * N_case` are
selected case and control counts, not numbers screened from the source
population.

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer, Eqs. 1.22, 6.1, and 6.3–6.5, pp. 26
and 324–326.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`qtl_threshold_chisq_mssn`](https://akilanthony.github.io/pawh/reference/qtl_threshold_chisq_mssn.md)
and
[`qtl_falconer_threshold_parameters`](https://akilanthony.github.io/pawh/reference/qtl_falconer_threshold_parameters.md).

## Examples

``` r
qtl_threshold_chisq_power(
  N_case = 126, alpha = 0.0001,
  qtl_var = 0.025, tau = 0.5, pd = 0.15,
  x_upper = 5, x_lower = 5, verbose = FALSE
)
```
