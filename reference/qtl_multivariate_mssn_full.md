# Multivariate Falconer Minimum Sample Size

Finds the minimum integer total sample size for Pillai MANOVA or the
minimum selected case-control sample for a joint-threshold genotype
chi-square test.

## Usage

``` r
qtl_multivariate_mssn_full(
  power,
  alpha,
  qtl_var,
  tau,
  pd,
  cor_matrix,
  test = c("pillai", "threshold_chisq"),
  x_upper = NULL,
  x_lower = NULL,
  k = 1,
  verbose = TRUE
)
```

## Arguments

- power:

  Target power in (0, 1).

- alpha:

  Significance level in (0, 1).

- qtl_var:

  Numeric vector of phenotype-specific QTL variances in (0, 1).

- tau:

  Numeric vector of phenotype-specific dominance/additivity ratios.

- pd:

  Shared increaser-allele frequency in (0, 1).

- cor_matrix:

  Positive-definite phenotype correlation matrix. Its order must match
  `qtl_var` and `tau`.

- test:

  Either `"pillai"` or `"threshold_chisq"`.

- x_upper, x_lower:

  Vectors of upper- and lower-tail percentages. A case satisfies every
  upper threshold (joint AND); a control satisfies every lower threshold
  (joint AND).

- k:

  Control/case ratio for the threshold chi-square design.

- verbose:

  Logical; whether to print a polished summary.

## Value

An object of class `"qtl_multivariate_mssn_full"`. Both modes include
`test`, target and achieved power, `alpha`, and the full `falconer`
model. The Pillai result includes integer total `N`,
`historical_fractional_mssn`, genotype counts, NCP, degrees of freedom,
critical value, and all `pillai` intermediates. The threshold result
includes selected case/control/total MSSNs, target and achieved NCPs,
internal `S`, thresholds and integration diagnostics, expected genotype
counts, sparse-cell diagnostics, and separately labelled expected
source-population screening counts.

## Details

Pillai MSSN is found by an integer search because its denominator
degrees of freedom depend on sample size. For threshold selection, the
statistical MSSN is kept separate from the expected population screening
burden. The Pillai branch recomputes its noncentral F distribution at
each candidate total `N`; the returned integer is the first design
attaining the target power. The historical fractional root is retained
only for comparison. The threshold branch uses direct
multivariate-normal rectangle integration from Eqs. 6.11–6.14 and
genotype chi-square Eq. 1.22; it solves for selected cases and controls
and reports source-population screening expectations separately. The two
branches therefore return different sample-size units.

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer, Chapter 6, Section 6.2, Eqs.
6.9–6.14, pp. 332–333; validation in Section 6.2.4, pp. 336–339.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

Gordon, D., Londono, D., Patel, P., Kim, W., Finch, S. J., & Heiman, G.
A. (2017). An analytic solution to computation of power and sample size
for genetic association studies under a pleiotropic mode of inheritance.
*Human Heredity*, 81(4), 194–209.
[doi:10.1159/000457135](https://doi.org/10.1159/000457135) .

Pillai, K. C. S. (1955). Some new test criteria in multivariate
analysis. *Annals of Mathematical Statistics*, 26(1), 117–121.
[doi:10.1214/aoms/1177728599](https://doi.org/10.1214/aoms/1177728599) .

Genz, A., & Bretz, F. (2009). *Computation of Multivariate Normal and t
Probabilities*. Springer.
[doi:10.1007/978-3-642-01689-9](https://doi.org/10.1007/978-3-642-01689-9)
.

## See also

[`qtl_multivariate_power_full`](https://akilanthony.github.io/pawh/reference/qtl_multivariate_power_full.md),
[`qtl_falconer_parameters`](https://akilanthony.github.io/pawh/reference/qtl_falconer_parameters.md),
and
[`qtl_threshold_chisq_mssn`](https://akilanthony.github.io/pawh/reference/qtl_threshold_chisq_mssn.md).

## Examples

``` r
cor_matrix <- matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
qtl_multivariate_mssn_full(
  power = 0.95, alpha = 5e-8, qtl_var = c(0.01, 0.005),
  tau = c(0, 0.5), pd = 0.25, cor_matrix = cor_matrix,
  test = "pillai", verbose = FALSE
)
```
