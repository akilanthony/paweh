# Specialized Case-Control Locus-Heterogeneity Functions

Convenience functions for case-control association calculations under
locus heterogeneity. These focused functions compute genotype chi-square
or genotype trend calculations one test at a time. The all-in-one
case-control functions remain the recommended interface when multiple
tests or modifiers are needed in one call.

## Usage

``` r
cc_chisq_mssn_locus_heterogeneity(
  power,
  alpha,
  g_case_assoc,
  g_ctrl,
  pi,
  k = 1,
  verbose = TRUE
)

cc_chisq_power_locus_heterogeneity(
  N_case,
  alpha,
  g_case_assoc,
  g_ctrl,
  pi,
  k = 1,
  verbose = TRUE
)

cc_trend_mssn_locus_heterogeneity(
  power,
  alpha,
  g_case_assoc,
  g_ctrl,
  pi,
  k = 1,
  w = c(0, 1, 2),
  verbose = TRUE
)

cc_trend_power_locus_heterogeneity(
  N_case,
  alpha,
  g_case_assoc,
  g_ctrl,
  pi,
  k = 1,
  w = c(0, 1, 2),
  verbose = TRUE
)
```

## Arguments

- power:

  Numeric in \\(0,1)\\. Desired target power for MSSN functions.

- alpha:

  Numeric in \\(0,1)\\. Significance level.

- g_case_assoc:

  Numeric vector of length 3. Case genotype frequencies under the
  associated-locus model, ordered as `c(g0, g1, g2)` and summing to 1.

- g_ctrl:

  Numeric vector of length 3. Control genotype frequencies, ordered as
  `c(g0, g1, g2)` and summing to 1.

- pi:

  Numeric in \\\[0,1\]\\. Locus-homogeneity fraction. The adjusted case
  genotype frequencies are \\\pi g\_{case,assoc} + (1 - \pi)
  g\_{ctrl}\\.

- k:

  Numeric \\\> 0\\. Control-to-case ratio \\N\_{ctrl}/N\_{case}\\.

- verbose:

  Logical. If `TRUE`, prints a formatted summary.

- N_case:

  Numeric \\\> 0\\. Number of cases for power functions.

- w:

  Numeric vector of length 3. Trend-test genotype scores. The three
  weights cannot all be equal.

## Value

A list with class matching the function name. MSSN functions include
target non-centrality parameter, sample sizes, internal `S` values, and
adjusted frequencies. Power functions include sample sizes,
non-centrality parameter, power, internal `S` values, and adjusted
frequencies. Trend functions additionally return the trend numerator and
denominator.

## Details

The associated case distribution is mixed with the control distribution
as \\\pi g\_{case,assoc} + (1-\pi)g\_{ctrl}\\, following textbook Eq.
2.16 (p. 88). Thus `pi` is the homogeneous/associated fraction and
\\1-\pi\\ is the heterogeneous fraction. The adjusted probabilities are
then passed to either the genotype chi-square NCP in Eq. 1.22 (p. 26) or
the trend-test NCP in Eq. 1.24 (p. 27). The trend construction under
locus heterogeneity is given more specifically by Eqs. 5.29a–b (pp.
287–288).

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

Armitage, P. (1955). Tests for linear trends in proportions and
frequencies. *Biometrics*, 11(3), 375–386.
[doi:10.2307/3001775](https://doi.org/10.2307/3001775) .

Slager, S. L., & Schaid, D. J. (2001). Case-control studies of genetic
markers: power and sample size approximations for Armitage's test for
trend. *Human Heredity*, 52(3), 149–153.
[doi:10.1159/000053370](https://doi.org/10.1159/000053370) .

## See also

[`cc_power`](https://akilanthony.github.io/paweh/reference/cc_power.md),
[`cc_mssn`](https://akilanthony.github.io/paweh/reference/cc_mssn.md),
and
[`case_control_genotype_misclassification`](https://akilanthony.github.io/paweh/reference/case_control_genotype_misclassification.md).

## Examples

``` r
cc_chisq_mssn_locus_heterogeneity(
  power = 0.8, alpha = 0.05,
  g_case_assoc = c(0.25, 0.50, 0.25),
  g_ctrl = c(0.36, 0.48, 0.16),
  pi = 0.8,
  verbose = FALSE
)

cc_chisq_power_locus_heterogeneity(
  N_case = 500, alpha = 0.05,
  g_case_assoc = c(0.25, 0.50, 0.25),
  g_ctrl = c(0.36, 0.48, 0.16), pi = 0.8, verbose = FALSE
)

cc_trend_mssn_locus_heterogeneity(
  power = 0.8, alpha = 0.05,
  g_case_assoc = c(0.25, 0.50, 0.25),
  g_ctrl = c(0.36, 0.48, 0.16), pi = 0.8, verbose = FALSE
)

cc_trend_power_locus_heterogeneity(
  N_case = 500, alpha = 0.05,
  g_case_assoc = c(0.25, 0.50, 0.25),
  g_ctrl = c(0.36, 0.48, 0.16),
  pi = 0.8,
  verbose = FALSE
)
```
