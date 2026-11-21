# Prospective Power for an LRTae Design

Computes prospective power for the likelihood-ratio test allowing for
error (LRTae). The function uses
[`lrtae_ncp()`](https://akilanthony.github.io/paweh/reference/lrtae_ncp.md)
and the corrected LRTae information matrix; it does not perform an
observed-data hypothesis test.

## Usage

``` r
lrtae_power(
  true_control_genotype_frequencies,
  true_case_genotype_frequencies,
  true_sampling_proportions,
  phenotype_misclassification_matrix,
  genotype_misclassification_matrix,
  n_fallible,
  n_pheno_validated,
  n_geno_validated,
  n_both_validated,
  alpha
)
```

## Arguments

- true_control_genotype_frequencies:

  Numeric vector of true genotype frequencies among controls.

- true_case_genotype_frequencies:

  Numeric vector of true genotype frequencies among cases, in the same
  category order as `true_control_genotype_frequencies`.

- true_sampling_proportions:

  Numeric vector of length two containing the true control and case
  sampling proportions, in that order.

- phenotype_misclassification_matrix:

  Numeric two-by-two matrix whose rows are true phenotype categories
  (control, case) and whose columns are observed phenotype categories
  (control, case).

- genotype_misclassification_matrix:

  Numeric square matrix whose rows are true genotype categories and
  whose columns are observed genotype categories, in the
  genotype-frequency order.

- n_fallible:

  Number sampled with fallible phenotype and genotype measurements only.

- n_pheno_validated:

  Number whose phenotype is validated while genotype is measured
  fallibly.

- n_geno_validated:

  Number whose genotype is validated while phenotype is measured
  fallibly.

- n_both_validated:

  Number whose phenotype and genotype are both validated.

- alpha:

  Significance level, strictly between zero and one.

## Value

An object of class `lrtae_power`, containing `power`, `ncp`, `df`,
`alpha`, `critical_value`, the four sampling-tier counts, and `N_total`.

## Details

For `k` genotype categories, the asymptotic reference distribution has
`k - 1` degrees of freedom. The rejection threshold is the upper `alpha`
quantile of the central chi-square distribution, and power is the
corresponding upper-tail probability from a noncentral chi-square
distribution with the NCP returned by
[`lrtae_ncp()`](https://akilanthony.github.io/paweh/reference/lrtae_ncp.md).

## Examples

``` r
lrtae_power(
  true_control_genotype_frequencies = c(0.81, 0.18, 0.01),
  true_case_genotype_frequencies = c(0.64, 0.32, 0.04),
  true_sampling_proportions = c(0.5, 0.5),
  phenotype_misclassification_matrix = matrix(
    c(0.95, 0.05, 0.05, 0.95), 2, byrow = TRUE
  ),
  genotype_misclassification_matrix = matrix(
    c(0.950, 0.050, 0.000,
      0.025, 0.950, 0.025,
      0.000, 0.050, 0.950),
    3, byrow = TRUE
  ),
  n_fallible = 700,
  n_pheno_validated = 100,
  n_geno_validated = 100,
  n_both_validated = 100,
  alpha = 0.05
)
#> $power
#> [1] 0.998997
#> 
#> $ncp
#> [1] 29.17305
#> 
#> $df
#> [1] 2
#> 
#> $alpha
#> [1] 0.05
#> 
#> $critical_value
#> [1] 5.991465
#> 
#> $n_fallible
#> [1] 700
#> 
#> $n_pheno_validated
#> [1] 100
#> 
#> $n_geno_validated
#> [1] 100
#> 
#> $n_both_validated
#> [1] 100
#> 
#> $N_total
#> [1] 1000
#> 
#> attr(,"class")
#> [1] "lrtae_power"
```
