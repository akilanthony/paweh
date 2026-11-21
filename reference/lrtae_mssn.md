# Minimum Total Sample Size for an LRTae Design

Finds the minimum integer total sample size required to attain a
requested prospective LRTae power while preserving a requested
allocation among the four measurement and validation tiers as closely as
integer counts permit. This is a prospective-design calculation and is
not an observed-data test.

## Usage

``` r
lrtae_mssn(
  true_control_genotype_frequencies,
  true_case_genotype_frequencies,
  true_sampling_proportions,
  phenotype_misclassification_matrix,
  genotype_misclassification_matrix,
  target_power,
  alpha,
  prop_fallible,
  prop_pheno_validated,
  prop_geno_validated,
  prop_both_validated
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

- target_power:

  Requested power, strictly greater than `alpha` and less than one.

- alpha:

  Significance level, strictly between zero and one.

- prop_fallible:

  Proportion allocated to fallible phenotype and genotype measurements
  only.

- prop_pheno_validated:

  Proportion allocated to phenotype validation with fallible genotype
  measurement.

- prop_geno_validated:

  Proportion allocated to genotype validation with fallible phenotype
  measurement.

- prop_both_validated:

  Proportion allocated to validation of both phenotype and genotype.

## Value

An object of class `lrtae_mssn`, containing `MSSN_total`,
`target_power`, `achieved_power`, `achieved_ncp`, `df`, `alpha`,
`critical_value`, requested allocation proportions, and the four
allocated integer tier counts.

## Details

The four allocation proportions must be nonnegative and sum to one. For
every candidate total `N`, integer counts are assigned by the
largest-remainder (Hamilton) rule: floor `N * proportion` for each tier,
then distribute the remaining individuals in descending order of
fractional remainder. Exact ties are resolved in the fixed order:
fallible, phenotype validated, genotype validated, both validated. Thus
allocated counts always sum exactly to `N`.

The search doubles an upper bound until it attains `target_power`,
applies integer binary search within the bracket, and explicitly
verifies the neighboring design at `N - 1`. Power at each candidate is
evaluated through
[`lrtae_power()`](https://akilanthony.github.io/paweh/reference/lrtae_power.md),
which uses the validated corrected
[`lrtae_ncp()`](https://akilanthony.github.io/paweh/reference/lrtae_ncp.md)
engine, `k - 1` degrees of freedom, and a noncentral chi-square
distribution.

## Examples

``` r
lrtae_mssn(
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
  target_power = 0.80,
  alpha = 0.05,
  prop_fallible = 0.7,
  prop_pheno_validated = 0.1,
  prop_geno_validated = 0.1,
  prop_both_validated = 0.1
)
#> $MSSN_total
#> [1] 331
#> 
#> $target_power
#> [1] 0.8
#> 
#> $achieved_power
#> [1] 0.8008275
#> 
#> $achieved_ncp
#> [1] 9.653685
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
#> $prop_fallible
#> [1] 0.7
#> 
#> $prop_pheno_validated
#> [1] 0.1
#> 
#> $prop_geno_validated
#> [1] 0.1
#> 
#> $prop_both_validated
#> [1] 0.1
#> 
#> $n_fallible
#> [1] 232
#> 
#> $n_pheno_validated
#> [1] 33
#> 
#> $n_geno_validated
#> [1] 33
#> 
#> $n_both_validated
#> [1] 33
#> 
#> attr(,"class")
#> [1] "lrtae_mssn"
```
