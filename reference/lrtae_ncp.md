# Noncentrality Parameter for Prospective LRTae Designs

Computes the analytic noncentrality parameter (NCP) for the
likelihood-ratio test allowing for error (LRTae) in a prospective
case-control design. The calculation incorporates fallible measurements
and validation subsamples for phenotype, genotype, or both. It does not
analyze observed study data.

## Usage

``` r
lrtae_ncp(
  true_control_genotype_frequencies,
  true_case_genotype_frequencies,
  true_sampling_proportions,
  phenotype_misclassification_matrix,
  genotype_misclassification_matrix,
  n_fallible,
  n_pheno_validated,
  n_geno_validated,
  n_both_validated
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

## Value

A single nonnegative numeric value: the LRTae NCP.

## References

Ji, F., Yang, Y., Haynes, C., Finch, S. J., & Gordon, D. (2005).
Computing asymptotic power and sample size for case-control genetic
association studies in the presence of phenotype and/or genotype
misclassification errors. *Statistical Applications in Genetics and
Molecular Biology*, 4(1), Article 37.
[doi:10.2202/1544-6115.1184](https://doi.org/10.2202/1544-6115.1184)

Barral, S., Haynes, C., Stone, M., & Gordon, D. (2006). LRTae: improving
statistical power for genetic association with case/control data when
phenotype and/or genotype misclassification errors are present. *BMC
Genetics*, 7, 24.
[doi:10.1186/1471-2156-7-24](https://doi.org/10.1186/1471-2156-7-24)

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)

## Examples

``` r
lrtae_ncp(
  true_control_genotype_frequencies = c(0.81, 0.18, 0.01),
  true_case_genotype_frequencies = c(0.64, 0.32, 0.04),
  true_sampling_proportions = c(0.5, 0.5),
  phenotype_misclassification_matrix = matrix(
    c(0.95, 0.05, 0.05, 0.95), nrow = 2, byrow = TRUE
  ),
  genotype_misclassification_matrix = matrix(
    c(0.950, 0.050, 0.000,
      0.025, 0.950, 0.025,
      0.000, 0.050, 0.950),
    nrow = 3, byrow = TRUE
  ),
  n_fallible = 700,
  n_pheno_validated = 100,
  n_geno_validated = 100,
  n_both_validated = 100
)
#> [1] 29.17305
```
