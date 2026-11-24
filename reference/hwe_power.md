# Plan power and sample size for a Pearson HWE test

Compute prospective power or the minimum whole-number sample size for
the ordinary, uncorrected Pearson Hardy-Weinberg equilibrium (HWE) test.
This is a noncentral chi-square planning approximation for independently
sampled diploid individuals at one biallelic autosomal locus.

## Usage

``` r
hwe_power(
  N,
  alpha = 0.05,
  genotype_probs = NULL,
  allele_frequency = NULL,
  misclassification_matrix = NULL
)

hwe_mssn(
  target_power = 0.8,
  alpha = 0.05,
  genotype_probs = NULL,
  allele_frequency = NULL,
  misclassification_matrix = NULL
)
```

## Arguments

- N:

  Number of analyzable individuals, a positive whole number at most
  `2^52`.

- alpha:

  Nominal significance level, strictly between 0 and 1.

- genotype_probs:

  Assumed observed probabilities for genotypes with 0, 1, and 2 copies
  of the same allele. Supply a finite, nonnegative vector summing to 1,
  or instead supply both `allele_frequency` and
  `misclassification_matrix`. Names, if given, must be `"0"`, `"1"`, and
  `"2"` in that order.

- allele_frequency:

  Frequency of the counted allele **before** error, in `[0, 1]`. It
  defines the true HWE probabilities.

- misclassification_matrix:

  Numeric 3 by 3 matrix with rows for true genotype and columns for
  observed genotype, both in 0, 1, 2 allele-copy order. Any row or
  column names must be `"0"`, `"1"`, `"2"` in that order. Entry
  `[i + 1, j + 1]` is the probability of calling `j` copies when there
  are truly `i` copies. Each row sums to 1. Multiplying the true
  probability row vector by this matrix gives the observed
  probabilities.

- target_power:

  Desired rejection probability, strictly between 0 and 1.

## Value

A list with the sample size `N`, achieved `power`, `alpha`, `df`,
`critical_value`, `lambda` (the noncentrality parameter),
`lambda_per_person`, `expected_counts`, `sparse_expected`,
`true_genotype_probs` (matrix scenario only), `observed_probs`, fitted
`allele_frequency`, `expected_probs`, and `F`. `hwe_mssn()` also returns
`target_power` and `previous_power` at the immediately smaller positive
size when available. An equilibrium alternative has power `alpha` at
every size in this approximation. For such an alternative and target
above `alpha`, `N` is `Inf` and achieved power and size-dependent
diagnostics are `NA`.

## Details

The fitted null uses the counted-allele frequency from the assumed
observed probabilities `g`. Write `h` for the HWE probabilities at that
frequency. The noncentrality parameter is
`lambda = N * sum((g - h)^2 / h) = N * F^2`, where `F = 1 - g[2] / h[2]`
uses the heterozygote probabilities. Power is the upper-tail probability
under this one-degree-of-freedom noncentral chi-square at the central
chi-square cutoff. No exact-test power or goodness of fit to a fixed
external frequency is calculated. The fitted distribution must contain
both alleles. Expected counts below five trigger a warning for finite
designs. They are a diagnostic for the approximation, rather than a
constraint in the sample-size search. Sample sizes above `2^52` are
outside the supported whole-number range.

## References

Rohlfs RV, Weir BS (2008). Distributions of Hardy-Weinberg equilibrium
test statistics. *Genetics* 180(3), 1609-1616.
[doi:10.1534/genetics.108.088005](https://doi.org/10.1534/genetics.108.088005)
.

## See also

[`hwe_test()`](https://akilanthony.github.io/paweh/reference/hwe_test.md)

## Examples

``` r
g <- c(0.35588, 0.48424, 0.15988)
hwe_power(1000, genotype_probs = g)$power
#> [1] 0.05591224
hwe_mssn(genotype_probs = g)$N
#> [1] 152612
```
