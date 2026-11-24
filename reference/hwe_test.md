# Test observed genotype counts for Hardy-Weinberg equilibrium

Test one biallelic autosomal locus using counts of diploid individuals.
The null is Hardy-Weinberg proportions with an unknown allele frequency.
Pearson inference estimates that frequency from the counts. Exact
inference conditions on the allele total. A large p-value does not
establish equilibrium.

## Usage

``` r
hwe_test(counts, method = c("exact", "chisq"))
```

## Arguments

- counts:

  Numeric vector of three nonnegative whole counts, in 0, 1, and 2
  copies of the same counted allele. If named, use `"0"`, `"1"`, and
  `"2"` in that order. The sum must be positive and no larger than
  `2^52`.

- method:

  One of `"exact"` (default) or `"chisq"` (uncorrected Pearson).

## Value

A list with `N`, `counts`, `observed_probs`, `allele_frequency`,
`expected_probs`, `expected_counts`, Pearson `statistic` and `df`, the
selected `method` and its `p_value`, the signed heterozygote deficit
`F = 1 - observed/expected`, and `sparse_expected` indicating any
expected count below five. Negative `F` indicates heterozygote excess.
The Pearson degrees of freedom describe `statistic`, even when the
selected p-value is exact.

## Details

The uncorrected Pearson statistic has an asymptotic chi-square reference
distribution with one degree of freedom: three cells minus one
total-count constraint and one estimated allele frequency. Expected
counts below five prompt a warning for Pearson inference. The two-sided
exact p-value conditions on the observed allele count and sums the
probabilities of tables at least as unlikely as the observed table,
including ties (Wigginton et al., 2005). Exact enumeration is limited to
one million conditional tables. At a monomorphic locus there is one
conditional table (exact p-value 1) and no regular Pearson test. The
function warns and returns `NA` for the Pearson statistic, degrees of
freedom, and `F`.

## References

Wigginton JE, Cutler DJ, Abecasis GR (2005). A note on exact tests of
Hardy-Weinberg equilibrium. *American Journal of Human Genetics* 76(5),
887-893. [doi:10.1086/429864](https://doi.org/10.1086/429864) .

## See also

[`hwe_power()`](https://akilanthony.github.io/paweh/reference/hwe_power.md),
[`hwe_mssn()`](https://akilanthony.github.io/paweh/reference/hwe_power.md)

## Examples

``` r
hwe_test(c(88, 10, 2))
#> $N
#> [1] 100
#> 
#> $counts
#>  0  1  2 
#> 88 10  2 
#> 
#> $observed_probs
#>    0    1    2 
#> 0.88 0.10 0.02 
#> 
#> $allele_frequency
#> [1] 0.07
#> 
#> $expected_probs
#>      0      1      2 
#> 0.8649 0.1302 0.0049 
#> 
#> $expected_counts
#>     0     1     2 
#> 86.49 13.02  0.49 
#> 
#> $statistic
#> [1] 5.380119
#> 
#> $df
#> [1] 1
#> 
#> $method
#> [1] "exact"
#> 
#> $p_value
#> [1] 0.06544427
#> 
#> $F
#> [1] 0.2319508
#> 
#> $sparse_expected
#> [1] TRUE
#> 
hwe_test(c(40, 30, 30), method = "chisq")
#> $N
#> [1] 100
#> 
#> $counts
#>  0  1  2 
#> 40 30 30 
#> 
#> $observed_probs
#>   0   1   2 
#> 0.4 0.3 0.3 
#> 
#> $allele_frequency
#> [1] 0.45
#> 
#> $expected_probs
#>      0      1      2 
#> 0.3025 0.4950 0.2025 
#> 
#> $expected_counts
#>     0     1     2 
#> 30.25 49.50 20.25 
#> 
#> $statistic
#> [1] 15.51882
#> 
#> $df
#> [1] 1
#> 
#> $method
#> [1] "chisq"
#> 
#> $p_value
#> [1] 8.16877e-05
#> 
#> $F
#> [1] 0.3939394
#> 
#> $sparse_expected
#> [1] FALSE
#> 
```
