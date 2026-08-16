# Power for One-Way ANOVA Under the Falconer Model

Calculates power for a one-way ANOVA comparing the means of the three
SNP genotype groups for one continuous quantitative trait.

## Usage

``` r
qtl_anova_power(
  N,
  alpha,
  qtl_var,
  tau,
  pd,
  count_method = c("rounded", "expected"),
  verbose = TRUE
)
```

## Arguments

- N:

  Integer total sample size greater than 3.

- alpha:

  Significance level in \\(0,1)\\.

- qtl_var:

  Variance in the standardized trait attributable to the QTL. Must be in
  \\(0,1)\\.

- tau:

  Finite dominance-to-additivity ratio \\\delta/a\\.

- pd:

  Frequency of the increaser allele. Must be in \\(0,1)\\.

- count_method:

  Genotype count method. `"rounded"` rounds expected genotype counts to
  the nearest integer before computing the ANOVA quantities.
  `"expected"` uses expected genotype counts directly without rounding.

- verbose:

  Logical. If `TRUE`, prints a concise summary.

## Value

An object of class `"qtl_anova_power"` containing power, sample size,
degrees of freedom, critical value, non-centrality parameter, genotype
counts, genotype means, residual variance, and the complete Falconer
parameter object.

## Details

The three genotype groups have the means returned by
[`qtl_falconer_parameters()`](https://akilanthony.github.io/pawh/reference/qtl_falconer_parameters.md)
and a common residual variance \\1-V\_{QTL}\\. Genotype group counts are
either rounded Hardy–Weinberg expectations or unrounded expected counts
according to `count_method`. The between-group effect gives the
noncentrality parameter in textbook Eq. 6.8, and power is the upper tail
of a noncentral F distribution with numerator degrees of freedom 2 and
denominator degrees of freedom \\N-3\\. Genotype counts follow Eq. 6.7
(Chapter 6, Section 6.1.1, p. 328).

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer, Chapter 6, Section 6.1.1, Eqs.
6.7–6.8, p. 328.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`qtl_anova_mssn`](https://akilanthony.github.io/pawh/reference/qtl_anova_mssn.md)
and
[`qtl_falconer_parameters`](https://akilanthony.github.io/pawh/reference/qtl_falconer_parameters.md).

## Examples

``` r
qtl_anova_power(
  N = 996, alpha = 0.0001, qtl_var = 0.025,
  tau = 0.5, pd = 0.15, verbose = FALSE
)
```
