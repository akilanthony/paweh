# Minimum Sample Size for One-Way ANOVA Under the Falconer Model

Finds the smallest sufficient integer total sample size for the
three-group one-way ANOVA. Both denominator degrees of freedom and the
non-centrality parameter are recalculated for every candidate sample
size.

## Usage

``` r
qtl_anova_mssn(
  power,
  alpha,
  qtl_var,
  tau,
  pd,
  count_method = c("rounded", "expected"),
  multiple_of_three = TRUE,
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

- count_method:

  Genotype count method. `"rounded"` rounds expected genotype counts to
  the nearest integer before computing the ANOVA quantities.
  `"expected"` uses expected genotype counts directly without rounding.

- multiple_of_three:

  Logical. If `TRUE`, searches sample sizes in multiples of three. If
  `FALSE`, searches every integer sample size.

- verbose:

  Logical. If `TRUE`, prints a concise summary.

## Value

An object of class `"qtl_anova_mssn"` containing the smallest sufficient
`N`, achieved power, degrees of freedom, non-centrality parameter,
genotype counts, and Falconer parameters.

## Details

For every candidate integer `N`, the genotype counts, denominator
degrees of freedom, noncentrality parameter, critical F value, and
achieved power are recomputed. With `count_method = "rounded"`,
candidates for which any genotype group has zero members are excluded.
The first candidate attaining `power` is returned; setting
`multiple_of_three = TRUE` restricts the search grid rather than
changing the ANOVA formula. The model and NCP follow textbook Eqs.
6.7–6.8 (p. 328).

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer, Chapter 6, Section 6.1.1, Eqs.
6.7–6.8, p. 328.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`qtl_anova_power`](https://akilanthony.github.io/paweh/reference/qtl_anova_power.md)
and
[`qtl_falconer_parameters`](https://akilanthony.github.io/paweh/reference/qtl_falconer_parameters.md).

## Examples

``` r
qtl_anova_mssn(
  power = 0.8, alpha = 0.0001, qtl_var = 0.025,
  tau = 0.5, pd = 0.15, verbose = FALSE
)
```
