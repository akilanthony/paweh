# Falconer Parameters for a Single Quantitative Trait

Calculates genotype-specific normal-mixture parameters for a
standardized continuous trait under the single-locus Falconer model.

## Usage

``` r
qtl_falconer_parameters(qtl_var, tau, pd, verbose = TRUE)
```

## Arguments

- qtl_var:

  Variance in the standardized trait attributable to the QTL. Must be in
  \\(0,1)\\.

- tau:

  Finite dominance-to-additivity ratio \\\delta/a\\.

- pd:

  Frequency of the increaser allele. Must be in \\(0,1)\\.

- verbose:

  Logical. If `TRUE`, prints a formatted parameter summary. Set to
  `FALSE` to suppress all console output.

## Value

An object of class `"qtl_falconer_parameters"`. The list contains
`qtl_var`, `tau`, `pd`, `p_plus`, additive effect `a`, dominance effect
`delta`, centering constant `m`, genotype means `mu`, mixing proportions
`pi`, `residual_variance`, `residual_sd`, `weighted_mean`, and
`total_variance`.

## Details

Genotypes are ordered as zero, one, and two copies of the increaser
allele. Their Hardy-Weinberg mixing proportions are \\((1-p_d)^2,
2p_d(1-p_d), p_d^2)\\. The population-centering constant makes the
genotype-weighted trait mean zero, and the residual variance is \\1 -
V\_{QTL}\\.

This is a classical, Falconer-style quantitative-genetic
parameterization as formulated for this study-design framework by
Gordon, Finch, and Kim (2020), Chapter 6, Section 6.1. The population
mixture is Eq. 6.1 (p. 324); the additive effect, dominance effect, and
centering formulas on p. 325 are unnumbered. The function name does not
assert that every exact implemented expression is a numbered equation in
Falconer and Mackay.

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer, Chapter 6, Section 6.1, pp. 324–325.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

Fisher, R. A. (1918). The correlation between relatives on the
supposition of Mendelian inheritance. *Transactions of the Royal Society
of Edinburgh*, 52, 399–433.
[doi:10.1017/S0080456800012163](https://doi.org/10.1017/S0080456800012163)
.

Lynch, M., & Walsh, B. (1998). *Genetics and Analysis of Quantitative
Traits*. Sinauer Associates.

Falconer, D. S., & Mackay, T. F. C. (1996). *Introduction to
Quantitative Genetics*, 4th ed. Longman.

## See also

[`qtl_falconer_threshold_parameters`](https://akilanthony.github.io/paweh/reference/qtl_falconer_threshold_parameters.md),
[`qtl_anova_power`](https://akilanthony.github.io/paweh/reference/qtl_anova_power.md),
and
[`qtl_anova_mssn`](https://akilanthony.github.io/paweh/reference/qtl_anova_mssn.md).

## Examples

``` r
qtl_falconer_parameters(qtl_var = 0.025, tau = 0.5, pd = 0.15)
#> 
#> ==========================================================================
#> Falconer Quantitative Trait Parameters
#> ==========================================================================
#> 
#> Input Parameters
#> QTL variance:                                        0.0250
#> Dominance/Additivity ratio (tau):                    0.5000
#> Increaser allele frequency:                          0.1500
#> --------------------------------------------------------------------------
#> 
#> Derived Parameters
#> Additive effect (a):                                 0.2280
#> Dominance effect (delta):                            0.1140
#> Centering mean (m):                                  0.1305
#> Residual variance:                                   0.9750
#> Residual SD:                                         0.9874
#> --------------------------------------------------------------------------
#> 
#> Genotype-Specific Quantities
#>                                  Genotype 0     Genotype 1     Genotype 2
#> Mixing proportion:                   0.7225         0.2550         0.0225
#> Trait mean:                         -0.0975         0.2445         0.3585
#> --------------------------------------------------------------------------
#> 
#> Validation
#> Weighted population mean:                          0.000000
#> Total variance:                                    1.000000
#> --------------------------------------------------------------------------
```
