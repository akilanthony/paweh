# Threshold-Selected Parameters Under the Falconer Model

Converts the genotype-specific Falconer trait distributions into
upper-tail case and lower-tail control penetrances and conditional
genotype frequencies.

## Usage

``` r
qtl_falconer_threshold_parameters(
  qtl_var,
  tau,
  pd,
  threshold_mode = c("percentile", "direct"),
  x_upper = NULL,
  x_lower = NULL,
  upper_threshold = NULL,
  lower_threshold = NULL,
  verbose = TRUE
)
```

## Arguments

- qtl_var:

  Variance in the standardized trait attributable to the QTL. Must be in
  \\(0,1)\\.

- tau:

  Finite dominance-to-additivity ratio \\\delta/a\\.

- pd:

  Frequency of the increaser allele. Must be in \\(0,1)\\.

- threshold_mode:

  Either `"percentile"` or `"direct"`.

- x_upper:

  Percentage selected from the upper population tail in percentile mode.

- x_lower:

  Percentage selected from the lower population tail in percentile mode.

- upper_threshold:

  Direct upper-tail threshold in direct mode.

- lower_threshold:

  Direct lower-tail threshold in direct mode.

- verbose:

  Logical. If `TRUE`, prints a formatted threshold-model summary. Set to
  `FALSE` to suppress all console output.

## Value

An object of class `"qtl_falconer_threshold_parameters"` containing
thresholds, genotype-specific penetrances, selected-population
prevalences, conditional genotype frequencies, and Falconer parameters.

## Details

Cases are selected above the upper threshold and controls below the
lower threshold. The middle portion is excluded, so the affected and
unaffected selection events are not complements. Percentile inputs refer
to standard normal population percentiles before conditioning on
genotype. The genotype-specific normal mixture follows Eq. 6.1 (pp.
324–325). Genotype-specific upper- and lower-tail probabilities, their
population prevalences, and conditional genotype frequencies follow Eqs.
6.3–6.5 (p. 326). Eq. 6.6 (p. 327) is a numerical worked example, not an
additional symbolic step. This is a classical/Falconer-style
parameterization as formulated for this design framework by Gordon,
Finch, and Kim (2020).

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer, Chapter 6, Section 6.1, Eqs. 6.3–6.5,
p. 326.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

Falconer, D. S., & Mackay, T. F. C. (1996). *Introduction to
Quantitative Genetics*, 4th ed. Longman.

## See also

[`qtl_falconer_parameters`](https://akilanthony.github.io/paweh/reference/qtl_falconer_parameters.md),
[`qtl_threshold_chisq_power`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_power.md),
and
[`qtl_threshold_chisq_mssn`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_mssn.md).

## Examples

``` r
qtl_falconer_threshold_parameters(
  qtl_var = 0.025, tau = 0.5, pd = 0.15,
  x_upper = 5, x_lower = 5
)
#> 
#> ==========================================================================
#> Falconer Threshold-Selected Trait Parameters
#> ==========================================================================
#> 
#> Input Parameters
#> QTL variance:                                        0.0250
#> Dominance/Additivity ratio (tau):                    0.5000
#> Increaser allele frequency:                          0.1500
#> --------------------------------------------------------------------------
#> 
#> Threshold Information
#> Threshold mode:                                  percentile
#> Upper percentile:                                     5.00%
#> Lower percentile:                                     5.00%
#> Upper threshold:                                    1.64485
#> Lower threshold:                                   -1.64485
#> --------------------------------------------------------------------------
#> 
#> Derived Falconer Parameters
#> QTL variance:                                        0.0250
#> Residual variance:                                   0.9750
#> Residual SD:                                         0.9874
#> Increaser allele frequency:                          0.1500
#> Dominance/Additivity ratio (tau):                    0.5000
#> Additive effect (a):                                 0.2280
#> Dominance effect (delta):                            0.1140
#> Centering mean (m):                                  0.1305
#> --------------------------------------------------------------------------
#> 
#> Genotype-Specific Quantities
#>                                  Genotype 0     Genotype 1     Genotype 2
#> Trait mean:                         -0.0975         0.2445         0.3585
#> Affected penetrance:                 0.0388         0.0781         0.0963
#> Unaffected penetrance:               0.0585         0.0278         0.0212
#> Affected genotype freq.:             0.5596         0.3972         0.0432
#> Unaffected genotype freq.:           0.8481         0.1424         0.0096
#> --------------------------------------------------------------------------
#> 
#> Selected Prevalences
#> Affected prevalence:                                0.05012
#> Unaffected prevalence:                              0.04988
#> --------------------------------------------------------------------------
```
