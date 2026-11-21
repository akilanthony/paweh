# Transmission Disequilibrium Test (TDT) Power from Expected Transmissions and Non-Transmissions

Computes the statistical power of the Transmission Disequilibrium Test
(TDT) given the expected number of transmissions (ET) and
non-transmissions (ENT) under a specified significance level. Implements
Eq. 1.25 (Chapter 1, Section 1.6.1.3, p. 27) of Gordon, Finch, and Kim
(2020).

## Usage

``` r
tdt_power_from_expected_counts(ET, ENT, alpha = 0.05)
```

## Arguments

- ET:

  Numeric. Expected number of transmissions.

- ENT:

  Numeric. Expected number of non-transmissions.

- alpha:

  Numeric. Significance level (default = 0.05).

## Value

A list containing:

- lambda:

  Non-centrality parameter.

- power:

  Computed power at the given alpha level.

- ET:

  Expected transmissions.

- ENT:

  Expected non-transmissions.

## Details

The function calculates the non-centrality parameter and statistical
power for the TDT using the chi-square distribution with 1 degree of
freedom. `ET` and `ENT` are expected counts, not probabilities, and must
refer to the same affected-trio design.

The non-centrality parameter is computed as: \$\$\lambda = \frac{(ET -
ENT)^2}{ET + ENT}\$\$

Power is then obtained as: \$\$1 - P(\chi^2\_{1,\lambda} \<
q\_{\chi^2\_{1,1-\alpha}})\$\$

## References

Spielman, R. S., McGinnis, R. E., & Ewens, W. J. (1993). Transmission
test for linkage disequilibrium: the insulin gene region and
insulin-dependent diabetes mellitus. *American Journal of Human
Genetics*, 52(3), 506–516. PMID: 8447318; PMCID: PMC1682161.

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`tdt_power`](https://akilanthony.github.io/paweh/reference/tdt_power.md)
and
[`tdt_power_from_model`](https://akilanthony.github.io/paweh/reference/tdt_power_from_model.md).

## Examples

``` r
# Example: compute power for ET = 140 and ENT = 100
tdt_power_from_expected_counts(ET = 140, ENT = 100, alpha = 0.05)
#> 
#> --- Transmission Disequilibrium Test (TDT) ---
#> Equation: 1.25  |  Input: Expected Transmissions and Non-Transmissions
#> -----------------------------------------------------------
#> Expected Transmissions (ET):             140.0000
#> Expected Non-Transmissions (ENT):        100.0000
#> Non-Centrality Parameter (lambda):         6.6667
#> Power at alpha = 0.05:                     0.7330
#> -----------------------------------------------------------
```
