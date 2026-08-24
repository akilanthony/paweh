# Transmission Disequilibrium Test (TDT) Power from Genetic Model Parameters

Computes the statistical power of the Transmission Disequilibrium Test
(TDT) using genetic model parameters such as allele frequency, relative
risks, disease prevalence, and the number of affected trios. Implements
the penetrance construction in Eqs. 1.6–1.7 (p. 13) and TDT Eq. 1.25 (p.
27) of Gordon, Finch, and Kim (2020).

## Usage

``` r
tdt_power_from_model(
  pd,
  N,
  delta_prime,
  f0 = NULL,
  f1 = NULL,
  f2 = NULL,
  prev = NULL,
  R1 = NULL,
  R2 = NULL,
  alpha = 0.05
)
```

## Arguments

- pd:

  Numeric. Frequency of the disease-associated allele.

- N:

  Numeric. Number of affected trios.

- delta_prime:

  Numeric. Linkage disequilibrium (LD) scale factor (default = 1).

- f0, f1, f2:

  Optional. Penetrances for genotypes with 0, 1, and 2 risk alleles. If
  not provided, they are computed internally using `prev`, `R1`, and
  `R2`.

- prev:

  Numeric. Disease prevalence.

- R1:

  Numeric. Relative risk for heterozygotes.

- R2:

  Numeric. Relative risk for homozygotes.

- alpha:

  Numeric. Significance level (default = 0.05).

## Value

A list containing:

- lambda:

  Non-centrality parameter.

- power:

  Computed power at the given significance level.

- ET:

  Expected transmissions.

- ENT:

  Expected non-transmissions.

- Penetrances:

  Vector of computed penetrances (f0, f1, f2).

## Details

When penetrances (`f0`, `f1`, `f2`) are not provided, they are derived
from the model parameters using: \$\$f0 = prev / Z, \quad f1 = R1 \* f0,
\quad f2 = R2 \* f0\$\$ where \$\$Z = (1 - p_d)^2 + 2 \* p_d \* (1 -
p_d) \* R1 + p_d^2 \* R2.\$\$

Expected transmission (\\ET\\) and non-transmission (\\ENT\\) counts are
computed based on the allele frequency and penetrance model, and the
power is derived from the non-central chi-square distribution with 1
degree of freedom. Here `N` is the number of affected-child trios.
`delta_prime` scales \\D=p_d(1-p_d)D'\\ in the implemented model; the
disease-locus and marker-locus assumptions should be considered when
interpreting `pd`.

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

[`tdt_power`](https://akilanthony.github.io/paweh/reference/tdt_power.md),
[`tdt_power_from_expected_counts`](https://akilanthony.github.io/paweh/reference/tdt_power_from_expected_counts.md),
and
[`tdt_mssn_from_model`](https://akilanthony.github.io/paweh/reference/tdt_mssn_from_model.md).

## Examples

``` r
# Example: model-based power computation
tdt_power_from_model(
  pd = 0.25, N = 10000, delta_prime = 1,
  prev = 0.05, R1 = 1, R2 = 1.1, alpha = 0.05
)
#> 
#> --- Transmission Disequilibrium Test (Model-Based) ---
#> Equation: 1.25  |  Inputs: Allele Frequency and Genetic Model Parameters
#> -----------------------------------------------------------
#> Allele Frequency (p_d):                    0.2500
#> Prevalence (phi1):                         0.0500
#> Relative Risks (R1,R2):                    1, 1.1
#> Number of Affected Trios (N):               10000
#> -----------------------------------------------------------
#> Non-Centrality Parameter (lambda):         1.1502
#> Power at alpha = 0.05:                     0.1886
#> Expected ET:                              3819.88  |  Expected ENT:    3726.71
#> -----------------------------------------------------------
```
