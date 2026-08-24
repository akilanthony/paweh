# TDT Minimum Sample Size Necessary from Genetic Model Parameters

Computes the required number of affected trios (\\N^\*\\) needed to
achieve a specified statistical power in the Transmission Disequilibrium
Test (TDT), given model parameters for allele frequency, relative risks,
disease prevalence, and heterogeneity. Implements the probability and
MSSN chain in Eqs. 5.33–5.34b (Chapter 5, Section 5.3.3, pp. 293–294) of
Gordon, Finch, and Kim (2020).

## Usage

``` r
tdt_mssn_from_model(
  power,
  alpha,
  df,
  pd,
  prev,
  R1,
  R2,
  delta_prime = 1,
  pi = 1
)
```

## Arguments

- power:

  Numeric. Desired power (e.g., 0.8).

- alpha:

  Numeric. Significance level (e.g., 0.05).

- df:

  Integer. Degrees of freedom (typically 1 for TDT).

- pd:

  Numeric. Frequency of the disease-associated allele.

- prev:

  Numeric. Disease prevalence.

- R1:

  Numeric. Relative risk for heterozygotes.

- R2:

  Numeric. Relative risk for homozygotes.

- delta_prime:

  Numeric. Linkage disequilibrium (LD) scale factor (default = 1).

- pi:

  Numeric in \\\[0,1\]\\. Linked/homogeneous trio fraction; `1` is
  complete homogeneity and \\1-\pi\\ is the heterogeneous fraction.

## Value

A list containing:

- lambda_star:

  Non-centrality parameter (\\\lambda^\*\\).

- gT_star:

  Expected transmission probability.

- gNT_star:

  Expected non-transmission probability.

- N_star:

  Required number of trios.

## Details

This function determines the non-centrality parameter (\\\lambda^\*\\)
via root-finding (`uniroot`) such that the test power equals the desired
level. It then computes the expected transmission (\\gT^\*\\) and
non-transmission (\\gNT^\*\\) probabilities, followed by the required
number of trios using: \$\$N^\* = \frac{\lambda^\*}{2} \frac{(gT^\* +
gNT^\*)}{(gT^\* - gNT^\*)^2}\$\$

The expected transmission and non-transmission components are calculated
under allele frequency and penetrance model assumptions in Eq. 5.33.
Their NCP and MSSN are Eqs. 5.34a–b.

## References

Chen, C., Yang, G., Buyske, S., Matise, T., Finch, S. J., & Gordon, D.
(2009). Transmission disequilibrium test power and sample size in the
presence of locus heterogeneity. *Statistical Applications in Genetics
and Molecular Biology*, 8, Article 44.
[doi:10.2202/1544-6115.1501](https://doi.org/10.2202/1544-6115.1501) .

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`tdt_mssn`](https://akilanthony.github.io/paweh/reference/tdt_mssn.md),
[`tdt_power_from_model`](https://akilanthony.github.io/paweh/reference/tdt_power_from_model.md),
and
[`tdt_expected_transmission_counts`](https://akilanthony.github.io/paweh/reference/tdt_expected_transmission_counts.md).

## Examples

``` r
# Example: compute required trios for 80% power at alpha = 0.05
tdt_mssn_from_model(
  power = 0.8, alpha = 0.05, df = 1,
  pd = 0.25, prev = 0.005, R1 = 2, R2 = 2,
  delta_prime = 1, pi = 1
)
#> 
#> --- Transmission Disequilibrium Test (Trios) ---
#> Equation: 5.34b  |  Computes Required Number of Trios (N_star)
#> -----------------------------------------------------------
#> Desired Power:                              0.800  |  Significance Level (alpha):  0.050
#> Allele Frequency (p_d):                     0.250  |  Prevalence (phi1):  0.005
#> Relative Risks (R1,R2):                      2, 2
#> Heterogeneity Parameter (pi):               1.000
#> -----------------------------------------------------------
#> Non-Centrality Parameter (lambda_star):     7.8488
#> Expected Transmission (gT_star):          0.26087
#> Expected Non-Transmission (gNT_star):     0.16304
#> Required Number of Trios (N_star):            174
#> -----------------------------------------------------------
```
