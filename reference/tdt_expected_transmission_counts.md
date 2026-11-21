# Expected Transmissions (ET\*) and Non-Transmissions (ENT\*) in the Transmission Disequilibrium Test (TDT)

Computes the expected number of transmissions (\\ET^\*\\) and
non-transmissions (\\ENT^\*\\) for a given number of affected trios
under heterogeneity and linkage disequilibrium models. Implements
Equations 5.31-5.32 from *Gordon et al.* (2020), *Heterogeneity in
Statistical Genetics*.

## Usage

``` r
tdt_expected_transmission_counts(
  N_star,
  pd,
  prev,
  R1,
  R2,
  delta_prime = 1,
  pi = 1,
  theta1 = NULL,
  digits = 6,
  verbose = TRUE
)
```

## Arguments

- N_star:

  Numeric. Number of affected trios.

- pd:

  Numeric. Frequency of the disease-associated allele.

- prev:

  Numeric. Disease prevalence (\\\phi_1\\).

- R1:

  Numeric. Relative risk for heterozygotes.

- R2:

  Numeric. Relative risk for homozygotes.

- delta_prime:

  Numeric. Linkage disequilibrium scale parameter (\\D'\\) (default =
  1).

- pi:

  Numeric in \\\[0,1\]\\. Linked/homogeneous trio fraction (default
  `1`); \\1-\pi\\ is the heterogeneous fraction.

- theta1:

  Numeric. Population allele frequency (defaults to `pd` if `NULL`).

- digits:

  Integer. Number of digits for printing (default = 6).

- verbose:

  Logical. If `TRUE`, prints intermediate quantities (default = TRUE).

## Value

A list containing:

- ET_star:

  Expected transmissions (\\ET^\*\\).

- ENT_star:

  Expected non-transmissions (\\ENT^\*\\).

- C, D:

  Contrast and LD terms.

- f0, f1, f2:

  Derived penetrance frequencies.

- pd, prev, R1, R2, pi, theta1:

  Input parameters.

## Details

The expected transmissions (\\ET^\*\\) and non-transmissions
(\\ENT^\*\\) are computed as:

\$\$ ET^\* = 2N^\* \left\[ \pi \left( \frac{D(p\_+ -
\theta_1)C}{\phi_1} + p_d p\_+ \right) + (1 - \pi)\left( \frac{D(p\_+ -
0.5)C}{\phi_1} + p_d p\_+ \right) \right\] \$\$

\$\$ ENT^\* = 2N^\* \left\[ \pi \left( \frac{D(\theta_1 -
p_d)C}{\phi_1} + p_d p\_+ \right) + (1 - \pi)\left( \frac{D(0.5 -
p_d)C}{\phi_1} + p_d p\_+ \right) \right\] \$\$

These are Eqs. 5.31–5.32 in Chapter 5, Section 5.3.3 (pp. 293–294). They
are expected counts over `N_star` affected trios. Dividing by \\2N^\*\\
gives the corresponding \\g_T^\*\\ and \\g\_{NT}^\*\\ probabilities. The
canonical `heter_rate` used by
[`tdt_power()`](https://akilanthony.github.io/paweh/reference/tdt_power.md)
and
[`tdt_mssn()`](https://akilanthony.github.io/paweh/reference/tdt_mssn.md)
equals \\1-\pi\\.

## References

Chen, C., Yang, G., Buyske, S., Matise, T., Finch, S. J., & Gordon, D.
(2009). Transmission disequilibrium test power and sample size in the
presence of locus heterogeneity. *Statistical Applications in Genetics
and Molecular Biology*, 8, Article 44.
[doi:10.2202/1544-6115.1501](https://doi.org/10.2202/1544-6115.1501) .
PMID: 19883370.

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`tdt_power`](https://akilanthony.github.io/paweh/reference/tdt_power.md),
[`tdt_mssn`](https://akilanthony.github.io/paweh/reference/tdt_mssn.md),
and
[`tdt_mssn_from_model`](https://akilanthony.github.io/paweh/reference/tdt_mssn_from_model.md).

## Examples

``` r
# Example 1: Homogeneous model (pi = 1)
tdt_expected_transmission_counts(
  N_star = 1000, pd = 0.25, prev = 0.005,
  R1 = 2, R2 = 2,
  delta_prime = 1, pi = 1
)
#> 
#> --- Transmission Disequilibrium Test: Expected ET* and ENT* ---
#> Implements Equations 5.31-5.32 (Heterogeneity model)
#> -----------------------------------------------------------
#> Number of Trios (N*):                        1000
#> Allele Frequency (p_d):                  0.250000
#> Prevalence (phi1):                       0.005000
#> LD Scale (delta_prime):                  1.000000
#> Heterogeneity Parameter (pi):            1.000000
#> Population Allele Freq (theta1):         0.250000
#> -----------------------------------------------------------
#> Contrast Term (C):                       0.002609
#> LD Term (D):                             0.187500
#> -----------------------------------------------------------
#> Expected Transmissions (ET*):          472.826087
#> Expected Non-Transmissions (ENT*):     375.000000
#> -----------------------------------------------------------

# Example 2: Heterogeneous model (pi = 0.7)
tdt_expected_transmission_counts(
  N_star = 1000, pd = 0.25, prev = 0.005,
  R1 = 2, R2 = 2,
  delta_prime = 1, pi = 0.7
)
#> 
#> --- Transmission Disequilibrium Test: Expected ET* and ENT* ---
#> Implements Equations 5.31-5.32 (Heterogeneity model)
#> -----------------------------------------------------------
#> Number of Trios (N*):                        1000
#> Allele Frequency (p_d):                  0.250000
#> Prevalence (phi1):                       0.005000
#> LD Scale (delta_prime):                  1.000000
#> Heterogeneity Parameter (pi):            0.700000
#> Population Allele Freq (theta1):         0.250000
#> -----------------------------------------------------------
#> Contrast Term (C):                       0.002609
#> LD Term (D):                             0.187500
#> -----------------------------------------------------------
#> Expected Transmissions (ET*):          458.152174
#> Expected Non-Transmissions (ENT*):     389.673913
#> -----------------------------------------------------------
```
