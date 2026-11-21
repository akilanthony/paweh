# Expected Transmission Probability (gT\*) in the Transmission Disequilibrium Test (TDT)

Computes the expected transmission probability (\\g_T^\*\\) for the
Transmission Disequilibrium Test (TDT), incorporating effects of linkage
disequilibrium (LD), relative risks, and phenotype misclassification.
Implements Equation 5.24 from *Gordon et al.* (2020), *Heterogeneity in
Statistical Genetics*.

## Usage

``` r
tdt_expected_transmission_probability(
  pd,
  prev,
  R1,
  R2,
  delta_prime = 1,
  pi01 = 0,
  theta1 = NULL,
  digits = 6,
  verbose = TRUE
)
```

## Arguments

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

- pi01:

  Numeric. Misclassification rate (default = 0).

- theta1:

  Numeric. Population allele frequency (defaults to `pd` if `NULL`).

- digits:

  Integer. Number of digits for printed results (default = 6).

- verbose:

  Logical. If `TRUE`, prints intermediate quantities (default = TRUE).

## Value

A list containing:

- gT_star:

  Expected transmission probability (\\g_T^\*\\).

- C, D:

  Contrast and LD terms.

- f0, f1, f2:

  Derived penetrance frequencies.

- pd, phi1, phi0, theta1:

  Input and derived parameters.

## Details

The function computes the penetrance-weighted genotype frequencies
\\f_0,f_1,f_2\\ and the contrast term \\C\\, representing
genotype-specific effects. Linkage disequilibrium is modeled as \\D = D'
p_d p\_+\\.

The expected transmission probability is given by:

\$\$ g_T^\* = \frac{ p_d p\_+ \phi_1 + D(p\_+ - \theta_1)C +
\pi\_{01}(p_d p\_+ \phi_0 + D(\theta_1 - p\_+)C) }{\phi_1 +
\pi\_{01}\phi_0} \$\$

When \\\pi\_{01}=0\\, this reduces to the standard non-misclassified
case. This is Eq. 5.24 in Chapter 5, Section 5.2.6 (p. 284). `gT_star`
is a probability, not an expected count.

## References

Buyske, S., Yang, G., Matise, T. C., & Gordon, D. (2009). When a case is
not a case: Effects of phenotype misclassification on power and sample
size requirements for the transmission disequilibrium test with affected
child trios. *Human Heredity*, 67(4), 287–292.
[doi:10.1159/000194981](https://doi.org/10.1159/000194981) . PMID:
19172087.

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`tdt_expected_nontransmission_probability`](https://akilanthony.github.io/paweh/reference/tdt_expected_nontransmission_probability.md),
[`tdt_expected_transmission_counts`](https://akilanthony.github.io/paweh/reference/tdt_expected_transmission_counts.md),
and
[`tdt_power`](https://akilanthony.github.io/paweh/reference/tdt_power.md).

## Examples

``` r
# Example 1 – no misclassification
tdt_expected_transmission_probability(
  pd = 0.25, prev = 0.005,
  R1 = 2, R2 = 2,
  delta_prime = 1, pi01 = 0
)
#> 
#> --- Transmission Disequilibrium Test: Expected Transmission (gT*) ---
#> Implements Equation 5.24 from Gordon et al. (2020)
#> -----------------------------------------------------------
#> Allele Frequency (p_d):                  0.250000
#> Prevalence (phi1):                       0.005000
#> LD Scale (delta_prime):                  1.000000
#> Misclassification Rate (pi01):           0.000000
#> Population Allele Freq (theta1):         0.250000
#> -----------------------------------------------------------
#> Contrast Term (C):                       0.002609
#> LD Term (D):                             0.187500
#> Expected Transmission (gT*):             0.236413
#> -----------------------------------------------------------

# Example 2 – 10% misclassification
tdt_expected_transmission_probability(
  pd = 0.25, prev = 0.005,
  R1 = 2, R2 = 2,
  delta_prime = 1, pi01 = 0.1
)
#> 
#> --- Transmission Disequilibrium Test: Expected Transmission (gT*) ---
#> Implements Equation 5.24 from Gordon et al. (2020)
#> -----------------------------------------------------------
#> Allele Frequency (p_d):                  0.250000
#> Prevalence (phi1):                       0.005000
#> LD Scale (delta_prime):                  1.000000
#> Misclassification Rate (pi01):           0.100000
#> Population Allele Freq (theta1):         0.250000
#> -----------------------------------------------------------
#> Contrast Term (C):                       0.002609
#> LD Term (D):                             0.187500
#> Expected Transmission (gT*):             0.189606
#> -----------------------------------------------------------
```
