# TDT Minimum Sample Size Necessary from an NCP under Phenotype Misclassification

Computes the expected transmissions (\\gT^\*\\) and non-transmissions
(\\gNT^\*\\) as well as the required number of trios (\\N^\*\\) for a
specified non-centrality parameter (\\\lambda^\*\\) under a
misclassification model. Implements Eq. 5.26 and Eqs. 5.27a–b (Chapter
5, Section 5.2.6, pp. 284–285) of Gordon, Finch, and Kim (2020). Eq.
5.28a–b is the book's numerical worked example.

## Usage

``` r
tdt_mssn_phenotype_misclassification_from_ncp(
  lambda_star,
  pd,
  prev,
  R1,
  R2,
  delta_prime = 1,
  pi01 = 0,
  digits = 5
)
```

## Arguments

- lambda_star:

  Numeric. Non-centrality parameter (\\\lambda^\*\\) derived from
  desired power (e.g., from
  [`tdt_mssn_from_model()`](https://akilanthony.github.io/pawh/reference/tdt_mssn_from_model.md)).

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

  Numeric. Misclassification rate for controls (default = 0).

- digits:

  Integer. Number of digits to round intermediate printed results
  (default = 5).

## Value

A list containing:

- gT_star:

  Expected transmission probability.

- gNT_star:

  Expected non-transmission probability.

- N_required:

  Required number of trios (\\N^\*\\).

- lambda_star:

  Non-centrality parameter.

- C, f0, f1, f2:

  Intermediate derived values used in Eq. 5.26.

## Details

This function adjusts the Transmission Disequilibrium Test (TDT) for
possible phenotype misclassification. When \\\pi\_{01} \> 0\\, a
fraction of unaffected individuals are incorrectly classified as
affected, inflating the apparent transmission probability. The expected
transmission (\\gT^\*\\) and non-transmission (\\gNT^\*\\) are
calculated as:

\$\$gT^\* = p_d p\_+ + \frac{D(p_T - \theta_1)C(1 - \pi\_{01})}{\phi_1 +
\pi\_{01}\phi_0}\$\$ \$\$gNT^\* = p_d p\_+ + \frac{D(p_A -
\theta_1)C(\pi\_{01} - 1)}{\phi_1 + \pi\_{01}\phi_0}\$\$

The required number of trios is then: \$\$N^\* = \frac{\lambda^\*
(gT^\* + gNT^\*)}{2 (gT^\* - gNT^\*)^2}\$\$

Setting \\\pi\_{01} = 0\\ reproduces the standard (non-misclassified)
TDT.

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

[`tdt_mssn`](https://akilanthony.github.io/pawh/reference/tdt_mssn.md),
[`tdt_expected_transmission_probability`](https://akilanthony.github.io/pawh/reference/tdt_expected_transmission_probability.md),
and
[`tdt_expected_nontransmission_probability`](https://akilanthony.github.io/pawh/reference/tdt_expected_nontransmission_probability.md).

## Examples

``` r
# Example: Compute N* with misclassification adjustment (pi01 = 0.1)
tdt_mssn_phenotype_misclassification_from_ncp(
  lambda_star = 7.8488,
  pd = 0.25, prev = 0.005,
  R1 = 2, R2 = 2,
  delta_prime = 1,
  pi01 = 0.1
)
#> 
#> --- Transmission Disequilibrium Test (Trios) with Misclassification ---
#> Implements Eq. 5.26 (gT_star, gNT_star) and Eq. 5.27b for N_star
#> -----------------------------------------------------------------------
#> Parameters
#> Misclassification Rate (pi01):              0.100
#> Allele Frequency (p_d):                     0.250  |  Prevalence (phi1):    0.005
#> Relative Risks (R1,R2):                      2, 2
#> LD scale c (delta_prime):                   1.000
#> -----------------------------------------------------------------------
#> Derived
#> p_plus (transmitting allele freq):        0.75000  |  phi0:  0.99500
#> f0:                                       0.00348  |  f1:  0.00696  |  f2:  0.00696
#> C:                                        0.00261  |  DpT:  0.14062  |  DpA:  0.04688
#> -----------------------------------------------------------------------
#> Expectations
#> gT_star:                                  0.19066
#> gNT_star:                                 0.18645
#> -----------------------------------------------------------------------
#> lambda_star (TDT):                        7.84880
#> Required Number of Trios (N_star):          83395
#> -----------------------------------------------------------------------
```
