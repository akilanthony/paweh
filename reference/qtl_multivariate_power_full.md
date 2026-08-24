# Multivariate Falconer Power

Computes prospective power for a one-way MANOVA using Pillai's trace or
for a genotype chi-square test after joint multivariate threshold
selection.

## Usage

``` r
qtl_multivariate_power_full(
  N = NULL,
  N_case = NULL,
  alpha,
  qtl_var,
  tau,
  pd,
  cor_matrix,
  test = c("pillai", "threshold_chisq"),
  x_upper = NULL,
  x_lower = NULL,
  k = 1,
  verbose = TRUE
)
```

## Arguments

- N:

  Integer total sample size for `test = "pillai"`.

- N_case:

  Selected case sample size for `test = "threshold_chisq"`.

- alpha:

  Significance level in (0, 1).

- qtl_var:

  Numeric vector of phenotype-specific QTL variances in (0, 1).

- tau:

  Numeric vector of phenotype-specific dominance/additivity ratios.

- pd:

  Shared increaser-allele frequency in (0, 1).

- cor_matrix:

  Positive-definite phenotype correlation matrix. Its order must match
  `qtl_var` and `tau`.

- test:

  Either `"pillai"` or `"threshold_chisq"`.

- x_upper, x_lower:

  Vectors of upper- and lower-tail percentages. A case satisfies every
  upper threshold (joint AND); a control satisfies every lower threshold
  (joint AND).

- k:

  Control/case ratio for the threshold chi-square design.

- verbose:

  Logical; whether to print a polished summary.

## Value

An object of class `"qtl_multivariate_power_full"`. Both modes contain
`test`, `alpha`, `power`, and `falconer`, the complete multivariate
mixture model. For `test = "pillai"`, additional components include
total `N`, genotype frequencies/counts, the NCP, numerator and
denominator degrees of freedom, critical value, and `pillai`, containing
contrasts, characteristic roots, Pillai trace, and matrix intermediates.
For `test = "threshold_chisq"`, components include selected `N_case`,
`N_control`, and `N_total`, `k`, two-df NCP, internal `S`, `thresholds`
(bounds, penetrances, prevalences, conditional genotype frequencies, and
integration diagnostics), expected genotype counts, and sparse-cell
diagnostics.

## Details

Genotype-specific means are stored in a matrix whose rows are phenotypes
and columns are genotypes 0, 1, and 2. All genotypes share the residual
covariance matrix formed by scaling `cor_matrix` by residual standard
deviations. Joint probabilities use
[`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
with the deterministic Miwa algorithm through 20 dimensions; higher
dimensions use reproducibly seeded Genz–Bretz integration and return its
error/status.

Pillai power follows Gordon et al.'s general matrix derivation: the
characteristic roots of Phi-star determine V-star and lambda = N s
V-star / (s - V-star). The null critical value and alternative power use
central and noncentral F distributions.

The two `test` branches are distinct. With `test = "pillai"`, `N` is a
total sample size and the function performs a three-genotype one-way
MANOVA using Pillai's trace. With `test = "threshold_chisq"`, `N_case`
is a selected case count; cases must exceed every upper threshold and
controls must fall below every lower threshold, after which their
conditional genotype probabilities are compared by a two-df chi-square
test. The latter is not a MANOVA.

The multivariate mixture and component densities correspond to textbook
Eqs. 6.9–6.10; rectangular probabilities to Eq. 6.11; and prevalence and
bound construction to Eqs. 6.12–6.14 (Chapter 6, Section 6.2, pp.
332–333). The threshold branch then uses genotype chi-square Eq. 1.22
(p. 26). Eq. 6.15 is not used by this direct-integration implementation.

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer, Chapter 6, Section 6.2, Eqs.
6.9–6.14, pp. 332–333; validation in Section 6.2.4, pp. 336–339.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

Gordon, D., Londono, D., Patel, P., Kim, W., Finch, S. J., & Heiman, G.
A. (2017). An analytic solution to computation of power and sample size
for genetic association studies under a pleiotropic mode of inheritance.
*Human Heredity*, 81(4), 194–209.
[doi:10.1159/000457135](https://doi.org/10.1159/000457135) .

Pillai, K. C. S. (1955). Some new test criteria in multivariate
analysis. *Annals of Mathematical Statistics*, 26(1), 117–121.
[doi:10.1214/aoms/1177728599](https://doi.org/10.1214/aoms/1177728599) .

Genz, A., & Bretz, F. (2009). *Computation of Multivariate Normal and t
Probabilities*. Springer.
[doi:10.1007/978-3-642-01689-9](https://doi.org/10.1007/978-3-642-01689-9)
.

## See also

[`qtl_multivariate_mssn_full`](https://akilanthony.github.io/paweh/reference/qtl_multivariate_mssn_full.md),
[`qtl_falconer_parameters`](https://akilanthony.github.io/paweh/reference/qtl_falconer_parameters.md),
and
[`qtl_threshold_chisq_power`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_power.md).

## Examples

``` r
cor_matrix <- matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
qtl_multivariate_power_full(
  N = 4514, alpha = 5e-8, qtl_var = c(0.01, 0.005),
  tau = c(0, 0.5), pd = 0.25, cor_matrix = cor_matrix,
  test = "pillai", verbose = FALSE
)
```
