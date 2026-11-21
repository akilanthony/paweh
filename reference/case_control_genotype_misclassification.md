# Specialized Case-Control Genotype Misclassification Functions

Convenience functions for genotype-only case-control chi-square
calculations under conditional model-based genotype frequencies and
genotype misclassification. These functions are narrower than
[`cc_mssn`](https://akilanthony.github.io/paweh/reference/cc_mssn.md)
and
[`cc_power`](https://akilanthony.github.io/paweh/reference/cc_power.md):
they compute only the genotype chi-square test for one specific
misclassification model.

## Usage

``` r
cc_chisq_mssn_genotype_misclassification_1p(
  power,
  alpha,
  prev,
  pd,
  R2,
  MOI = c("M", "D", "Rec"),
  k = 1,
  e = 0,
  verbose = TRUE
)

cc_chisq_power_genotype_misclassification_1p(
  N_case,
  alpha,
  prev,
  pd,
  R2,
  MOI = c("M", "D", "Rec"),
  k = 1,
  e = 0,
  verbose = TRUE
)

cc_chisq_mssn_genotype_misclassification_2p(
  power,
  alpha,
  prev,
  pd,
  R2,
  MOI = c("M", "D", "Rec"),
  k = 1,
  e1 = 0,
  e2 = 0,
  verbose = TRUE
)

cc_chisq_power_genotype_misclassification_2p(
  N_case,
  alpha,
  prev,
  pd,
  R2,
  MOI = c("M", "D", "Rec"),
  k = 1,
  e1 = 0,
  e2 = 0,
  verbose = TRUE
)

cc_chisq_mssn_genotype_misclassification_3p(
  power,
  alpha,
  prev,
  pd,
  R2,
  MOI = c("M", "D", "Rec"),
  k = 1,
  e01 = 0,
  e02 = 0,
  e03 = 0,
  verbose = TRUE
)

cc_chisq_power_genotype_misclassification_3p(
  N_case,
  alpha,
  prev,
  pd,
  R2,
  MOI = c("M", "D", "Rec"),
  k = 1,
  e01 = 0,
  e02 = 0,
  e03 = 0,
  verbose = TRUE
)

cc_chisq_mssn_genotype_misclassification_differential_3p(
  power,
  alpha,
  prev,
  pd,
  R2,
  MOI = c("M", "D", "Rec"),
  k = 1,
  case_e01 = 0,
  case_e02 = 0,
  case_e03 = 0,
  ctrl_e01 = 0,
  ctrl_e02 = 0,
  ctrl_e03 = 0,
  verbose = TRUE
)

cc_chisq_power_genotype_misclassification_differential_3p(
  N_case,
  alpha,
  prev,
  pd,
  R2,
  MOI = c("M", "D", "Rec"),
  k = 1,
  case_e01 = 0,
  case_e02 = 0,
  case_e03 = 0,
  ctrl_e01 = 0,
  ctrl_e02 = 0,
  ctrl_e03 = 0,
  verbose = TRUE
)
```

## Arguments

- power:

  Numeric in \\(0,1)\\. Desired target power for MSSN functions.

- alpha:

  Numeric in \\(0,1)\\. Significance level.

- prev:

  Numeric in \\(0,1)\\. Disease prevalence.

- pd:

  Numeric in \\(0,1)\\. Disease-allele frequency.

- R2:

  Numeric \\\> 0\\. Homozygote relative risk.

- MOI:

  Character. Mode of inheritance: `"M"`, `"D"`, or `"Rec"`.

- k:

  Numeric \\\> 0\\. Control-to-case ratio \\N\_{ctrl}/N\_{case}\\.

- e:

  Numeric in \\\[0,0.5\]\\. One-parameter symmetric error probability
  assigned to each adjacent off-diagonal call. The matrix diagonal is
  \\1-2e\\; textbook Eq. 2.5 uses a total error parameter
  \\\epsilon=2e\\.

- verbose:

  Logical. If `TRUE`, prints a formatted summary.

- N_case:

  Numeric \\\> 0\\. Number of cases for power functions.

- e1, e2:

  Numeric. Two-parameter genotype misclassification rates.

- e01, e02, e03:

  Numeric. Three-parameter non-differential genotype misclassification
  rates.

- case_e01, case_e02, case_e03:

  Numeric. Case-specific three-parameter genotype misclassification
  rates for differential misclassification.

- ctrl_e01, ctrl_e02, ctrl_e03:

  Numeric. Control-specific three-parameter genotype misclassification
  rates for differential misclassification.

## Value

A list with class matching the function name. MSSN functions include
target non-centrality parameter, case/control sample sizes, model
parameters, penetrances, misclassification matrices, internal `S`
values, and true and observed genotype frequencies. Power functions
include case/control sample sizes, non-centrality parameter, power,
model parameters, penetrances, misclassification matrices, internal `S`
values, and true and observed genotype frequencies.

## Details

These are specialized convenience functions retained for users who want
a focused genotype-only calculation. The all-in-one case-control
functions are recommended when combining model-free inputs, locus
heterogeneity, trend tests, or multiple modifier choices in one call.

Conditional genotype probabilities are constructed from prevalence,
allele frequency, and relative risks using Chapter 1, Section 1.4.2 and
Eqs. 1.6–1.7 (p. 13) of Gordon, Finch, and Kim (2020). The genotype test
uses the two-degree-of-freedom NCP in Eq. 1.22 (p. 26). The error
matrices are the one-, two-, and three-parameter families in Eqs. 2.5,
2.6, and 2.7 (pp. 57–58). Douglas et al. (2002) and Sobel et al. (2002)
support the corresponding error-model families; Gordon et al. (2002)
supports their use in case-control power and sample-size calculations.

For the 1p functions, package `e` is the probability in each adjacent
off-diagonal cell, and the diagonal is \\1-2e\\. Textbook Eq. 2.5
instead writes each off-diagonal as \\\epsilon/2\\ and the diagonal as
\\1-\epsilon\\; therefore \\\epsilon=2e\\.

Differential 3p functions use separate Eq. 2.7-family matrices for cases
and controls, as discussed in Chapter 2, Section 2.5.2 (pp. 61–69).
Their power is nominal asymptotic power, and their MSSN is nominal MSSN,
evaluated with the usual chi-square critical value. Different case and
control error mechanisms can distort the null distribution and inflate
type I error. These functions do not independently recalibrate the null
distribution for arbitrary differential genotyping error.

The returned objects keep internal effect-size components such as `S`
and, for differential misclassification, per-genotype components, for
validation and debugging.

## References

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

Gordon, D., Finch, S. J., Nothnagel, M., & Ott, J. (2002). Power and
sample size calculations for case-control genetic association tests when
errors are present: application to single nucleotide polymorphisms.
*Human Heredity*, 54(1), 22–33.
[doi:10.1159/000066696](https://doi.org/10.1159/000066696) .

Douglas, J. A., Skol, A. D., & Boehnke, M. (2002). Probability of
detection of genotyping errors and mutations as inheritance
inconsistencies in nuclear-family data. *American Journal of Human
Genetics*, 70(2), 487–495.
[doi:10.1086/338919](https://doi.org/10.1086/338919) .

Sobel, E., Papp, J. C., & Lange, K. (2002). Detection and integration of
genotyping errors in statistical genetics. *American Journal of Human
Genetics*, 70(2), 496–508.
[doi:10.1086/338920](https://doi.org/10.1086/338920) .

Ahn, K., Gordon, D., & Finch, S. J. (2009). Increase of rejection rate
in case-control studies with the differential genotyping error rates.
*Statistical Applications in Genetics and Molecular Biology*, 8(1),
Article 25.
[doi:10.2202/1544-6115.1429](https://doi.org/10.2202/1544-6115.1429) .

## See also

[`cc_power`](https://akilanthony.github.io/paweh/reference/cc_power.md),
[`cc_mssn`](https://akilanthony.github.io/paweh/reference/cc_mssn.md),
[`case_control_locus_heterogeneity`](https://akilanthony.github.io/paweh/reference/case_control_locus_heterogeneity.md),
and
[`case_control_phenotype_misclassification`](https://akilanthony.github.io/paweh/reference/case_control_phenotype_misclassification.md).

## Examples

``` r
cc_chisq_mssn_genotype_misclassification_1p(
  power = 0.8, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
  MOI = "D", e = 0.02, verbose = FALSE
)

cc_chisq_power_genotype_misclassification_1p(
  N_case = 500, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
  MOI = "D", e = 0.02, verbose = FALSE
)

cc_chisq_mssn_genotype_misclassification_2p(
  power = 0.8, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
  MOI = "D", e1 = 0.02, e2 = 0.01, verbose = FALSE
)

cc_chisq_power_genotype_misclassification_2p(
  N_case = 500, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
  MOI = "D", e1 = 0.02, e2 = 0.01, verbose = FALSE
)

cc_chisq_mssn_genotype_misclassification_3p(
  power = 0.8, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
  MOI = "D", e01 = 0.02, e02 = 0.01, e03 = 0.005,
  verbose = FALSE
)

cc_chisq_power_genotype_misclassification_3p(
  N_case = 500, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
  MOI = "D", e01 = 0.02, e02 = 0.01, e03 = 0.005,
  verbose = FALSE
)

cc_chisq_mssn_genotype_misclassification_differential_3p(
  power = 0.8, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
  MOI = "D",
  case_e01 = 0.02, case_e02 = 0.01, case_e03 = 0.005,
  ctrl_e01 = 0.01, ctrl_e02 = 0.005, ctrl_e03 = 0.002,
  verbose = FALSE
)

cc_chisq_power_genotype_misclassification_differential_3p(
  N_case = 500, alpha = 0.05, prev = 0.1, pd = 0.3, R2 = 1.8,
  MOI = "D",
  case_e01 = 0.02, case_e02 = 0.01, case_e03 = 0.005,
  ctrl_e01 = 0.01, ctrl_e02 = 0.005, ctrl_e03 = 0.002,
  verbose = FALSE
)
```
