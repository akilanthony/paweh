# Specialized Case-Control Phenotype Misclassification Functions

Convenience functions for genotype chi-square case-control calculations
in the presence of phenotype misclassification only. These focused
functions are narrower than
[`cc_mssn`](https://akilanthony.github.io/pawh/reference/cc_mssn.md) and
[`cc_power`](https://akilanthony.github.io/pawh/reference/cc_power.md):
they use true affected and true unaffected genotype frequencies supplied
directly by the user and do not apply locus heterogeneity or genotype
misclassification.

## Usage

``` r
cc_chisq_power_phenotype_misclassification(
  N_case,
  alpha = 0.05,
  g_aff,
  g_unaff,
  prev,
  theta = 0,
  phi = 0,
  k = 1
)

cc_chisq_mssn_phenotype_misclassification(
  target_power = 0.8,
  alpha = 0.05,
  g_aff,
  g_unaff,
  prev,
  theta = 0,
  phi = 0,
  k = 1
)
```

## Arguments

- N_case:

  Numeric \\\> 0\\. Number of cases for power calculations.

- alpha:

  Numeric in \\(0,1)\\. Significance level.

- g_aff:

  Numeric vector of true affected genotype frequencies, ordered
  consistently across affected and unaffected groups and summing to 1.

- g_unaff:

  Numeric vector of true unaffected genotype frequencies, ordered
  consistently with `g_aff` and summing to 1.

- prev:

  Numeric in \\(0,1)\\. Disease prevalence.

- theta:

  Numeric in \\\[0,1)\\. Probability that a truly affected individual is
  classified as a control, `Pr(affected -> control)`.

- phi:

  Numeric in \\\[0,1)\\. Probability that a truly unaffected individual
  is classified as a case, `Pr(unaffected -> case)`.

- k:

  Numeric \\\> 0\\. Control-to-case ratio \\N\_{ctrl}/N\_{case}\\.

- target_power:

  Numeric in \\(0,1)\\. Desired target power for MSSN calculations.

## Value

`cc_chisq_power_phenotype_misclassification()` returns a list containing
sample sizes, phenotype-error parameters, observed case/control genotype
frequencies, internal `S`, non-centrality parameter `lambda`, and power.
`cc_chisq_mssn_phenotype_misclassification()` returns a list containing
target power, phenotype-error parameters, observed case/control genotype
frequencies, internal `S`, target non-centrality parameter
`lambda_star`, and case/control/total MSSN.

## Details

Phenotype misclassification is applied by mixing the true affected
genotype distribution `g_aff` and true unaffected genotype distribution
`g_unaff` using disease prevalence and the phenotype-error
probabilities. Here `theta = Pr(affected -> control)` and
`phi = Pr(unaffected -> case)`. These specialized functions apply only
phenotype misclassification; use the full case-control functions for
combined phenotype misclassification, locus heterogeneity, genotype
misclassification, or trend tests. The observed case distribution is a
prevalence-weighted mixture of true affected and unaffected genotype
probabilities conditional on being classified as a case; the control
distribution is constructed analogously.

## References

Edwards, B. J., Haynes, C., Levenstien, M. A., Finch, S. J., & Gordon,
D. (2005). Power and sample size calculations in the presence of
phenotype errors for case/control genetic association studies. *BMC
Genetics*, 6, 18.
[doi:10.1186/1471-2156-6-18](https://doi.org/10.1186/1471-2156-6-18) .
PMID: 15819990.

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## See also

[`cc_power`](https://akilanthony.github.io/pawh/reference/cc_power.md),
[`cc_mssn`](https://akilanthony.github.io/pawh/reference/cc_mssn.md),
and
[`case_control_genotype_misclassification`](https://akilanthony.github.io/pawh/reference/case_control_genotype_misclassification.md).

## Examples

``` r
g_aff <- c((1 - 0.05)^2, 2 * 0.05 * (1 - 0.05), 0.05^2)
g_unaff <- c((1 - 0.15)^2, 2 * 0.15 * (1 - 0.15), 0.15^2)

cc_chisq_power_phenotype_misclassification(
  N_case = 250, alpha = 0.01,
  g_aff = g_aff, g_unaff = g_unaff,
  prev = 0.05, theta = 0, phi = 0.01
)
#> $N_case
#> [1] 250
#> 
#> $N_ctrl
#> [1] 250
#> 
#> $alpha
#> [1] 0.01
#> 
#> $df
#> [1] 2
#> 
#> $prev
#> [1] 0.05
#> 
#> $theta
#> [1] 0
#> 
#> $phi
#> [1] 0.01
#> 
#> $g_case_obs
#> [1] 0.873760504 0.120546218 0.005693277
#> 
#> $g_ctrl_obs
#> [1] 0.7225 0.2550 0.0225
#> 
#> $S
#> [1] 0.07248965
#> 
#> $lambda
#> [1] 18.12241
#> 
#> $power
#> [1] 0.9134913
#> 

cc_chisq_mssn_phenotype_misclassification(
  target_power = 0.8, alpha = 0.01,
  g_aff = g_aff, g_unaff = g_unaff,
  prev = 0.05, theta = 0, phi = 0.01
)
#> $target_power
#> [1] 0.8
#> 
#> $alpha
#> [1] 0.01
#> 
#> $df
#> [1] 2
#> 
#> $prev
#> [1] 0.05
#> 
#> $theta
#> [1] 0
#> 
#> $phi
#> [1] 0.01
#> 
#> $g_case_obs
#> [1] 0.873760504 0.120546218 0.005693277
#> 
#> $g_ctrl_obs
#> [1] 0.7225 0.2550 0.0225
#> 
#> $S
#> [1] 0.07248965
#> 
#> $lambda_star
#> [1] 13.8807
#> 
#> $N_case
#> [1] 192
#> 
#> $N_ctrl
#> [1] 192
#> 
#> $N_total
#> [1] 384
#> 
```
