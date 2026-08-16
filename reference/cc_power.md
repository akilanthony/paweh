# Case-Control Power for Conditional Genotype Frequencies

Computes power for fixed case-control sample sizes using conditional
genotype frequencies. The function supports model-based genotype
frequencies, model-free genotype frequencies, optional locus
heterogeneity, optional phenotype misclassification, optional genotype
misclassification, genotype chi-square tests and genotype trend tests.

## Usage

``` r
cc_power(
  N_case,
  alpha,
  input_mode = c("model_based", "model_free"),
  prev = NULL,
  pd = NULL,
  R2 = NULL,
  MOI = c("M", "D", "Rec"),
  g1 = NULL,
  g0 = NULL,
  locus_het = FALSE,
  pi = 1,
  pheno_misclass = FALSE,
  theta = 0,
  phi = 0,
  k = 1,
  w = c(0, 1, 2),
  geno_misclass = c("none", "1p", "2p", "3p", "diff3p"),
  e = 0,
  e1 = 0,
  e2 = 0,
  e01 = 0,
  e02 = 0,
  e03 = 0,
  case_e01 = 0,
  case_e02 = 0,
  case_e03 = 0,
  ctrl_e01 = 0,
  ctrl_e02 = 0,
  ctrl_e03 = 0,
  diff_source = c("explicit", "case", "ctrl"),
  diff_multiplier = 1,
  verbose = TRUE
)
```

## Arguments

- N_case:

  Numeric \\\> 0\\. Number of cases.

- alpha:

  Numeric in \\(0,1)\\. Significance level.

- input_mode:

  Character. One of `"model_based"` or `"model_free"`. See Details.

- prev:

  Numeric in \\(0,1)\\. Disease prevalence for
  `input_mode = "model_based"`.

- pd:

  Numeric in \\(0,1)\\. Disease-allele frequency for
  `input_mode = "model_based"`.

- R2:

  Numeric \\\> 0\\. Homozygote relative risk for
  `input_mode = "model_based"`.

- MOI:

  Character. Mode of inheritance for model-based frequencies: `"M"` for
  multiplicative, `"D"` for dominant, or `"Rec"` for recessive.

- g1, g0:

  Numeric vectors of length 3 for `input_mode = "model_free"`. `g1`
  gives case genotype frequencies and `g0` gives control genotype
  frequencies, each ordered as `c(g0, g1, g2)` and summing to 1.

- locus_het:

  Logical. If `TRUE`, applies locus heterogeneity to the case genotype
  frequencies before phenotype and genotype misclassification.

- pi:

  Numeric in \\\[0,1\]\\. Locus-homogeneity fraction used when
  `locus_het = TRUE`.

- pheno_misclass:

  Logical. If `TRUE`, applies phenotype misclassification before
  genotype misclassification.

- theta:

  Numeric in \\\[0,1)\\. Probability that a truly affected individual is
  classified as a control.

- phi:

  Numeric in \\\[0,1)\\. Probability that a truly unaffected individual
  is classified as a case.

- k:

  Numeric \\\> 0\\. Control-to-case sample size ratio \\N\_{ctrl} /
  N\_{case}\\.

- w:

  Numeric vector of length 3. Genotype trend-test scores. The three
  weights cannot all be equal.

- geno_misclass:

  Character. Genotype misclassification model: `"none"`, `"1p"`, `"2p"`,
  `"3p"`, or `"diff3p"`.

- e:

  Numeric in \\\[0,0.5\]\\. For `geno_misclass = "1p"`, the probability
  assigned to each adjacent off-diagonal genotype call. The textbook Eq.
  2.5 parameter satisfies \\\epsilon=2e\\; see Details.

- e1, e2:

  Numeric. Two-parameter misclassification rates for
  `geno_misclass = "2p"`.

- e01, e02, e03:

  Numeric. Non-differential three-parameter misclassification rates for
  `geno_misclass = "3p"`.

- case_e01, case_e02, case_e03:

  Numeric. Case-specific three-parameter misclassification rates for
  `geno_misclass = "diff3p"`.

- ctrl_e01, ctrl_e02, ctrl_e03:

  Numeric. Control-specific three-parameter misclassification rates for
  `geno_misclass = "diff3p"`.

- diff_source:

  Character. For `geno_misclass = "diff3p"`, one of `"explicit"`,
  `"case"`, or `"ctrl"`. See Details.

- diff_multiplier:

  Numeric \\\ge 0\\. Multiplier used when `diff_source = "case"` or
  `diff_source = "ctrl"`.

- verbose:

  Logical. If `TRUE`, prints a clean formatted summary.

## Value

An object of class `"cc_power"`, containing:

- alpha, N_case, N_ctrl, N_total:

  Significance level and numbers of case, control, and total
  individuals.

- input_mode, k, w, locus_het:

  Input mode, control-to-case ratio, trend scores, and
  locus-heterogeneity settings.

- errors:

  Genotype-error model and matrices, plus phenotype-error settings and
  intermediate frequencies.

- model_info:

  Model-based penetrances and risk-model inputs, or model-free
  identifying information.

- tests\$genotypes, tests\$trend:

  Test label, degrees of freedom, NCP `lambda`, internal `S`, and power.
  The trend result also contains its numerator and denominator.

- freqs:

  Baseline, post-heterogeneity (true), post-phenotype-error, and final
  observed case and control genotype-probability vectors.

## Details

The workflow is:

1.  Construct baseline conditional genotype frequencies for cases and
    controls.

2.  Optionally apply locus heterogeneity to cases as \\g\_{case,true} =
    \pi g\_{case,base} + (1 - \pi) g\_{ctrl,base}\\.

3.  Optionally apply phenotype misclassification, where
    `theta = Pr(affected -> control)` and
    `phi = Pr(unaffected -> case)`.

4.  Optionally apply genotype misclassification matrices to the
    resulting case and control genotype frequencies.

5.  Compute genotype chi-square and genotype trend-test non-centrality
    parameters and powers from the observed genotype frequencies.

With `input_mode = "model_based"`, conditional case and control genotype
frequencies are derived from `prev`, `pd`, `R2`, and `MOI` using Chapter
1, Section 1.4.2, Eqs. 1.6–1.7 (p. 13) of Gordon, Finch, and Kim (2020).
With `input_mode = "model_free"`, the user supplies `g1` and `g0`
directly; when phenotype misclassification is enabled, `g1` is treated
as the true affected genotype distribution and `g0` is treated as the
true unaffected genotype distribution.

Phenotype misclassification requires `prev` in both input modes because
disease prevalence is used to mix the true affected and unaffected
genotype distributions into observed case and control genotype
distributions. It is applied after optional locus heterogeneity and
before optional genotype misclassification.

The genotype misclassification models are: `"none"` for identity
matrices, `"1p"` for one symmetric error rate, `"2p"` for adjacent
homozygote/heterozygote and heterozygote error rates, `"3p"` for a
non-differential three-parameter matrix, and `"diff3p"` for separate
case and control three-parameter matrices. These matrices correspond to
Eqs. 2.5–2.7 (pp. 57–58). For `"1p"`, package `e` is the probability in
each adjacent off-diagonal cell and the diagonal is \\1-2e\\; textbook
Eq. 2.5 uses off-diagonal \\\epsilon/2\\, so \\\epsilon=2e\\. The
genotype and trend tests use Eqs. 1.22 and 1.24 (pp. 26–27), with the
trend statistic defined by Eqs. 1.20–1.21 (p. 24). Locus heterogeneity
follows Eq. 2.16 (p. 88), where `pi` is the homogeneous fraction; its
trend-test form is given in Eqs. 5.29a–b (pp. 287–288).

For `geno_misclass = "diff3p"`, `diff_source = "explicit"` uses the case
and control error parameters exactly as supplied. With
`diff_source = "case"`, control parameters are computed by multiplying
the case parameters by `diff_multiplier`. With `diff_source = "ctrl"`,
case parameters are computed by multiplying the control parameters by
`diff_multiplier`.

With `geno_misclass = "diff3p"`, returned power is nominal asymptotic
power evaluated with the usual chi-square critical value. Different case
and control error mechanisms can distort the null distribution and
inflate type I error. This function does not independently recalibrate
the null distribution for arbitrary differential genotyping error.

Internal effect-size components `S` are retained in the returned object
for validation and debugging, but are not printed in the clean verbose
output.

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

Armitage, P. (1955). Tests for linear trends in proportions and
frequencies. *Biometrics*, 11(3), 375–386.
[doi:10.2307/3001775](https://doi.org/10.2307/3001775) .

Slager, S. L., & Schaid, D. J. (2001). Case-control studies of genetic
markers: power and sample size approximations for Armitage's test for
trend. *Human Heredity*, 52(3), 149–153.
[doi:10.1159/000053370](https://doi.org/10.1159/000053370) .

Edwards, B. J., Haynes, C., Levenstien, M. A., Finch, S. J., & Gordon,
D. (2005). Power and sample size calculations in the presence of
phenotype errors for case/control genetic association studies. *BMC
Genetics*, 6, 18.
[doi:10.1186/1471-2156-6-18](https://doi.org/10.1186/1471-2156-6-18) .

## See also

[`cc_mssn`](https://akilanthony.github.io/pawh/reference/cc_mssn.md),
[`case_control_genotype_misclassification`](https://akilanthony.github.io/pawh/reference/case_control_genotype_misclassification.md),
[`case_control_locus_heterogeneity`](https://akilanthony.github.io/pawh/reference/case_control_locus_heterogeneity.md),
and
[`case_control_phenotype_misclassification`](https://akilanthony.github.io/pawh/reference/case_control_phenotype_misclassification.md).

## Examples

``` r
cc_power(
  N_case = 500, alpha = 0.05,
  input_mode = "model_free",
  g1 = c(0.25, 0.50, 0.25),
  g0 = c(0.36, 0.48, 0.16),
  geno_misclass = "none",
  verbose = FALSE
)

cc_power(
  N_case = 500, alpha = 0.05,
  input_mode = "model_free",
  g1 = c(0.25, 0.50, 0.25),
  g0 = c(0.36, 0.48, 0.16),
  locus_het = TRUE, pi = 0.8,
  verbose = FALSE
)

cc_power(
  N_case = 500, alpha = 0.05,
  input_mode = "model_free",
  g1 = c(0.25, 0.50, 0.25),
  g0 = c(0.36, 0.48, 0.16),
  locus_het = TRUE, pi = 0.8,
  geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005,
  verbose = FALSE
)

cc_power(
  N_case = 500, alpha = 0.05,
  input_mode = "model_based",
  prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M",
  pheno_misclass = TRUE, theta = 0.05, phi = 0.01,
  verbose = FALSE
)

cc_power(
  N_case = 500, alpha = 0.05,
  input_mode = "model_based",
  prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M",
  locus_het = TRUE, pi = 0.8,
  geno_misclass = "3p", e01 = 0.02, e02 = 0.01, e03 = 0.005,
  verbose = FALSE
)

cc_power(
  N_case = 500, alpha = 0.05,
  input_mode = "model_free",
  g1 = c(0.25, 0.50, 0.25),
  g0 = c(0.36, 0.48, 0.16),
  geno_misclass = "diff3p",
  diff_source = "explicit",
  case_e01 = 0.02, case_e02 = 0.01, case_e03 = 0.005,
  ctrl_e01 = 0.01, ctrl_e02 = 0.005, ctrl_e03 = 0.002,
  verbose = FALSE
)

cc_power(
  N_case = 500, alpha = 0.05,
  input_mode = "model_free",
  g1 = c(0.25, 0.50, 0.25),
  g0 = c(0.36, 0.48, 0.16),
  geno_misclass = "diff3p",
  diff_source = "case", diff_multiplier = 0.5,
  case_e01 = 0.02, case_e02 = 0.01, case_e03 = 0.005,
  verbose = FALSE
)
```
