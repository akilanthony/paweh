# Analytic MSSN for a Case-Control Sequencing Study

Computes the minimum sample size necessary (MSSN) for a model-based
case-control sequencing trend design. It evaluates the same independent
sequencing-only, phenotype, heterogeneity, and post-call genotype-error
scenarios as
[`cc_ngs_power`](https://akilanthony.github.io/paweh/reference/cc_ngs_power.md).

## Usage

``` r
cc_ngs_mssn(
  power,
  alpha,
  prev,
  pd,
  R2,
  coverage = NULL,
  seq_error = NULL,
  MOI = c("M", "D", "Rec"),
  k = 1,
  verbose = TRUE,
  locus_het = FALSE,
  pi = 1,
  case_coverage = coverage,
  ctrl_coverage = coverage,
  case_seq_error = seq_error,
  ctrl_seq_error = seq_error,
  pheno_misclass = FALSE,
  theta = 0,
  phi = 0,
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
  diff_multiplier = 1
)
```

## Arguments

- power:

  Numeric in \\(0,1)\\. Requested power.

- alpha:

  Numeric in \\(0,1)\\. Significance level.

- prev:

  Numeric in \\(0,1)\\. Disease prevalence.

- pd:

  Numeric in \\(0,1)\\. Disease-allele frequency.

- R2:

  Numeric \\\> 0\\. Homozygote relative risk.

- coverage:

  Common positive integer fixed sequencing depth. May be omitted when
  both group-specific depths are supplied.

- seq_error:

  Symmetric per-read sequencing-error probability in \\\[0,0.5)\\.
  Common shorthand; may be omitted when both group-specific errors are
  supplied.

- MOI:

  Character mode of inheritance: `"M"` for multiplicative, `"D"` for
  dominant, or `"Rec"` for recessive.

- k:

  Numeric \\\> 0\\. Planned control-to-case sample-size ratio.

- verbose:

  Logical. If `TRUE`, print a concise result summary.

- locus_het:

  Logical. If `TRUE`, apply the canonical PAWEH case-control
  locus-heterogeneity mixture before sequencing observation.

- pi:

  Numeric in \\\[0,1\]\\. Locus-homogeneity fraction used when
  `locus_het = TRUE`. One retains the original associated-case
  distribution; zero makes the case distribution equal to controls. When
  `locus_het = FALSE`, `pi` must remain at its default value of 1.

- case_coverage, ctrl_coverage:

  Case/control fixed depths, respectively. Each defaults to `coverage`;
  an explicit value overrides it.

- case_seq_error, ctrl_seq_error:

  Case/control per-read error probabilities. Each defaults to
  `seq_error`; an explicit value overrides it. Effective depths and
  errors obey the common parameter constraints.

- pheno_misclass:

  Logical. If `TRUE`, apply the ordinary PAWEH case-control
  phenotype-misclassification model before sequencing.

- theta:

  Numeric in \\\[0,1)\\. Probability that a truly affected individual is
  classified as a control.

- phi:

  Numeric in \\\[0,1)\\. Probability that a truly unaffected individual
  is classified as a case.

- geno_misclass:

  Character genotype-misclassification model: `"none"`, `"1p"`, `"2p"`,
  `"3p"`, or `"diff3p"`. The selected model is applied after sequencing
  calls.

- e:

  Numeric in \\\[0,0.5\]\\. One-parameter symmetric error rate.

- e1, e2:

  Numeric two-parameter homozygote-to-heterozygote and
  heterozygote-to-homozygote error rates.

- e01, e02, e03:

  Numeric non-differential three-parameter error rates.

- case_e01, case_e02, case_e03:

  Case-specific three-parameter error rates.

- ctrl_e01, ctrl_e02, ctrl_e03:

  Control-specific three-parameter error rates.

- diff_source:

  Character source for differential three-parameter errors:
  `"explicit"`, `"case"`, or `"ctrl"`.

- diff_multiplier:

  Nonnegative multiplier used to derive the other group's rates when
  `diff_source` is `"case"` or `"ctrl"`.

## Value

Invisibly, an object of class `"cc_ngs_mssn"` containing a
scenario-first `scenarios` list. Each scenario records target NCP,
`D`/`Q`, continuous and integer sample sizes, achieved NCP and power,
rounding metadata, staged frequencies, and sequencing/error matrices.
Legacy fields follow the `compatibility_scenario` convention documented
for
[`cc_ngs_power`](https://akilanthony.github.io/paweh/reference/cc_ngs_power.md).

## Details

Locus heterogeneity uses \\g\_{case,H}=\pi
g\_{case}+(1-\pi)g\_{control}\\ in its independent branch. The
prevalence-weighted phenotype transformation from
[`cc_mssn`](https://akilanthony.github.io/paweh/reference/cc_mssn.md) is
applied directly to the common biological baseline in its branch. Case
and control sequencing settings are then applied to the corresponding
pre-sequencing groups. Conventional genotype misclassification is
applied after sequencing only in its own branch. See
[`cc_ngs_power`](https://akilanthony.github.io/paweh/reference/cc_ngs_power.md)
for the scenario paths and matrix convention. When \\\pi=0\\ and
sequencing mechanisms are equal, no finite MSSN exists for target power
greater than `alpha` if the observed contrast is zero. Differential
mechanisms may create an observed contrast even under a biological null.
MSSN is a nominal asymptotic calculation using the existing chi-square
critical value, without separate null-distribution or Type I error
recalibration.

The function numerically inverts the one-degree-of-freedom noncentral
chi-square distribution only to obtain the target NCP. It then solves
the Ahn/Chapman-Nam trend-test sample-size equation analytically.
Planned cases are the ceiling of the continuous requirement and planned
controls are `ceiling(k * MSSN_case)`, following the PAWEH convention.
Achieved NCP and power are recomputed using these actual integer sample
sizes, with a local boundary adjustment if required to ensure target
attainment and integer minimality.

Scores are `c(0,1,2)` for `"M"`, `c(0,1,1)` for `"D"`, and `c(0,0,1)`
for `"Rec"`. Finite depth can cause genotype-call uncertainty even when
`seq_error = 0`, because a true heterozygote can yield reads from only
one allele.

This is an analytic study-design calculation, not LTTae,NGS, a raw-read
latent-genotype likelihood or EM method, downstream association testing,
or simulation.

## References

Ahn, K., Haynes, C., Kim, W., St. Fleur, R., Gordon, D., & Finch, S. J.
(2007). The effects of SNP genotyping errors on the power of the
Cochran-Armitage linear trend test for case/control association studies.
*Annals of Human Genetics*, 71, 249–261.
[doi:10.1111/j.1469-1809.2006.00318.x](https://doi.org/10.1111/j.1469-1809.2006.00318.x)
.

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## Examples

``` r
cc_ngs_mssn(
  power = 0.80, alpha = 0.05,
  prev = 0.05, pd = 0.30, R2 = 1.8,
  coverage = 20, seq_error = 0.01,
  MOI = "M", verbose = FALSE
)
```
