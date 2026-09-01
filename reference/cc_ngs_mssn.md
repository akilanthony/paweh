# Analytic MSSN for a Case-Control Sequencing Study

Computes the minimum sample size necessary (MSSN) for a model-based
case-control sequencing trend design. It uses the same fixed-depth,
symmetric sequencing-error model and deterministic maximum-likelihood
genotype calls as
[`cc_ngs_power`](https://akilanthony.github.io/paweh/reference/cc_ngs_power.md).

## Usage

``` r
cc_ngs_mssn(
  power,
  alpha,
  prev,
  pd,
  R2,
  coverage,
  seq_error,
  MOI = c("M", "D", "Rec"),
  k = 1,
  verbose = TRUE,
  locus_het = FALSE,
  pi = 1
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

  Positive integer sequencing depth.

- seq_error:

  Symmetric per-read sequencing-error probability in \\\[0,0.5)\\.

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

## Value

Invisibly, an object of class `"cc_ngs_mssn"` containing the target and
achieved power and NCP, continuous and integer sample sizes, model and
sequencing inputs, true and called genotype frequencies, and the
true-to-called transition matrix.

## Details

Locus heterogeneity is applied to true case genotype probabilities as
\\g\_{case,H}=\pi g\_{case}+(1-\pi)g\_{control}\\, using the same
parameterization as ordinary PAWEH case-control design. Sequencing
observation follows this mixture. Under the current nondifferential
error model, applying the common transition matrix before or after
forming the mixture is algebraically equivalent; this need not hold for
future differential case/control sequencing error. When \\\pi=0\\, no
finite MSSN exists for target power greater than `alpha` because the
trend contrast is zero.

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
