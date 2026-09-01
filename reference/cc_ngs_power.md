# Analytic Power for a Case-Control Sequencing Study

Computes prospective asymptotic power for a model-based case-control
sequencing design. The calculation constructs true case and control
genotype probabilities, applies a fixed-depth symmetric sequencing-error
model with deterministic maximum-likelihood genotype calls, and
evaluates the Ahn/Chapman-Nam Cochran-Armitage trend-test noncentrality
parameter on the resulting called-genotype probabilities.

## Usage

``` r
cc_ngs_power(
  N_case,
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

- N_case:

  Numeric \\\> 0\\. Number of cases.

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

  Numeric \\\> 0\\. Control-to-case sample-size ratio
  \\N\_{ctrl}/N\_{case}\\.

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

Invisibly, an object of class `"cc_ngs_power"` containing the design
inputs, trend scores, NCP and power, model information, true and called
genotype frequencies, and the true-to-called transition matrix.

## Details

Trend scores are selected from `MOI`: `"M"` uses `c(0,1,2)`, `"D"` uses
`c(0,1,1)`, and `"Rec"` uses `c(0,0,1)`. Power is the upper-tail
probability beyond the central one-degree-of-freedom chi-square critical
value under a noncentral chi-square distribution with the calculated
NCP.

Locus heterogeneity uses the same parameterization as ordinary PAWEH
case-control design: \$\$g\_{case,H} = \pi g\_{case} +
(1-\pi)g\_{control},\$\$ with the control distribution unchanged. This
biological mixture is applied to true genotype probabilities before
sequencing observation. Because the current sequencing-error model is
nondifferential, the same transition matrix is applied to cases and
controls, so mixing and sequencing commute by matrix linearity. This
identity does not extend automatically to future differential
case/control sequencing-error models. At \\\pi=0\\, the case-control
contrast and NCP are zero, and asymptotic power equals `alpha`.

This is an analytic study-design calculation applied to
sequencing-derived called genotypes. It is not a raw-read likelihood or
EM analysis, performs no downstream association testing, and uses no
simulation. Even when `seq_error = 0`, finite depth can cause call
uncertainty because a true heterozygote can yield reads from only one
allele.

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
cc_ngs_power(
  N_case = 1000, alpha = 0.05,
  prev = 0.05, pd = 0.30, R2 = 1.8,
  coverage = 20, seq_error = 0.01,
  MOI = "M", verbose = FALSE
)
```
