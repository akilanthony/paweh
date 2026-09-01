# Analytic MSSN for a TDT1-NGS Sequencing Study

Computes the prospective minimum sample size necessary (MSSN), measured
in complete father-mother-affected-child trios, for the same
single-variant TDT1-NGS design as
[`tdt_ngs_power()`](https://akilanthony.github.io/paweh/reference/tdt_ngs_power.md).

## Usage

``` r
tdt_ngs_mssn(power, pd, R1, coverage, seq_error, alpha = 0.05, verbose = TRUE)
```

## Arguments

- power:

  A single finite target power strictly between 0 and 1.

- pd:

  A single finite disease/risk-allele frequency strictly between 0
  and 1. The supplied allele labeling is retained.

- R1:

  A single finite positive heterozygote genotype relative risk under the
  multiplicative model. The homozygote relative risk is \\R_2 = R_1^2\\.

- coverage:

  A single finite integer greater than or equal to 2. This is the equal
  fixed read depth for the father, mother, and affected child.

- seq_error:

  A single finite symmetric per-read sequencing-error probability in
  \\\[0,0.5)\\. Internally, the directional error parameters are
  evaluated at \\\epsilon_0 = \epsilon_1 =\\ `seq_error`.

- alpha:

  A single finite significance level in \\(0,1)\\. Defaults to 0.05.

- verbose:

  Logical scalar. If `TRUE`, print a concise result summary.

## Value

Invisibly, an object of class `"tdt_ngs_mssn"` containing the target,
continuous and integer trio requirements, achieved power and NCP,
per-trio NCP coefficient, multiplicative-model parameters, efficient
information, the 11 by 11 information matrix, numerical diagnostics, and
model metadata.

## Details

The calculation uses Kim's Appendix B TDT1-NGS noncentrality parameter
\$\$\lambda = N \delta^2 I\_{eff},\$\$ where \\N\\ is the number of
complete trios, \\I\_{eff}\\ is the nuisance-adjusted per-trio
information under the null, and under the multiplicative model \\\delta
= \log(R_1)\\. The target one-degree-of- freedom NCP, \\\lambda\_\*\\,
is obtained by numerical inversion of the noncentral chi-square power
function. Trio MSSN is then solved analytically: \$\$N\_{continuous} =
\frac{\lambda\_\*}{\log(R_1)^2 I\_{eff}}.\$\$ The planned MSSN is the
ceiling of this quantity, with achieved power and the immediately
smaller design checked at the integer boundary.

If target `power` is no greater than `alpha`, the target NCP is zero and
the minimum supported design is one trio. If `R1 = 1` and target power
exceeds alpha, no finite MSSN exists because the transmission effect is
zero.

The method uses equal fixed coverage and symmetric public sequencing
error. Coverage 1 is unsupported because efficient information is not
identifiable under the implemented 11-parameter nuisance model. This
prospective calculation uses raw sequencing read-count information and
latent trio genotype states; it performs no simulation, genotype
calling, or EM fitting.

## References

Kim, W. (2015). Transmission disequilibrium tests based on read counts
for low-coverage next-generation sequence data. *Human Heredity*, 80(1),
36–49. [doi:10.1159/000434645](https://doi.org/10.1159/000434645) .

## See also

[`tdt_ngs_power`](https://akilanthony.github.io/paweh/reference/tdt_ngs_power.md)

## Examples

``` r
tdt_ngs_mssn(
  power = 0.80, pd = 0.325, R1 = 1.2,
  coverage = 12, seq_error = 0.005,
  alpha = 5e-8, verbose = FALSE
)
```
