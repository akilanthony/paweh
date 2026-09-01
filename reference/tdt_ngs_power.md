# Analytic Power for a TDT1-NGS Sequencing Study

Computes prospective analytic power for a single-variant TDT1-NGS study
of complete father-mother-affected-child trios. The implementation uses
the published latent-state, sequencing-read-count likelihood rather than
hard genotype calls.

## Usage

``` r
tdt_ngs_power(N, pd, R1, coverage, seq_error, alpha = 0.05, verbose = TRUE)
```

## Arguments

- N:

  A single finite integer greater than or equal to 1. Number of complete
  father-mother-affected-child trios.

- pd:

  A single finite disease/risk-allele frequency strictly between 0
  and 1. The null information is reflection-symmetric in `pd` and
  `1 - pd`; the supplied allele labeling is retained.

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

Invisibly, an object of class `"tdt_ngs_power"` containing the design
inputs, power and NCP, multiplicative-model parameters, efficient
information, the 11 by 11 information matrix, compact numerical
diagnostics, and model metadata.

## Details

TDT1-NGS is evaluated for one biallelic variant under Hardy-Weinberg
parental genotype frequencies, random mating, a multiplicative disease
model, equal fixed coverage, and symmetric public sequencing error. With
\\t = R_1/(1+R_1)\\, the transmission parameter is \\\delta =
\log\\t/(1-t)\\ = \log(R_1)\\ and \\R_2 = R_1^2\\.

Kim's Appendix B noncentrality parameter is \$\$\lambda = N \delta^2
I\_{eff},\$\$ where \\I\_{eff}\\ is the nuisance-adjusted per-trio
information evaluated under the null from raw sequencing read-count
probabilities and the 15 latent Mendelian trio states. Power is the
upper-tail probability beyond the central one-degree-of-freedom
chi-square critical value under a noncentral chi-square distribution
with NCP \\\lambda\\.

Coverage 1 is unsupported under the full published nuisance model. Its
eight observable read-count triples provide at most seven independent
probability dimensions for an 11-parameter information model, so
efficient information is not identifiable without changing the model or
using a generalized inverse.

This prospective calculation performs no simulation or EM fitting. It
does not implement TDT2-NGS, unequal member-specific coverage, locus
heterogeneity, phenotype misclassification, or multi-locus testing.

## References

Kim, W. (2015). Transmission disequilibrium tests based on read counts
for low-coverage next-generation sequence data. *Human Heredity*, 80(1),
36–49. [doi:10.1159/000434645](https://doi.org/10.1159/000434645) .

## Examples

``` r
tdt_ngs_power(
  N = 5000, pd = 0.325, R1 = 1.2,
  coverage = 12, seq_error = 0.005,
  alpha = 5e-8, verbose = FALSE
)
```
