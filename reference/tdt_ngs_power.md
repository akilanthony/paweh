# Analytic Power for a TDT1-NGS Sequencing Study

Computes prospective analytic power for a single-variant TDT1-NGS study
of complete father-mother-affected-child trios. The implementation uses
the published latent-state, sequencing-read-count likelihood rather than
hard genotype calls.

## Usage

``` r
tdt_ngs_power(
  N,
  pd,
  R1,
  coverage = NULL,
  seq_error = NULL,
  alpha = 0.05,
  verbose = TRUE,
  father_coverage = coverage,
  mother_coverage = coverage,
  child_coverage = coverage,
  epsilon0 = seq_error,
  epsilon1 = seq_error,
  pheno_misclass = FALSE,
  prev = NULL,
  pi01 = 0,
  locus_het = FALSE,
  heter_rate = 0
)
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

  A single finite integer greater than or equal to 2. Common fixed-depth
  shorthand for the father, mother, and affected child. It may be
  omitted when all three member-specific depths are supplied.

- seq_error:

  A single finite symmetric per-read sequencing-error probability in
  \\\[0,0.5)\\. Internally, the directional error parameters default to
  \\\epsilon_0 = \epsilon_1 =\\ `seq_error`. It may be omitted when both
  directional parameters are supplied.

- alpha:

  A single finite significance level in \\(0,1)\\. Defaults to 0.05.

- verbose:

  Logical scalar. If `TRUE`, print a concise result summary.

- father_coverage, mother_coverage, child_coverage:

  Optional fixed depths for each trio member. An explicitly supplied
  member depth overrides `coverage`; otherwise that member inherits
  `coverage`.

- epsilon0:

  Directional probability that a reference-allele read is observed as
  the alternative allele. Defaults to `seq_error`.

- epsilon1:

  Directional probability that an alternative-allele read is observed as
  the reference allele. Defaults to `seq_error`.

- pheno_misclass:

  Logical scalar. If `TRUE`, add a phenotype- misclassification
  sensitivity scenario using the ordinary-TDT NCP bridge.

- prev:

  Disease prevalence in \\(0,1)\\. Required when either phenotype
  misclassification or locus heterogeneity is requested. It is used only
  by the ordinary-TDT bridges and does not alter the raw-read
  likelihood.

- pi01:

  Probability that an unaffected individual is misclassified or
  ascertained as affected, in \\\[0,1)\\. It must be zero unless
  `pheno_misclass = TRUE`.

- locus_het:

  Logical scalar. If `TRUE`, add an independent locus- heterogeneity
  sensitivity scenario using the ordinary-TDT NCP bridge.

- heter_rate:

  Heterogeneous trio fraction in \\\[0,1\]\\. The ordinary- TDT
  linked/homogeneous fraction is \\\pi=1-\\`heter_rate`. It must be zero
  unless `locus_het = TRUE`.

## Value

Invisibly, an object of class `"tdt_ngs_power"` containing the design
inputs, power and NCP, multiplicative-model parameters, efficient
information, the 11 by 11 information matrix, compact numerical
diagnostics, effective member-specific coverage and directional-error
metadata, model metadata, and a `scenarios` list that always contains
`sequencing_only` and, when requested, `phenotype_misclassification` and
`heterogeneity`. Modifier scenarios record their ordinary-TDT bridge
inputs, coefficients, attenuation factor, and confirmation that the
raw-read likelihood was not modified. Legacy `coverage` and `seq_error`
fields retain the supplied common shorthand values.

## Details

TDT1-NGS is evaluated for one biallelic variant under Hardy-Weinberg
parental genotype frequencies, random mating, a multiplicative disease
model, fixed member-specific coverage, and common directional sequencing
error. With \\t = R_1/(1+R_1)\\, the transmission parameter is \\\delta
= \log\\t/(1-t)\\ = \log(R_1)\\ and \\R_2 = R_1^2\\.

Kim's Appendix B noncentrality parameter is \$\$\lambda = N \delta^2
I\_{eff},\$\$ where \\I\_{eff}\\ is the nuisance-adjusted per-trio
information evaluated under the null from raw sequencing read-count
probabilities and the 15 latent Mendelian trio states. Power is the
upper-tail probability beyond the central one-degree-of-freedom
chi-square critical value under a noncentral chi-square distribution
with NCP \\\lambda\\.

When requested, phenotype misclassification is a sensitivity bridge
rather than a phenotype-aware sequencing likelihood. Let
\\c=2(g_T-g\_{NT})^2/(g_T+g\_{NT})\\ be the ordinary-TDT per-trio NCP
coefficient. The bridge computes
\$\$r\_{pheno}=c\_{pheno}/c\_{noerror}\$\$ from the ordinary-TDT
expected transmission and nontransmission probabilities, using the fixed
assumptions \\R_2=R_1^2\\, \\\delta'=1\\, and \\\theta_1=p_d\\. It then
reports \\\lambda\_{pheno}=\lambda\_{NGS}r\_{pheno}\\. Kim's raw-read
likelihood, information matrix, and efficient information remain
unchanged.

Locus heterogeneity uses the analogous independent bridge
\$\$r\_{het}=c\_{heterogeneity}/c\_{noerror},\qquad
\lambda\_{het}=\lambda\_{NGS}r\_{het}.\$\$ The ordinary-TDT expected
counts are evaluated at one trio using \\\pi=1\\ for the homogeneous
reference and \\\pi=1-\\`heter_rate` for the requested sensitivity
scenario. It uses the same fixed \\R_2\\, \\\delta'\\, and \\\theta_1\\
assumptions as the phenotype bridge. If both modifiers are requested,
they are reported separately against the same sequencing-only baseline;
their factors are never multiplied and no combined result is calculated.

The existing directional read model is \$\$q_G=\epsilon_0 +
(1-\epsilon_0-\epsilon_1)G/2,\$\$ for genotype \\G\in\\0,1,2\\\\. Thus
\\q_0=\epsilon_0\\, \\q_2=1-\epsilon_1\\, and
\\q_1=(1+\epsilon_0-\epsilon_1)/2\\; the heterozygote probability is
one-half in the symmetric case. The effective errors must be nonnegative
and sum to less than one.

Coverage 1 is unsupported under the full published nuisance model. Its
eight observable read-count triples provide at most seven independent
probability dimensions for an 11-parameter information model, so
efficient information is not identifiable without changing the model or
using a generalized inverse.

This prospective calculation performs no simulation or EM fitting. It
does not implement TDT2-NGS, conventional genotype
misclassification/TDTae, a joint modifier likelihood, or multi-locus
testing. Coverage is fixed for each member rather than random or
sample-specific.

## References

Kim, W. (2015). Transmission disequilibrium tests based on read counts
for low-coverage next-generation sequence data. *Human Heredity*, 80(1),
36–49. [doi:10.1159/000434645](https://doi.org/10.1159/000434645) .

Buyske, S., Yang, G., Matise, T. C., & Gordon, D. (2009). When a case is
not a case: Effects of phenotype misclassification on power and sample
size requirements for the transmission disequilibrium test with affected
child trios. *Human Heredity*, 67(4), 287–292.
[doi:10.1159/000194981](https://doi.org/10.1159/000194981) .

Chen, C., Yang, G., Buyske, S., Matise, T., Finch, S. J., & Gordon, D.
(2009). Transmission disequilibrium test power and sample size in the
presence of locus heterogeneity. *Statistical Applications in Genetics
and Molecular Biology*, 8, Article 44.
[doi:10.2202/1544-6115.1501](https://doi.org/10.2202/1544-6115.1501) .

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

## Examples

``` r
tdt_ngs_power(
  N = 5000, pd = 0.325, R1 = 1.2,
  coverage = 12, seq_error = 0.005,
  alpha = 5e-8, verbose = FALSE
)
tdt_ngs_power(
  N = 5000, pd = 0.325, R1 = 1.2, coverage = 12, seq_error = 0.005,
  locus_het = TRUE, prev = 0.01, heter_rate = 0.25,
  alpha = 5e-8, verbose = FALSE
)
tdt_ngs_power(
  N = 5000, pd = 0.325, R1 = 1.2, coverage = 12, seq_error = 0.005,
  pheno_misclass = TRUE, prev = 0.01, pi01 = 0.05,
  alpha = 5e-8, verbose = FALSE
)
```
