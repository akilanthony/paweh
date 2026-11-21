# Analytic MSSN for a TDT1-NGS Sequencing Study

Computes the prospective minimum sample size necessary (MSSN), measured
in complete father-mother-affected-child trios, for the same
single-variant TDT1-NGS design as
[`tdt_ngs_power()`](https://akilanthony.github.io/paweh/reference/tdt_ngs_power.md).

## Usage

``` r
tdt_ngs_mssn(
  power,
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

- power:

  A single finite target power strictly between 0 and 1.

- pd:

  A single finite disease/risk-allele frequency strictly between 0
  and 1. The supplied allele labeling is retained.

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

Invisibly, an object of class `"tdt_ngs_mssn"` containing the target,
continuous and integer trio requirements, achieved power and NCP,
per-trio NCP coefficient, multiplicative-model parameters, efficient
information, the 11 by 11 information matrix, numerical diagnostics,
effective sequencing-design metadata, model metadata, and scenario-level
base and independently adjusted coefficients and sample sizes. Legacy
`coverage` and `seq_error` fields retain supplied shorthand.

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
smaller design checked at the integer boundary. For the optional
phenotype scenario, the ordinary-TDT attenuation factor \\r\_{pheno}\\
described in
[`tdt_ngs_power`](https://akilanthony.github.io/paweh/reference/tdt_ngs_power.md)
multiplies the per-trio coefficient before inversion:
\$\$N\_{continuous,pheno}= \frac{\lambda\_\*}{\log(R_1)^2
I\_{eff}r\_{pheno}}.\$\$ The bridge fixes \\R_2=R_1^2\\, \\\delta'=1\\,
and \\\theta_1=p_d\\. It does not modify Kim's likelihood or information
matrix and is not a TDTae or simulation model.

For locus heterogeneity, \\r\_{het}\\ is the ratio of the ordinary-TDT
heterogeneous and homogeneous per-trio NCP coefficients, with
\\\pi=1-\\`heter_rate`. The adjusted coefficient is \\\log(R_1)^2
I\_{eff}r\_{het}\\ before inversion. Phenotype and heterogeneity
requests remain separate scenarios; no product of their attenuation
factors is used.

If target `power` is no greater than `alpha`, the target NCP is zero and
the minimum supported design is one trio. If `R1 = 1` and target power
exceeds alpha, no finite MSSN exists because the transmission effect is
zero.

Each member uses its resolved fixed coverage, and the directional read
model is \\q_G=\epsilon_0+(1-\epsilon_0-\epsilon_1)G/2\\. See
[`tdt_ngs_power`](https://akilanthony.github.io/paweh/reference/tdt_ngs_power.md)
for parameter precedence and interpretation. Coverage 1 is unsupported
for every member because efficient information is not identifiable under
the implemented 11-parameter nuisance model. This prospective
calculation uses raw sequencing read-count information and latent trio
genotype states; it performs no simulation, genotype calling, or EM
fitting.

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

## See also

[`tdt_ngs_power`](https://akilanthony.github.io/paweh/reference/tdt_ngs_power.md)

## Examples

``` r
tdt_ngs_mssn(
  power = 0.80, pd = 0.325, R1 = 1.2,
  coverage = 12, seq_error = 0.005,
  alpha = 5e-8, verbose = FALSE
)
tdt_ngs_mssn(
  power = 0.80, pd = 0.325, R1 = 1.2, coverage = 12, seq_error = 0.005,
  locus_het = TRUE, prev = 0.01, heter_rate = 0.25,
  alpha = 5e-8, verbose = FALSE
)
tdt_ngs_mssn(
  power = 0.80, pd = 0.325, R1 = 1.2, coverage = 12, seq_error = 0.005,
  pheno_misclass = TRUE, prev = 0.01, pi01 = 0.05,
  alpha = 5e-8, verbose = FALSE
)
```
