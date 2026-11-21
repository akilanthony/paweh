# Analytic Power for a Case-Control Sequencing Study

Computes prospective asymptotic power for a model-based case-control
sequencing design. Sequencing is the common reference mechanism.
Optional phenotype misclassification, locus heterogeneity, and
conventional genotype misclassification are evaluated as independent
sensitivity scenarios from the same biological baseline and sequencing
design.

## Usage

``` r
cc_ngs_power(
  N_case,
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

Invisibly, an object of class `"cc_ngs_power"`. The `scenarios` list
always contains `sequencing_only` and each explicitly requested modifier
scenario. Every scenario records staged frequencies, sequencing and
error matrices, Ahn components `D` and `Q`, NCP, and power. Legacy
top-level result fields mirror sequencing only with zero or multiple
optional modifiers and the requested modifier with exactly one;
`compatibility_scenario` identifies it. The frequency list distinguishes
post-sequencing calls from final post-genotype-misclassification values.
The `sequencing` list records the four effective group parameters and
`case_transition_matrix`/`ctrl_transition_matrix`. Legacy
`coverage`/`seq_error` retain supplied common inputs (or NULL);
`transition_matrix` aliases the case matrix, shared when mechanisms are
equal. Use the two explicit matrices for differential designs.

## Details

Trend scores are selected from `MOI`: `"M"` uses `c(0,1,2)`, `"D"` uses
`c(0,1,1)`, and `"Rec"` uses `c(0,0,1)`. Power is the upper-tail
probability beyond the central one-degree-of-freedom chi-square critical
value under a noncentral chi-square distribution with the calculated
NCP.

Locus heterogeneity uses the same parameterization as ordinary PAWEH
case-control design: \$\$g\_{case,H} = \pi g\_{case} +
(1-\pi)g\_{control},\$\$ with the control distribution unchanged. The
scenario paths are: sequencing only, baseline to sequencing to test;
phenotype, baseline to phenotype misclassification to sequencing to
test; heterogeneity, baseline to locus heterogeneity to sequencing to
test; and genotype, baseline to sequencing to conventional genotype
misclassification to test. Supplying multiple modifiers returns multiple
independent scenarios; modifiers are not combined. The same effective
case/control coverage and sequencing-error settings are used in every
scenario. Phenotype misclassification uses the same prevalence-weighted
transformation as
[`cc_power`](https://akilanthony.github.io/paweh/reference/cc_power.md),
with `theta = Pr(affected -> control)` and
`phi = Pr(unaffected -> case)`. The resulting observed-case and
observed-control distributions enter sequencing. Separate row-true,
column-called matrices are then constructed using each observed group's
effective depth and error, and applied as `t(M_case) %*% g_case` and
`t(M_ctrl) %*% g_control`. The ordinary PAWEH 1p, 2p, 3p, or
differential 3p genotype-error matrices are applied only in their
independent scenario, after sequencing calls, with the same row-current,
column-observed convention. Only equal sequencing mechanisms permit
commuting a shared biological mixture and sequencing. At \\\pi=0\\,
equal mechanisms with no phenotype-induced difference give zero NCP and
power equal to `alpha`.

Differential sequencing or genotype-misclassification mechanisms can
create observed case/control differences even under a biological null.
Results remain nominal asymptotic calculations using the existing
chi-square critical value; the null distribution and Type I error are
not separately recalibrated here.

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
