# Plot TDT MSSN Sensitivity to Locus Heterogeneity

Plot required number of trios vs heterogeneity rate (1 - pi), holding
misclassification fixed.

## Usage

``` r
plot_tdt_mssn_locus_heterogeneity(
  pd,
  prev,
  R1,
  R2,
  alpha = 0.05,
  delta_prime = 1,
  target_power,
  heter_seq = seq(0, 0.5, by = 0.05),
  misclass_fixed = 0,
  title = "Required trios vs heterogeneity (1 - pi)"
)
```

## Arguments

- pd:

  Disease allele frequency.

- prev:

  Disease prevalence.

- R1:

  Genotype relative risk for heterozygotes.

- R2:

  Genotype relative risk for risk homozygotes.

- alpha:

  Significance level.

- delta_prime:

  Linkage disequilibrium scaling parameter.

- target_power:

  Desired statistical power.

- heter_seq:

  Numeric vector. Sequence of heterogeneity rates.

- misclass_fixed:

  Numeric. Misclassification rate held fixed (default 0). (Other params
  same meaning as in plot_tdt_power_phenotype_misclassification.)

- title:

  Plot title.

## Value

A `ggplot` object, returned invisibly after it is printed.

## Details

The x-axis is the heterogeneous trio fraction \\1-\pi\\; the
homogeneous/linked fraction is \\\pi\\. The y-axis is MSSN, the required
number of affected trios returned by
[`tdt_mssn()`](https://akilanthony.github.io/pawh/reference/tdt_mssn.md)
for fixed `target_power`. `misclass_fixed` and genetic-model parameters
remain fixed.

## See also

[`tdt_mssn`](https://akilanthony.github.io/pawh/reference/tdt_mssn.md),
[`plot_tdt_mssn_phenotype_misclassification`](https://akilanthony.github.io/pawh/reference/plot_tdt_mssn_phenotype_misclassification.md),
and
[`plot_tdt_mssn`](https://akilanthony.github.io/pawh/reference/plot_tdt_mssn.md).

## Examples

``` r
# Required sample size sensitivity to locus heterogeneity
plot_tdt_mssn_locus_heterogeneity(
  target_power = 0.80,
  pd   = 0.30,
  prev = 0.05,
  R1   = 1.5,
  R2   = 2.25,
  heter_seq      = c(0, 0.1, 0.2),
  misclass_fixed = 0
)

```
