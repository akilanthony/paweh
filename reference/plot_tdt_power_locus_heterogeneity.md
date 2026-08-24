# Plot TDT Power Sensitivity to Locus Heterogeneity

Plot TDT power vs locus heterogeneity rate (1 - pi), holding
misclassification fixed.

## Usage

``` r
plot_tdt_power_locus_heterogeneity(
  pd,
  prev,
  R1,
  R2,
  alpha = 0.05,
  delta_prime = 1,
  N,
  heter_seq = seq(0, 0.5, by = 0.05),
  misclass_fixed = 0,
  title = "TDT power vs heterogeneity (1 - pi)"
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

- N:

  Number of affected trios.

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
homogeneous/linked fraction is \\\pi\\. The y-axis is
heterogeneity-scenario power from
[`tdt_power()`](https://akilanthony.github.io/paweh/reference/tdt_power.md).
`N`, `misclass_fixed`, and all genetic-model parameters remain fixed
while `heter_seq` is swept.

## See also

[`tdt_power`](https://akilanthony.github.io/paweh/reference/tdt_power.md),
[`plot_tdt_power_phenotype_misclassification`](https://akilanthony.github.io/paweh/reference/plot_tdt_power_phenotype_misclassification.md),
and
[`plot_tdt_power`](https://akilanthony.github.io/paweh/reference/plot_tdt_power.md).

## Examples

``` r
# Power sensitivity to locus heterogeneity
plot_tdt_power_locus_heterogeneity(
  N    = 600,
  pd   = 0.30,
  prev = 0.05,
  R1   = 1.5,
  R2   = 2.25,
  heter_seq      = c(0, 0.1, 0.2),
  misclass_fixed = 0
)

```
