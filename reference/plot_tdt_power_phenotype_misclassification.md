# Plot TDT Power Sensitivity to Phenotype Misclassification

Plot TDT power vs phenotype misclassification rate (pi01), holding
heterogeneity fixed.

## Usage

``` r
plot_tdt_power_phenotype_misclassification(
  pd,
  prev,
  R1,
  R2,
  alpha = 0.05,
  delta_prime = 1,
  N,
  misclass_seq = seq(0, 0.15, by = 0.01),
  heter_fixed = 0,
  title = "TDT power vs misclassification (pi01)"
)
```

## Arguments

- pd:

  Numeric in (0,1). High-risk allele frequency at the marker.

- prev:

  Numeric in (0,1). Disease prevalence (phi1).

- R1, R2:

  Numeric. Genotype relative risks.

- alpha:

  Numeric. Significance level (default 0.05).

- delta_prime:

  Numeric. LD scale parameter (default 1).

- N:

  Integer. Number of affected trios.

- misclass_seq:

  Numeric vector. Sequence of misclassification rates.

- heter_fixed:

  Numeric. Heterogeneity rate held fixed (default 0).

- title:

  Character. Plot title.

## Value

A `ggplot` object, returned invisibly after it is printed.

## Details

The x-axis is phenotype misclassification probability \\\pi\_{01}\\ and
the y-axis is the corresponding power from
[`tdt_power()`](https://akilanthony.github.io/pawh/reference/tdt_power.md).
`N` is a fixed number of affected trios; `heter_fixed` and all
genetic-model parameters remain fixed while `misclass_seq` is swept.

## See also

[`tdt_power`](https://akilanthony.github.io/pawh/reference/tdt_power.md),
[`plot_tdt_power_locus_heterogeneity`](https://akilanthony.github.io/pawh/reference/plot_tdt_power_locus_heterogeneity.md),
and
[`plot_tdt_power`](https://akilanthony.github.io/pawh/reference/plot_tdt_power.md).

## Examples

``` r
# Power sensitivity to phenotype misclassification
plot_tdt_power_phenotype_misclassification(
  N    = 600,
  pd   = 0.30,
  prev = 0.05,
  R1   = 1.5,
  R2   = 2.25,
  misclass_seq = c(0, 0.03, 0.06),
  heter_fixed  = 0
)

```
