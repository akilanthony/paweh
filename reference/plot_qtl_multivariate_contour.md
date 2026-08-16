# Plot a Two-Phenotype Falconer Mixture Surface

Visualizes the population mixture of three genotype-specific bivariate
normal distributions for exactly two quantitative phenotypes. It can
show either the mixture density or the mixture lower-tail CDF, together
with the joint threshold regions used by the Chapter 6.2
selected-sampling design.

## Usage

``` r
plot_qtl_multivariate_contour(
  qtl_var,
  tau,
  pd,
  cor_matrix,
  x_upper = NULL,
  x_lower = NULL,
  surface = c("density", "cdf"),
  show_thresholds = TRUE,
  show_means = TRUE,
  show_labels = TRUE,
  grid_n = 150L,
  title = NULL,
  return_data = FALSE
)
```

## Arguments

- qtl_var:

  Numeric vector of two phenotype-specific QTL variances.

- tau:

  Numeric vector of two phenotype-specific dominance/additivity ratios.

- pd:

  Shared increaser-allele frequency in `(0, 1)`.

- cor_matrix:

  A positive-definite `2 x 2` phenotype correlation matrix.

- x_upper, x_lower:

  Optional length-two vectors of upper- and lower-tail selection
  percentages. Supply both or neither.

- surface:

  Either `"density"` or `"cdf"`.

- show_thresholds:

  Logical; display joint threshold regions when thresholds are supplied.

- show_means:

  Logical; mark the three genotype-specific mean vectors.

- show_labels:

  Logical; label the affected and unaffected joint regions.

- grid_n:

  Integer number of grid points along each phenotype axis.

- title:

  Optional plot title.

- return_data:

  Logical; return the surface grid instead of the ggplot. The validated
  model and threshold details are retained as attributes.

## Value

A `ggplot` object, or the plotting grid when `return_data = TRUE`.
Returned data retain the validated model and threshold details as
attributes.

## Details

In density mode, each grid value is the genotype-weighted mixture
`sum(pi[j] * f[j](y1, y2))`, where each `f[j]` is a bivariate-normal
density. In CDF mode, it is the lower-tail mixture probability
`sum(pi[j] * P(Y1 <= y1, Y2 <= y2 | G = j))`. Density and CDF surfaces
therefore represent different mathematical quantities.

When thresholds are supplied, affected subjects occupy only the joint
upper-right region `Y1 >= TU1 AND Y2 >= TU2`. Unaffected subjects occupy
only the joint lower-left region `Y1 <= TL1 AND Y2 <= TL2`. Subjects in
all other regions are not selected. The statistical backend can support
more traits where documented; this conventional contour visualization is
intentionally restricted to two.

## See also

[`qtl_multivariate_power_full()`](https://akilanthony.github.io/pawh/reference/qtl_multivariate_power_full.md),
[`qtl_multivariate_mssn_full()`](https://akilanthony.github.io/pawh/reference/qtl_multivariate_mssn_full.md),
and
[`plot_qtl_multivariate_surface3d()`](https://akilanthony.github.io/pawh/reference/plot_qtl_multivariate_surface3d.md).

## Examples

``` r
cor_matrix <- matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
plot_qtl_multivariate_contour(
  qtl_var = c(0.95, 0.92), tau = c(0, 0.5), pd = 0.5,
  cor_matrix = cor_matrix,
  x_upper = c(10, 10), x_lower = c(15, 15),
  surface = "cdf", grid_n = 30
)

```
