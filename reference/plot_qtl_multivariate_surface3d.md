# Plot an Interactive 3D Two-Phenotype Falconer Surface

Creates an interactive Plotly surface for the population mixture of
three genotype-specific bivariate-normal distributions. The plot
displays exactly two quantitative phenotypes and complements
[`plot_qtl_multivariate_contour()`](https://akilanthony.github.io/pawh/reference/plot_qtl_multivariate_contour.md),
which provides a static 2D view.

## Usage

``` r
plot_qtl_multivariate_surface3d(
  qtl_var,
  tau,
  pd,
  cor_matrix,
  x_upper = NULL,
  x_lower = NULL,
  surface = c("density", "cdf"),
  show_means = TRUE,
  show_thresholds = TRUE,
  show_labels = TRUE,
  grid_n = 80L,
  z_scale = c("raw", "normalized"),
  title = NULL
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

- show_means:

  Logical; show the three genotype-specific mean vectors at their
  corresponding mixture-surface heights.

- show_thresholds:

  Logical; draw the two joint threshold rectangles on the `z = 0` base
  plane when thresholds are supplied.

- show_labels:

  Logical; label the affected and unaffected rectangles.

- grid_n:

  Integer number of grid points along each phenotype axis.

- z_scale:

  Either `"raw"` for the calculated surface values or `"normalized"` to
  linearly rescale the displayed heights to `[0, 1]`.

- title:

  Optional plot title.

## Value

An object inheriting from `"plotly"` and `"htmlwidget"`. Validated model
quantities, raw and displayed surface data, thresholds, and mean-marker
data are attached as attributes.

## Details

Density mode plots the genotype-weighted mixture
`sum(pi[j] * f[j](y1, y2))`, where each `f[j]` is a bivariate-normal
density. CDF mode plots the lower-tail mixture
`sum(pi[j] * P(Y1 <= y1, Y2 <= y2 | G = j))`. Thus, the density surface
has genotype-related peaks, whereas the CDF surface is cumulative and
monotone. Normalizing the z axis changes only its displayed scale; the
raw calculated values remain available in the plot's `plot_data`
attribute.

Affected subjects are defined only by the joint upper-right region
`Y1 >= TU1 AND Y2 >= TU2`. Unaffected subjects are defined only by the
joint lower-left region `Y1 <= TL1 AND Y2 <= TL2`. Other phenotype
combinations are not selected. The statistical backend can handle more
traits where documented, but a surface over phenotype space is
restricted here to two.

This function requires the optional package plotly.

## See also

[`plot_qtl_multivariate_contour()`](https://akilanthony.github.io/pawh/reference/plot_qtl_multivariate_contour.md),
[`qtl_multivariate_power_full()`](https://akilanthony.github.io/pawh/reference/qtl_multivariate_power_full.md),
and
[`qtl_multivariate_mssn_full()`](https://akilanthony.github.io/pawh/reference/qtl_multivariate_mssn_full.md).

## Examples

``` r
if (requireNamespace("plotly", quietly = TRUE)) {
  cor_matrix <- matrix(c(1, 0.15, 0.15, 1), 2, byrow = TRUE)
  p <- plot_qtl_multivariate_surface3d(
    qtl_var = c(0.1, 0.08), tau = c(0, 0.5), pd = 0.3,
    cor_matrix = cor_matrix,
    x_upper = c(10, 10), x_lower = c(15, 15),
    surface = "density", grid_n = 20
  )
}
```
