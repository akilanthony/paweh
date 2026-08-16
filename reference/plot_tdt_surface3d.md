# Plot a Generalized 3D TDT Sensitivity Surface

Builds an interactive sensitivity surface for model-based TDT power or
minimum sample size. Any two supported model parameters can be swept
while the remaining parameters are held fixed. Every grid point is
evaluated by the canonical
[`tdt_power()`](https://akilanthony.github.io/pawh/reference/tdt_power.md)
or
[`tdt_mssn()`](https://akilanthony.github.io/pawh/reference/tdt_mssn.md)
backend.

## Usage

``` r
plot_tdt_surface3d(
  metric = c("power", "mssn"),
  scenario = c("misclassification", "heterogeneity", "no_error"),
  x = "pd",
  y = "misclass_rate",
  x_values = NULL,
  y_values = NULL,
  N = 600,
  target_power = 0.8,
  pd = 0.3,
  prev = 0.05,
  R1 = 1.5,
  R2 = 2.25,
  alpha = 0.05,
  delta_prime = 1,
  misclass_rate = 0.01,
  heter_rate = 0.1,
  ceiling_N = TRUE,
  title = NULL
)
```

## Arguments

- metric:

  Character. Surface height: `"power"` (default) or `"mssn"`.

- scenario:

  Character. One of `"misclassification"` (default), `"heterogeneity"`,
  or `"no_error"`.

- x, y:

  Character scalars naming two distinct swept parameters. Supported axes
  are `"pd"`, `"prev"`, `"R1"`, `"R2"`, `"alpha"`, `"delta_prime"`,
  `"misclass_rate"`, and `"heter_rate"`.

- x_values, y_values:

  Numeric grid values for `x` and `y`. Each must contain at least two
  distinct, finite values. If `NULL`, documented parameter-specific
  defaults are used; see Details.

- N:

  Numeric scalar greater than zero. Number of affected trios, used when
  `metric = "power"`.

- target_power:

  Numeric scalar in `(0, 1)`, used when `metric = "mssn"`.

- pd:

  Numeric in `(0, 1)`. Fixed risk-allele frequency when `pd` is not an
  axis.

- prev:

  Numeric in `(0, 1)`. Fixed disease prevalence when `prev` is not an
  axis.

- R1, R2:

  Positive numeric scalars. Fixed genotype relative risks when they are
  not axes.

- alpha:

  Numeric in `(0, 1)`. Fixed significance level when `alpha` is not an
  axis.

- delta_prime:

  Numeric in `[0, 1]`. Fixed positive linkage-disequilibrium scale when
  `delta_prime` is not an axis.

- misclass_rate:

  Numeric in `[0, 1)`. Fixed phenotype misclassification rate when it is
  not an axis.

- heter_rate:

  Numeric in `[0, 1)`. Fixed locus-heterogeneity rate when it is not an
  axis.

- ceiling_N:

  Logical scalar. For MSSN surfaces, whether to plot the integer ceiling
  of the canonical sample-size result (default `TRUE`). The unrounded
  result remains available in the attached surface data.

- title:

  Optional character scalar. If `NULL`, an informative title is
  generated from `metric` and `scenario`.

## Value

An object inheriting from `"plotly"` and `"htmlwidget"`. The Cartesian
grid is attached as `attr(x, "surface_data")`, with dynamically named
axis columns, `metric_value` (the plotted height), and
`raw_metric_value`. The plotted z matrix is attached as
`attr(x, "surface_matrix")`, and calculation settings as
`attr(x, "surface_spec")`.

## Details

The supported axes and their default grids are:
`pd = seq(0.10, 0.50, length.out = 20)`,
`prev = seq(0.01, 0.10, length.out = 20)`,
`R1 = seq(1.1, 2, length.out = 20)`,
`R2 = seq(1.2, 3, length.out = 20)`,
`alpha = seq(0.01, 0.10, length.out = 20)`,
`delta_prime = seq(0.25, 1, length.out = 20)`,
`misclass_rate = seq(0, 0.20, length.out = 20)`, and
`heter_rate = seq(0, 0.50, length.out = 20)`. The package's implemented
LD convention is `D = delta_prime * pd * (1 - pd)`, where `delta_prime`
is the proportion of maximum positive disequilibrium: `0` represents
linkage equilibrium and `1` represents maximum positive LD under the
model assumptions. Negative LD would require a different,
allele-frequency-dependent normalization and is not represented by this
parameterization.

The surface is deliberately model-based. `N` and `target_power` define
the design calculation and therefore remain fixed rather than becoming
axes. The model-free inputs `ET`, `ENT`, and `n_trios`, as well as
nonnumeric controls such as `verbose`, are outside this visualization
interface. MSSN is measured in affected-child trios. `heter_rate` is the
heterogeneous fraction, `1 - pi`.

The canonical TDT backends report separate scenarios, not a
combined-error scenario. Consequently, `misclass_rate` is an active axis
only for `scenario = "misclassification"`, and `heter_rate` is active
only for `scenario = "heterogeneity"`. An inactive modifier axis is
rejected rather than silently producing a flat, misleading surface.

This function requires the optional `plotly` package. The plotting layer
introduces no independent statistical approximation: it only constructs
the grid, delegates calculations to
[`tdt_power()`](https://akilanthony.github.io/pawh/reference/tdt_power.md)
or
[`tdt_mssn()`](https://akilanthony.github.io/pawh/reference/tdt_mssn.md),
extracts the selected scenario, and renders the result.

## See also

[`tdt_power()`](https://akilanthony.github.io/pawh/reference/tdt_power.md),
[`tdt_mssn()`](https://akilanthony.github.io/pawh/reference/tdt_mssn.md),
[`plot_tdt_power()`](https://akilanthony.github.io/pawh/reference/plot_tdt_power.md),
[`plot_tdt_mssn()`](https://akilanthony.github.io/pawh/reference/plot_tdt_mssn.md)

## Examples

``` r
if (requireNamespace("plotly", quietly = TRUE)) {
  power_surface <- plot_tdt_surface3d(
    metric = "power", scenario = "misclassification",
    x = "pd", y = "misclass_rate",
    x_values = c(0.2, 0.3, 0.4),
    y_values = c(0, 0.03, 0.06),
    N = 600
  )

  mssn_surface <- plot_tdt_surface3d(
    metric = "mssn", scenario = "heterogeneity",
    x = "prev", y = "heter_rate",
    x_values = c(0.02, 0.05, 0.08),
    y_values = c(0, 0.1, 0.2),
    target_power = 0.8
  )
}
```
