# Plot One-Way ANOVA Minimum Sample Size Under the Falconer Model

Sweeps one model or design parameter and repeatedly calls
[`qtl_anova_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_anova_mssn.md).

## Usage

``` r
plot_qtl_anova_mssn(
  x_var,
  x_values,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  return_data = FALSE,
  ...
)
```

## Arguments

- x_var:

  Parameter to vary. One of `"power"`, `"alpha"`, `"qtl_var"`, `"tau"`,
  or `"pd"`.

- x_values:

  Numeric vector of at least two finite x-axis values.

- title:

  Optional plot-title override.

- x_label:

  Optional x-axis-label override.

- y_label:

  Optional y-axis-label override.

- return_data:

  Logical. If `TRUE`, returns plotting data rather than a ggplot object.

- ...:

  Fixed arguments passed to
  [`qtl_anova_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_anova_mssn.md).

## Value

A `ggplot` object, or a data frame with the swept parameter and
`required_N` if `return_data = TRUE`.

## Details

The x-axis is `x_values` for the selected `x_var`; the y-axis is the
minimum total sample size returned by
[`qtl_anova_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_anova_mssn.md).
Target power and other arguments in `...` remain fixed unless selected
as `x_var`.

## See also

[`qtl_anova_mssn`](https://akilanthony.github.io/paweh/reference/qtl_anova_mssn.md),
[`plot_qtl_anova_power`](https://akilanthony.github.io/paweh/reference/plot_qtl_anova_power.md).

## Examples

``` r
plot_qtl_anova_mssn(
  x_var = "qtl_var", x_values = c(0.015, 0.025, 0.05),
  power = 0.8, alpha = 0.0001, tau = 0.5, pd = 0.15
)

```
