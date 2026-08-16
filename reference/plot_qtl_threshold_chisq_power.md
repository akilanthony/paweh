# Plot Threshold-Selected Genotype Chi-Square Power

Sweeps one model, selection, or design parameter and repeatedly calls
[`qtl_threshold_chisq_power()`](https://akilanthony.github.io/pawh/reference/qtl_threshold_chisq_power.md).

## Usage

``` r
plot_qtl_threshold_chisq_power(
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

  Parameter to vary. One of `"N_case"`, `"alpha"`, `"qtl_var"`, `"tau"`,
  `"pd"`, `"x_upper"`, `"x_lower"`, or `"k"`.

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
  [`qtl_threshold_chisq_power()`](https://akilanthony.github.io/pawh/reference/qtl_threshold_chisq_power.md).

## Value

A `ggplot` object, or a data frame with the swept parameter and power if
`return_data = TRUE`.

## Details

The x-axis is `x_values` for the selected model, selection, or design
parameter. The y-axis is power for the two-df genotype chi-square test
after upper-tail case and lower-tail control selection. Arguments in
`...` remain fixed while `x_var` is swept.

## See also

[`qtl_threshold_chisq_power`](https://akilanthony.github.io/pawh/reference/qtl_threshold_chisq_power.md),
[`plot_qtl_threshold_chisq_mssn`](https://akilanthony.github.io/pawh/reference/plot_qtl_threshold_chisq_mssn.md).

## Examples

``` r
plot_qtl_threshold_chisq_power(
  x_var = "N_case", x_values = c(75, 100, 125, 150),
  alpha = 0.0001, qtl_var = 0.025, tau = 0.5, pd = 0.15,
  x_upper = 5, x_lower = 5, k = 1
)

```
