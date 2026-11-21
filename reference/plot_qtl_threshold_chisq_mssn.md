# Plot Threshold-Selected Genotype Chi-Square Minimum Sample Size

Sweeps one model, selection, or design parameter and repeatedly calls
[`qtl_threshold_chisq_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_mssn.md).

## Usage

``` r
plot_qtl_threshold_chisq_mssn(
  x_var,
  x_values,
  sample_size = c("total", "case", "control"),
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
  `"pd"`, `"x_upper"`, `"x_lower"`, or `"k"`.

- x_values:

  Numeric vector of at least two finite x-axis values.

- sample_size:

  Selected-sample result to plot: `"total"`, `"case"`, or `"control"`.
  Population screening counts are deliberately excluded because they are
  not statistical MSSN values.

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
  [`qtl_threshold_chisq_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_mssn.md).

## Value

A `ggplot` object, or a data frame with the swept parameter and the
selected MSSN result if `return_data = TRUE`.

## Details

The x-axis is `x_values` for the selected model, selection, or design
parameter. The y-axis is the selected case, control, or total MSSN
requested by `sample_size`; it is not a source-population screening
count. Other arguments remain fixed while `x_var` is swept.

## See also

[`qtl_threshold_chisq_mssn`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_mssn.md),
[`plot_qtl_threshold_chisq_power`](https://akilanthony.github.io/paweh/reference/plot_qtl_threshold_chisq_power.md).

## Examples

``` r
plot_qtl_threshold_chisq_mssn(
  x_var = "qtl_var", x_values = c(0.015, 0.025, 0.05),
  power = 0.8, alpha = 0.0001, tau = 0.5, pd = 0.15,
  x_upper = 5, x_lower = 5, k = 1
)

```
