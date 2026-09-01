# Plot Case-Control Minimum Sample Size Necessary

Sweeps one x-axis parameter and repeatedly calls
[`cc_mssn()`](https://akilanthony.github.io/paweh/reference/cc_mssn.md)
to plot the minimum sample size necessary (MSSN) for case-control tests.

## Usage

``` r
plot_cc_mssn(
  x_var,
  x_values,
  test = c("genotypes", "trend"),
  input_mode = c("model_based", "model_free"),
  sample_size = c("total", "case", "control"),
  compare_tests = FALSE,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  return_data = FALSE,
  ...
)
```

## Arguments

- x_var:

  Character. Parameter to vary on the x-axis.

- x_values:

  Numeric vector of x-axis values.

- test:

  One of `"genotypes"` or `"trend"`.

- input_mode:

  One of `"model_based"` or `"model_free"`.

- sample_size:

  One of `"total"`, `"case"`, or `"control"`.

- compare_tests:

  Logical. If TRUE, plot genotype and trend tests together and ignore
  `test`.

- title:

  Optional character title override.

- x_label:

  Optional x-axis label override.

- y_label:

  Optional y-axis label override.

- return_data:

  Logical. If TRUE, return the data frame instead of a ggplot.

- ...:

  Arguments passed to
  [`cc_mssn()`](https://akilanthony.github.io/paweh/reference/cc_mssn.md).

## Value

A ggplot object, or a data frame if `return_data = TRUE`.

## Details

The supported `test` values are `"genotypes"` and `"trend"`.
`sample_size` selects whether to plot required cases, controls, or total
sample size. `compare_tests = TRUE` plots both tests together.

Supported `x_var` values are the same heterogeneity and
misclassification variables documented for
[`plot_cc_power()`](https://akilanthony.github.io/paweh/reference/plot_cc_power.md),
plus `"power"` for target power. Fixed sample-size variables such as
`"N_case"` are not valid because sample size is the MSSN output.

Multiplier behavior matches
[`plot_cc_power()`](https://akilanthony.github.io/paweh/reference/plot_cc_power.md):
`pheno_error_multiplier` multiplies baseline `theta_base` and
`phi_base`; `geno_error_multiplier` multiplies the corresponding
baseline genotype-error parameters for the selected genotype
misclassification model. All arguments in `...` remain fixed while
`x_var` is swept. The y-axis is the required number of selected cases,
controls, or total individuals returned by
[`cc_mssn()`](https://akilanthony.github.io/paweh/reference/cc_mssn.md),
according to `sample_size`. For locus-heterogeneity sweeps, an exact
`pi = 0` design has no finite MSSN when target power exceeds `alpha`.
That scientifically structural boundary is retained with MSSN `NA`,
`finite_mssn = FALSE`, and `status = "no finite MSSN"`; other errors are
propagated unchanged.

## See also

[`cc_mssn`](https://akilanthony.github.io/paweh/reference/cc_mssn.md),
[`plot_cc_power`](https://akilanthony.github.io/paweh/reference/plot_cc_power.md).

## Examples

``` r
plot_cc_mssn(
  x_var = "geno_error_multiplier",
  x_values = c(0, 1, 2),
  test = "trend",
  input_mode = "model_based",
  power = 0.80,
  alpha = 0.05,
  prev = 0.10,
  pd = 0.25,
  R2 = 2,
  MOI = "M",
  geno_misclass = "3p",
  e01_base = 0.02,
  e02_base = 0.01,
  e03_base = 0.005,
  k = 1
)

```
