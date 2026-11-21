# Construct the paweh Shiny Dashboard

Creates the modular `bslib` dashboard for `paweh` study design. The
Case-Control and TDT / Family workspaces provide canonical power,
minimum-sample-size, sensitivity, and study-specific visual results. The
Quantitative Trait workspace collects the canonical function inputs and
displays that function's verbose console output. Dashboard calculations
delegate to the same canonical package functions as the programmatic R
interface.

## Usage

``` r
paweh_app()
```

## Value

A Shiny application object inheriting from `shiny.appobj`.

## Details

`paweh_app()` returns the application object without launching it. This
makes the app safe to construct in package code, tests, and deployment
tooling. Launch it explicitly with
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Examples

``` r
app <- paweh_app()
inherits(app, "shiny.appobj")
#> [1] TRUE
if (interactive()) {
  shiny::runApp(app)
}
```
