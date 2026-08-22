# Construct the pawh Shiny Dashboard

Creates the modular `bslib` dashboard for `pawh` study design. The
Case-Control and TDT / Family workspaces provide canonical power,
minimum-sample-size, sensitivity, and study-specific visual results. The
concise results are complemented by collapsed calculation details and
reproducible canonical R calls. The Quantitative Trait workspace is
clearly marked as forthcoming.

## Usage

``` r
pawh_app()
```

## Value

A Shiny application object inheriting from `shiny.appobj`.

## Details

`pawh_app()` returns the application object without launching it. This
makes the app safe to construct in package code, tests, and deployment
tooling. Launch it explicitly with
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html).

## Examples

``` r
app <- pawh_app()
inherits(app, "shiny.appobj")
#> [1] TRUE
if (interactive()) {
  shiny::runApp(app)
}
```
