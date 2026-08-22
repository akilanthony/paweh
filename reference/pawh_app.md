# Construct the pawh Shiny Dashboard

Creates the modular `bslib` dashboard for `pawh` study design. The
current application provides the landing page, study navigation, and
transparent placeholders for future Case-Control, TDT / Family, and
Quantitative Trait workflows. It does not yet perform dashboard
calculations.

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
