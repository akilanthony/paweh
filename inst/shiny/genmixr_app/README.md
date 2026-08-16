# genmixr alpha-test Shiny app

This is a stepwise Shiny interface for exploring `genmixr`
case-control and TDT power/sample-size workflows.

## Run locally

After installing the package, run:

```r
shiny::runApp(system.file("shiny/genmixr_app", package = "genmixr"))
```

During package development, you can also run the app from the repository root:

```r
shiny::runApp("inst/shiny/genmixr_app")
```

## Workflow

The app starts with two large cards:

- TDT
- Case-Control

The TDT path then offers Power, Sample Size, and Plots workflows.

The Case-Control path first asks for model-based or model-free input and a
checklist of tests, then offers Power, Sample Size, and Plots workflows.

## Backend functions

- TDT Power: `tdt_power()`
- TDT Sample Size: `tdt_mssn()`
- TDT Plots: `plot_tdt_power()` and `plot_tdt_mssn()`
- Case-Control Power: `cc_power()`
- Case-Control Sample Size: `cc_mssn()`
- Case-Control Plots: `plot_cc_power()` and `plot_cc_mssn()`

No deployment credentials, tokens, or secrets are included.
