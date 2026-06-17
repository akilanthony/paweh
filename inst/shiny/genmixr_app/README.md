# genmixr alpha-test Shiny app

This is a draft stepwise Shiny interface for exploring `genmixr`
case-control and TDT power/sample-size workflows.

Draft note: equation/page references will be finalized after textbook
verification.

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

- TDT Power: `tdt_power_full()`
- TDT Sample Size: `tdt_required_trios_full()`
- TDT Plots: `tdt_plot_power()` and `tdt_plot_mssn()`
- Case-Control Power: `cc_power_conditional_full()`
- Case-Control Sample Size: `cc_mssn_conditional_full()`
- Case-Control Plots: `cc_plot_power()` and `cc_plot_mssn()`

No deployment credentials, tokens, or secrets are included.
