# genmixr draft Shiny app

This is a draft Shiny interface for exploring selected `genmixr` case-control
and TDT power/sample-size functions.

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

## Current scope

- Case-control tab: calls `cc_power_conditional_full()` and
  `cc_mssn_conditional_full()`.
- TDT tab: calls `tdt_power_from_ET_ENT()` and `tdt_required_trios()`.
- TDT Plots tab: calls `tdt_plot_power_misclassification()` with sliders for
  misclassification-rate sensitivity plots and a simple PNG download.
- Results panels show captured clean console output, key summary tables,
  observed genotype frequencies, and allele frequencies when available.

## TODOs

- Finalize textbook equation and page references.
- Refine UI labels and layout after manual review.
- Add additional TDT functions if needed.
- Add additional TDT plotting functions if useful.
- Decide whether any app-specific automated tests should be added later.

No deployment credentials, tokens, or secrets are included.
