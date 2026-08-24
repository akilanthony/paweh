#' paweh: Power and Sample Size Analysis for Genetic Studies
#'
#' `paweh` provides prospective power and minimum sample size necessary (MSSN)
#' calculations for genetic and genomic association study design. Many standard
#' calculations use homogeneous, error-free assumptions. The package implements
#' published statistical-genetics methods for evaluating selected departures
#' from those assumptions, particularly locus heterogeneity and supported forms
#' of genotype and phenotype misclassification.
#'
#' @section Study-design families:
#' The principal interfaces are [cc_power()] and [cc_mssn()] for case-control
#' designs; [tdt_power()] and [tdt_mssn()] for transmission disequilibrium test
#' designs with affected-child trios; and [qtl_anova_power()],
#' [qtl_anova_mssn()], [qtl_threshold_chisq_power()],
#' [qtl_threshold_chisq_mssn()], [qtl_multivariate_power_full()], and
#' [qtl_multivariate_mssn_full()] for quantitative-trait designs. Calculations
#' are primarily parameter and model based and support prospective planning, not
#' downstream association testing.
#'
#' @section Visualization and interactive exploration:
#' High-level plotting functions include [plot_cc_power()], [plot_cc_mssn()],
#' [plot_tdt_power()], [plot_tdt_mssn()], and the quantitative-trait plotting
#' interfaces documented with the corresponding analysis functions.
#' [paweh_app()] returns an optional Shiny dashboard that uses the same canonical
#' calculations for interactive study-design exploration.
#'
#' @section Vignettes:
#' Start with `vignette("paweh-01-getting-started")`. Study-specific guidance is
#' available in `vignette("paweh-02-case-control-study-design")`,
#' `vignette("paweh-03-tdt-study-design")`, and
#' `vignette("paweh-04-quantitative-trait-study-design")`. The dashboard user
#' guide is `vignette("paweh-05-interactive-dashboard")`.
#'
#' @name paweh-package
#' @aliases paweh
#' @docType package
#' @keywords internal
"_PACKAGE"
