# paweh: Power and Sample Size Analysis for Genetic Studies

`paweh` provides prospective power and minimum sample size necessary
(MSSN) calculations for genetic and genomic association study design.
Many standard calculations use homogeneous, error-free assumptions. The
package implements published statistical-genetics methods for evaluating
selected departures from those assumptions, particularly locus
heterogeneity and supported forms of genotype and phenotype
misclassification.

## Study-design families

The principal interfaces are
[`cc_power()`](https://akilanthony.github.io/paweh/reference/cc_power.md)
and
[`cc_mssn()`](https://akilanthony.github.io/paweh/reference/cc_mssn.md)
for case-control designs;
[`tdt_power()`](https://akilanthony.github.io/paweh/reference/tdt_power.md)
and
[`tdt_mssn()`](https://akilanthony.github.io/paweh/reference/tdt_mssn.md)
for transmission disequilibrium test designs with affected-child trios;
and
[`qtl_anova_power()`](https://akilanthony.github.io/paweh/reference/qtl_anova_power.md),
[`qtl_anova_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_anova_mssn.md),
[`qtl_threshold_chisq_power()`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_power.md),
[`qtl_threshold_chisq_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_mssn.md),
[`qtl_multivariate_power_full()`](https://akilanthony.github.io/paweh/reference/qtl_multivariate_power_full.md),
and
[`qtl_multivariate_mssn_full()`](https://akilanthony.github.io/paweh/reference/qtl_multivariate_mssn_full.md)
for quantitative-trait designs. Calculations are primarily parameter and
model based and support prospective planning, not downstream association
testing.

## Visualization and interactive exploration

High-level plotting functions include
[`plot_cc_power()`](https://akilanthony.github.io/paweh/reference/plot_cc_power.md),
[`plot_cc_mssn()`](https://akilanthony.github.io/paweh/reference/plot_cc_mssn.md),
[`plot_tdt_power()`](https://akilanthony.github.io/paweh/reference/plot_tdt_power.md),
[`plot_tdt_mssn()`](https://akilanthony.github.io/paweh/reference/plot_tdt_mssn.md),
and the quantitative-trait plotting interfaces documented with the
corresponding analysis functions.
[`paweh_app()`](https://akilanthony.github.io/paweh/reference/paweh_app.md)
returns an optional Shiny dashboard that uses the same canonical
calculations for interactive study-design exploration.

## Vignettes

Start with
[`vignette("paweh-01-getting-started")`](https://akilanthony.github.io/paweh/articles/paweh-01-getting-started.md).
Study-specific guidance is available in
[`vignette("paweh-02-case-control-study-design")`](https://akilanthony.github.io/paweh/articles/paweh-02-case-control-study-design.md),
[`vignette("paweh-03-tdt-study-design")`](https://akilanthony.github.io/paweh/articles/paweh-03-tdt-study-design.md),
and
[`vignette("paweh-04-quantitative-trait-study-design")`](https://akilanthony.github.io/paweh/articles/paweh-04-quantitative-trait-study-design.md).
The dashboard user guide is
[`vignette("paweh-05-interactive-dashboard")`](https://akilanthony.github.io/paweh/articles/paweh-05-interactive-dashboard.md).

## See also

Useful links:

- <https://github.com/akilanthony/paweh>

- <https://akilanthony.github.io/paweh/>

- Report bugs at <https://github.com/akilanthony/paweh/issues>

## Author

**Maintainer**: Akil Anthony <akilanthony19@gmail.com>
([ORCID](https://orcid.org/0009-0008-8813-8560))

Authors:

- Akil Anthony <akilanthony19@gmail.com>
  ([ORCID](https://orcid.org/0009-0008-8813-8560))
