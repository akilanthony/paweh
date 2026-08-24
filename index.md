# pawh

`pawh` is an R package for prospective power and minimum sample size
necessary (MSSN) calculations in genetic association studies. It helps
researchers evaluate study operating characteristics before data
collection or genotyping: given a proposed sample size, `pawh`
calculates power; given a target power, it calculates the required
sample size.

The package supports case-control association studies, affected-child
trio designs using the Transmission Disequilibrium Test (TDT), and
continuous, extreme-phenotype, and multivariate quantitative-trait
designs. Depending on the design, calculations can account for genotype
misclassification, phenotype misclassification, and locus heterogeneity.
Canonical case-control and TDT interfaces support model-based and
model-free workflows where those input modes are applicable.

## What `pawh` supports

- **Case-control association:** two-degree-of-freedom genotype
  chi-square and one-degree-of-freedom genotype trend designs, with
  supported models for genotype or phenotype misclassification and locus
  heterogeneity.
- **Affected-child trios:** TDT power and required-trio calculations
  under no-error, phenotype-misclassification, and locus-heterogeneity
  scenarios.
- **Quantitative traits:** one-way genotype-group ANOVA for continuously
  measured traits, genotype chi-square designs for selected phenotype
  tails, and multivariate designs using Pillai MANOVA or joint threshold
  selection.

Model-based calculations express assumptions through quantities such as
prevalence, disease-allele frequency, genotype relative risk,
inheritance mode, or quantitative-trait variance. Where model-free input
is supported, users can instead supply conditional genotype
probabilities or expected transmission and non-transmission counts
derived from pilot data, prior work, or external calculations.
Model-free input does not eliminate assumptions; it places them in the
supplied conditional quantities.

Most functions return structured objects containing both the primary
design answer and intermediate quantities useful for validation.
Examples below extract only the result needed for interpretation;
consult the reference documentation when diagnostic components are
required.

## Documentation

The complete documentation is available at
<https://akilanthony.github.io/pawh/>.

Five articles provide the main entry points:

- [Getting
  Started](https://akilanthony.github.io/pawh/articles/pawh-01-getting-started.html)
- [Case-Control Study
  Design](https://akilanthony.github.io/pawh/articles/pawh-02-case-control-study-design.html)
- [TDT Study
  Design](https://akilanthony.github.io/pawh/articles/pawh-03-tdt-study-design.html)
- [Quantitative-Trait Study
  Design](https://akilanthony.github.io/pawh/articles/pawh-04-quantitative-trait-study-design.html)
- [Interactive
  Dashboard](https://akilanthony.github.io/pawh/articles/pawh-05-interactive-dashboard.html)

Function-level reference documentation is also available from the
[pkgdown reference
index](https://akilanthony.github.io/pawh/reference/).

## Installation

`pawh` is currently under development and is not yet available through
Bioconductor. Install the development version from GitHub:

``` r

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("akilanthony/pawh")
```

## Interactive dashboard

`pawh` includes an interactive Shiny interface for:

- Case-Control studies
- TDT / Family studies
- Quantitative Trait studies

The dashboard uses the same canonical package functions as the
programmatic R interface and provides power, minimum-sample-size,
sensitivity, visualization, advanced calculation details, and
reproducible R calls. Construct and launch the dashboard explicitly
with:

``` r

app <- pawh_app()
shiny::runApp(app)
```

[`pawh_app()`](https://akilanthony.github.io/pawh/reference/pawh_app.md)
returns a Shiny application object and does not launch a browser or
server on its own.

See the [Interactive
Dashboard](https://akilanthony.github.io/pawh/articles/pawh-05-interactive-dashboard.html)
vignette for a guided walkthrough of all three study-design workspaces.

## Quick start

The canonical case-control interface can calculate genotype-test power
for a fixed design. This example assumes 500 cases, an equal number of
controls, 5% disease prevalence, a disease-allele frequency of 0.30, and
a genotype relative risk of 1.8 under the package’s multiplicative
model:

``` r

library(pawh)

power_result <- cc_power(
    N_case = 500,
    alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05,
    pd = 0.30,
    R2 = 1.8,
    MOI = "M",
    k = 1,
    verbose = FALSE
)

power_result$tests$genotypes$power
#> [1] 0.8366557
```

The inverse calculation asks how many subjects are required for 80%
power under the same assumptions:

``` r

mssn_result <- cc_mssn(
    power = 0.80,
    alpha = 0.05,
    input_mode = "model_based",
    prev = 0.05,
    pd = 0.30,
    R2 = 1.8,
    MOI = "M",
    k = 1,
    verbose = FALSE
)

c(
    cases = mssn_result$tests$genotypes$MSSN_case,
    controls = mssn_result$tests$genotypes$MSSN_ctrl,
    total = mssn_result$tests$genotypes$MSSN_total
)
#>    cases controls    total
#>      457      457      914
```

These are prospective operating characteristics under the stated model,
not estimates fitted from subject-level genomic data. See the
case-control article for model-free inputs and supported error and
heterogeneity scenarios.

Power and MSSN are inverse design questions, but integer recruitment
targets may not reproduce a requested power exactly because required
sample sizes are rounded up. Results also inherit all stated
distributional, prevalence, frequency, relative-risk,
linkage-disequilibrium, and error assumptions. The study-design articles
document these assumptions and identify settings that are outside the
current implementation.

## Study-design workflows

| Design | Primary functions | Documentation |
|----|----|----|
| Case-control association | [`cc_power()`](https://akilanthony.github.io/pawh/reference/cc_power.md), [`cc_mssn()`](https://akilanthony.github.io/pawh/reference/cc_mssn.md) | [Case-Control Study Design](https://akilanthony.github.io/pawh/articles/pawh-02-case-control-study-design.html) |
| Affected-child trios / TDT | [`tdt_power()`](https://akilanthony.github.io/pawh/reference/tdt_power.md), [`tdt_mssn()`](https://akilanthony.github.io/pawh/reference/tdt_mssn.md) | [TDT Study Design](https://akilanthony.github.io/pawh/articles/pawh-03-tdt-study-design.html) |
| Continuous quantitative trait | [`qtl_anova_power()`](https://akilanthony.github.io/pawh/reference/qtl_anova_power.md), [`qtl_anova_mssn()`](https://akilanthony.github.io/pawh/reference/qtl_anova_mssn.md) | [Quantitative-Trait Study Design](https://akilanthony.github.io/pawh/articles/pawh-04-quantitative-trait-study-design.html) |
| Extreme quantitative trait | [`qtl_threshold_chisq_power()`](https://akilanthony.github.io/pawh/reference/qtl_threshold_chisq_power.md), [`qtl_threshold_chisq_mssn()`](https://akilanthony.github.io/pawh/reference/qtl_threshold_chisq_mssn.md) | [Quantitative-Trait Study Design](https://akilanthony.github.io/pawh/articles/pawh-04-quantitative-trait-study-design.html) |
| Multivariate quantitative traits | [`qtl_multivariate_power_full()`](https://akilanthony.github.io/pawh/reference/qtl_multivariate_power_full.md), [`qtl_multivariate_mssn_full()`](https://akilanthony.github.io/pawh/reference/qtl_multivariate_mssn_full.md) | [Quantitative-Trait Study Design](https://akilanthony.github.io/pawh/articles/pawh-04-quantitative-trait-study-design.html) |

## Published-example validation

`pawh` includes literature-backed regression tests and worked examples
for case-control genotype and phenotype misclassification, TDT phenotype
misclassification and locus heterogeneity, and multivariate
quantitative-trait power and sample size. The articles distinguish exact
or near-exact published-example reproductions from papers used only as
methodological context or motivating applications.

## Scope

`pawh` is a prospective statistical study-design package. It does not
perform GWAS association testing, sequence processing, variant
annotation, or general genomic data management. Instead, it is intended
to complement downstream genetic-analysis workflows by helping
researchers examine power, sample size, and sensitivity to supported
design assumptions before collecting data or performing genotyping.

## Citation

A formal package citation will be added as `pawh` approaches its initial
Bioconductor release. In the meantime, please consult the methodological
references in the package documentation and study-design articles.

## License

`pawh` is available under the MIT License.
