# paweh

`paweh` is statistical-genetics study-design software for prospective
power and minimum sample size necessary (MSSN) calculations. Given a
proposed sample size, it calculates power; given a target power, it
calculates the sample size needed under the specified design
assumptions. Its calculations are primarily parameter and model based
and are intended for planning before data collection or genotyping, not
for downstream association testing.

The package supports case-control association studies, affected-child
trio designs using the Transmission Disequilibrium Test (TDT), and
continuous, extreme-phenotype, and multivariate quantitative-trait
designs. Depending on the design, calculations can account for genotype
misclassification, phenotype misclassification, and locus heterogeneity.
Canonical case-control and TDT interfaces support model-based and
model-free workflows where those input modes are applicable.

Genetic study-design calculations are commonly formulated under
homogeneous, error-free assumptions. Selected departures from those
assumptions can materially affect power and MSSN while being difficult
to incorporate analytically. `paweh` provides accessible, reproducible
implementations of published methods for examining supported
heterogeneity and misclassification mechanisms alongside conventional
power and sample-size calculations.

## What `paweh` supports

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
- **Sequencing-based designs:** analytic power, MSSN, and sensitivity
  plots for fixed equal-depth case-control and TDT1-NGS studies with
  symmetric per-read sequencing error. Case-control NGS supports locus
  heterogeneity.

For case-control NGS, `paweh` constructs a sequencing-to-called-genotype
transition matrix and then applies the published Ahn / Chapman–Nam
trend-test framework. TDT1-NGS implements the Kim (2015) raw-read
likelihood and information framework under the documented fixed
equal-depth multiplicative model. In the matching Kim AIS design, the
package reproduces the reported requirements of 654 complete trios at 4x
and 416 at 25x.

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
<https://akilanthony.github.io/paweh/>.

Six articles provide the main entry points:

- [Getting
  Started](https://akilanthony.github.io/paweh/articles/paweh-01-getting-started.html)
- [Case-Control Study
  Design](https://akilanthony.github.io/paweh/articles/paweh-02-case-control-study-design.html)
- [TDT Study
  Design](https://akilanthony.github.io/paweh/articles/paweh-03-tdt-study-design.html)
- [Quantitative-Trait Study
  Design](https://akilanthony.github.io/paweh/articles/paweh-04-quantitative-trait-study-design.html)
- [Interactive
  Dashboard](https://akilanthony.github.io/paweh/articles/paweh-05-interactive-dashboard.html)
- [Sequencing-Based Genetic Study
  Design](https://akilanthony.github.io/paweh/articles/paweh-06-sequencing-study-design.html)

Function-level reference documentation is also available from the
[pkgdown reference
index](https://akilanthony.github.io/paweh/reference/).

## Installation

`paweh` is currently under development and is not yet available through
Bioconductor. Install the development version from GitHub:

``` r

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("akilanthony/paweh")
```

## Interactive dashboard

`paweh` includes an interactive Shiny interface for:

- Case-Control studies
- TDT / Family studies
- Quantitative Trait studies

The dashboard uses the same canonical package functions as the
programmatic R interface and provides power, minimum-sample-size,
sensitivity, visualization, advanced calculation details, and
reproducible R calls. Construct and launch the dashboard explicitly
with:

``` r

app <- paweh_app()
shiny::runApp(app)
```

[`paweh_app()`](https://akilanthony.github.io/paweh/reference/paweh_app.md)
returns a Shiny application object and does not launch a browser or
server on its own.

See the [Interactive
Dashboard](https://akilanthony.github.io/paweh/articles/paweh-05-interactive-dashboard.html)
vignette for a guided walkthrough of all three study-design workspaces.

## Quick start

The canonical case-control interface can calculate genotype-test power
for a fixed design. This example assumes 500 cases, an equal number of
controls, 5% disease prevalence, a disease-allele frequency of 0.30, and
a genotype relative risk of 1.8 under the package’s multiplicative
model:

``` r

library(paweh)

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
| Case-control association | [`cc_power()`](https://akilanthony.github.io/paweh/reference/cc_power.md), [`cc_mssn()`](https://akilanthony.github.io/paweh/reference/cc_mssn.md) | [Case-Control Study Design](https://akilanthony.github.io/paweh/articles/paweh-02-case-control-study-design.html) |
| Affected-child trios / TDT | [`tdt_power()`](https://akilanthony.github.io/paweh/reference/tdt_power.md), [`tdt_mssn()`](https://akilanthony.github.io/paweh/reference/tdt_mssn.md) | [TDT Study Design](https://akilanthony.github.io/paweh/articles/paweh-03-tdt-study-design.html) |
| Continuous quantitative trait | [`qtl_anova_power()`](https://akilanthony.github.io/paweh/reference/qtl_anova_power.md), [`qtl_anova_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_anova_mssn.md) | [Quantitative-Trait Study Design](https://akilanthony.github.io/paweh/articles/paweh-04-quantitative-trait-study-design.html) |
| Extreme quantitative trait | [`qtl_threshold_chisq_power()`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_power.md), [`qtl_threshold_chisq_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_mssn.md) | [Quantitative-Trait Study Design](https://akilanthony.github.io/paweh/articles/paweh-04-quantitative-trait-study-design.html) |
| Multivariate quantitative traits | [`qtl_multivariate_power_full()`](https://akilanthony.github.io/paweh/reference/qtl_multivariate_power_full.md), [`qtl_multivariate_mssn_full()`](https://akilanthony.github.io/paweh/reference/qtl_multivariate_mssn_full.md) | [Quantitative-Trait Study Design](https://akilanthony.github.io/paweh/articles/paweh-04-quantitative-trait-study-design.html) |
| Case-control NGS | [`cc_ngs_power()`](https://akilanthony.github.io/paweh/reference/cc_ngs_power.md), [`cc_ngs_mssn()`](https://akilanthony.github.io/paweh/reference/cc_ngs_mssn.md) | [Sequencing-Based Genetic Study Design](https://akilanthony.github.io/paweh/articles/paweh-06-sequencing-study-design.html) |
| TDT1-NGS | [`tdt_ngs_power()`](https://akilanthony.github.io/paweh/reference/tdt_ngs_power.md), [`tdt_ngs_mssn()`](https://akilanthony.github.io/paweh/reference/tdt_ngs_mssn.md) | [Sequencing-Based Genetic Study Design](https://akilanthony.github.io/paweh/articles/paweh-06-sequencing-study-design.html) |

## Published-example validation

`paweh` includes literature-backed regression tests and worked examples
for case-control genotype and phenotype misclassification, TDT phenotype
misclassification and locus heterogeneity, and multivariate
quantitative-trait power and sample size. The articles distinguish exact
or near-exact published-example reproductions from papers used only as
methodological context or motivating applications.

## Scope

`paweh` is a prospective statistical study-design package. It does not
perform GWAS association testing, sequence processing, variant
annotation, or general genomic data management. Instead, it is intended
to complement downstream genetic-analysis workflows by helping
researchers examine power, sample size, and sensitivity to supported
design assumptions before collecting data or performing genotyping.

## Citation

A formal package citation will be added as `paweh` approaches its
initial Bioconductor release. In the meantime, please consult the
methodological references in the package documentation and study-design
articles.

## License

`paweh` is available under the MIT License.
