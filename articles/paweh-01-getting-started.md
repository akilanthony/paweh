# Getting started with paweh

## What problem does `paweh` solve?

`paweh` supports prospective power and minimum-sample-size-necessary
(MSSN) calculations for genetic association study designs. It is
intended for the stage at which a researcher is deciding what data to
collect, how many subjects or trios are needed, and how sensitive that
design may be to model assumptions.

The relevant design inputs differ across studies. They can include
sample size, allele frequency, genotype relative risk, genotype or
phenotype misclassification, locus heterogeneity, family-based
transmission, extreme-phenotype sampling, and correlation among
quantitative traits. These quantities are assumptions or design targets,
not estimates produced from raw genomic data by `paweh`.

## Power and MSSN answer inverse questions

- **Power:** given a sample size and model, what probability of
  rejection is expected under the specified alternative?
- **MSSN:** given a target power and model, how many cases, controls,
  trios, or total individuals are required?

For example, suppose pilot work supplies conditional genotype
probabilities for cases and controls. A compact power calculation is:

``` r

g_case <- c(0.25, 0.50, 0.25)
g_control <- c(0.36, 0.48, 0.16)

starter_power <- cc_power(
  N_case = 500,
  alpha = 0.05,
  input_mode = "model_free",
  g1 = g_case,
  g0 = g_control,
  k = 1,
  verbose = FALSE
)

starter_power$tests$genotypes$power
#> [1] 0.9852151
```

The corresponding
[`cc_mssn()`](https://akilanthony.github.io/paweh/reference/cc_mssn.md)
calculation reverses the design question by replacing `N_case` with a
target `power`.

## Three study-design families

### Case-control association

Cases and controls are compared using a genotype chi-square or trend
test. The canonical
[`cc_power()`](https://akilanthony.github.io/paweh/reference/cc_power.md)
and
[`cc_mssn()`](https://akilanthony.github.io/paweh/reference/cc_mssn.md)
interfaces support model-based and model-free inputs, along with
compatible heterogeneity and error modifiers.

### Affected-child trios

The transmission disequilibrium test (TDT) compares transmitted and
non-transmitted alleles in affected-child trios.
[`tdt_power()`](https://akilanthony.github.io/paweh/reference/tdt_power.md)
and
[`tdt_mssn()`](https://akilanthony.github.io/paweh/reference/tdt_mssn.md)
report power or required trio counts under distinct no-error,
phenotype-misclassification, and locus-heterogeneity scenarios.

### Quantitative traits

For continuous traits, `paweh` supports one-way genotype-group ANOVA,
threshold-selected genotype chi-square designs, and multivariate designs
using Pillai MANOVA or joint threshold selection.

## Model-based and model-free inputs

These modes locate the scientific assumptions in different places.

For **case-control** studies, model-based input uses prevalence,
disease-allele frequency, genotype relative risk, and mode of
inheritance. Model-free input uses supplied conditional genotype
probabilities: `g1` for cases or affected subjects and `g0` for controls
or unaffected subjects.

For the **TDT**, model-based input uses `pd`, `prev`, `R1`, `R2`, and
`delta_prime`. Model-free input starts from expected transmission and
non-transmission counts (`ET` and `ENT`), with the associated number of
trios needed for MSSN calculations.

Model-free input is useful when these quantities come from pilot data,
prior publications, external calculations, or direct design assumptions.
It does not remove assumptions; it moves them into the supplied
conditional quantities.

## Where should I start?

| Scientific design | Canonical starting functions | Detailed vignette |
|----|----|----|
| Case-control association | [`cc_power()`](https://akilanthony.github.io/paweh/reference/cc_power.md), [`cc_mssn()`](https://akilanthony.github.io/paweh/reference/cc_mssn.md) | Case-Control Study Design |
| Affected-child trios | [`tdt_power()`](https://akilanthony.github.io/paweh/reference/tdt_power.md), [`tdt_mssn()`](https://akilanthony.github.io/paweh/reference/tdt_mssn.md) | TDT Study Design |
| Continuous single-trait QTL | [`qtl_anova_power()`](https://akilanthony.github.io/paweh/reference/qtl_anova_power.md), [`qtl_anova_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_anova_mssn.md) | Quantitative-Trait Study Design |
| Extreme-trait QTL | [`qtl_threshold_chisq_power()`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_power.md), [`qtl_threshold_chisq_mssn()`](https://akilanthony.github.io/paweh/reference/qtl_threshold_chisq_mssn.md) | Quantitative-Trait Study Design |
| Multivariate quantitative traits | [`qtl_multivariate_power_full()`](https://akilanthony.github.io/paweh/reference/qtl_multivariate_power_full.md), [`qtl_multivariate_mssn_full()`](https://akilanthony.github.io/paweh/reference/qtl_multivariate_mssn_full.md) | Quantitative-Trait Study Design |

## Why Bioconductor?

`paweh` is being prepared for Bioconductor because it addresses
prospective design questions in genetic association research. These
calculations can complement downstream genomic analysis by examining
operating characteristics before data generation or collection. The
package does not currently claim sequencing-file, genomic-range, or
`SummarizedExperiment` integration.

## Published validation philosophy

The remaining vignettes distinguish three uses of literature:

- **Published-example reproduction:** a reported analytic result is
  evaluated with a matching `paweh` design.
- **Published methodological context:** a paper explains an important
  risk or method, but its exact analysis is not claimed to be
  reproduced.
- **Motivating published application:** an application motivates a
  design question while using a different statistical procedure.

The study-design vignettes reproduce validated examples across
case-control, TDT, and multivariate quantitative-trait settings and then
extend those fixed designs with sensitivity calculations.

## Installation

Once `paweh` is available through Bioconductor, the standard
installation path will be:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("paweh")
```

No installation is performed while this vignette is built.

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] paweh_0.0.0.9000 BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        jsonlite_2.0.0      dplyr_1.2.1        
#>  [4] compiler_4.6.1      BiocManager_1.30.27 tidyselect_1.2.1   
#>  [7] jquerylib_0.1.4     systemfonts_1.3.2   scales_1.4.0       
#> [10] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
#> [13] ggplot2_4.0.3       R6_2.6.1            generics_0.1.4     
#> [16] knitr_1.51          htmlwidgets_1.6.4   tibble_3.3.1       
#> [19] bookdown_0.47       desc_1.4.3          bslib_0.12.0       
#> [22] pillar_1.11.1       RColorBrewer_1.1-3  rlang_1.3.0        
#> [25] cachem_1.1.0        xfun_0.60           fs_2.1.0           
#> [28] sass_0.4.10         S7_0.2.2            otel_0.2.0         
#> [31] cli_3.6.6           pkgdown_2.2.1       magrittr_2.0.5     
#> [34] digest_0.6.39       grid_4.6.1          mvtnorm_1.4-2      
#> [37] lifecycle_1.0.5     vctrs_0.7.3         evaluate_1.0.5     
#> [40] glue_1.8.1          farver_2.1.2        ragg_1.5.2         
#> [43] rmarkdown_2.31      pkgconfig_2.0.3     tools_4.6.1        
#> [46] htmltools_0.5.9
```
