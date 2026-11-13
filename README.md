# paweh: power and sample size for genetic association studies

If you can recruit 500 people with a disease and 500 without it, how likely is
your study to detect a difference at a particular genetic variant? `paweh`
calculates **power**: the chance that the planned test will detect an
association of the size you assume. It also answers the recruitment question
in reverse: the **minimum sample size necessary (MSSN)** to reach a chosen
power under those same assumptions.

People with the same diagnosis can have different underlying genetic causes
(*locus heterogeneity*). A lab can misread a genotype, and a diagnosis can be
recorded incorrectly. Each can change the difference a test sees between
groups. Use `paweh` before data collection to compare a baseline plan with
the departures that your chosen design supports.

## Methodological foundation

Gordon D, Finch SJ, Kim W. [*Heterogeneity in Statistical Genetics: How to
Assess, Address, and Account for Mixtures in Association Studies*](https://doi.org/10.1007/978-3-030-61121-7).
Springer; 2020. Derek Gordon is a coauthor of both the book and `paweh`.
The book provides a foundation for understanding mixtures and measurement
error. The package also draws on other published methods; each study-design
article identifies the sources for its calculations.

## Install and start here

The package is available from GitHub and is not yet on CRAN. Choose **one**
installation. For the stable package:

```r
# Run once if needed: install.packages("remotes")
remotes::install_github("akilanthony/paweh", ref = "main")
```

For the experimental `develop` branch **with the revised tutorials**, use the
following instead. Building the articles needs Pandoc (included with RStudio);
`dependencies = TRUE` installs the suggested R packages for the examples.

```r
remotes::install_github(
  "akilanthony/paweh", ref = "develop", dependencies = TRUE,
  build_vignettes = TRUE
)
```

Start with [Getting Started](vignettes/paweh-01-getting-started.Rmd), or open
the rendered tutorial after the second installation with
`vignette("paweh-01-getting-started", package = "paweh")`.

## Quick start: is the recruitment budget enough?

Suppose you plan to compare unrelated people with a disease (**cases**) and
without it (**controls**) at one prespecified variant. First choose a test:
here the genotype chi-square test asks whether the distribution of the three
genotypes differs between groups. These are illustrative assumptions, not
estimates for a particular disease.

| Argument | What it means in this study |
|---|---|
| `N_case = 500`, `k = 1` | Recruit 500 cases. The control-to-case ratio is one, so recruit 500 controls too (1,000 people in total). |
| `alpha = 0.05` | Use a 5% false-positive threshold for this one test under no association. A genome-wide scan needs a threshold accounting for many tests. |
| `prev = 0.05` | Assume 5% of the source population has the disease. This is different from the deliberately balanced recruitment split. |
| `pd = 0.30` | Assume 30% of allele copies in the population are the disease-associated version of the variant. |
| `R2 = 1.8` | People with two copies have 1.8 times the disease risk of people with zero copies. This is a relative risk, not an odds ratio. |
| `MOI = "M"` | Multiplicative inheritance: each additional copy multiplies risk by the same factor, here `sqrt(1.8)`. |
| `power = 0.80` | For the inverse calculation, aim to detect the assumed association in 80% of repeated studies. |

Give `cc_power()` the recruitment budget. With `input_mode = "model_based"`,
the function derives the expected genotype frequencies from the disease
model; the default assumes no measurement error or heterogeneity.

```r
library(paweh)

fixed <- cc_power(
  N_case = 500, k = 1, alpha = 0.05, input_mode = "model_based",
  prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M", verbose = FALSE
)
fixed$tests$genotypes$power
#> [1] 0.8366557
```

The budget gives **about 84% power**: if these assumptions hold, about 84 of
100 such studies would detect the association. To plan for 80% instead, give
`cc_mssn()` the target power and keep the test and assumptions fixed:

```r
needed <- cc_mssn(
  power = 0.80, k = 1, alpha = 0.05, input_mode = "model_based",
  prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M", verbose = FALSE
)
g <- needed$tests$genotypes
c(cases = g$MSSN_case, controls = g$MSSN_ctrl, total = g$MSSN_total)
#>    cases controls    total
#>      457      457      914
```

Plan for **457 cases plus 457 controls**, or **914 people**. The function
rounds up to whole participants. If the expected contrast is smaller or
measurement is less reliable, rerun the design with those assumptions before
committing to a recruitment target.

## Choose your study design

| Design and guide | Planning question | Start with |
|---|---|---|
| [Case-control](vignettes/paweh-02-case-control-study-design.Rmd) | How many cases and controls would reveal a difference in genotype frequencies? | `cc_power()`, `cc_mssn()` |
| [Transmission disequilibrium test (TDT)](vignettes/paweh-03-tdt-study-design.Rmd) | In families with an affected child, is one allele passed on more often than expected? How many complete trios are needed? | `tdt_power()`, `tdt_mssn()` |
| [Quantitative traits](vignettes/paweh-04-quantitative-trait-study-design.Rmd) | How many people are needed to detect a difference in a measurement, such as height, across genotypes? The guide also covers sampling the extremes and studying multiple traits. | `qtl_anova_power()`, `qtl_anova_mssn()` |
| [Sequencing](vignettes/paweh-06-sequencing-study-design.Rmd) | At a specified read depth and error rate, how many cases and controls or complete trios are needed? | `cc_ngs_power()`, `cc_ngs_mssn()`; `tdt_ngs_power()`, `tdt_ngs_mssn()` |
| [Validation and measurement error (LRTae)](vignettes/paweh-07-lrtae-study-design.Rmd) | What power can a case-control study achieve when a subset receives gold-standard measurements? | `lrtae_power()`, `lrtae_mssn()` |

Prefer to explore assumptions interactively? The [dashboard guide](vignettes/paweh-05-interactive-dashboard.Rmd)
covers case-control, trio and quantitative-trait studies. Launch it with
`shiny::runApp(paweh::paweh_app())`. Use the
[function reference](https://akilanthony.github.io/paweh/reference/) for
arguments and return values, or `?paweh::cc_power` for help matching the
version you installed.

## Related software and scope

[Quanto](https://keck.usc.edu/biostatistics/software/) calculates power and
sample size for genetic association studies, including gene–environment and
gene–gene interactions. Michigan’s [CaTS](https://csg.sph.umich.edu/abecasis/cats/)
addresses large association studies, including two-stage genome-wide designs;
its related [GAS calculator](https://csg.sph.umich.edu/abecasis/cats/gas_power_calculator/)
calculates power for one-stage studies.

Choose the tool that models your sampling plan and intended test. `paweh`
focuses on prospective calculations, including selected heterogeneity,
misclassification, sequencing and validation models. Error options differ by
design; for example, ordinary case-control calculations evaluate separate
error scenarios. The package does not process sequencing data or run
genome-wide association tests on participants.

## Citation and license

Run `citation("paweh")` for the software citation and the background book.
Also cite the method used in your analysis: for example,
[Gordon et al. (2002)](https://doi.org/10.1159/000066696) for case-control
genotyping error or [Edwards et al. (2005)](https://doi.org/10.1186/1471-2156-6-18)
for errors in disease classification. The package uses the MIT license.
