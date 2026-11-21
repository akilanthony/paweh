# Power and Sample Size for Genetic Association Studies

paweh

Plan genetic association studies before collecting data. Estimate power
or the sample size you need under assumptions about genetic
heterogeneity, genotype error, and phenotype misclassification.

[Get
Started](https://akilanthony.github.io/paweh/articles/paweh-01-getting-started.md)
[Study Designs](#study-designs) [Methods &
Sources](https://akilanthony.github.io/paweh/articles/paweh-08-methods-and-sources.md)
[Function
Reference](https://akilanthony.github.io/paweh/reference/index.md)

## What can paweh calculate?

Both calculations start with a proposed effect and a planned statistical
test. In this [worked
example](https://akilanthony.github.io/paweh/articles/paweh-01-getting-started.md),
the two calculations use the same case-control assumptions and genotype
chi-square test.

### How much power does my study have?

98.52% power

With 500 cases and 500 controls, this is the chance the test detects the
assumed association.

### How many participants do I need?

482 participants

For 80% power, recruit 241 cases and 241 controls. In a family study,
plan for complete trios instead.

## Choose your study design

Start with the people you can recruit and the measurements you expect to
collect.

- ### [Case-Control Studies](https://akilanthony.github.io/paweh/articles/paweh-02-case-control-study-design.md)

  Plan how many people with and without a condition to recruit for a
  comparison of genotypes.

- ### [TDT / Family Studies](https://akilanthony.github.io/paweh/articles/paweh-03-tdt-study-design.md)

  Plan a study of affected children and both biological parents using
  informative allele transmissions.

- ### [Quantitative Traits](https://akilanthony.github.io/paweh/articles/paweh-04-quantitative-trait-study-design.md)

  Plan a study of measurements such as height. Compare a full sample,
  selected extremes, or several traits.

- ### [Sequencing-Based Studies](https://akilanthony.github.io/paweh/articles/paweh-06-sequencing-study-design.md)

  Compare how sequencing depth and read errors affect recruitment for
  case-control studies and trios.

- ### [Validation-Based / LRTae Studies](https://akilanthony.github.io/paweh/articles/paweh-07-lrtae-study-design.md)

  Plan how many people to recruit when only some receive a reference
  diagnosis or a confirmatory genotype assay.

- ### [Interactive Dashboard](https://akilanthony.github.io/paweh/articles/paweh-05-interactive-dashboard.md)

  Explore case-control, trio, and quantitative-trait designs in a Shiny
  interface.

[Browse all
articles](https://akilanthony.github.io/paweh/articles/index.md)

## Install and get started

paweh is currently available from GitHub. The package has yet to be
released on CRAN. Install the stable `main` branch with the command
below.

    remotes::install_github("akilanthony/paweh", ref = "main")

Install `remotes` first if needed. For the current development package
and all eight tutorials, see the `develop` installation in [Getting
Started](https://akilanthony.github.io/paweh/articles/paweh-01-getting-started.md).

## Methodological foundation

[![Springer cover of Heterogeneity in Statistical Genetics by Derek
Gordon, Stephen J. Finch, and Wonkuk
Kim](reference/figures/textbook-cover.jpg)](https://doi.org/10.1007/978-3-030-61121-7)

Cover artwork by Springer. Used with permission.

Heterogeneity in Statistical Genetics: How to Assess, Address, and
Account for Mixtures in Association Studies

Derek Gordon, Stephen J. Finch, and Wonkuk Kim. Springer, 2020.

The book explains how mixtures and measurement errors can change a
genetic study’s power. Derek Gordon is a coauthor of both the book and
paweh. The package also draws on other publications for individual
methods.

[Read Methods &
Sources](https://akilanthony.github.io/paweh/articles/paweh-08-methods-and-sources.md)

## Capabilities and limitations

paweh supports prospective study design using specified effect sizes,
sampling plans, and measurement assumptions. Its sequencing calculations
use planned read depths and error rates. Use separate analysis tools to
process raw sequencing files or run downstream GWAS association tests.

[Cite paweh](https://akilanthony.github.io/paweh/authors.html#citation)
