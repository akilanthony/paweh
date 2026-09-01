# Sequencing-Based Genetic Study Design with PAWEH

## Overview

`paweh` supports prospective analytic power and
minimum-sample-size-necessary (MSSN) calculations for two scientifically
different sequencing designs. The inputs describe a planned experiment;
these functions do not process FASTQ, BAM, or VCF files and do not
perform downstream association testing.

For case-control sequencing designs, the pathway is:

> true genotype -\> sequencing read model -\> deterministic genotype
> call -\> 3 by 3 transition matrix -\> Ahn/Chapman–Nam trend NCP -\>
> power or MSSN.

For affected-child sequencing trios, the pathway is:

> latent Mendelian trio state -\> raw read-count likelihood -\>
> efficient information -\> Kim TDT1-NGS NCP -\> power or MSSN.

The case-control calculation therefore propagates uncertainty in called
genotypes. TDT1-NGS instead retains raw-read likelihood information and
does not hard-call the three genotypes. These are not interchangeable
methods.

## Two sequencing workflows in `paweh`

| Feature | Case-control NGS | TDT1-NGS |
|----|----|----|
| Sampling unit | Cases and controls | Complete father-mother-affected-child trios |
| Sequencing information | Called-genotype probabilities | Raw read-count likelihood |
| Association kernel | Ahn/Chapman–Nam trend NCP | Kim efficient-information NCP |
| Public power function | [`cc_ngs_power()`](https://akilanthony.github.io/paweh/reference/cc_ngs_power.md) | [`tdt_ngs_power()`](https://akilanthony.github.io/paweh/reference/tdt_ngs_power.md) |
| Public MSSN function | [`cc_ngs_mssn()`](https://akilanthony.github.io/paweh/reference/cc_ngs_mssn.md) | [`tdt_ngs_mssn()`](https://akilanthony.github.io/paweh/reference/tdt_ngs_mssn.md) |
| Locus heterogeneity | Available | Not currently implemented |

Both designs currently use equal fixed coverage and a symmetric public
per-read error parameter. Neither uses simulation.

## Exact published TDT1-NGS replication: the Kim AIS example

Kim (2015) used adolescent idiopathic scoliosis (AIS) to illustrate a
low-coverage trio design. The motivating variant, rs1400180, had a
reported global minor-allele frequency near 0.349; the design
calculation used 0.35. The published assumptions were HWE, symmetric
read error 0.025, genome-wide `alpha = 5e-8`, a multiplicative effect
reported as OR at least 2, and target power 0.90. In the TDT1-NGS
parameterization this effect is `R1 = 2`, with `R2 = R1^2`.

``` r

kim_4x <- tdt_ngs_mssn(
  power = 0.90,
  pd = 0.35,
  R1 = 2,
  coverage = 4,
  seq_error = 0.025,
  alpha = 5e-8,
  verbose = FALSE
)

kim_25x <- tdt_ngs_mssn(
  power = 0.90,
  pd = 0.35,
  R1 = 2,
  coverage = 25,
  seq_error = 0.025,
  alpha = 5e-8,
  verbose = FALSE
)

data.frame(
  coverage = c("4x", "25x"),
  published_trios = c(654, 416),
  paweh_trios = c(kim_4x$MSSN_trios, kim_25x$MSSN_trios)
)
#>   coverage published_trios paweh_trios
#> 1       4x             654         654
#> 2      25x             416         416
```

Under Kim’s published HWE design with MAF 0.35, symmetric read error
0.025, genome-wide alpha $`5\times10^{-8}`$, and multiplicative effect
`R1 = 2`, PAWEH exactly reproduces the reported integer requirements of
654 complete trios at 4x and 416 at 25x for 90% power. This statement
applies to this specific published design and is not a claim about every
model in the paper.

Kim’s full Table 3 pairs three target powers with three effect sizes.
All six integer results reproduce exactly:

``` r

kim_table3 <- data.frame(
  coverage = rep(c(4, 25), times = 3),
  power = rep(c(0.90, 0.70, 0.10), each = 2),
  R1 = rep(c(2.0, 1.8, 1.5), each = 2),
  published_trios = c(654, 416, 716, 456, 733, 466)
)

kim_table3$paweh_trios <- c(
  kim_4x$MSSN_trios,
  kim_25x$MSSN_trios,
  vapply(3:6, function(i) {
    tdt_ngs_mssn(
      power = kim_table3$power[i], pd = 0.35, R1 = kim_table3$R1[i],
      coverage = kim_table3$coverage[i], seq_error = 0.025,
      alpha = 5e-8, verbose = FALSE
    )$MSSN_trios
  }, numeric(1))
)
kim_table3
#>   coverage power  R1 published_trios paweh_trios
#> 1        4   0.9 2.0             654         654
#> 2       25   0.9 2.0             416         416
#> 3        4   0.7 1.8             716         716
#> 4       25   0.7 1.8             456         456
#> 5        4   0.1 1.5             733         733
#> 6       25   0.1 1.5             466         466
```

## TDT1-NGS power and MSSN

For Kim’s one-degree-of-freedom multiplicative model, the NCP is

``` math
\lambda=N\{\log(R_1)\}^2 I_{eff},
```

where $`N`$ is the number of complete trios and $`I_{eff}`$ is the
nuisance-adjusted per-trio information obtained from the raw-read
likelihood under the null. The public interface does not require users
to manipulate the full information matrix.

``` r

tdt_seq_power <- tdt_ngs_power(
  N = 5000,
  pd = 0.325,
  R1 = 1.2,
  coverage = 12,
  seq_error = 0.005,
  alpha = 5e-8,
  verbose = FALSE
)

tdt_seq_mssn <- tdt_ngs_mssn(
  power = 0.80,
  pd = 0.325,
  R1 = 1.2,
  coverage = 12,
  seq_error = 0.005,
  alpha = 5e-8,
  verbose = FALSE
)

data.frame(
  quantity = c("Power at 5,000 trios", "Required trios", "Total individuals"),
  value = c(
    tdt_seq_power$power,
    tdt_seq_mssn$MSSN_trios,
    tdt_seq_mssn$total_individuals
  )
)
#>               quantity        value
#> 1 Power at 5,000 trios 7.072653e-01
#> 2       Required trios 5.507000e+03
#> 3    Total individuals 1.652100e+04
```

The MSSN calculation first finds the target one-degree-of-freedom NCP
and then solves analytically for $`N`$. It is not an iterative
simulation over possible trio counts.

## TDT sequencing sensitivity

A small coverage-by-error grid helps identify whether a proposed design
is in a steep low-depth region or near a higher-depth plateau.

``` r

plot_ngs_power(
  design = "tdt",
  coverage = c(2, 4, 8, 12, 20),
  seq_error = c(0.005, 0.01),
  target_power = 0.80,
  N = 5000, pd = 0.325, R1 = 1.2, alpha = 5e-8
)
```

![TDT1-NGS power across equal fixed coverage for two symmetric per-read
error
assumptions.](paweh-06-sequencing-study-design_files/figure-html/tdt-ngs-power-plot-1.png)

TDT1-NGS power across equal fixed coverage for two symmetric per-read
error assumptions.

``` r

plot_ngs_mssn(
  design = "tdt",
  coverage = c(2, 4, 8, 12, 20),
  seq_error = c(0.005, 0.01),
  power = 0.80, pd = 0.325, R1 = 1.2, alpha = 5e-8
)
```

![Required complete trios across equal fixed coverage for two symmetric
per-read error
assumptions.](paweh-06-sequencing-study-design_files/figure-html/tdt-ngs-mssn-plot-1.png)

Required complete trios across equal fixed coverage for two symmetric
per-read error assumptions.

Low-to-moderate depth gains can be large, whereas gains at higher depth
may plateau. The effect of sequencing error remains setting-dependent
and should not be assumed negligible. The current fixed-depth TDT1-NGS
implementation does not support a design in which every member of every
trio is modeled at exactly 1x coverage under the full nuisance
parameterization. This is a specific identifiability boundary of the
implemented model, not a claim that 1x sequencing is universally
unusable.

## Published case-control framework: the Ahn example

Ahn et al. (2007) studied ordinary genotype misclassification, not
sequencing depth. The paper supplies genotype-error probabilities
directly. PAWEH implements Ahn Equation (1), its analytic MSSN
rearrangement, and the published trend scores T001 = `(0,0,1)`, T011 =
`(0,1,1)`, and T012 = `(0,1,2)`.

The stronger numerical example uses prevalence 0.005, marker and disease
allele frequencies 0.30, $`D'=0.5`$, and relative risks 2 and 3. The
following conditional frequencies were independently reconstructed from
Ahn’s genetic model during package validation. Adjacent-genotype errors
of 5% correspond to the PAWEH three-parameter matrix with
`e01 = e02 = 0.05` and `e03 = 0`.

``` r

ahn_case <- c(0.398125, 0.472500, 0.129375)
ahn_control <- c(0.490461683417085, 0.419736180904523, 0.089802135678392)

ahn_no_error <- cc_power(
  N_case = 1000, alpha = 0.001,
  input_mode = "model_free",
  g1 = ahn_case, g0 = ahn_control,
  k = 1, w = c(0, 1, 2),
  geno_misclass = "none", verbose = FALSE
)

ahn_adjacent_error <- cc_power(
  N_case = 1000, alpha = 0.001,
  input_mode = "model_free",
  g1 = ahn_case, g0 = ahn_control,
  k = 1, w = c(0, 1, 2),
  geno_misclass = "3p", e01 = 0.05, e02 = 0.05, e03 = 0,
  verbose = FALSE
)

data.frame(
  scenario = c("No error", "No error", "5% adjacent error", "5% adjacent error"),
  test = c("T012", "Genotype chi-square", "T012", "Genotype chi-square"),
  published_approximate_power = c(0.88, 0.81, 0.79, 0.70),
  paweh_power = c(
    ahn_no_error$tests$trend$power,
    ahn_no_error$tests$genotypes$power,
    ahn_adjacent_error$tests$trend$power,
    ahn_adjacent_error$tests$genotypes$power
  )
)
#>            scenario                test published_approximate_power paweh_power
#> 1          No error                T012                        0.88   0.8742681
#> 2          No error Genotype chi-square                        0.81   0.8060845
#> 3 5% adjacent error                T012                        0.79   0.7875922
#> 4 5% adjacent error Genotype chi-square                        0.70   0.6992310
```

This is a **Published-method replication**, not a CC-NGS
sequencing-depth replication. PAWEH implements Ahn’s Equation (1)
exactly and closely reproduces the paper’s ordinary
genotype-misclassification power example.

## The PAWEH CC-NGS sequencing bridge

Ahn takes $`P(\text{observed genotype}\mid\text{true genotype})`$ as an
input. PAWEH adds an upstream analytic bridge from sequencing
assumptions to those probabilities. For true genotype $`g\in\{0,1,2\}`$,
alternate-read count $`X`$, and fixed depth $`v`$,

``` math
X\mid G=g,V=v\sim\operatorname{Binomial}(v,q_g),\qquad
(q_0,q_1,q_2)=(\epsilon,0.5,1-\epsilon).
```

A deterministic maximum-likelihood rule calls genotype $`k`$. Each
transition probability

``` math
E_{gk}(v,\epsilon)=P(\text{called genotype}=k\mid\text{true genotype}=g)
```

is obtained by finite summation over all possible read counts. Rows of
$`E`$ are true genotypes and columns are called genotypes. This
coverage-to-matrix bridge is a PAWEH composition; it was not proposed in
the Ahn paper.

``` r

ngs_genotype_error_matrix(coverage = 4, seq_error = 0.005)
#>            called_0  called_1     called_2
#> true_0 9.801495e-01 0.0198505 6.250000e-10
#> true_1 6.250000e-02 0.8750000 6.250000e-02
#> true_2 6.250000e-10 0.0198505 9.801495e-01
```

## CC-NGS power and MSSN

The public case-control interfaces construct true case and control
genotype probabilities from the disease model, apply the sequencing
transition matrix, and then use the Ahn/Chapman–Nam trend NCP.

``` r

cc_seq_power <- cc_ngs_power(
  N_case = 1000,
  alpha = 5e-8,
  prev = 0.05,
  pd = 0.30,
  R2 = 1.8,
  coverage = 4,
  seq_error = 0.005,
  MOI = "M",
  k = 1,
  verbose = FALSE
)

cc_seq_mssn <- cc_ngs_mssn(
  power = 0.80,
  alpha = 5e-8,
  prev = 0.05,
  pd = 0.30,
  R2 = 1.8,
  coverage = 4,
  seq_error = 0.005,
  MOI = "M",
  k = 1,
  verbose = FALSE
)

data.frame(
  quantity = c(
    "Power at 1,000 cases", "Required cases", "Required controls", "Total MSSN"
  ),
  value = c(
    cc_seq_power$power,
    cc_seq_mssn$MSSN_case,
    cc_seq_mssn$MSSN_ctrl,
    cc_seq_mssn$MSSN_total
  )
)
#>               quantity        value
#> 1 Power at 1,000 cases    0.1197042
#> 2       Required cases 2168.0000000
#> 3    Required controls 2168.0000000
#> 4           Total MSSN 4336.0000000
```

Sequencing uncertainty is upstream of the published trend NCP. The MSSN
calculation remains analytic after inversion of the target NCP.

## Why zero per-read error is not perfect at finite depth

Setting `seq_error = 0` eliminates incorrect base reads, but it does not
make finite-depth heterozygote calls perfect. For a true heterozygote,

``` math
X\sim\operatorname{Binomial}(v,0.5).
```

At low $`v`$, every observed read can come from the same allele. Such an
observation can be more likely under a homozygous genotype, creating
call uncertainty even with no per-read sequencing error. This
probability decreases with depth but is nonzero at every finite depth.

``` r

list(
  coverage_2 = ngs_genotype_error_matrix(2, 0),
  coverage_10 = ngs_genotype_error_matrix(10, 0)
)
#> $coverage_2
#>        called_0 called_1 called_2
#> true_0     1.00      0.0     0.00
#> true_1     0.25      0.5     0.25
#> true_2     0.00      0.0     1.00
#> 
#> $coverage_10
#>            called_0  called_1     called_2
#> true_0 1.0000000000 0.0000000 0.0000000000
#> true_1 0.0009765625 0.9980469 0.0009765625
#> true_2 0.0000000000 0.0000000 1.0000000000
```

## CC sequencing sensitivity

``` r

plot_ngs_power(
  design = "cc",
  coverage = c(2, 4, 8, 12, 20),
  seq_error = c(0.005, 0.01),
  target_power = 0.80,
  N_case = 1000, alpha = 5e-8,
  prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M", k = 1
)
```

![Case-control sequencing trend power across equal fixed coverage for
two symmetric per-read error
assumptions.](paweh-06-sequencing-study-design_files/figure-html/cc-ngs-power-plot-1.png)

Case-control sequencing trend power across equal fixed coverage for two
symmetric per-read error assumptions.

``` r

plot_ngs_mssn(
  design = "cc",
  coverage = c(2, 4, 8, 12, 20),
  seq_error = c(0.005, 0.01),
  power = 0.80, alpha = 5e-8,
  prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M", k = 1
)
```

![Required case sample size across equal fixed coverage for two
symmetric per-read error
assumptions.](paweh-06-sequencing-study-design_files/figure-html/cc-ngs-mssn-plot-1.png)

Required case sample size across equal fixed coverage for two symmetric
per-read error assumptions.

These are deterministic scenario curves, not uncertainty bands. A robust
design should be examined over scientifically plausible depth and error
values rather than only at a single optimistic point.

## CC locus heterogeneity

Case-control sequencing designs can dilute the associated-case
distribution with the control distribution:

``` math
g_{case,H}=\pi g_{case}+(1-\pi)g_{control}.
```

Here `pi = 1` retains the full locus-associated case distribution and
`pi = 0` removes the locus-specific contrast. The master switch is
strict: when `locus_het = FALSE`, `pi` must remain 1.

``` r

cc_heterogeneous <- cc_ngs_power(
  N_case = 1000, alpha = 0.05,
  prev = 0.05, pd = 0.30, R2 = 1.8,
  coverage = 4, seq_error = 0.005,
  MOI = "M", k = 1,
  locus_het = TRUE, pi = 0.75,
  verbose = FALSE
)

data.frame(
  pi = cc_heterogeneous$locus_het$pi,
  lambda = cc_heterogeneous$lambda,
  power = cc_heterogeneous$power
)
#>     pi   lambda     power
#> 1 0.75 10.41674 0.8975178
```

## The exact `pi = 0` boundary

At `locus_het = TRUE, pi = 0`, case and control genotype distributions
are identical. Therefore the NCP is zero and power equals `alpha`. If
requested power exceeds `alpha`, no finite sample size exists.

The plotting-data interface preserves this scientific boundary instead
of inventing a very large finite MSSN:

``` r

pi_boundary <- plot_ngs_mssn(
  design = "cc",
  coverage = 4,
  seq_error = 0.005,
  return_data = TRUE,
  power = 0.80, alpha = 0.05,
  prev = 0.05, pd = 0.30, R2 = 1.8, MOI = "M", k = 1,
  locus_het = TRUE, pi = c(1, 0.75, 0.50, 0)
)

pi_boundary[, c("pi", "MSSN_case", "finite_mssn", "status")]
#>     pi MSSN_case finite_mssn         status
#> 1 1.00       430        TRUE         finite
#> 2 0.75       754        TRUE         finite
#> 3 0.50      1671        TRUE         finite
#> 4 0.00        NA       FALSE no finite MSSN
```

The `pi = 0` row remains visible with `MSSN_case = NA`,
`finite_mssn = FALSE`, and status `"no finite MSSN"`.

## Fixed-depth assumptions and limitations

Current sequencing APIs assume one equal fixed coverage value for the
relevant study members. They do not yet integrate over a realistic
empirical or budget-generated coverage distribution in which a locus may
be 0x in some samples, 1x or 2x in others, and 4x or greater elsewhere.
Variable/BGE coverage-distribution modeling is deferred.

| Assumption | CC-NGS | TDT1-NGS |
|----|----|----|
| Locus | Single biallelic | Single biallelic |
| Sampling | Cases and controls | Complete affected-child trios |
| Genetic model | M, dominant, or recessive trend scores | Multiplicative |
| Parental HWE/random mating | Not applicable | Assumed by current public design |
| Coverage | Equal fixed depth | Equal fixed depth |
| Public read error | Symmetric | Symmetric |
| Sequencing analysis | Deterministic genotype calls | Raw-read latent-state likelihood |
| Case/control observation error | Nondifferential | Not applicable |
| Locus heterogeneity | Optional canonical mixture | Not implemented |
| Simulation | None | None |

Deferred capabilities include variable/BGE coverage, simulation
validation, cost-depth optimization, TDT sequencing heterogeneity, and a
sequencing Shiny module. The existing dashboard covers the
non-sequencing study-design modules.

## A practical design workflow

1.  Choose an unrelated case-control or affected-child trio design.
2.  Specify allele-frequency and effect assumptions from defensible
    prior data.
3.  Choose candidate fixed coverage and per-read error scenarios.
4.  Calculate power at the available sample size.
5.  Calculate MSSN for the target power.
6.  Plot coverage and sequencing-error sensitivity.
7.  For case-control designs, stress-test plausible locus heterogeneity.
8.  Report the model, fixed-depth assumption, sampling unit, error
    assumption, and sensitivity range with the design recommendation.

## References

Ahn K, Haynes C, Kim W, St Fleur R, Gordon D, Finch SJ. The effects of
SNP genotyping errors on the power of the Cochran-Armitage linear trend
test for case/control association studies. *Annals of Human Genetics*.
2007;71(2):249–261. <https://doi.org/10.1111/j.1469-1809.2006.00318.x>.

Chapman DG, Nam JM. Asymptotic power of chi square tests for linear
trends in proportions. *Biometrics*. 1968;24(2):315–327.
<https://doi.org/10.2307/2528035>.

Gordon D, Finch SJ, Kim W. *Heterogeneity in Statistical Genetics: How
to Assess, Address, and Account for Mixtures in Association Studies*.
Springer; 2020. <https://doi.org/10.1007/978-3-030-61121-7>.

Kim W. Transmission disequilibrium tests based on read counts for
low-coverage next-generation sequence data. *Human Heredity*.
2015;80(1):36–49. <https://doi.org/10.1159/000434645>.

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
#> [1] paweh_0.99.0     BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        jsonlite_2.0.0      dplyr_1.2.1        
#>  [4] compiler_4.6.1      BiocManager_1.30.27 tidyselect_1.2.1   
#>  [7] jquerylib_0.1.4     systemfonts_1.3.2   scales_1.4.0       
#> [10] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
#> [13] ggplot2_4.0.3       R6_2.6.1            labeling_0.4.3     
#> [16] generics_0.1.4      knitr_1.51          htmlwidgets_1.6.4  
#> [19] tibble_3.3.1        bookdown_0.48       desc_1.4.3         
#> [22] bslib_0.12.0        pillar_1.11.1       RColorBrewer_1.1-3 
#> [25] rlang_1.3.0         cachem_1.1.0        xfun_0.60          
#> [28] fs_2.1.0            sass_0.4.10         S7_0.2.2           
#> [31] otel_0.2.0          cli_3.6.6           withr_3.0.3        
#> [34] pkgdown_2.2.1       magrittr_2.0.5      digest_0.6.39      
#> [37] grid_4.6.1          mvtnorm_1.4-2       lifecycle_1.0.5    
#> [40] vctrs_0.7.3         evaluate_1.0.5      glue_1.8.1         
#> [43] farver_2.1.2        ragg_1.5.2          rmarkdown_2.31     
#> [46] pkgconfig_2.0.3     tools_4.6.1         htmltools_0.5.9
```
