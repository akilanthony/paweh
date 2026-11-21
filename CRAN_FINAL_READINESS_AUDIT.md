# paweh final CRAN-readiness audit

**21 November 2026 · Phase A5.1 · NO GO**

## Executive summary

The audited release candidate is `447cd8cdc57fa388176d5c8b48c34ba3c341cd71` on `develop`, package version `0.99.0`. It builds, installs, and runs the eight vignettes on the primary local setup. Independent calculations support the main formulas and many examples. Three issues block submission in the current form: the LRTae search can overstate an attainable minimum enrollment, rounded QTL ANOVA can report power for an internally inconsistent enrollment, and redistribution rights for the Springer cover in the source package remain unverified. Resolve the scientific issues or withdraw the affected claims or options before release, and clear or remove the artwork before packaging. The TDT publication discrepancies require qualified documentation. The evidence does not support changing the TDT formulas to fit the printed counts.

The audit patch adds this report and exclusion rules and corrects TDT help. It leaves statistical code execution, benchmarks, public function signatures, the version, and the deployed site untouched. CI and portability results apply to the audited source plus the help corrections.

## Release candidate and repository scope

| Ref | Verified during source synchronization |
| --- | --- |
| Local and remote `develop` before synchronization | Local `447cd8cdc57fa388176d5c8b48c34ba3c341cd71`; remote `3e0d7b652562274af92d387697543b987672d8db` |
| Local and remote `develop` after Task 1 | `447cd8cdc57fa388176d5c8b48c34ba3c341cd71` |
| `main` | `7a5ed7bdfcab9f9f63e7e260204058bfcadf9cce` |
| `gh-pages` / deployed Pages build | `d41430c813c3a8738b16f2f448e304eea6803c8d` |

The comparison `origin/main...447cd8c` covered 33 files, with 4,761 additions and 1,703 deletions. The changes concern documentation and citation, site sources and assets, contributor files, and the case-control and QTL Shiny interfaces and their tests. The statistical backends and NAMESPACE show no difference in that comparison. The older `CRAN_PREPARATION_AUDIT.md` records earlier decisions, including the Bioconductor withdrawal and formatting conversion. Use this report for the eight-article candidate. Keep the older report as project history. The earlier documentation ZIP and original LRTae workbook were unavailable for comparison.

## CRAN metadata, dependencies, and structure

Reviewed the current [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html), [CRAN Submission Checklist](https://cran.r-project.org/web/packages/submission_checklist.html), and [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html) on 21 November 2026. The checks below separate repository requirements from scientific release decisions.

| ID | Severity | Files and evidence | Action | Blocks submission |
| --- | --- | --- | --- | --- |
| C1 | RESOLVED | `DESCRIPTION` has an informative title in title case and an English description without a package-name opening. It quotes 'shiny' as a package name, gives the book in CRAN DOI syntax, and credits it alongside method-specific sources. `Authors@R` lists two authors, a creator with email, and the creator's ORCID. Repository, site, and issue URLs identify the expected destinations. The UTF-8 declaration matches the rendered documentation. No Bioconductor metadata remains. | Keep the metadata and verify the maintainer email at submission. Use the release date or omit the optional date when staging the release. | No |
| C2 | RESOLVED | `inst/CITATION` derives the software citation from metadata and gives the textbook as a separate background citation. README, NEWS, and Methods & Sources describe prospective planning. | Preserve method-specific citations in help and articles. | No |
| C3 | IMPORTANT | `DESCRIPTION` advertises testthat ≥ 3.0.0. Multiple tests call `testthat::local_mocked_bindings()`, introduced in 3.1.7 and declared stable in 3.2.0. Current locally installed versions pass these calls. [testthat changelog](https://testthat.r-lib.org/news/index.html) No explicit R floor is declared. bslib's modern navbar/sidebar UI also needs a verified minimum. | Test a floor of 3.2.0, which supports the stable API, and check the rest of the suite against that floor. Check bslib's APIs and establish a defensible supported R floor with a dependency matrix. Avoid inferring dependency minimums from one populated library. | No, recommended before release |
| C4 | RESOLVED | `NAMESPACE` registers four S3 print methods for sequencing results. All 58 exports are documented and indexed. No package S4 classes, native code, compiled package library, or external system dependency is declared. | Keep build tooling in CI. Pandoc is needed for vignette builds and a TeX installation for the optional PDF manual. | No |
| C5 | RESOLVED | `R/falconer_multivariate.R` scopes the Genz–Bretz seed and `R/plot_falconer_distribution.R` scopes seeded simulation with withr. `R/` has help URLs and a local Shiny interface. Default tested analyses, examples, and vignettes operate on local or generated data. Package examples ran in about 11 seconds in aggregate across the follow-up checks. | Keep test scratch files under `tempdir()`. Check figure sizes and routine runtimes again on the final tarball. | No |
| C6 | IMPORTANT | `tests/testthat/` ran in 112 seconds on the standalone run. That does not establish CRAN timing on all builders. Trio raw-read likelihood storage grows with the product of the three depth supports. | Profile the slowest cases on Linux and Windows. Clarify practical depth limits as in S9 and keep execution within CRAN resource guidance. | No, unless target-platform checks exceed limits |
| C7 | RESOLVED | Full source archive is 2,799,602 bytes. Installed vignette directory is about 2.7 MiB including examples and sources. The main development files and older audit are outside the archive. The PDF manual was built separately. | Keep generated site output, original papers, workbooks, and report artifacts outside the tarball. See L1 for the cover, which is currently inside it. | No for size or structure |

### Dependency inventory

| Declaration | Classification and use |
| --- | --- |
| `bslib` | Imports, required installation dependency for the optional Shiny UI. Used for layout and navigation. Check the minimum UI API versions before setting a bound. |
| `ggplot2` | Imports, runtime plotting for planned-design visualizations and dashboards. |
| `mvtnorm` | Imports, numerical multivariate normal probabilities in the QTL calculations. |
| `rlang` | Imports, `.data` for tidy evaluation in plotting code. |
| `shiny` | Imports, required installation dependency for the optional `paweh_app()` interface and its modules and tests. |
| `withr` | Imports, runtime preservation of RNG state in integration and simulation. Also used directly in tests. |
| `knitr` | Suggests and VignetteBuilder, executes and builds the eight vignettes. |
| `rmarkdown` | Suggests, HTML vignette formatting and build tooling. |
| `plotly` | Suggests, optional interactive TDT surfaces and dashboard, with conditional test coverage. |
| `testthat` | Suggests, testing only, edition 3. The declared version bound needs review under C3. |
| `stats` | R distribution, explicitly imported functions for probability and numerical routines. |
| `pkgdown`, Pandoc, TeX, HTML Tidy, devtools, roxygen2 | Site, vignette, manual, validation, or maintainer tooling. None are package runtime dependencies. |

All declared packages have identified callers. No dependency was removed. No aggregate package-specific thread-pool setting was required in the tested environment. Retain checks of CPU usage and any threaded numerical libraries on supported platforms.

## Scientific findings

### Case-control association

**S1 · RESOLVED · Scientific · Does not block.** Independently derived `cc_power()` and `cc_mssn()`'s genotype Pearson comparison from the pooled mixture and the trend comparison from the two group score means and pooled score variance. The weighting uses the planned case-control allocation. Error matrices transform the true genotype frequencies to the observed distribution, and the locus-mixture model replaces a fraction of the case distribution by the control distribution. The Gordon, Edwards, and Ahn comparisons below support the stated interpretation in their tested regimes. Preserve the distinction between genotype and scored trend alternatives. [Gordon et al.](https://doi.org/10.1159/000066696), [Edwards et al.](https://doi.org/10.1186/1471-2156-6-18), [Ahn et al.](https://doi.org/10.1111/j.1469-1809.2006.00318.x)

**S2 · IMPORTANT · Scientific/API · Does not block if the supported domain is made explicit.** In `R/case_control.R`, the combined interface rejects a valid genotype-only alternative. Cases `c(.375,.25,.375)` and controls `c(.25,.5,.25)` have equal mean scores at `0,1,2`. The genotype frequency difference gives `0.9151752` genotype-test power at 100 participants per arm and alpha 0.05 in the independent calculation. Both combined entry points currently stop with “Trend numerator is 0; implies no weighted mean difference.” Decide how to represent the trend-null scenario in the combined return, including the absence of a finite trend MSSN, and test both endpoints. Review zero pooled cells as a separate choice about the supported degrees of freedom. The calculation also deserves a dedicated check of discrete allocation when the requested recruitment ratio does not produce an integer control count.

### TDT and affected-child trios

**S3 · RESOLVED · Scientific · Does not block.** Independent parental-haplotype enumeration, Mendelian transmissions, and affected-child ascertainment reproduced the expected transmitted and non-transmitted contrasts for the audited models. The central threshold uses 1 df and `lambda = (ET - ENT)^2 / (ET + ENT)`. Expected counts use two parental transmission opportunities per complete trio. Model-free inputs require ET and ENT for the stated trio count and need the additional allele-frequency or prevalence inputs only when required by the modifier. Keep canonical linked calculations, legacy helper recombination arguments, and separate modifier scenarios explicit. The phenotype-error and locus-heterogeneity sensitivity calculations are distinct from the model-based adjacent-genotype call-error branch. [Textbook §§1.6.1.4, 5.2.6 and 5.3.3](https://doi.org/10.1007/978-3-030-61121-7)

**S3a · RESOLVED · Scientific limitation · Does not block.** In the separate adjacent-genotype error branch (`R/tdt_genotype_misclassification.R`), independent finite-state calculation over parental genotype counts and Mendelian child probabilities reproduced the current model-based illustration under complete LD. At 600 retained trios with allele frequency .30, relative risks 1.5 and 2.25, and adjacent error .01, the independently derived ET and ENT are 324.9530 and 223.0449, with NCP 18.95129 and retention .986267. Applying the same error rule to the null states gives nominal-threshold rejection .057528 at nominal .05. Re-estimating the nominal cutoff was outside the implemented sensitivity model. Keep both conditioning on retained calls and the nominal null rejection rate visible in interpretation. [Textbook §5.2.5](https://doi.org/10.1007/978-3-030-61121-7)

**S3b · IMPORTANT · Scientific input validation · Does not block if the intended parameter domain is explicit.** The model-based TDT can return a result after warning that implied penetrances exceed 1. With `pd=.01, prev=.9, R1=1, R2=20`, the normalized homozygote penetrance is approximately 17.97; the genotype-error branch warns and returns power .09764 at 100 retained trios and adjacent error .01. Treat infeasible penetrances as invalid model input. Decide whether to reject such inputs uniformly across planning functions or require a clearly labeled extrapolation outside the probability model. Include tests around the upper bound.

**S4 · IMPORTANT · Scientific provenance · Does not block if reported as an unresolved absolute reproduction.** Buyske et al. use prevalence 0.01 and marker/disease-allele frequencies of 0.5 under the illustrated multiplicative model. The printed example gives homozygote relative risk 2.5, so heterozygote risk is its square root. At no false positives versus a control-to-case false-positive probability of 0.05, with true cases retained, independent enumeration gives `ET/N, ENT/N` of `0.612574113, 0.387425887` and `0.517974018, 0.482025982`. At alpha 0.00001 and target 0.8, the inverse noncentral chi-square gives `lambda*=27.6549211398`. The continuous sizes are 545.551009 and 21,400.40953 complete trios. Ceilings give 546 and 21,401. The continuous ratio is 39.227 versus the paper's 39.245. Integer rounding cannot explain the absolute gap of about 27%. The original p. 290 Eq. (5) omits the penetrance contrast in the unmodified transmitted term. At the stated parameters, reading that term literally would produce a per-parent probability of 12.75; the corrected textbook Eqs. 5.24–5.27 and independent enumeration support the current contrast.

Changing the homozygote risk to 2.25, equivalently the heterozygote risk to 1.5, produces `0.6, 0.4` and `0.515966387, 0.484033613`, yielding 691.373 and 27,120.591 before rounding. This is a strong numerical reconciliation. It does not establish the authors' actual inputs. Confirm with author code or contemporaneous calculation records if available; until then label the printed values and exact modern results separately. [Buyske et al.](https://doi.org/10.1159/000194981), [textbook §5.2.6](https://doi.org/10.1007/978-3-030-61121-7)

**S5 · IMPORTANT · Scientific provenance · Does not block with a qualification.** Chen et al. Fig. 2 uses disease and marker frequency 0.25, homozygote risk 1.5 with its square root for the heterozygote, strong linkage/LD for linked trios, and a recombination fraction of 0.5 for unlinked trios. At 0% and 50% locus heterogeneity, respectively, `ET/N, ENT/N` are `0.434846923, 0.355051026` and `0.414897949, 0.375`. The exact modern NCP gives fractional sizes 3,430.6964 and 13,722.7856, yielding 3,431 and 13,723 with ceiling. Nearest rounding of those same fractions alone gives 3,431 and 13,723 as well. With nearest-integer reporting and a shared historical NCP in `[27.6517689429, 27.6523303194]` (about 0.011% below the modern root), all four printed results can be matched if Buyske used the 2.25 risk above. For example, 27.652 gives `691.300, 27117.727, 3430.334, 13721.336`. This identifies a plausible precision and reporting mechanism for Chen. The original software settings remain unverified, and the surviving TDT-HET spreadsheet lacks the relevant alpha's pretabulated root. Keep the precise solver; qualify the discrepancy and seek the original calculation before stating its exact cause. [Chen et al.](https://doi.org/10.2202/1544-6115.1501)

### Quantitative traits

**S6 · RESOLVED · Scientific · Does not block.** Independent genotype-weighted centering and variance agree with the Falconer parameterization in `R/falconer_quantitative.R`. At the textbook §6.1.3 inputs, the parameter 0.5 is the dominance ratio `delta/a`, and the independently derived means are `-0.097462244, 0.244510542, 0.358501471` with residual SD `sqrt(.975)`. The first printed mean of `-0.0985` does not follow from the equations; the other two round as printed. Keep the equation-derived values and note the discrepancy. The expected screening burden for 126 selected per tail is `2513.738` for the upper tail and `2526.190` for the lower tail. In a common screening cohort, powers and quotients of selected count by eligibility do not guarantee enough people in both tails. A screening assurance calculation needs the joint recruitment distribution. [Textbook Ch. 6](https://doi.org/10.1007/978-3-030-61121-7), [Li et al.](https://doi.org/10.1038/s41598-019-51790-w)

**S7 · BLOCKER · Scientific · Blocks the current rounded-design claim.** `qtl_anova_power()` and the rounded design search in `qtl_anova_mssn()` use independently rounded genotype counts with nominal `N` for residual degrees of freedom. For `N=6`, `pd=.30`, dominance ratio `.5`, and QTL variance `.1`, the counts are `(3,3,1)`, totaling 7, while the reported residual df is 3. With those fixed counts the df is 4, and weighting and re-centering the means gives NCP `0.8072268821` instead of `0.8240783548`. This is an internal accounting inconsistency. Use integer allocation summing to `N`, count-weighted centering, and consistent residual df, then revisit whether the rounded-design search guarantees minimality. When a genotype group is absent, specify what design is estimable; current code rejects a zero group in the rounded calculation. Alternatively withhold the rounded design option and claims pending validation. The 996-person continuous example reproduces under its documented assumptions. The nominal and actual counts can be inspected with:

```r
check <- paweh::qtl_anova_power(N = 6, alpha = .05, qtl_var = .1,
  tau = .5, pd = .30, count_method = "rounded", verbose = FALSE)
c(N = check$N, count_total = sum(check$genotype_counts), df2 = check$df2)
```

**S8 · IMPORTANT · Scientific limitation · Does not block with clear disclosure.** Independently verified the weighted between-genotype covariance and the Pillai approximation in `R/falconer_multivariate.R`: at 325 people power is 0.7984717 and at 326 it is 0.8014862. The HWE two-copy probability is `0.05^2 = .0025`. Expected two-copy count is 0.815 and the probability of no two-copy participants with random enrollment is `(1-.0025)^326 = .4421879`. The 326 value is defensible as a reproduced planning approximation, with no assurance of an estimable three-group sample. Numerical integration of the bivariate normal mixture also independently gives 73 individuals per tail and high/low eligibility 0.0332161/0.0333035. Miwa handles suitable low dimensions and the Genz–Bretz integration uses a locally scoped seed. Make rare-cell and screening assumptions prominent across multivariate results, and add occupancy or simulation sensitivity separately if needed. [Pillai](https://doi.org/10.1214/aoms/1177728599), [Gordon et al.](https://doi.org/10.1159/000457135), [Miwa et al.](https://doi.org/10.1111/1467-9868.00382), [Genz and Bretz](https://doi.org/10.1007/978-3-642-01689-9)

### Sequencing

**S9 · IMPORTANT · Scientific scope and resource use · Does not block with clear disclosure and reasonable guards.** For `ngs_genotype_error_matrix()` and `cc_ngs_*()`, independently derived the three binomial emission distributions and integrated them over called genotypes. At depth 2 and zero symmetric base error, the transition rows are `(1,0,0), (.25,.5,.25), (0,0,1)` under the implemented call rule. The case-control trend pools the called-genotype frequencies for the null variance in the Ahn / Chapman–Nam calculation. This is a two-stage hard-call design. The trio TDT1-NGS calculation instead integrates raw-read likelihood over feasible genotypes in affected-child trios. Independent finite differences over nine mating combinations and 15 feasible genotype states gave efficient information `0.162788103730` versus `0.162788103689` for depths 6/5/4, allele frequency .325, and directional errors .01/.02. Fixed role-specific depths and constant conditional read error are design assumptions. Depth-1 and limiting-error scenarios can lose identifiability; the implementation screens unsupported boundaries rather than treating them as estimable. Trio locus heterogeneity and phenotype error scale the raw-read baseline with ordinary-TDT NCP ratios under the documented parameter bridge. These ratios are package adaptations.

Differential observation can induce apparent signal under a biological null. At `N_case=5000, pd=.30, prev=.05, R2=1, coverage=10, case_coverage=1, case_seq_error=.1, ctrl_seq_error=.01, alpha=.05`, `cc_ngs_power()` predicts nominal-threshold rejection `0.9977661` (NCP 23.06934). Surface this distinction at the entry points and in the dashboard. Profile or guard large depths: a single double matrix over 201 cubed read outcomes and 15 states uses about 0.97 GB before intermediates. [Ahn et al.](https://doi.org/10.1111/j.1469-1809.2006.00318.x), [Chapman and Nam](https://doi.org/10.2307/2528037), [Kim](https://doi.org/10.1159/000434645), [Moskvina et al.](https://doi.org/10.1159/000092553)

**S10 · RESOLVED · Reproduction scope · Does not block.** All six Kim (2015) regression values recomputed, as shown in the table. Direct published evidence was confirmed for 654 and 416 from the abstract. Treat the four remaining numbers as verified fixtures until corresponding full published rows can be inspected. Chapman–Nam DOI `10.2307/2528037` is consistent in current source and rendered site.

### LRTae validation designs

**S11 · RESOLVED · Scientific implementation and provenance · Does not block.** In `R/lrtae.R`, matrix rows index true states and columns recorded states. Fallible, phenotype-validated, genotype-validated, and both-validated observations retain their corresponding truth/recorded margins. The model uses known observation matrices, conditional independence of the phenotype and genotype recording mechanisms given truth, a common genotype distribution under the null, and nuisance coordinates for baseline genotype frequencies and true-group sampling. An independently differentiated observed-cell likelihood and nuisance projection recover the historical-example NCP `22.408269183092` versus `22.408269182990` in the package with tier counts 248/135/72/45. Independent integer evaluation confirms the example's minimum 216 at target .8 (power `.802186466`; predecessor `.799841124`). The package checks rank and estimability before using generalized inverses on singular systems. The historical workbook was unavailable. The repository carries its settings as a historical demonstration. The true-case-to-recorded-control rate there is .04, versus .05 in the related book example. The book's `24.1626` is an observed statistic. Maintain separate provenance for these quantities and request the original workbook before calling its figures an independently verified publication reproduction. [Ji et al.](https://pubmed.ncbi.nlm.nih.gov/16646856/), [Barral et al.](https://doi.org/10.1186/1471-2156-7-24), [textbook §4.2.1.4](https://doi.org/10.1007/978-3-030-61121-7)

**S12 · BLOCKER · Scientific · Blocks the claimed global minimum.** Hamilton largest-remainder allocation can move a participant out of a validation tier when total enrollment rises. The exponential bracket and binary search in `lrtae_mssn()` treat power as monotone and verify only the immediate predecessor. This fails with controls `(.75,.20,.05)`, cases `(.18,.32,.50)`, true group shares `(.5,.5)`, phenotype matrix diagonal .55 and off-diagonal .45, genotype matrix diagonal .40 and off-diagonal .30, tier shares `(3,2,3,1)/9`, alpha .05 and target .756. Enumerating Hamilton sizes in order gives the first crossing at **103** (tiers `34,23,34,12`, power `.7724066`); 104 drops to `.7386427` (tiers `35,23,35,11`); the package returns **105** (power `.7735136`). At target `.144516920589074` the same design first crosses at 13, drops at 14 and returns 15. Both counterexamples survive checking one predecessor. Use an exact bounded integer search with precomputed per-tier information or a search with a proof that every smaller size is excluded. Keep Hamilton and its tie-break fixed when validating the change. Add the full witness and singular scenarios as regression cases.

```r
library(paweh)
inputs <- list(
  c(.75, .20, .05), c(.18, .32, .50), c(.5, .5),
  matrix(c(.55, .45, .45, .55), 2, byrow = TRUE),
  matrix(.30, 3, 3) + diag(.10, 3)
)
tiers <- c("fallible", "pheno_validated", "geno_validated", "both_validated")
shares <- c(3, 2, 3, 1) / 9
fit <- do.call(lrtae_mssn, c(inputs, list(alpha = .05, target_power = .756),
  as.list(setNames(shares, paste0("prop_", tiers)))))
power_at <- function(n) {
  ideal <- n * shares
  counts <- floor(ideal)
  top <- head(order(-(ideal - counts), seq_along(shares)), n - sum(counts))
  counts[top] <- counts[top] + 1L
  do.call(lrtae_power, c(inputs, list(alpha = .05),
    as.list(setNames(counts, paste0("n_", tiers)))))$power
}
fit$MSSN_total
which(vapply(1:fit$MSSN_total, power_at, numeric(1)) >= .756)[1]
vapply(103:105, power_at, numeric(1))
```

## Published-result reproduction

The primary checks used independent counts, moments, quadrature, or likelihood derivatives where stated. Regression fixtures and paper/title or abstract checks provide different levels of evidence.

| Source and example | Published or historical value | Independent calculation or tested output | Assessment |
| --- | --- | --- | --- |
| Gordon et al. (2002), Table 1b / Fig. 1b | 386 and 477 | 386.07774 and 476.54765 before ceiling; 387 and 477 whole participants per planned group | Approximate, one participant difference in baseline |
| Edwards et al. (2005), Fig. 2 | About 91%, 33%, 76%, 11% | 91.349%, 33.205%, 76.339%, 10.983% | Agreement at published precision |
| Ahn et al. ordinary case-control trend and genotype illustrations, p. 257 | About .88/.81 and .79/.70 | .874268/.806084 and .787592/.699231 | Approximate numerical agreement |
| Buyske et al. (2009) | 691 and 27,118 trios | 546 and 21,401 with stated risk, plus numerical reconciliation in S4–S5 | Relative increase reproduced; absolute figures unresolved |
| Chen et al. (2009), Fig. 2 | 3,430 and 13,721 trios | 3,431 and 13,723 | Near agreement with the same model; precision explanation plausible |
| Textbook §6.1.3 | Means −.0985, .2445, .3585 | −.097462, .244511, .358501 from the equations | Discrepant lower-genotype mean |
| Textbook continuous and extreme QTL examples | 996 continuous participants; 126 per tail in the audited extreme scenario | Power .8000478 at 996; 126 per tail independently | Reproduced under the stated approximations |
| Gordon et al. multivariate QTL | 326 for full-population design | .7984717 at 325 and .8014862 at 326 | Minimum reproduced under the Pillai approximation |
| Bivariate threshold scenario | Tutorial 73 per tail | Independent normal-mixture quadrature gives 73 per tail | Model output verified; no separate paper result asserted |
| Kim (2015), three stored comparisons | 654/416 directly checked in the abstract; additional stored pairs 716/456 and 733/466 | 654/416, 716/456, 733/466 | First pair checked against publication; later pairs rerun against fixtures |
| Ji/Barral and historical workbook example | Published method and a later historical scenario | NCP 22.40826918 and MSSN 216 with repository inputs | Independent implementation agreement; no verified exact publication benchmark for those inputs |

## Documentation consistency

All eight articles were built from source and re-executed in package checks. Reviewed visible return values, explanatory text and figures, equations and method links, and article/reference navigation. The following actions reflect findings from that review.

| Article | Finding, evidence, and action | Severity / blocks |
| --- | --- | --- |
| 01 Getting Started | The power and MSSN walkthrough builds and distinguishes likely detection from planned recruitment. Keep installation tied to the development branch until a CRAN release exists. Review which worked example is highlighted if any policy decision limits rounded or validation designs. | RESOLVED / No |
| 02 Case-Control | Code demonstrates inputs and returned quantities. Clarify zero trend signal, fixed scores, and nominal rejection under differential error. After S2, demonstrate the genotype-only alternative. | IMPORTANT / See S2 |
| 03 TDT | Trios and transmitted/non-transmitted counts are explained. Add the evidenced S4–S5 precision and parameterization qualification and distinguish complete trios from individuals in every result. Keep the conditional population/ascertainment assumptions and separate sensitivity scenarios visible. | IMPORTANT / See S4–S5 |
| 04 Quantitative Traits | Textbook lower mean differs from the equations, and the 326-person example risks over-interpreting a rare group. Keep the expected-count versus realized-design distinction near the tables and resolve rounded-count behavior before teaching that design as fixed enrollment. | BLOCKER / See S7; S8 is important |
| 05 Interactive Dashboard | Module examples and figures build. Controls and model summaries reflect the code inspected. Check labels and help once S2, S7, and S9 decisions are settled, and confirm the current snapshots in Windows and Linux browser smoke tests before a release. | MINOR / No |
| 06 Sequencing | Distinguishes the hard-call and raw-read models and uses the correct Chapman–Nam DOI. Give nominal type I implications and depth limits closer to the inputs, and separate the first Kim pair from the four rerun fixtures in benchmark labels. | IMPORTANT / See S9 |
| 07 Validation-Based / LRTae | Four validation tiers, orientation, and the observed-versus-expected distinction are documented. Qualify the workbook provenance and replace the global minimum claim only after S12 is addressed. | BLOCKER / See S12 |
| 08 Methods & Sources | Paper and textbook origins are distinguished across the five families. Update the evidence captions for the edge cases and historical examples using the reproduction table. | IMPORTANT / Scientific claims follow S4–S12 |

**D1 · RESOLVED · Documentation · Does not block.** Corrected `theta1` for the three low-level TDT helpers. It is a recombination fraction in `[0, .5]`. The legacy default equals `pd`; zero gives the canonical linked setting. Corrected the legacy power helper's result names (including `Power`) and the pointer to §1.6.1.4, Eq. 1.25. The canonical TDT help now uses the same section pointer. The existing verbose legacy message still labels `theta1` as a population allele frequency. Fix that message in a separate reviewable runtime change. This audit regenerated the six affected help pages from roxygen.

## Licensing and public website

**L1 · BLOCKER · CRAN redistribution · Blocks submission.** `man/figures/textbook-cover.jpg` (102,878 bytes) is displayed on the site, referenced by README and included in the real source archive. The image and book identification link to Springer, and the site supplies alt text. Author permission alone does not establish permission from the artwork rights holder for website distribution, CRAN redistribution, mirrors, and perpetual archiving. No publisher permission covering those uses was available for review. Obtain written authorization that permits the planned distribution and record the asset-specific attribution/licensing terms, or remove or replace the cover with material cleared for the intended uses. Test README and homepage after the resolution. Review existing public copy as part of that decision. [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html), [Springer Nature permissions](https://support.springernature.com/en/support/solutions/articles/6000085109-permissions-requests)

**L2 · MINOR · Project attribution · Does not block on current evidence.** The R `LICENSE` template names Akil Anthony as copyright holder and `LICENSE.md` gives the MIT text. Confirm contributor ownership and any distinct asset terms before tagging. The text and screenshots generated for the package have identified project origins; the cover has a separate rights review in L1.

**W1 · RESOLVED · Site and workflow · Does not block.** Rebuilt the committed `447cd8c` source with pkgdown 2.2.1. A structural check found 80 HTML pages and checked 2,240 local anchors, 44 image references, 650 script paths, and 275 linked resources, including local fragments. All eight articles appear in the index, and each of the 58 exports appears once in the reference index and resolves to a page or an alias redirect. Verified the home example, textbook asset, Methods & Sources entry and Chapman–Nam link. The live homepage returned HTTP 200 with the redesigned worked example. The published Pages revision remains `d41430c813c3a8738b16f2f448e304eea6803c8d`. The website is built by ordinary pushes or pull requests unless CI is explicitly skipped. Both the build-artifact step and deployment job require `workflow_dispatch` on `refs/heads/develop`; the deploy job also needs a successful build and `PAWEH_PAGES_TOKEN`, with force disabled. The audit report is outside the site build. [Live site](https://akilanthony.github.io/paweh/), [GitHub dispatch documentation](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow)

## Technical validation

The main environment was macOS Tahoe 26.6.2 on Apple silicon (aarch64), UTF-8, with Pandoc 3.8.3 and TeX Live 2024. R 4.6.1 used bslib 0.12.0, ggplot2 4.0.3, mvtnorm 1.4.2, rlang 1.3.0, shiny 1.14.0, withr 3.0.3, knitr 1.51, rmarkdown 2.31, plotly 4.12.1 and testthat 3.3.2. The older local R 4.4.1 environment included bslib 0.9.0, shiny 1.13.0, mvtnorm 1.3.3, ggplot2 4.0.2 and testthat 3.3.1.

| Platform and R | Operation | Result and exact diagnostic |
| --- | --- | --- |
| macOS arm64, R 4.6.1 | `devtools::document()`, standalone `devtools::test()`, source build and staged `R CMD INSTALL` | Help regenerated. Standalone tests: **10,324 pass, 0 fail, 0 warn, 0 skip** in 112.0 s. Full source build produced all eight HTML vignettes plus their sources. Installation passed and loaded from temporary and final locations. |
| macOS arm64, R 4.6.1 | `devtools::check()` invoking full `R CMD check --as-cran`, remote incoming on, system Tidy | **0 errors, 0 warnings, 2 notes.** `New submission` and `Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.` Tests, examples, PDF manual (421,228 bytes), and vignette re-build passed. |
| macOS arm64, R 4.6.1 | Same archive, `R CMD check --as-cran --timings`, Tidy 5.8.0, remote incoming disabled | **0 errors, 0 warnings, 1 note.** `Unable to verify current time. To disable remote verification, set environment variable _R_CHECK_SYSTEM_CLOCK_ to a false value.` The HTML manual, PDF manual, tests and rebuilt vignettes passed. |
| macOS arm64, R 4.4.1 | Same archive, `R CMD check --as-cran --timings`, Tidy 5.8.0, remote incoming disabled | **0 errors, 0 warnings, 2 notes.** `unable to verify current time`, and `Author field differs from that derived from Authors@R`. The newer build formats the URL with an `ORCID:` prefix, which older R omits during derivation. Repository index lookups emitted network-access warnings in this run, outside the final check-warning count; locally installed dependencies resolved. PDF/HTML manuals, tests and rebuilt vignettes passed. Rebuild the archive with the target R version when investigating ORCID-format differences. |
| macOS arm64, current R 4.6.2 | Same archive, `R CMD check --as-cran --timings`, Tidy 5.8.0, remote incoming on | First run: **0 errors, 0 warnings, 2 notes.** `New submission` and `Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.` The packaged R was run from a temporary extraction of the official CRAN download, which matched CRAN's published checksum. |
| macOS arm64, current R 4.6.2 | Repeat with Pandoc and Tidy 5.8.0 on PATH, remote incoming on | **0 errors, 0 warnings, 1 note.** `New submission`. README/NEWS validation, all 10,324 assertions, all eight vignette rebuilds, PDF and HTML manual validation, and clock verification passed. |
| R-devel; Windows, Linux and Intel macOS | Full CRAN check matrix | Unrun. No external R-hub or Win-builder job was started. Run Win-builder release/devel and add Linux release/devel and Intel macOS checks where practical on the final reviewed candidate. |

The full checked tarball's SHA-256 is `d60f0ba2be1b4b6e75c473b38b4375e35a5240ee79eeb7b375705a5c9d035507`. A build with the audit file present and vignette building disabled confirmed its exclusion and the continued inclusion of the cover. Relative paths, package loading, and native-encoding checks passed in the full local checks. A separate C-locale smoke test loaded the installed package, found 58 exports, and executed the case-control power example. Cross-platform floating-point, file-system, and locale behavior remains open. Parsed expressions in both touched R files exactly match those in `447cd8c`, establishing that the audit changed comments only. The new-submission note is expected for an initial submission. Modern Tidy removes the system-Tidy note. [CRAN current R sources](https://cran.r-project.org/sources.html), [HTML Tidy 5.8.0](https://github.com/htacg/tidy-html5/releases/tag/5.8.0)

## Prioritized release plan

1. **Resolve the three blockers.** Review S12's search proof and regression case, S7's genotype-count and NCP accounting, and L1's explicit rights evidence or cover replacement. Confirm that the resulting public claims match the methods that ship.
2. **Complete the pre-release scientific and portability work.** Document the Buyske and Chen evidence and remaining provenance uncertainty, rare genotype occupancy, and the null-rejection example. Resolve or scope the valid genotype-only case and the handling of infeasible penetrances. Set and test dependency floors. Run exact candidate checks on current R and R-devel and at least one additional major platform. Examine package and example timings and any notes in those runs.
3. **Optional follow-up.** Obtain original author scripts and workbooks, add finite-sample QTL and screening assurance sensitivity, and improve sequencing memory use. Fix the legacy `theta1` verbose label and review screenshot aging after interface changes.

### Version and Git strategy

Recommend **1.0.0** for the initial CRAN release once the reviewed API and assumptions are ready to support. `0.99.0` is valid and remains the development version in this audit. Finish distinct, reviewable fixes on `develop`. Approve a candidate commit and run all final checks on its release tarball. After separate approval, bring the verified commit to `main` via PR, confirm the merge tree matches the tested tree, and create an annotated `v1.0.0` tag on the verified release commit. Build from that tag and record the tarball SHA-256. Submit the source tarball only after the maintainer signs off. After acceptance, update `develop` from the released source and increment to a development suffix such as `1.0.0.9000`; maintain NEWS and tagged release history.

### Final CRAN submission checklist

- [ ] Resolve S7, S12, and L1 and verify the approved fixes.
- [ ] Resolve or explicitly document the S2, S4, S5, S8, S9 and S3b scientific limits, metadata/version changes, and dependency minimums.
- [ ] Run current R and R-devel candidate checks across platforms, including Windows or Linux, using the submitted source tarball where practical. [Win-builder](https://win-builder.r-project.org/)
- [ ] Run tests, examples, source build, eight vignette builds, PDF and current-Tidy validation, and review every surviving check note.
- [ ] Recheck redistribution rights, file and documentation sizes, and source exclusions. Record model assumptions and review the citation, titles, URLs, NEWS, and contributor metadata.
- [ ] With explicit approval, merge and tag the reviewed candidate, verify the tarball checksum, and submit as a new package. Describe how each final check note was investigated in submission comments.

## Exact audit files

- `CRAN_FINAL_READINESS_AUDIT.md`, `.Rbuildignore`, `.github/workflows/pkgdown.yaml`, and `pkgdown/README.md` add the audit and ensure it remains a repository document.
- `R/tdt_conditional.R` and `R/tdt_functions.R` correct roxygen prose. `man/tdt_expected_nontransmission_probability.Rd`, `man/tdt_expected_transmission_counts.Rd`, `man/tdt_expected_transmission_probability.Rd`, `man/tdt_power.Rd`, `man/tdt_power_from_expected_counts.Rd`, and `man/tdt_power_from_model.Rd` are regenerated help.
