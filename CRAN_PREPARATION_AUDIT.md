# CRAN preparation audit: Phase A1

Audit date: 13 November 2026. Scope: minimal CRAN compatibility cleanup,
dependency review, and an educational documentation plan. The statistical
algorithms and the substantive vignette content are unchanged.

## 1. Provenance and branch preservation

- Fetched personal `origin`, verified as `https://github.com/akilanthony/paweh.git`.
- The clean post-interoperability-removal baseline is
  `ce0aa46ef84c82e1bfbcc3fd83e22a7a500d7166` on `bioconductor-revision`.
  Local and remote copies agreed; the working tree was clean.
- Existing `develop` was at `b466e867a59b8a7e98798ca68495558b3dc4916d`.
  Its history diverged: four develop-side commits and 23 baseline-side commits.
  Resetting or force-pushing would have lost development history.
- Merge `c75c4f4` preserves both histories and the develop-only project,
  license, instruction, pkgdown, and deployment files. The only merge conflict
  was README's six-versus-seven article count and LRTae links; the validated
  seven-article wording was retained.
- After the merge, `R/`, `man/`, `tests/`, NAMESPACE, DESCRIPTION and all seven
  vignettes matched the post-removal baseline. A1 subsequently changes only
  the files listed below. Stable `main` is unchanged.
- The existing workflow publishes to `gh-pages` on develop pushes. Commit
  messages use `[skip ci]` for this delivery so pushing only `develop` does
  not implicitly publish a website. No site deployment or submission is part
  of A1.

## 2. Official guidance used

Reviewed the current official sources, rather than treating check success as
a guarantee of CRAN acceptance:

- [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html),
  revision 6985: publication quality, dependency availability, resource use,
  network behavior, and author-year references with identifiers.
- [CRAN submission checklist](https://cran.r-project.org/web/packages/submission_checklist.html):
  informative metadata, DOI syntax, `--as-cran`, and additional platform checks.
- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html),
  sections 1.1.1, 1.1.3, 1.3, 1.4 and 1.9: metadata, dependency declarations,
  checking, executable vignettes, and citation infrastructure.
- [rmarkdown's html_vignette reference](https://rmarkdown.rstudio.com/docs/reference/html_vignette.html):
  a lightweight format intended for package vignettes, including CRAN use.
- [Publisher's textbook record](https://link.springer.com/book/10.1007/978-3-030-61121-7):
  Gordon, Finch and Kim, Springer, 2020.

These references support the compatibility decisions below. They do not
independently validate the package's statistical implementations.

## 3. Exact A1 file changes

| File | Change |
|---|---|
| `DESCRIPTION` | Remove `biocViews` and BiocStyle; shorten Title; quote the external software name 'shiny'; add a scoped author-year textbook reference with `<doi:...>` syntax. |
| `.Rbuildignore` | Exclude this developer audit from the source package. Existing develop-only exclusions are preserved. |
| `inst/CITATION` (new) | Generate the software citation from package metadata and provide a separately scoped background book entry. |
| `README.md` | Add `citation("paweh")` and the book reference with a method-specific attribution caveat. |
| `NEWS.md` | Record the compatibility and citation changes without changing the version. |
| `vignettes/paweh-01-getting-started.Rmd` | Change output header only. |
| `vignettes/paweh-02-case-control-study-design.Rmd` | Change output header only. |
| `vignettes/paweh-03-tdt-study-design.Rmd` | Change output header only. |
| `vignettes/paweh-04-quantitative-trait-study-design.Rmd` | Change output header only. |
| `vignettes/paweh-05-interactive-dashboard.Rmd` | Change output header only. |
| `vignettes/paweh-06-sequencing-study-design.Rmd` | Change output header only. |
| `vignettes/paweh-07-lrtae-study-design.Rmd` | Change output header only. |
| `CRAN_PREPARATION_AUDIT.md` (new) | This audit and future revision plan. |

The vignette format is now `rmarkdown::html_vignette`, with a table of contents
and numbered sections. Floating Bootstrap navigation is removed. The bodies,
code, equations, tables, captions, image inputs, and figure settings are
byte-identical to `ce0aa46`. No examples or mathematical explanations were
rewritten. NAMESPACE was regenerated with roxygen2 and is unchanged. No R
implementation, test, or man-page change is introduced by A1.

BiocStyle itself is not inherently prohibited by CRAN; this change removes an
unneeded Bioconductor formatting dependency and uses a lighter conventional
format. The removed interoperability functions, helper, tests and Applications
vignette have not been restored.

## 4. DESCRIPTION and citation assessment

The old Title was informative but long enough to be truncated in package
listings. The new Title is **Power and Sample Size for Genetic Association
Studies**: title case, plain text, no trailing period or repeated package
name. Heterogeneity and misclassification remain explicit in Description.

Description remains one paragraph of complete sentences. It states supported
designs and the prospective scope, expands MSSN, and uses 'shiny' as an
external software name. Terminology such as locus heterogeneity and phenotype
misclassification still warrants the future glossary proposed below, but A1
does not make a substantive explanatory rewrite of the metadata.

The new reference is explicitly one source among others:
Gordon, Finch and Kim (2020) `<doi:10.1007/978-3-030-61121-7>`.
The citation does not attribute every implemented method to the textbook.
`inst/CITATION` follows the official `citation(auto = meta)` plus `bibentry()`
pattern. It neither invents a package paper nor calls `packageDescription()`
while the citation file is being read. Method-specific papers in the existing
help pages and vignettes remain necessary.

Version `0.99.0` is unchanged. That version format is valid for CRAN; changing
it is a release decision, not a necessary A1 compatibility fix. Author roles,
maintainer address, ORCID and MIT licensing are preserved. Their substantive
accuracy and permission to distribute third-party material require the
maintainer's eventual release review.

## 5. Dependency inventory

This inventory covers every declared direct dependency before and after A1.
Classification describes the current package structure, not a proposed API
redesign. Namespaced calls and imported symbols were inspected in R code,
examples, tests and vignettes.

| Dependency | Declaration / classification | Evidence and decision |
|---|---|---|
| bslib | Imports; required runtime for the dashboard | Themes, cards, sidebars and navigation in `R/shiny-*.R`. Retain; making the UI optional at installation would need a separate loading/API audit. |
| ggplot2 | Imports; required runtime | Public plotting functions, dashboard plots and imported plot constructors. Retain. |
| mvtnorm | Imports; required runtime | Multivariate tail probabilities in `R/falconer_multivariate.R`, including Miwa/Genz-Bretz integration; regression tests also use it. Retain. |
| rlang | Imports; required runtime | Imported `.data` pronoun in plot aesthetics. Retain even though a simple `rlang::` search alone would miss the import. |
| shiny | Imports; required runtime for the dashboard | Namespaced UI/server/application objects and module tests. Retain; the dashboard is optional to use, but currently not an optional installation dependency. |
| withr | Imports; required runtime | Locally seeded integration and genotype-distribution plotting; tests verify preservation of the caller's random state. Retain. |
| plotly | Suggests; optional runtime | `plot_tdt_surface3d()` and the optional dashboard TDT surface. Availability is checked; relevant examples/tests are guarded. Retain. |
| knitr | Suggests and VignetteBuilder; documentation-only | Executes R Markdown and includes dashboard images. Retain. |
| rmarkdown | Suggests; documentation-only | Existing vignette engine and now the output format. Retain. |
| testthat (>= 3.0.0) | Suggests; testing-only | Test harness, expectations and `local_mocked_bindings()`. Retain. The tested local version is recorded below; verify that the declared lower bound covers these helpers in a later minimum-version check. |
| BiocStyle | Former Suggests; unused after format conversion | Its only source uses were the seven output headers. Removed after replacing all seven. |

Base/recommended R facilities such as stats, utils, graphics and grDevices are
also used. Existing stats imports remain in NAMESPACE. No new external
runtime package is required for the citation file.

Development tools are distinct from package dependencies: roxygen2 and
devtools generate/check the package; pkgdown and r-lib/actions build the
website; remotes appears only in non-executed installation instructions.
These uses do not justify adding those tools to Imports. No unused remaining
declared dependency was identified. SummarizedExperiment, S4Vectors and an
explicit methods dependency remain absent.

## 6. Professor-authored reference draft

Inspected `/Users/akilanthony19/Downloads/paweh-revised-docs.zip` as reference
material only. It contains README, DESCRIPTION and seven vignette drafts.
Copies were read outside the checkout; no archive file replaced current
source and no code from it was executed.

Useful ideas for the later educational phase include introductory clinical
questions, a genotype/allele/prevalence glossary, parameter-meaning tables,
annotated first calculations, clear sample units, and explanations of why
an error mechanism changes the design. The trio draft supplies a helpful
plain-language account of detected versus Mendelian-consistent genotype
errors. The QTL draft explains variance explained, dominance and screening.

Do not restore the draft's Bioconductor submission/install promises,
SummarizedExperiment/VCF walkthrough, dependency additions, or container
integration claims. Its Chapter 5 references need verification against the
book and current help pages. Its broad claims that errors always reduce power
or increase sample size also need qualification: differential mechanisms can
inflate nominal signal without valid Type I error control. Preserve the
current distinctions between independent scenarios, method reproduction,
methodological context, and motivating applications.

## 7. Proposed README improvements (later phase)

1. Add a short reader-oriented choice: cases/controls, complete affected-child
   trios, continuous trait, extreme selection, or validation subsampling.
   Keep NGS and LRTae routes visible rather than embedding them only in prose.
2. Explain the first model's `pd`, `prev`, `R2`, `MOI`, `alpha` and `k` next to
   the code. State that the 457-case/457-control result is conditional on those
   assumptions, and distinguish case counts from total recruitment.
3. Add an interpretation sentence after the displayed 0.8366557 power rather
   than relying on readers to convert a bare numeric output to a decision.
4. Introduce or link NCP, HWE, locus heterogeneity, genotype/phenotype error,
   ANOVA, Pillai and LD on first use. Expand AIS and distinguish risk from an
   odds ratio when discussing the Kim example.
5. The GitHub installation command currently installs default `main`, which
   is valid but does not include this develop work. Later distinguish stable
   GitHub installation from `ref = "develop"`. Add CRAN installation only
   after the package is actually available there; do not promise acceptance.
6. Audit the claim that the website is the complete/current documentation
   against the deployed revision. Link directly to `citation("paweh")` and
   method-specific references without treating the background book as a
   universal source.

Acceptance: a reader can identify a function, interpret one numerical design,
and install the intended version without reading the source or textbook.

## 8. Educational revision plan, in vignette order

### 01. Getting Started

Evidence: `tiny-example` prints power without interpreting its value. The
opening lists relative risk, heterogeneity, extreme sampling and correlation
without introductory definitions. The five-family list includes sequencing
and LRTae but the novice selection path remains mostly function names.

1. Introduce a short planning question and a shared glossary using the useful
   draft terminology, checked for scientific scope. Explain genotype order
   0/1/2, allele frequency versus disease prevalence, and risk relative to the
   zero-copy genotype.
2. Explain the source and meaning of `g_case` and `g_control` before
   `tiny-example`; state group-specific frequencies sum to one. Interpret the
   displayed power as repeated-study detection probability under fixed inputs.
3. Provide one short inverse MSSN call and interpret cases, controls and total
   people, including integer rounding. Retain separate prospective and
   association-testing workflows.
4. Explain where model-free summaries can come from using numeric vectors,
   without resurrecting container wrappers or file readers.
5. Add a conceptual map to the seven articles and link the applicable
   methodological references. Installation currently uses GitHub correctly;
   clarify release/branch choice when the release plan is settled.

Acceptance: the first complete example explains every input and output and
requires no knowledge of the specialist vignettes.

### 02. Case-Control Study Design

Evidence: `basic-cc` provides a result table without an explicit recruitment
recommendation. `gordon-reproduction` uses the acronym DSB without expanding
its paper-specific meaning. `diff_error_design$tests$genotypes$power` and
`independent-cc-scenarios` leave readers to interpret a scalar or nested list.

1. Define genotype chi-square versus trend tests, degrees of freedom and
   score vector `w` before comparing results. Explain the multiplicative
   relation between `R1` and `R2`, and the control-to-case meaning of `k`.
2. Add a small input table for the first design and a numerical interpretation
   after it. Distinguish model assumptions from fitted/pilot estimates.
3. Expand DSB from the original paper after checking its exact terminology.
   Explain the error matrix and `e2 = epsilon / 2` mapping in words while
   preserving the 387-versus-approximately-386 rounding discussion and the
   477-case error reproduction.
4. Explain phenotype false-positive/false-negative directions before Edwards;
   retain the rare-disease prevalence explanation and add a sentence comparing
   the displayed powers.
5. Explain `pi` versus `1 - pi` beside the heterogeneity plot and turn the
   independent-scenario output into a narrated comparison in the later phase.
6. Retain the differential-error warning. A larger nominal power is not proof
   of improved valid detection when the null distribution is distorted.
7. Add precise citations for the introductory genotype/trend methods and a
   scoped textbook link where helpful; keep the existing primary papers.

Acceptance: a reader can distinguish the error directions, choose a test,
identify the recruitment unit and explain each sensitivity scenario.

### 03. TDT Study Design

Evidence: `basic-tdt` jumps to `R1`, `R2`, `delta_prime`, rates and NCP.
`expected-counts` presents gT/gNT and counts with little explanation of their
denominators. The heading and several descriptions say only "Chapter 5".

1. Define complete trio, informative transmission, heterozygote, population
   stratification and positive linkage disequilibrium. Explain why one trio
   represents three people without changing the package's recruitment unit.
2. Explain the first inputs and interpret the three rows of `basic-tdt`.
   Introduce NCP as a measure of expected test separation before the equation.
3. Clarify the normalization of gT/gNT and their conversion to expected counts;
   do not imply they must sum to one. Distinguish `N` from model-free
   `n_trios`, and interpret `model-free-power` equivalence.
4. Replace vague Chapter 5 prose with the named Gordon/Finch/Kim book and a
   verified section/equation/page reference. The draft points to section
   5.2.5; verify it rather than mechanically importing it.
   The existing adjacent-error paragraph also writes `\le` outside math
   delimiters; correct that inline typesetting in the educational phase.
5. Introduce adjacent errors, Mendelian-inconsistent removal and retained
   consistent trios before the matrix. Explain why this differs from ordinary
   case-control error and from sequencing errors.
6. Preserve the Buyske/Chen reproduced values, separate mechanisms and
   percent-versus-fold explanation. Tie each curve back to a recruitment
   decision and keep the optional 3D example explicitly interactive.

Acceptance: readers can trace assumptions to transmission quantities and trio
counts and understand what the genotype-error model does and excludes.

### 04. Quantitative-Trait Study Design

Evidence: the introduction invokes a Falconer-style model before defining a
QTL in everyday language. `anova-design` prints an unlabeled-by-context table;
`gordon_mv_power$power` is a bare scalar. Pillai, pleiotropy, Miwa and
Genz-Bretz are advanced first-use terms. Joint threshold selection is explained
mathematically but has no complete worked power calculation.

1. Introduce a familiar continuous phenotype and explain ANOVA, QTL and
   pleiotropy. Use the draft's three-design distinction as an outline, not a
   source overwrite.
2. Explain `qtl_var = 0.025` as 2.5% of modeled variance, the scale and sign of
   `tau`, and residual versus total correlation. Cross-check definitions
   against current help/code before simplifying them.
3. Interpret the genotype-distribution figure in terms of overlap. Explain
   deterministic bins separately from simulated observations.
4. Explain `count_method`, `multiple_of_three`, the chosen alpha and why `N`
   is total people. Add a numerical conclusion after ANOVA and its curve.
5. Preserve the distinction between selected sample size and screening burden
   and explain upper/lower percentages in `threshold-parameters`.
6. Explain Pillai's test and the 326-person Gordon reproduction in plain
   language; interpret the returned achieved power and rounding.
7. In the later phase add a small executable joint-threshold example and an
   interpretation of AND selection, mixed regions and screening effort.
8. Identify an exact source for the Falconer parameterization; the current
   reference list covers other components but has no full Falconer reference.
   Verify any new bibliography entry rather than guessing it. Keep numerical
   integration details in a clearly marked advanced subsection.

Acceptance: readers can select a QTL design and distinguish modeled variance,
selected recruitment and source-population screening.

### 05. Interactive Dashboard

Evidence: the screenshots and navigation instructions are useful, but the
Case-Control `cc-example` table has no explicit numerical interpretation.
"Genotype misclassification (Chapter 5)" refers to a control label without
naming the book locally. No scientific reference section is present.

1. Retain screenshots, the four-step workflow and frozen-result behavior.
   Add a single end-to-end example connecting a GUI value to its R argument,
   output and practical interpretation.
2. Explain proportions versus GUI percentages, case/control allocation, R1/R2
   and D-prime through short tool-to-article links. Do not change the UI in
   this documentation phase.
3. Identify the textbook behind the exact Chapter 5 control label and link to
   the verified TDT section. Preserve all incompatibility and scenario rules.
4. Explain the units of the QTL console results and what a reader should take
   from a minimum-sample-size result. Keep the current limitation that QTL
   does not have the richer Case-Control/TDT tabs.
5. Add local scientific references or precise links to the method articles.
   Do not imply the dashboard implements sequencing or LRTae workflows.
6. Check screenshot freshness and accessibility after any later UI change;
   keep `runApp()` and optional 3D interactions unevaluated in static builds.

Acceptance: users can reproduce one GUI design in R and interpret its units
without treating advanced diagnostics as a second statistical method.

### 06. Sequencing-Based Study Design

Evidence: the opening immediately uses Ahn/Chapman-Nam and efficient-information
NCPs. HWE, MAF, OR and BGE need clearer first-use explanation. Several chunks
print implementation lists (`ordinary_tdt_bridge`, `$sequencing`, phenotype
errors, and genotype-error matrices) with no small decision-oriented summary.

1. Start with depth versus sample-count planning and explain coverage,
   per-read error, genotype calls and the raw-read/called-genotype distinction.
   Preserve the two-workflow comparison and state sampling units early.
2. Provide concise parameter tables, using the draft as a reference, for
   case-control and trio designs. Explain why CC uses `R2` and TDT1-NGS uses
   `R1`, and explain directional read errors before their formula.
3. Expand HWE, MAF and OR before the Kim AIS example. Retain the exact six
   reproduced integer results and qualify the relationship between the
   reported effect and this model's relative-risk parameter.
4. Explain how 654/416 trios translate to complete families and people, and
   interpret `tdt-ngs-workflow` and `cc-ngs-workflow` outputs numerically.
5. Put attenuation-bridge equations and detailed lists in an advanced section
   in the later phase; first state the design consequence and which fields
   matter. Preserve the distinction from a new error-aware likelihood.
6. Narrate a transition-matrix row and explain finite-depth heterozygote
   ambiguity at zero per-read error. Maintain the independent scenario paths
   and differential-error calibration caveat.
7. Interpret each coverage curve, the `pi = 0` no-finite-MSSN boundary and the
   one-read trio identifiability boundary. Define BGE if the acronym remains;
   do not advertise deferred features as available.
8. Preserve exact attribution: the CC sequencing bridge is the package's
   composition with the Ahn method, not an example claimed to be published by
   Ahn. Keep method-specific citations adjacent to the relevant workflow.

Acceptance: readers understand the effect of a depth/error assumption and do
not confuse called genotypes, raw reads, trio counts and participant counts.

### 07. LRTae Power and Sample-Size Design

Evidence: the tier and matrix-orientation tables are strong. The example's
`p_control`, `q` and matrices come from a historical workbook whose provenance
is not independently accessible in the vignette. The NCP explanation advances
quickly to efficient information and a Schur complement. Power and MSSN lists
are printed without a concise numerical recruitment conclusion.

1. Expand the historical method name using the primary source and explain
   gold-standard validation versus statistical validation of software. Retain
   the distinction from fitting an observed-data likelihood or running LTTae.
2. Explain every input vector and matrix, category order, sampling fractions
   and why the four tiers are mutually exclusive. Label workbook values as
   deterministic design assumptions, not patient data or universal rates.
3. Document a distributable provenance/reference for the workbook if possible;
   otherwise use a self-contained synthetic rationale in the later rewrite
   without changing validated numeric examples silently.
4. Present one numerical power and recruitment conclusion before advanced
   information-matrix algebra. Define nuisance parameters, identifiability
   and Schur complement only where that deeper explanation is needed.
5. Explain Hamilton allocation with one numerical rounding example and show
   that tier counts sum to the total. Distinguish its MSSN allocation rule
   from a global cost-optimal design.
6. Interpret the four-allocation comparison and 80% reference line; explain
   that singular designs and validation requirements are assumptions of this
   model, not universal rules for every misclassification method.
7. Preserve Ji/Barral and textbook references, matching each to prospective
   information calculations versus the historical association-fitting method.

Acceptance: readers can define a validation allocation, interpret total/tier
recruitment, and state which measurement and identifiability assumptions apply.

## 9. Website and project-infrastructure audit

The newer baseline did not contain website source configuration, but the
existing develop branch did. Preserving that branch restored `_pkgdown.yml`
and `.github/workflows/pkgdown.yaml`; these are excluded from the CRAN tarball.

Findings and recommendations, intentionally not a website redesign:

1. `_pkgdown.yml` still names removed `plot_qtl_multivariate_contour` and
   `plot_qtl_multivariate_surface3d` references. Remove those obsolete entries
   in a dedicated website correction; do not restore removed statistical or
   plotting functionality to satisfy navigation.
2. Its article list contains six entries and omits the seventh LRTae article.
   Its reference categories also omit `lrtae_ncp`, `lrtae_power` and
   `lrtae_mssn`. Bring navigation into agreement with the current seven
   vignettes and exports before the next site publication.
3. Recommend a short Start Here route followed by study-family groups, with
   sequencing and validation studies visible. Keep lower-level statistical
   helpers separate from canonical entry points. Do not change order/content
   during A1.
4. Fetched `origin/gh-pages` is at `057d9d3`, a deployment identified as coming
   from `d2e3a6b`, before interoperability withdrawal. A later website phase
   should rebuild from an explicitly selected, verified source and check for
   obsolete installation/submission wording and broken links. A1 does not
   publish or claim to have updated the live site.
5. The workflow runs on develop pushes, pull requests and manual dispatch;
   it deploys non-PR builds to `gh-pages`. Keep build validation distinct from
   publication, make the source branch explicit and consider an approval or
   manual release gate in a later deployment decision. A1 uses `[skip ci]` to
   honor the requested develop-only delivery, using GitHub's documented
   [workflow-skip behavior](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/skip-workflow-runs).
   This does not represent a successful remote CI run; validation reported
   here was performed locally.
6. The workflow supplies pkgdown as a build tool; it is not a package runtime
   dependency. The current `_pkgdown.yml` and `.github` remain unchanged.
7. Restored `AGENTS.md` contains outdated statements that TDT genotype error
   is unimplemented. These are development instructions, not user-facing
   method documentation; reconcile them with the validated current code in a
   later instructions-only update. Do not change code to match stale notes.

## 10. Validation record

Validation runs use local macOS arm64, R 4.6.1, and the RStudio-bundled Pandoc.
This is not an R-devel, Windows, Linux, or CRAN-farm result.

| Validation | Actual result |
|---|---|
| `devtools::document()` | Passed with roxygen2 8.1.0; NAMESPACE and man pages unchanged. |
| `devtools::test()` | 476 test blocks; 10,297 passing assertions; zero failures, errors, warnings or skips. testthat 3.3.2. |
| `R CMD build` | Passed, building all seven vignettes with knitr 1.51 and rmarkdown 2.31. |
| `R CMD INSTALL --library=<temporary library> paweh_0.99.0.tar.gz` | Passed, including load tests from temporary and final locations. |
| Installed `citation("paweh")` | Two entries: metadata-derived software citation and the separately scoped book entry. |
| Installed vignette registration | All seven source/HTML articles present. |
| `R CMD check --as-cran` through `devtools::check()` | Completed with 0 errors, 0 warnings and 2 notes. Remote incoming checks and forced Suggests were enabled; `NOT_CRAN=false`; manual generation was enabled. |
| Installed tests inside `--as-cran` | 10,301 passing assertions, zero failures, warnings or skips. This is a separate installed-package run from the load-all test run. |
| Vignette rebuild during `--as-cran` | Passed for all seven articles. |
| PDF manual | Built successfully. |
| HTML-manual validation | Not run successfully: the check skipped validation because the available HTML Tidy is too old. |
| Browser rendering review | All seven pages loaded with contents navigation, preserved tables, 18 loaded figures/screenshots and 38 rendered mathematical expressions; no broken images detected. Standard MathJax loading requires network access. |
| Source preservation | Every vignette body and all R implementations, tests, man pages and exports are byte-identical to the validated post-removal baseline. Only vignette format headers differ. |
| Source-package contents | CITATION is included; the audit, AGENTS, RStudio project, development license copy, pkgdown configuration and GitHub workflow are excluded. No BiocStyle dependency or `biocViews` remains. |

Exact check invocation used `devtools::check(document=FALSE, manual=TRUE,
cran=TRUE, remote=TRUE, incoming=TRUE, force_suggests=TRUE, args="--as-cran",
env_vars=c(NOT_CRAN="false"), error_on="never")`. The emitted command was
`R CMD check --as-cran`; no `--no-manual` flag was used.

The two notes were:

1. **CRAN incoming feasibility:** maintainer identity and "New submission".
   This is the expected first-submission notice, not a claim of submission or
   acceptance.
2. **HTML version of manual:** "Skipping checking HTML validation: 'tidy'
   doesn't look like recent enough HTML Tidy." The machine provides
   `/usr/bin/tidy`, released 31 October 2006 (Apple build 13670). Install or
   select a current HTML Tidy and repeat that validation before submission.
   The browser vignette review does not substitute for this skipped check.

An initial sandboxed attempt stopped with DNS failure while accessing CRAN's
package index; it did not complete the check. It also encountered the user's
globally configured Bioconductor repository entries. The successful rerun
used authorized network access and `options(repos=c(CRAN=
"https://cloud.r-project.org"))` in that check process only. No global R
repository settings were changed. The successful run reported no repository
access errors.

These results establish local compatibility, not universal CRAN readiness.
R 4.6.1 is the available local version; a current R-devel/current-release
cross-platform check has not been run. Logs, the structured test results,
reference-draft copies and rendered validation artifacts are outside the
checkout under `/private/tmp/paweh-cran-a1-validation`.

## 11. Outstanding release decisions

- Complete the educational plan in a separate phase, preserving published
  reproductions and the validated equations.
- Resolve the website reference/navigation issues before publishing the next
  site; this is independent of the CRAN source package.
- Run current R-devel and additional-platform checks before an actual
  submission; local checks do not replace them.
- Complete HTML-manual validation with a current HTML Tidy. Retain an
  explanation of the new-submission note in the eventual submission record.
- Verify minimum supported versions for runtime APIs and test helpers;
  testing against current installed packages does not prove the oldest
  declared versions work.
- Review author approval, licensing, screenshot redistribution, contact
  details, release version and manuscript/package citation strategy before
  submission. No new package-paper citation has been invented.
- Review any actual errors, warnings, notes or environmental limitations
  recorded in the final validation section before deciding readiness.

No CRAN submission, main-branch merge, Bioconductor interaction, or substantive
statistical/documentation redesign is performed in A1.
