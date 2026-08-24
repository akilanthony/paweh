# Family-Based (TDT) Power

Computes power for the transmission disequilibrium test (TDT) at a fixed
number of affected trios under three scenarios: (i) no error, (ii)
phenotype misclassification only, and (iii) locus heterogeneity only.
`input_mode` lets the transmission probabilities come either from a
genetic model (`"model_based"`, the default) or directly from
user-supplied expected transmission and non-transmission counts
(`"model_free"`).

## Usage

``` r
tdt_power(
  N,
  input_mode = c("model_based", "model_free"),
  pd = NULL,
  prev = NULL,
  R1 = NULL,
  R2 = NULL,
  alpha = 0.05,
  delta_prime = 1,
  misclass_rate = 0,
  heter_rate = 0,
  ET = NULL,
  ENT = NULL,
  verbose = TRUE
)
```

## Arguments

- N:

  Numeric \\\> 0\\. Number of affected trios.

- input_mode:

  Character. One of `"model_based"` (default) or `"model_free"`. See
  Details.

- pd:

  Numeric in (0,1). Frequency of the disease/high-risk allele at the
  marker locus. Required for `input_mode = "model_based"`. Optional for
  `"model_free"`: required only if `heter_rate` or `misclass_rate` is
  non-zero, and if omitted in that case it is solved from `ET`/`ENT`
  (see Details).

- prev:

  Numeric in (0,1). Disease prevalence (\\\phi_1\\). Required for
  `input_mode = "model_based"`. Required for `"model_free"` only if
  `misclass_rate` is non-zero.

- R1:

  Numeric \\\> 0\\. Genotype relative risk for heterozygotes. Required
  for `input_mode = "model_based"`; unused for `"model_free"`.

- R2:

  Numeric \\\> 0\\. Genotype relative risk for homozygotes. Required for
  `input_mode = "model_based"`; unused for `"model_free"`.

- alpha:

  Numeric in (0,1). Significance level for the TDT (default `0.05`).

- delta_prime:

  Numeric. Linkage disequilibrium scale parameter \\D'\\ (default `1`).
  Used for `input_mode = "model_based"` only.

- misclass_rate:

  Numeric in \\\[0,1)\\. Phenotype misclassification rate for controls
  (\\\pi\_{01}\\). A value of `0` corresponds to no misclassification
  and is the default.

- heter_rate:

  Numeric in \\\[0,1)\\. Proportion of trios whose affection status is
  *not* due to the locus of interest (\\1 - \pi\\). A value of `0`
  corresponds to complete homogeneity and is the default.

- ET, ENT:

  Numeric \\\ge 0\\. Expected transmission and non-transmission counts
  for `input_mode = "model_free"`, assumed to be accumulated over the
  same `N` trios that power is computed for.

- verbose:

  Logical. If `TRUE` (default), prints a formatted summary of no-error
  power and any requested modifier-specific power.

## Value

An object of class `"tdt_power"`, containing:

- alpha, N, input_mode:

  Significance level, affected-trio count, and input mode.

- lambda, power:

  Named `no_error`, `misclassification`, and `heterogeneity` NCP and
  power values.

- power_loss:

  Absolute power loss for each modifier relative to the no-error
  scenario.

- gT_star, gNT_star:

  Scenario-specific transmission and non-transmission probabilities.

- ET, ENT:

  Scenario-specific expected transmission and non-transmission counts
  for the supplied `N` trios.

- model_parameters:

  Supplied and derived genetic-model, LD, misclassification, and
  heterogeneity parameters.

## Details

With `input_mode = "model_based"`, penetrances \\f_0, f_1, f_2\\ are
derived from `prev`, `R1`, `R2`, and `pd`, and the expected transmission
and non-transmission probabilities under each of the three scenarios are
computed from the standard TDT formulation in Eq. 1.25 (Chapter 1,
Section 1.6.1.3, p. 27), the phenotype-misclassification probabilities
in Eqs. 5.24–5.27 (Section 5.2.6, pp. 284–285), and the
locus-heterogeneity construction in Eqs. 5.30–5.34a (Section 5.3.3, pp.
293–294) of Gordon, Finch, and Kim (2020). Eq. 5.28 is a numerical
worked example, not the general symbolic formula. Phenotype
misclassification and locus heterogeneity are reported as separate
sensitivity scenarios; supplying both rates does not create a
combined-error scenario.

`N` is the number of affected-child trios, with both parents genotyped.
`ET` and `ENT` are expected counts accumulated over trios, whereas
\\g_T\\ and \\g\_{NT}\\ are per-parental-allele transmission and
non-transmission probabilities. `heter_rate` is the heterogeneous
fraction; in lower-level formulas written with the linked/homogeneous
fraction \\\pi\\; equivalently, the heterogeneous fraction is \\1-\pi\\.

With `input_mode = "model_free"`, the no-error scenario uses \\g_T = ET
/ (2N)\\ and \\g\_{NT} = ENT / (2N)\\ directly. The misclassification
and heterogeneity scenarios are then computed from the following
identities, which are algebraically equivalent to the the corresponding
model-based identities. Let \\A = g_T - g\_{NT}\\ (no-error) and \\p\_+
= 1 - p_d\\: \$\$\text{heterogeneity: } g_T = p_d p\_+ + A (p\_+ - 0.5
\times \text{heter\\rate}), \quad g\_{NT} = p_d p\_+ + A (-p_d + 0.5
\times \text{heter\\rate}),\$\$ \$\$\text{misclassification: } m =
\frac{\phi_1 (1 - \pi\_{01})} {\phi_1 + \pi\_{01} (1 - \phi_1)}, \quad
g_T = p_d p\_+ + p\_+ A m, \quad g\_{NT} = p_d p\_+ - p_d A m.\$\$ A
scenario is only computed from these identities when its rate is
non-zero; at a rate of exactly `0` the scenario reuses the no-error
\\g_T\\/\\g\_{NT}\\ directly, so `pd` is not required unless at least
one rate is non-zero, and `prev` is not required unless `misclass_rate`
is non-zero.

If `pd` is needed but not supplied, it is solved from `ET` and `ENT` via
\\2 p_d^2 - 2 p_d (1 - A) + (g_T + g\_{NT} - A) = 0\\. This quadratic
has two roots summing to \\1 - A\\; the root in \\(0, 0.5)\\ is used,
with a message reporting the derived value. An error is raised if no
such unique root exists – supplying `pd` directly is preferred.

## References

Spielman, R. S., McGinnis, R. E., & Ewens, W. J. (1993). Transmission
test for linkage disequilibrium: the insulin gene region and
insulin-dependent diabetes mellitus. *American Journal of Human
Genetics*, 52(3), 506–516. PMID: 8447318; PMCID: PMC1682161.

Gordon, D., Finch, S. J., & Kim, W. (2020). *Heterogeneity in
Statistical Genetics: How to Assess, Address, and Account for Mixtures
in Association Studies*. Springer.
[doi:10.1007/978-3-030-61121-7](https://doi.org/10.1007/978-3-030-61121-7)
.

Buyske, S., Yang, G., Matise, T. C., & Gordon, D. (2009). When a case is
not a case: Effects of phenotype misclassification on power and sample
size requirements for the transmission disequilibrium test with affected
child trios. *Human Heredity*, 67(4), 287–292.
[doi:10.1159/000194981](https://doi.org/10.1159/000194981) .

Chen, C., Yang, G., Buyske, S., Matise, T., Finch, S. J., & Gordon, D.
(2009). Transmission disequilibrium test power and sample size in the
presence of locus heterogeneity. *Statistical Applications in Genetics
and Molecular Biology*, 8, Article 44.
[doi:10.2202/1544-6115.1501](https://doi.org/10.2202/1544-6115.1501) .

## See also

[`tdt_mssn`](https://akilanthony.github.io/paweh/reference/tdt_mssn.md)
for the MSSN counterpart.

## Examples

``` r
# Model-based input
tdt_power(
  N = 600, input_mode = "model_based",
  pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
  misclass_rate = 0.01, heter_rate = 0.10,
  verbose = FALSE
)$power$no_error
#> [1] 0.9967475

# model_free: supply expected transmissions/non-transmissions directly
tdt_power(
  N = 600, input_mode = "model_free",
  ET = 140, ENT = 100,
  pd = 0.30, prev = 0.05,
  misclass_rate = 0.01, heter_rate = 0.10,
  verbose = FALSE
)$power$no_error
#> [1] 0.73304
```
