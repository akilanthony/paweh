# Family-Based (TDT) Minimum Sample Size Necessary

Computes the minimum number of affected trios required to achieve a
specified power for the transmission disequilibrium test (TDT) under
three scenarios: (i) no error, (ii) phenotype misclassification only,
and (iii) locus heterogeneity only. `input_mode` lets the transmission
probabilities come either from a genetic model (`"model_based"`, the
default) or directly from user-supplied expected transmission and
non-transmission counts (`"model_free"`).

## Usage

``` r
tdt_mssn(
  target_power,
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
  n_trios = NULL,
  verbose = TRUE
)
```

## Arguments

- target_power:

  Numeric in (0,1). Desired power for the TDT.

- input_mode:

  Character. One of `"model_based"` (default) or `"model_free"`. See
  Details.

- pd:

  Numeric in (0,1). Frequency of the disease (high-risk) allele at the
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
  (\\\pi\_{01}\\) in the misclassification scenario. Zero means no
  misclassification and is the default.

- heter_rate:

  Numeric in \\\[0,1)\\. Proportion of trios whose affection status is
  not due to the locus of interest (\\1 - \pi\\) in the heterogeneity
  scenario. Zero means complete homogeneity and is the default.

- ET, ENT:

  Numeric \\\ge 0\\. Expected transmission and non-transmission counts
  for `input_mode = "model_free"`.

- n_trios:

  Numeric \\\> 0\\. Number of affected trios that the supplied `ET` and
  `ENT` correspond to. Required for `input_mode = "model_free"` (unlike
  the power function, there is no `N` argument here to reuse, since this
  function solves for the required number of trios).

- verbose:

  Logical. If `TRUE` (default), prints a formatted summary of the
  no-error required number of trios and any requested modifier-specific
  required counts and percentage increases.

## Value

An object of class `"tdt_mssn"`, containing:

- alpha, target_power, lambda_star, input_mode:

  Design targets, target one-df NCP, and input mode.

- N:

  Named `no_error`, `misclassification`, and `heterogeneity` required
  affected-trio counts.

- percent_increase:

  Modifier-specific conventional percentage inflation,
  \\100(N\_{adjusted}/N\_{no-error}-1)\\, relative to the no-error
  required count. Values are `NA` when the no-error required count is
  not finite and positive.

- power_at_N_no_error, power_loss_at_N_no_error:

  Power and power loss obtained if the no-error design size is used
  under each modifier.

- gT_star, gNT_star:

  Scenario-specific transmission and non-transmission probabilities used
  in the MSSN formulas.

- model_parameters:

  Supplied and derived genetic-model, LD, misclassification,
  heterogeneity, and model-free count-scale information.

## Details

With `input_mode = "model_based"`, penetrances are derived from `prev`,
`R1`, `R2`, and `pd`, and \$\$N^\* = \frac{\lambda^\* (g_T^\* +
g\_{NT}^\*)}{2 (g_T^\* - g\_{NT}^\*)^2}\$\$ is computed for each of the
three scenarios. The no-error calculation follows Eq. 1.25 (p. 27);
phenotype misclassification follows Eqs. 5.24–5.27 (pp. 284–285); and
locus heterogeneity follows Eqs. 5.30–5.34b (pp. 293–294). Eq. 5.28 is a
numerical example. Here `heter_rate` is the heterogeneous fraction,
equal to \\1-\pi\\ when lower-level formulas use homogeneous fraction
\\\pi\\. Phenotype misclassification and locus heterogeneity are
reported as separate sensitivity scenarios; supplying both rates does
not create a combined-error scenario.

With `input_mode = "model_free"`, the no-error scenario uses \\g_T = ET
/ (2\\n\_{trios})\\ and \\g\_{NT} = ENT / (2\\n\_{trios})\\. The
misclassification and heterogeneity scenarios then use the same
closed-form identities as
[`tdt_power`](https://akilanthony.github.io/pawh/reference/tdt_power.md)
(see its Details for the formulas and the `pd`-solving fallback),
applied to this no-error \\g_T\\/\\g\_{NT}\\ pair.

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

[`tdt_power`](https://akilanthony.github.io/pawh/reference/tdt_power.md)
for the power counterpart.

## Examples

``` r
# Model-based input
tdt_mssn(
  target_power = 0.80, input_mode = "model_based",
  pd = 0.30, prev = 0.05, R1 = 1.5, R2 = 2.25,
  verbose = FALSE
)$N$no_error
#> [1] 214.9086

# model_free: supply expected transmissions/non-transmissions directly
tdt_mssn(
  target_power = 0.80, input_mode = "model_free",
  ET = 140, ENT = 100, n_trios = 120,
  pd = 0.30, prev = 0.05,
  verbose = FALSE
)$N$no_error
#> [1] 141.279
```
