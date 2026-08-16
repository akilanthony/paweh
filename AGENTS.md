# AGENTS.md

This is an R package for statistical genetics power and sample size calculations.

General instructions:
- Preserve existing function naming conventions.
- Follow the existing TDT module style for roxygen2 documentation, examples, return objects, and testthat tests.
- Do not modify existing TDT functions unless explicitly asked.
- Do not change mathematical formulas unless there is a clear bug. If a possible bug is found, explain it before editing.
- Keep user-facing console output clean.
- Do not print internal S components in clean output.
- Keep internal S values in returned objects for validation and debugging.
- Use roxygen2 for documentation.
- Use testthat for tests.
- Run devtools::document(), devtools::test(), and devtools::check() after changes.
- Never include secrets, tokens, credentials, or personal access tokens in source files, documentation, examples, tests, commits, or comments.

# Case-control equation map

Please use this file to connect the case-control functions to the textbook.

Functions:
- cc_mssn()
- cc_power()

Tests implemented:
- Genotype chi-square test of independence
- Genotype trend test

Model modifiers:
- Locus heterogeneity: g_case_after = pi * g_case_before + (1 - pi) * g_ctrl
- Genotype misclassification:
  - none
  - 1p
  - 2p
  - 3p
  - diff3p

Textbook mapping:
- Genotype chi-square equation: TODO add textbook equation number/page
- Trend test equation: TODO add textbook equation number/page
- Locus heterogeneity equation or section: TODO add textbook equation number/page
- Genotype misclassification equation or section: TODO add textbook equation number/page

Do not invent citations. If the exact citation is missing, leave a TODO.

# Family-based (TDT) equation map

Canonical TDT framework (R/tdt_conditional.R):
- tdt_power() is the top-level power interface.
- tdt_mssn() is the top-level MSSN interface.

Their validated model-based formulas and three-scenario reporting (no error /
misclassification / heterogeneity) are preserved. Both also support
input_mode = c("model_based", "model_free"):

- "model_based" (default) uses genetic-model inputs.
- "model_free" lets a user who already has expected transmission and
  non-transmission counts (ET, ENT) supply them directly instead of
  prev/R1/R2. gT = ET / (2 * n_trios), gNT = ENT / (2 * n_trios), where
  n_trios is N itself for the power function (ET/ENT are for the same N
  trios power is computed for) and a separate n_trios argument for the MSSN
  function (ET/ENT are absolute counts that embed a sample size, which need
  not equal the N* being solved for).

heter_rate and misclass_rate work in both modes. In model_free mode they use
closed-form identities in terms of gT, gNT, pd, and prev alone (verified
algebraically equivalent to the calc_gTgNT_heter()/calc_gTgNT_misclass()
identities used by the model-based implementation). Let A = gT - gNT
(no-error) and
p_plus = 1 - pd:
  heterogeneity: gT = pd*p_plus + A*(p_plus - 0.5*heter_rate),
                 gNT = pd*p_plus + A*(-pd + 0.5*heter_rate)
  misclassification: m = prev*(1 - misclass_rate) / (prev + misclass_rate*(1 - prev)),
                      gT = pd*p_plus + p_plus*A*m, gNT = pd*p_plus - pd*A*m
A scenario is only computed from these identities when its rate is non-zero;
at a rate of exactly 0 the scenario reuses the no-error gT/gNT directly, so
pd is not required unless at least one rate is non-zero, and prev is not
required unless misclass_rate is non-zero. If pd is needed but not supplied,
it is solved from 2*pd^2 - 2*pd*(1-A) + (gT + gNT - A) = 0, taking the root
in (0, 0.5) with a message; supplying pd directly is preferred.

R1/R2 are never used in model_free mode.

Tests implemented:
- Standard transmission disequilibrium test (df = 1)

Model modifiers:
- Phenotype misclassification and locus heterogeneity: implemented, exactly
  as in tdt_power()/tdt_mssn() for model_based, and via
  the closed-form identities above for model_free.
- Genotype misclassification: not yet implemented in this framework.

Textbook mapping (Gordon, Finch, and Nothnagel 2020):
- Penetrances from prevalence and genotype relative risks: Eq. 1.6 and Eq. 1.7,
  pp. 13.
- TDT non-centrality parameter, lambda = (ET - ENT)^2 / (ET + ENT), with
  ET = 2N[pd*p_plus + delta*p_plus*C/phi1],
  ENT = 2N[pd*p_plus - delta*pd*C/phi1],
  C = pd*f2 + (1 - 2*pd)*f1 - p_plus*f0, delta = delta_prime*pd*p_plus:
  Eq. 1.25, p. 27.
- Null and alternative distributions (central and non-central chi-square,
  1 df): Table 1.5, p. 25.
- Two-locus haplotype frequencies and the disequilibrium parameter D: Table 1.2
  and Sect. 1.4.2, pp. 13-14.
- TDT statistic definition and transmitted/non-transmitted counts: Sect. 1.6.1.3
  and Table 1.4, p. 24.
- Phenotype misclassification for the TDT, gT*/gNT* under pi01: Eqs. 5.26-5.28,
  Sect. 5.2.6 (as implemented in calc_gTgNT_misclass() / tdt_expected_transmission_probability()).
- Locus heterogeneity for the TDT, gT*/gNT* mixture and NCP: Eq. 5.34,
  Sect. 5.3.3, Eqs. 5.30-5.34b, pp. 293-294. The NCP under locus heterogeneity
  is due to Chen, Yang, Buyske, Matise, Finch, and Gordon (2009), Statistical
  Applications in Genetics and Molecular Biology, 8(44), building on Deng and
  Chen (2001), Genetical Research, 78(3), 289-302.

Reserved for future modifier stages:
- Genotype misclassification for TDT, including type I error inflation:
  Sect. 5.2.5, pp. 277 and 281.

Note on style: tdt_power() and tdt_mssn() each define local copies of
calc_gTgNT_misclass(), calc_gTgNT_heter(), and the model_free pd-solving
helper rather than sharing file-level helpers.

Do not invent citations. If the exact citation is missing, leave a TODO.

## Next steps / open questions

- Genotype misclassification for the TDT (Sect. 5.2.5, pp. 277 and 281) is
  not implemented anywhere in the package yet, for either input_mode. Unlike
  phenotype misclassification and locus heterogeneity, there is no existing
  calc_gTgNT_*()-style formula to copy from R/tdt_functions.R -- this is new
  work, not a port, and needs its own derivation and citation before coding.
- The specialized 2D and transitional 3D TDT plotting functions call the
  canonical tdt_power()/tdt_mssn() backends. A future phase may consolidate
  the four 3D functions into a generalized surface API.
