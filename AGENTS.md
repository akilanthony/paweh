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
- cc_mssn_conditional_full()
- cc_power_conditional_full()

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

Generalized conditional framework (R/tdt_conditional.R):
- tdt_power_conditional_full()
- tdt_mssn_conditional_full()

These mirror the case-control generalized functions. Input modes are
"model_based" (prev, pd, R1/R2 or MOI, delta_prime, n_trios) and "model_free"
(ET, ENT; plus n_trios for MSSN, because ET and ENT are absolute counts that
embed a sample size).

Tests implemented:
- Standard transmission disequilibrium test (df = 1)

Model modifiers:
- Locus heterogeneity (stage 2 of .tdt_transmission_pipeline(), implemented):
  gT/gNT = pi * baseline + (1 - pi) * null_term, where pi is the probability a
  trio is linked to the disease locus (pi = 1 is full homogeneity and
  reproduces the stage-1 baseline exactly). Enabled via locus_het = TRUE and
  pi in tdt_power_conditional_full() and tdt_mssn_conditional_full(). Requires
  input_mode = "model_based", because the null term depends on p_plus, p_d, C,
  and delta, which are not recoverable from user-supplied ET/ENT alone.
- Stages 3-4 of .tdt_transmission_pipeline() are reserved for phenotype
  misclassification and genotype misclassification, and currently pass the
  stage-2 gT/gNT through unchanged while recording disabled modifier slots in
  the returned object.

Textbook mapping (Gordon, Finch, and Kim 2020):
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
- Locus heterogeneity for the TDT, gT*/gNT* mixture and NCP: Sect. 5.3.3,
  Eqs. 5.30-5.34b, pp. 293-294. The NCP under locus heterogeneity is due to
  Chen, Yang, Buyske, Matise, Finch, and Gordon (2009), Statistical
  Applications in Genetics and Molecular Biology, 8(44), building on Deng and
  Chen (2001), Genetical Research, 78(3), 289-302.

Reserved for future modifier stages:
- Phenotype misclassification for TDT: Sect. 5.2.6, p. 284.
- Genotype misclassification for TDT, including type I error inflation:
  Sect. 5.2.5, pp. 277 and 281.

Note on style: unlike case_control.R, which defines its helper functions inside
each exported function, R/tdt_conditional.R uses file-level internal helpers
prefixed with ".tdt_" and marked @noRd. This avoids duplicating the pipeline
between the power and MSSN entry points.

Do not invent citations. If the exact citation is missing, leave a TODO.
