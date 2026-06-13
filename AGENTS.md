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
- Optional allelic chi-square test
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
- Allelic chi-square equation: TODO add textbook equation number/page
- Trend test equation: TODO add textbook equation number/page
- Locus heterogeneity equation or section: TODO add textbook equation number/page
- Genotype misclassification equation or section: TODO add textbook equation number/page

Do not invent citations. If the exact citation is missing, leave a TODO.
