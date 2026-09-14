# Internal helpers for the Chapter 5 single-parameter adjacent-genotype TDT
# error model in Gordon, Finch, and Kim (2020), Section 5.2.5. These helpers are
# deliberately separate from both the case-control genotype-error models and
# the TDT-NGS raw-read likelihood.

.tdt_genotype_misclassification_matrix <- function(e) {
  if (!is.numeric(e) || length(e) != 1L || !is.finite(e) || e < 0 || e > 0.5) {
    stop("genotype_misclassification_rate must be a single finite number in [0, 0.5].")
  }

  M <- matrix(c(
    1 - e, e, 0,
    e, 1 - 2 * e, e,
    0, e, 1 - e
  ), nrow = 3L, byrow = TRUE,
  dimnames = list(true = 0:2, observed = 0:2))

  .validate_genotype_misclassification_matrix(M)
}

.tdt_all_trio_states <- function() {
  states <- expand.grid(
    father = 0:2,
    mother = 0:2,
    child = 0:2,
    KEEP.OUT.ATTRS = FALSE
  )
  states$label <- paste0(states$father, states$mother, states$child)
  states
}

.tdt_mendelian_consistent <- function(father, mother, child) {
  if (length(father) != length(mother) || length(father) != length(child) ||
      any(!is.finite(c(father, mother, child))) ||
      any(!(c(father, mother, child) %in% 0:2))) {
    stop("father, mother, and child must be equal-length genotype vectors in {0, 1, 2}.")
  }

  vapply(seq_along(father), function(i) {
    father_gametes <- if (father[[i]] == 0) 0 else if (father[[i]] == 2) 1 else 0:1
    mother_gametes <- if (mother[[i]] == 0) 0 else if (mother[[i]] == 2) 1 else 0:1
    child[[i]] %in% unique(as.vector(outer(father_gametes, mother_gametes, `+`)))
  }, logical(1))
}

.tdt_genotype_transition_kernel <- function(e) {
  true_states <- .tdt_ngs_trio_states()[, c("father", "mother", "child")]
  observed_states <- .tdt_all_trio_states()
  M <- .tdt_genotype_misclassification_matrix(e)

  kernel <- matrix(
    0,
    nrow = nrow(true_states),
    ncol = nrow(observed_states),
    dimnames = list(
      paste0("T", apply(true_states, 1L, paste0, collapse = "")),
      paste0("O", observed_states$label)
    )
  )
  for (i in seq_len(nrow(true_states))) {
    for (j in seq_len(nrow(observed_states))) {
      kernel[i, j] <-
        M[true_states$father[[i]] + 1L, observed_states$father[[j]] + 1L] *
        M[true_states$mother[[i]] + 1L, observed_states$mother[[j]] + 1L] *
        M[true_states$child[[i]] + 1L, observed_states$child[[j]] + 1L]
    }
  }

  if (any(!is.finite(kernel)) || any(kernel < 0) ||
      any(abs(rowSums(kernel) - 1) > 1e-12)) {
    stop("Internal TDT genotype-misclassification transition probabilities are invalid.")
  }
  kernel
}

# Construct the affected-child marker-trio distribution implied by the same
# two-locus LD and penetrance model used by the ordinary TDT formulas. The
# disease and marker alleles have the package's common frequency pd, and
# D = delta_prime * pd * (1 - pd). Enumerating ordered parental haplotypes and
# transmissions reproduces the ordinary TDT ET/ENT formulas at e = 0.
.tdt_genotype_true_trio_probabilities <- function(pd, prev, R1, R2,
                                                  delta_prime) {
  scalars <- c(pd = pd, prev = prev, R1 = R1, R2 = R2,
               delta_prime = delta_prime)
  if (any(!is.finite(scalars)) || pd <= 0 || pd >= 1 || prev <= 0 || prev >= 1 ||
      R1 <= 0 || R2 <= 0) {
    stop("The genotype-misclassification branch requires finite model parameters with pd and prev in (0, 1) and R1 and R2 > 0.")
  }

  qd <- 1 - pd
  D <- delta_prime * pd * qd
  haplotypes <- data.frame(
    disease = c(1L, 1L, 0L, 0L),
    marker = c(1L, 0L, 1L, 0L),
    frequency = c(pd^2 + D, pd * qd - D, qd * pd - D, qd^2 + D)
  )
  tolerance <- 1e-14
  if (any(haplotypes$frequency < -tolerance) ||
      abs(sum(haplotypes$frequency) - 1) > 1e-12) {
    stop("delta_prime implies invalid two-locus haplotype probabilities for the genotype-misclassification branch.")
  }
  haplotypes$frequency[haplotypes$frequency < 0] <- 0

  Z <- qd^2 + 2 * pd * qd * R1 + pd^2 * R2
  penetrance <- c(prev / Z, R1 * prev / Z, R2 * prev / Z)
  if (any(penetrance < 0)) {
    stop("The genotype-misclassification branch produced invalid penetrances.")
  }
  if (any(penetrance > 1)) {
    warning("Computed penetrances f1 or f2 exceed 1; check prev/R1/R2 inputs.")
  }

  states <- .tdt_ngs_trio_states()
  state_labels <- paste0(states$father, states$mother, states$child)
  probability <- numeric(nrow(states))
  names(probability) <- state_labels

  for (father_h1 in seq_len(nrow(haplotypes))) {
    for (father_h2 in seq_len(nrow(haplotypes))) {
      for (mother_h1 in seq_len(nrow(haplotypes))) {
        for (mother_h2 in seq_len(nrow(haplotypes))) {
          parental_probability <- prod(haplotypes$frequency[
            c(father_h1, father_h2, mother_h1, mother_h2)
          ])
          if (parental_probability == 0) next

          father_marker <- haplotypes$marker[[father_h1]] +
            haplotypes$marker[[father_h2]]
          mother_marker <- haplotypes$marker[[mother_h1]] +
            haplotypes$marker[[mother_h2]]

          for (father_transmission in c(father_h1, father_h2)) {
            for (mother_transmission in c(mother_h1, mother_h2)) {
              child_marker <- haplotypes$marker[[father_transmission]] +
                haplotypes$marker[[mother_transmission]]
              child_disease <- haplotypes$disease[[father_transmission]] +
                haplotypes$disease[[mother_transmission]]
              label <- paste0(father_marker, mother_marker, child_marker)
              probability[[label]] <- probability[[label]] +
                parental_probability * 0.25 * penetrance[[child_disease + 1L]]
            }
          }
        }
      }
    }
  }

  affected_probability <- sum(probability)
  if (!is.finite(affected_probability) || affected_probability <= 0) {
    stop("The genotype-misclassification true-trio distribution has zero or invalid mass.")
  }
  probability <- probability / affected_probability
  probability <- unname(probability[state_labels])
  names(probability) <- paste0("T", state_labels)

  if (any(!is.finite(probability)) || any(probability < 0) ||
      abs(sum(probability) - 1) > 1e-12) {
    stop("Internal true-trio probabilities are invalid.")
  }
  probability
}

.tdt_genotype_error_expected_counts <- function(true_probability, e, N) {
  if (!is.numeric(true_probability) || length(true_probability) != 15L ||
      any(!is.finite(true_probability)) || any(true_probability < 0) ||
      abs(sum(true_probability) - 1) > 1e-10) {
    stop("true_probability must contain 15 finite nonnegative probabilities summing to 1.")
  }
  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0) {
    stop("N must be a single finite number greater than zero.")
  }

  observed_states <- .tdt_all_trio_states()
  observed_probability <- as.numeric(
    true_probability %*% .tdt_genotype_transition_kernel(e)
  )
  names(observed_probability) <- paste0("O", observed_states$label)
  if (any(!is.finite(observed_probability)) || any(observed_probability < 0) ||
      abs(sum(observed_probability) - 1) > 1e-10) {
    stop("Internal observed-trio probabilities are invalid.")
  }

  is_mc <- .tdt_mendelian_consistent(
    observed_states$father, observed_states$mother, observed_states$child
  )
  retained_probability <- sum(observed_probability[is_mc])
  if (!is.finite(retained_probability) || retained_probability <= 0) {
    stop("No Mendelian-consistent observed trios remain after genotype misclassification.")
  }
  retained <- observed_probability[is_mc] / retained_probability
  names(retained) <- observed_states$label[is_mc]

  ET <- N * sum(retained[c("011", "101", "111", "112", "112", "122", "212")])
  ENT <- N * sum(retained[c("010", "100", "111", "110", "110", "121", "211")])
  denominator <- ET + ENT
  lambda <- if (denominator == 0) 0 else (ET - ENT)^2 / denominator

  list(
    ET = unname(ET),
    ENT = unname(ENT),
    lambda = unname(lambda),
    retained_probability = unname(retained_probability),
    mendelian_inconsistent_probability = unname(1 - retained_probability),
    observed_probability = observed_probability,
    retained_probability_by_state = retained
  )
}

.tdt_genotype_error_model_counts <- function(pd, prev, R1, R2,
                                              delta_prime, e, N) {
  true_probability <- .tdt_genotype_true_trio_probabilities(
    pd = pd, prev = prev, R1 = R1, R2 = R2,
    delta_prime = delta_prime
  )
  out <- .tdt_genotype_error_expected_counts(true_probability, e, N)
  out$true_probability <- true_probability
  out
}

.tdt_genotype_error_null_diagnostic <- function(pd, e, N = 5000,
                                                alpha = 5e-8) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single finite number in (0, 1).")
  }
  mu <- .tdt_ngs_hwe_mating_freqs(pd)
  true_probability <- .tdt_ngs_state_probabilities(mu, t = 0.5)
  counts <- .tdt_genotype_error_expected_counts(true_probability, e, N)
  critical <- stats::qchisq(1 - alpha, df = 1)
  actual_alpha <- stats::pchisq(
    critical, df = 1, ncp = counts$lambda, lower.tail = FALSE
  )

  list(
    ET0 = counts$ET,
    ENT0 = counts$ENT,
    lambda_null = counts$lambda,
    actual_alpha = actual_alpha,
    alpha_inflation_ratio = actual_alpha / alpha,
    log10_alpha_inflation_ratio = log10(actual_alpha / alpha),
    retained_probability = counts$retained_probability,
    mendelian_inconsistent_probability =
      counts$mendelian_inconsistent_probability
  )
}

.tdt_validate_effect <- function(effect, effect_missing, misclass_rate,
                                 heter_rate,
                                 genotype_misclassification_rate,
                                 input_mode) {
  if (effect_missing) {
    if (!is.null(genotype_misclassification_rate)) {
      stop("Set effect = 'genotype_misclassification' when supplying genotype_misclassification_rate.")
    }
    return("legacy_scenarios")
  }

  if (effect == "none" && (misclass_rate != 0 || heter_rate != 0 ||
                            !is.null(genotype_misclassification_rate))) {
    stop("effect = 'none' cannot be combined with a nonzero error or heterogeneity parameter.")
  }
  if (effect == "phenotype_misclassification" && heter_rate != 0) {
    stop("Phenotype misclassification and locus heterogeneity cannot be combined.")
  }
  if (effect == "locus_heterogeneity" && misclass_rate != 0) {
    stop("Locus heterogeneity and phenotype misclassification cannot be combined.")
  }
  if (effect == "genotype_misclassification") {
    if (misclass_rate != 0 || heter_rate != 0) {
      stop("Genotype misclassification cannot be combined with phenotype misclassification or locus heterogeneity.")
    }
    if (is.null(genotype_misclassification_rate)) {
      stop("genotype_misclassification_rate is required when effect = 'genotype_misclassification'.")
    }
    .tdt_genotype_misclassification_matrix(genotype_misclassification_rate)
    if (input_mode != "model_based") {
      stop("The genotype-misclassification effect requires input_mode = 'model_based'; ET and ENT alone do not identify the 15-state trio distribution.")
    }
  } else if (!is.null(genotype_misclassification_rate)) {
    stop("genotype_misclassification_rate is used only when effect = 'genotype_misclassification'.")
  }
  effect
}
