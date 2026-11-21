# Sequencing-Derived Genotype Error Matrix

Constructs the analytic genotype transition matrix induced by fixed
sequencing coverage, symmetric per-read sequencing error, and
deterministic maximum-likelihood genotype calling.

## Usage

``` r
ngs_genotype_error_matrix(coverage, seq_error)
```

## Arguments

- coverage:

  A single finite integer greater than or equal to 1 giving sequencing
  coverage (read depth).

- seq_error:

  A single finite numeric sequencing-error probability in `[0, 0.5)`.

## Value

A numeric `3 x 3` matrix. Rows `true_0`, `true_1`, and `true_2` identify
the true genotype; columns `called_0`, `called_1`, and `called_2`
identify the maximum-likelihood call.

## Details

For true genotype \\G \in \\0,1,2\\\\, the alternate-read count follows
\$\$X \mid G,V,\epsilon \sim \mathrm{Binomial}(V,q_G),\$\$ where
\\q_0=\epsilon\\, \\q_1=0.5\\, and \\q_2=1-\epsilon\\. For every
possible alternate-read count, the called genotype maximizes its
binomial likelihood among the three genotype models.

Numerical likelihood ties within `1e-12` are resolved symmetrically:
below an alternate-read fraction of one half, the lower tied genotype is
selected; above one half, the higher tied genotype is selected; at
exactly one half, the tied genotype closest to 1 is selected.

Rows correspond to true genotypes and columns to called genotypes. Thus,
if `E` is the returned matrix and `g_true` is a length-three vector of
true genotype probabilities, called-genotype probabilities are obtained
as `as.numeric(t(E) %*% g_true)`. Every row sums to one within numerical
precision.

At finite coverage, even with zero sequencing error, a true heterozygote
can produce all-reference or all-alternate reads. Consequently the
finite-depth matrix need not be exactly the identity, especially at very
low coverage.

## Examples

``` r
E <- ngs_genotype_error_matrix(coverage = 10, seq_error = 0.01)
g_true <- c(0.70, 0.25, 0.05)
as.numeric(t(E) %*% g_true)
#> [1] 0.69969921 0.24782856 0.05247224
```
