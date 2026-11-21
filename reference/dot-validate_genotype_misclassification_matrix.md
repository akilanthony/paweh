# Validate a True-to-Observed Genotype Matrix

Validate a True-to-Observed Genotype Matrix

## Usage

``` r
.validate_genotype_misclassification_matrix(M, tolerance = 1e-10)
```

## Arguments

- M:

  Numeric three-by-three matrix with true genotypes in rows and observed
  genotypes in columns, both ordered 0, 1, 2.

- tolerance:

  Tolerance for row sums, not for detecting zero entries.

## Value

The validated matrix, retaining dimensions and dimnames.

## Details

An exactly all-zero row is replaced by its corresponding identity row,
assuming correct classification for that genotype. One warning is issued
per repaired matrix. Small nonzero entries are never treated as zero,
and other malformed rows are rejected without normalization.
