#!/usr/bin/env Rscript
# One-off cross-check of the numpy TMM / log2-CPM in `sead-sexdiff export`
# against edgeR. Reads the files written by `run_export(debug_dir=...)`:
#   <debug_dir>/<cohort>_<omic>_counts.csv   donors x filtered features
#   <debug_dir>/<cohort>_<omic>_tmm.csv      donor, lib_size, factor
#   <debug_dir>/<cohort>_<omic>_logcpm.csv   donors x filtered features
#
# Usage (inside the Nix shell):
#   Rscript scripts/check_tmm_edger.R <debug_dir> [cohort] [omic ...]
# Exits non-zero if TMM factors differ by more than 1e-8 or log-CPM by more
# than 1e-10 (absolute).

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: check_tmm_edger.R <debug_dir> [cohort] [omic ...]")
debug_dir <- args[1]
cohort <- if (length(args) >= 2) args[2] else "cohort_30"
omics <- if (length(args) >= 3) args[3:length(args)] else c("rna", "atac")

tol_factor <- 1e-8
tol_logcpm <- 1e-10

# scan() instead of read.csv(): the ATAC matrices have ~130k columns
read_matrix <- function(path) {
  con <- file(path, "r")
  header <- strsplit(readLines(con, n = 1), ",", fixed = TRUE)[[1]]
  vals <- scan(con, what = character(), sep = ",", quiet = TRUE)
  close(con)
  m <- matrix(vals, ncol = length(header), byrow = TRUE)
  out <- matrix(as.numeric(m[, -1]), nrow = nrow(m))
  dimnames(out) <- list(m[, 1], header[-1])
  out
}

ok <- TRUE
for (omic in omics) {
  stem <- file.path(debug_dir, paste0(cohort, "_", omic))
  counts <- t(read_matrix(paste0(stem, "_counts.csv")))  # features x donors
  py_tmm <- read.csv(paste0(stem, "_tmm.csv"), row.names = 1)
  py_lcpm <- t(read_matrix(paste0(stem, "_logcpm.csv")))

  y <- DGEList(counts = counts)
  y <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 2)

  d_lib <- max(abs(y$samples$lib.size - py_tmm[colnames(counts), "lib_size"]))
  d_factor <- max(abs(y$samples$norm.factors - py_tmm[colnames(counts), "factor"]))
  d_lcpm <- max(abs(lcpm - py_lcpm[rownames(lcpm), colnames(lcpm)]))
  pass <- d_lib == 0 && d_factor <= tol_factor && d_lcpm <= tol_logcpm
  ok <- ok && pass
  cat(sprintf(
    "%s %s: %d features x %d donors; max |d lib.size| = %g, max |d factor| = %.3g, max |d log-CPM| = %.3g -> %s\n",
    cohort, omic, nrow(counts), ncol(counts), d_lib, d_factor, d_lcpm,
    if (pass) "PASS" else "FAIL"
  ))
}
cat(sprintf("edgeR %s\n", as.character(packageVersion("edgeR"))))
if (!ok) quit(status = 1)
