#!/usr/bin/env Rscript
# Cross-check of the numpy cyclic loess (trend normalization) of the
# preprocessing-sensitivity export against limma. Reads the files written by
# `_export_cohort(debug_dir=...)` for a trend-normalized omic:
#   <debug_dir>/<cohort>_<omic>_logcpm.csv       donors x filtered features (step-5 log2-CPM)
#   <debug_dir>/<cohort>_<omic>_normalized.csv   donors x filtered features (numpy cyclic loess)
#   <debug_dir>/<cohort>_<omic>_loess.csv        rule, span, iterations, robustness_iterations, n_features
# and runs limma::normalizeCyclicLoess(t(logcpm), span = <recorded span>,
# adaptive.span = FALSE, iterations = <recorded>, method = "fast") on the
# same log2-CPM.
#
# Usage (inside the Nix shell):
#   Rscript scripts/check_cyclic_loess_limma.R <debug_dir> [cohort] [omic ...]
# Exits non-zero if the normalized values differ by more than 1e-6 (absolute,
# log2-CPM), or if the recorded span is not limma's adaptive span for the
# number of features.

suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: check_cyclic_loess_limma.R <debug_dir> [cohort] [omic ...]")
debug_dir <- args[1]
cohort <- if (length(args) >= 2) args[2] else "cohort_30"
omics <- if (length(args) >= 3) args[3:length(args)] else c("rna", "atac")

tol <- 1e-6

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
  settings <- read.csv(paste0(stem, "_loess.csv"))
  if (settings$rule != "cyclic_loess") stop(sprintf("%s: rule %s, expected cyclic_loess", stem, settings$rule))
  lcpm <- t(read_matrix(paste0(stem, "_logcpm.csv")))      # features x donors
  py_norm <- t(read_matrix(paste0(stem, "_normalized.csv")))
  if (!identical(dimnames(lcpm), dimnames(py_norm))) stop(sprintf("%s: logcpm and normalized differ in features or donors", stem))

  adaptive <- chooseLowessSpan(nrow(lcpm), small.n = 50, min.span = 0.3, power = 1/3)
  span_ok <- settings$n_features == nrow(lcpm) && abs(settings$span - adaptive) <= 1e-12
  t0 <- proc.time()[["elapsed"]]
  r_norm <- normalizeCyclicLoess(
    lcpm, span = settings$span, adaptive.span = FALSE, iterations = settings$iterations, method = "fast"
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  d <- abs(r_norm - py_norm)
  d_max <- max(d)
  where <- arrayInd(which.max(d), dim(d))
  pass <- span_ok && d_max <= tol
  ok <- ok && pass
  cat(sprintf(
    paste0(
      "%s %s: %d features x %d donors; span %.17g (limma adaptive %.17g, %s), %d iterations; ",
      "max |d normalized| = %.3g (feature %s, donor %s), median %.3g, max |limma correction| = %.3g; %.1f s -> %s\n"
    ),
    cohort, omic, nrow(lcpm), ncol(lcpm), settings$span, adaptive, if (span_ok) "equal" else "DIFFERENT",
    settings$iterations, d_max, rownames(d)[where[1]], colnames(d)[where[2]], median(d), max(abs(r_norm - lcpm)),
    elapsed, if (pass) "PASS" else "FAIL"
  ))
}
cat(sprintf("limma %s, tolerance %g\n", as.character(packageVersion("limma")), tol))
if (!ok) quit(status = 1)
