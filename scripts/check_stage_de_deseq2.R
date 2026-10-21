#!/usr/bin/env Rscript
# Cross-check of one pydeseq2 fit of `sead-sexdiff stage-de` against R's
# DESeq2. Reads the files written by `run_stage_de(debug_dir=...)`:
#   <debug_dir>/<cohort>_<omic>_counts.csv                 donors x tested features (raw counts)
#   <debug_dir>/<cohort>_<omic>_design.csv                 donor_id, sex, stage, stage_num, the omic's covariates
#   <debug_dir>/<cohort>_<omic>_<model>_trend.csv          the dispersion trend pydeseq2 used (type, a0, a1, mean)
#   <debug_dir>/<cohort>_<omic>_<model>_dispersions.csv    pydeseq2's per-feature dispersions
#   <debug_dir>/<cohort>_<omic>_<model>_size_factors.csv   pydeseq2's size factors
# and the step's tables <stage_de_dir>/<cohort>/<omic>/<model>__<coefficient>.csv.
#
# Two comparisons per tested coefficient, with sex (reference Female) and stage
# (reference 0) as factors and the covariates entered on their recorded scale
# in the design file's column order:
#   full     DESeq() on the same counts with the dispersion trend type pydeseq2
#            used (pydeseq2 switches a non-converged parametric trend to a mean
#            trend, R's DESeq2 to a local one, so the type is passed
#            explicitly), Wald test, no LFC prior, results(alpha = 0.05). This
#            is task 3.4's criterion: log2FoldChange and the Wald statistic
#            within 1e-3 relative (+ 1e-6 absolute), and the same features at
#            padj <= 0.05.
#   plug-in  R's nbinomWaldTest() with pydeseq2's size factors and final
#            dispersions, so the GLM, Wald test and independent filtering are
#            compared without the dispersion estimation; same tolerance.
# The tolerance is evaluated on the features whose fitted means are all at
# least 0.5 and whose pydeseq2 LFC fit converged: DESeq2 floors mu at 0.5
# (minmu) in the Wald standard errors and pydeseq2 does not, and a feature
# whose pydeseq2 IRLS did not converge (`lfc_converged = 0` in the dispersion
# file) was fitted by its fallback optimizer, so both differ by construction.
# They are counted and reported apart, and they still enter the comparison of
# the significant sets and of the non-finite values.
# Exit status: 0 when both comparisons pass; 2 when only the full comparison's
# tolerance fails (the plug-in passes and the significant sets of the full fit
# are equal: the difference lies in the dispersion estimation); 3 when the
# dispersion trend types differ (the plug-in is still run and reported); 4
# when a comparison fails otherwise; 1 for an R error.
#
# With [out_dir], writes per tested coefficient <out_dir>/<cohort>_<omic>_<model>__<coefficient>_deseq2.csv
# (feature; log2FoldChange, stat, padj of the full and the plug-in fit;
# DESeq2's betaConv of each) to locate the differences.
#
# Usage (inside the Nix shell):
#   Rscript scripts/check_stage_de_deseq2.R <debug_dir> <stage_de_dir> [cohort] [omic] [model] [out_dir]

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: check_stage_de_deseq2.R <debug_dir> <stage_de_dir> [cohort] [omic] [model] [out_dir]")
debug_dir <- args[1]
stage_de_dir <- args[2]
cohort <- if (length(args) >= 3) args[3] else "cohort_30"
omic <- if (length(args) >= 4) args[4] else "rna"
model <- if (length(args) >= 5) args[5] else "factor_interaction"
out_dir <- if (length(args) >= 6) args[6] else NA

tol_rel <- 1e-3
tol_abs <- 1e-6
alpha <- 0.05

# STAGE_DE_MODELS: design terms and, per tested coefficient, the step's label
# and DESeq2's resultsNames() entry
models <- list(
  factor_additive = list(
    terms = c("sex", "stage"),
    tested = c(stage_1_vs_0 = "stage_1_vs_0", stage_2_vs_0 = "stage_2_vs_0")
  ),
  factor_interaction = list(
    terms = c("sex", "stage", "sex:stage"),
    tested = c(male_x_stage_1 = "sexMale.stage1", male_x_stage_2 = "sexMale.stage2")
  ),
  trend_additive = list(
    terms = c("sex", "stage_num"),
    tested = c(stage_num = "stage_num")
  ),
  trend_interaction = list(
    terms = c("sex", "stage_num", "sex:stage_num"),
    tested = c(male_x_stage_num = "sexMale.stage_num")
  )
)
if (!(model %in% names(models))) stop(sprintf("unknown model %s", model))
spec <- models[[model]]

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

stem <- file.path(debug_dir, paste0(cohort, "_", omic))
counts <- t(read_matrix(paste0(stem, "_counts.csv")))  # features x donors
if (max(counts) > .Machine$integer.max) stop("counts exceed the integer range")
storage.mode(counts) <- "integer"
design <- read.csv(paste0(stem, "_design.csv"), colClasses = c(donor_id = "character", sex = "character", stage = "character"))
rownames(design) <- design$donor_id
if (!identical(rownames(design), colnames(counts))) stop("the counts and the design list different donors")
design$sex <- factor(design$sex, levels = c("Female", "Male"))
design$stage <- factor(design$stage, levels = c("0", "1", "2"))
covariates <- setdiff(colnames(design), c("donor_id", "sex", "stage", "stage_num"))
formula <- as.formula(paste("~", paste(c(spec$terms, covariates), collapse = " + ")))

py_trend <- read.csv(paste0(stem, "_", model, "_trend.csv"))
fit_type <- as.character(py_trend$type[1])
if (!(fit_type %in% c("parametric", "mean"))) stop(sprintf("unexpected pydeseq2 trend type %s", fit_type))
py_disp <- read.csv(paste0(stem, "_", model, "_dispersions.csv"), colClasses = c(feature = "character"), na.strings = "")
py_sf <- read.csv(paste0(stem, "_", model, "_size_factors.csv"), colClasses = c(donor_id = "character"))
if (!identical(py_disp$feature, rownames(counts))) stop("the dispersion file lists other features than the counts")
if (!identical(py_sf$donor_id, colnames(counts))) stop("the size-factor file lists other donors than the counts")

t0 <- proc.time()[["elapsed"]]
dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = formula)
notes <- character()
dds <- withCallingHandlers(
  DESeq(dds, fitType = fit_type, quiet = TRUE),
  message = function(m) {
    notes <<- c(notes, trimws(conditionMessage(m)))
    invokeRestart("muffleMessage")
  }
)
elapsed <- proc.time()[["elapsed"]] - t0
r_fit_type <- attr(dispersionFunction(dds), "fitType")
missing <- setdiff(spec$tested, resultsNames(dds))
if (length(missing)) stop(sprintf("resultsNames lacks %s (has %s)", paste(missing, collapse = ", "), paste(resultsNames(dds), collapse = ", ")))
not_converged <- sum(!mcols(dds)$betaConv, na.rm = TRUE)
sf_rel <- max(abs(sizeFactors(dds) - py_sf$size_factor) / py_sf$size_factor)
disp_rel <- abs(dispersions(dds) - py_disp$dispersion) / py_disp$dispersion
cat(sprintf(
  paste0(
    "%s %s %s: %s; %d features x %d donors; DESeq2 %s; %.1f s; dispersion trend pydeseq2 %s, DESeq2 %s; rows not ",
    "converged in beta: %d; size factors max rel diff %.3g; final dispersions rel diff median %.3g, 99%% %.3g, max %.3g\n"
  ),
  cohort, omic, model, deparse(formula, width.cutoff = 500L), nrow(counts), ncol(counts),
  as.character(packageVersion("DESeq2")), elapsed, fit_type, r_fit_type, not_converged, sf_rel,
  median(disp_rel, na.rm = TRUE), unname(quantile(disp_rel, 0.99, na.rm = TRUE)), max(disp_rel, na.rm = TRUE)
))
for (n in unique(notes)) if (nchar(n)) cat(sprintf("  DESeq2 note: %s\n", n))
trend_same <- identical(r_fit_type, fit_type)
if (!trend_same) cat(sprintf("  trend types differ (pydeseq2 %s, DESeq2 %s): the full comparison is not made\n", fit_type, r_fit_type))
r_trend <- dispersionFunction(dds)
if (trend_same && fit_type == "parametric") {
  cat(sprintf("  trend coefficients a0 %.6g / %.6g, a1 %.6g / %.6g (pydeseq2 / DESeq2)\n",
              py_trend$a0[1], attr(r_trend, "coefficients")[["asymptDisp"]],
              py_trend$a1[1], attr(r_trend, "coefficients")[["extraPois"]]))
} else if (trend_same && fit_type == "mean") {
  cat(sprintf("  mean dispersion %.6g / %.6g (pydeseq2 / DESeq2)\n", py_trend$mean_dispersion[1], attr(r_trend, "mean")))
}

dds_plug <- dds
sizeFactors(dds_plug) <- setNames(py_sf$size_factor, py_sf$donor_id)
dispersions(dds_plug) <- py_disp$dispersion
dds_plug <- nbinomWaldTest(dds_plug, quiet = TRUE)
cat(sprintf("  plug-in rows not converged in beta: %d\n", sum(!mcols(dds_plug)$betaConv, na.rm = TRUE)))

py_lfc_nc <- !is.na(py_disp$lfc_converged) & py_disp$lfc_converged == 0
apart_of <- function(d) (apply(assays(d)[["mu"]], 1, min) < 0.5) | py_lfc_nc

compare <- function(py, r, low_mu) {
  same_finite <- identical(is.finite(py), is.finite(r))
  keep <- is.finite(py) & is.finite(r) & !low_mu
  low <- is.finite(py) & is.finite(r) & low_mu
  d <- abs(py[keep] - r[keep])
  bound <- tol_rel * abs(r[keep]) + tol_abs
  list(
    max_abs = if (any(keep)) max(d) else 0, n_fail = sum(d > bound),
    low_max_abs = if (any(low)) max(abs(py[low] - r[low])) else 0,
    pass = all(d <= bound) && same_finite, same_finite = same_finite
  )
}

describe <- function(label, res, py, low_mu) {
  lfc <- compare(py$log2FoldChange, res$log2FoldChange, low_mu)
  stat <- compare(py$stat, res$stat, low_mu)
  sig_py <- py$feature[!is.na(py$padj) & py$padj <= alpha]
  sig_r <- rownames(res)[!is.na(res$padj) & res$padj <= alpha]
  only_py <- setdiff(sig_py, sig_r)
  only_r <- setdiff(sig_r, sig_py)
  na_same <- identical(is.na(py$padj), is.na(res$padj))
  sets_equal <- length(only_py) == 0 && length(only_r) == 0
  list(
    line = sprintf(
      paste0(
        "log2FC max abs %.3g (%d beyond tolerance); stat max abs %.3g (%d beyond tolerance); %d features with ",
        "fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC %.3g, stat %.3g); padj <= %.2g: ",
        "%d pydeseq2, %d DESeq2, %d only ",
        "pydeseq2, %d only DESeq2; padj NA: %d pydeseq2, %d DESeq2 (%s features)"
      ),
      lfc$max_abs, lfc$n_fail, stat$max_abs, stat$n_fail, sum(low_mu), lfc$low_max_abs, stat$low_max_abs, alpha,
      length(sig_py), length(sig_r), length(only_py), length(only_r), sum(is.na(py$padj)), sum(is.na(res$padj)),
      if (na_same) "same" else "different"
    ),
    tol_pass = lfc$pass && stat$pass,
    sets_equal = sets_equal,
    only = sprintf("only pydeseq2: %s; only DESeq2: %s", paste(head(only_py, 10), collapse = ", "),
                   paste(head(only_r, 10), collapse = ", "))
  )
}

full_ok <- TRUE
plug_ok <- TRUE
sets_ok <- TRUE
for (label in names(spec$tested)) {
  py <- read.csv(
    file.path(stage_de_dir, cohort, omic, sprintf("%s__%s.csv", model, label)),
    colClasses = c(feature = "character"), na.strings = ""
  )
  res_full <- results(dds, name = spec$tested[[label]], alpha = alpha)
  res_plug <- results(dds_plug, name = spec$tested[[label]], alpha = alpha)
  if (!identical(py$feature, rownames(res_full))) stop(sprintf("%s: the pydeseq2 table lists other features", label))
  full <- describe(label, res_full, py, apart_of(dds))
  plug <- describe(label, res_plug, py, apart_of(dds_plug))
  full_ok <- full_ok && trend_same && full$tol_pass && full$sets_equal
  plug_ok <- plug_ok && plug$tol_pass && plug$sets_equal
  sets_ok <- sets_ok && full$sets_equal
  cat(sprintf("  %s (%s)\n", label, spec$tested[[label]]))
  if (trend_same) {
    cat(sprintf("    full:    %s -> %s\n", full$line, if (full$tol_pass && full$sets_equal) "PASS" else "FAIL"))
    if (!full$sets_equal) cat(sprintf("      %s\n", full$only))
  }
  cat(sprintf("    plug-in: %s -> %s\n", plug$line, if (plug$tol_pass && plug$sets_equal) "PASS" else "FAIL"))
  if (!plug$sets_equal) cat(sprintf("      %s\n", plug$only))
  if (!is.na(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(data.frame(
      feature = rownames(res_full),
      full_log2FoldChange = res_full$log2FoldChange, full_stat = res_full$stat, full_padj = res_full$padj,
      full_betaConv = mcols(dds)$betaConv,
      plug_log2FoldChange = res_plug$log2FoldChange, plug_stat = res_plug$stat, plug_padj = res_plug$padj,
      plug_betaConv = mcols(dds_plug)$betaConv
    ), file.path(out_dir, sprintf("%s_%s_%s__%s_deseq2.csv", cohort, omic, model, label)), row.names = FALSE)
  }
}
cat(sprintf(
  "tolerance %g relative + %g absolute on features with fitted mu >= 0.5 and a converged pydeseq2 LFC fit (%d of %d with a non-converged one); full %s, plug-in %s\n", tol_rel, tol_abs, sum(py_lfc_nc), length(py_lfc_nc),
  if (!trend_same) "not made (trend types differ)" else if (full_ok) "PASS" else "FAIL", if (plug_ok) "PASS" else "FAIL"
))
if (!trend_same) quit(status = 3)
if (full_ok && plug_ok) quit(status = 0)
if (plug_ok && sets_ok) quit(status = 2)
quit(status = 4)
