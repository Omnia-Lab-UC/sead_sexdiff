#!/usr/bin/env Rscript
# limma-voom side of `sead-sexdiff stage-de-limma` (data.py: stage_de_limma_fit,
# run_stage_de_limma). It exchanges files with Python, so one run can be
# repeated by hand from a debug directory.
#
# A model directory holds what Python writes:
#   counts.csv         features x donors, integer counts (export and donor order)
#   export_tmm.csv     donor_id, factor: the export manifest's TMM factors
#   design_matrix.csv  donors x design columns (%.17g): the formulaic matrix
#                      that pydeseq2 fits in `stage-de`, with its column names
#   design_table.csv   donors x design variables (sex, stage, stage_num, covariates)
#   settings.json      model, omic, formula, factor levels, tested and joint
#                      columns, alpha, the TMM tolerance and STAGE_DE_LIMMA_METHOD
# and, after a fit, what this script writes back:
#   tmm.csv, coef__<label>.csv, joint__<label>.csv, sample_weights.csv, fit.json
#
# Usage (inside the Nix shell; one BLAS thread):
#   Rscript --vanilla stage_de_limma.R <model_dir>                    fit one model
#   Rscript --vanilla stage_de_limma.R --check <omic_dir> <model_dir>...
#       read <omic_dir>/counts.csv and export_tmm.csv once, check the donor ids,
#       the TMM factors and every model's design; write <omic_dir>/check.json
#   Rscript --vanilla stage_de_limma.R --versions                     R, limma, edgeR
# A failed check or fit stops with a message naming the omic and model (exit 1).

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})

fail <- function(...) stop(..., call. = FALSE)

versions <- function() {
  list(
    R = paste(R.version$major, R.version$minor, sep = "."),
    limma = as.character(packageVersion("limma")),
    edgeR = as.character(packageVersion("edgeR")),
    jsonlite = as.character(packageVersion("jsonlite"))
  )
}

write_json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, digits = I(17), null = "null", na = "null", pretty = TRUE), path)
}

# 17 significant digits, so every double round-trips through Python's float()
fmt <- function(x) {
  out <- formatC(x, digits = 17, format = "g", width = 1)
  out[is.na(x)] <- "NA"
  out
}

write_table <- function(columns, path) {
  cells <- lapply(columns, function(x) if (is.character(x)) x else fmt(x))
  lines <- do.call(paste, c(cells, sep = ","))
  writeLines(c(paste(names(columns), collapse = ","), lines), path)
}

# scan() rather than read.csv(): the ATAC matrices have ~130k rows
read_counts <- function(path, where) {
  con <- file(path, "r")
  on.exit(close(con))
  header <- strsplit(readLines(con, n = 1), ",", fixed = TRUE)[[1]]
  cols <- scan(con, what = c(list(""), rep(list(0), length(header) - 1)), sep = ",", quiet = TRUE)
  m <- do.call(cbind, cols[-1])
  dimnames(m) <- list(cols[[1]], header[-1])
  if (anyNA(m) || any(m < 0) || any(m != round(m))) fail(where, ": ", path, " holds a value that is not a non-negative integer")
  if (anyDuplicated(rownames(m))) fail(where, ": ", path, " repeats a feature name")
  m
}

read_settings <- function(dir) jsonlite::fromJSON(file.path(dir, "settings.json"), simplifyVector = FALSE)

# The installed defaults must be the predeclared settings (STAGE_DE_LIMMA_METHOD)
check_defaults <- function(method, where) {
  default_of <- function(fun, arg) {
    d <- formals(fun)[[arg]]
    if (is.null(d)) NULL else eval(d)
  }
  same <- function(a, b) {
    if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
    if (is.numeric(a) && is.numeric(b)) return(length(a) == length(b) && all(a == b))
    identical(a, b)
  }
  checks <- list()
  add <- function(name, fun, arg, value) {
    d <- default_of(fun, arg)
    if (arg == "method" && is.character(d) && length(d) > 1) d <- d[1]
    ok <- same(d, value)
    checks[[name]] <<- ok
    if (!ok) {
      fail(where, ": the installed default of ", name, " is ", deparse(d), " but STAGE_DE_LIMMA_METHOD records ",
           deparse(value), "; the method was predeclared on limma 3.68.4 / edgeR 4.10.1 defaults")
    }
  }
  n <- method$normalization
  for (arg in c("logratioTrim", "sumTrim", "doWeighting", "Acutoff")) {
    add(paste0("calcNormFactors.", arg), edgeR:::calcNormFactors.DGEList, arg, unlist(n[[arg]]))
  }
  add("calcNormFactors.method", edgeR:::calcNormFactors.DGEList, "method", n$method)
  add("calcNormFactors.refColumn", edgeR:::calcNormFactors.DGEList, "refColumn", NULL)
  w <- method$weights
  for (arg in c("normalize.method", "span", "adaptive.span", "var.design", "var.group", "method", "maxiter", "tol")) {
    add(paste0("voomWithQualityWeights.", arg), voomWithQualityWeights, arg, if (is.null(w[[arg]])) NULL else unlist(w[[arg]]))
  }
  add("arrayWeights.prior.n", arrayWeights, "prior.n", w$prior.n)
  add("lmFit.method", lmFit, "method", method$fit$method)
  e <- method$moderation
  for (arg in c("proportion", "stdev.coef.lim", "trend", "robust", "legacy")) {
    add(paste0("eBayes.", arg), eBayes, arg, if (is.null(e[[arg]])) NULL else unlist(e[[arg]]))
  }
  checks
}

check_ids <- function(counts, export_tmm, design_matrix, design_table, where) {
  donors <- colnames(counts)
  for (item in list(
    list("export_tmm.csv", export_tmm$donor_id),
    list("design_matrix.csv", if (is.null(design_matrix)) donors else rownames(design_matrix)),
    list("design_table.csv", if (is.null(design_table)) donors else rownames(design_table))
  )) {
    if (!identical(as.character(item[[2]]), donors)) {
      fail(where, ": the donors of ", item[[1]], " are not the count columns in the same order")
    }
  }
  donors
}

tmm_check <- function(counts, export_tmm, settings, where) {
  n <- settings$method$normalization
  y <- DGEList(counts = counts, lib.size = colSums(counts))
  # edgeR 4 renamed calcNormFactors to normLibSizes; the alias passes every
  # argument through and only adds a message
  y <- suppressMessages(calcNormFactors(y, method = n$method, refColumn = NULL, logratioTrim = n$logratioTrim,
                                        sumTrim = n$sumTrim, doWeighting = n$doWeighting, Acutoff = n$Acutoff))
  diff <- abs(y$samples$norm.factors - export_tmm$factor)
  max_diff <- max(diff)
  if (!(max_diff <= settings$tmm_atol)) {
    fail(where, ": the TMM factors computed on the tested counts differ from the export manifest's by up to ",
         format(max_diff, digits = 3), " (absolute; tolerance ", settings$tmm_atol, ", donor ",
         colnames(counts)[which.max(diff)], "); the counts changed since the export; rerun `sead-sexdiff export` ",
         "(then `sead-sexdiff analyze`)")
  }
  list(y = y, table = list(donor_id = colnames(counts), lib_size = y$samples$lib.size,
                           factor = y$samples$norm.factors, export_factor = export_tmm$factor, abs_diff = diff),
       result = list(max_abs_diff = max_diff, atol = settings$tmm_atol, passed = TRUE,
                     lib_size = "column sums of the tested counts", n_features = nrow(counts)))
}

read_design <- function(dir, where) {
  dm <- read.csv(file.path(dir, "design_matrix.csv"), check.names = FALSE, colClasses = "character")
  X <- sapply(dm[-1], as.numeric)
  X <- matrix(X, nrow = nrow(dm), dimnames = list(dm[[1]], names(dm)[-1]))
  if (anyNA(X)) fail(where, ": design_matrix.csv holds a value that is not a number")
  dt <- read.csv(file.path(dir, "design_table.csv"), check.names = FALSE, colClasses = "character")
  rownames(dt) <- dt[[1]]
  dt <- dt[-1]
  list(X = X, table = dt)
}

# formulaic's names → R's: Intercept → (Intercept), var[T.level] → varlevel
r_name <- function(x) {
  x <- gsub("([A-Za-z_][A-Za-z0-9_.]*)\\[T\\.([^]]+)\\]", "\\1\\2", x)
  x[x == "Intercept"] <- "(Intercept)"
  x
}

design_check <- function(design, settings, where) {
  X <- design$X
  dt <- design$table
  data <- data.frame(row.names = rownames(dt))
  for (v in names(dt)) {
    if (v %in% names(settings$levels)) {
      levels <- unlist(settings$levels[[v]])
      bad <- setdiff(unique(dt[[v]]), levels)
      if (length(bad)) fail(where, ": design_table.csv has ", v, " value(s) ", paste(bad, collapse = ", "), " outside the levels")
      data[[v]] <- factor(dt[[v]], levels = levels)
    } else {
      data[[v]] <- as.numeric(dt[[v]])
    }
  }
  R <- model.matrix(as.formula(settings$formula), data = data)
  mapped <- r_name(colnames(X))
  if (!setequal(mapped, colnames(R)) || length(mapped) != ncol(R)) {
    fail(where, ": the design rebuilt by model.matrix (", paste(colnames(R), collapse = ", "),
         ") does not have the columns of design_matrix.csv under the name map (", paste(mapped, collapse = ", "),
         "); formulaic and R disagree on the coding of ", settings$model)
  }
  diffs <- vapply(seq_along(mapped), function(j) max(abs(X[, j] - R[, mapped[j]])), numeric(1))
  if (any(diffs != 0)) {
    fail(where, ": the design rebuilt by model.matrix differs from design_matrix.csv in column(s) ",
         paste(colnames(X)[diffs != 0], collapse = ", "), " of ", settings$model,
         " (the treatment coding or the levels disagree)")
  }
  list(passed = TRUE, r_columns = colnames(R), name_map = as.list(setNames(mapped, colnames(X))),
       order_equal = identical(mapped, colnames(R)), max_abs_diff = max(diffs), rank = qr(X)$rank,
       n_columns = ncol(X))
}

exact_f <- function(efit, v, X, columns, limma_f, alpha) {
  idx <- match(columns, colnames(X))
  b <- efit$coefficients[, columns, drop = FALSE]
  W <- v$weights
  G <- nrow(b)
  fx <- numeric(G)
  su_diff <- 0
  for (g in seq_len(G)) {
    q <- qr(X * sqrt(W[g, ]))
    V <- matrix(0, ncol(X), ncol(X))
    V[q$pivot, q$pivot] <- chol2inv(qr.R(q))
    Vb <- V[idx, idx]
    su_diff <- max(su_diff, abs(sqrt(diag(Vb)) / efit$stdev.unscaled[g, columns] - 1))
    fx[g] <- drop(crossprod(b[g, ], solve(Vb, b[g, ]))) / (length(idx) * efit$s2.post[g])
  }
  df2 <- efit$df.prior + efit$df.residual
  p <- if (df2[1] > 1e6) pchisq(length(idx) * fx, length(idx), lower.tail = FALSE) else pf(fx, length(idx), df2, lower.tail = FALSE)
  adj <- p.adjust(p, method = "BH")
  rel <- abs(fx - limma_f$F) / limma_f$F
  changed <- which((adj <= alpha) != (limma_f$adj <= alpha))
  list(
    definition = paste(
      "per feature, b' V^-1 b / (2 s2.post) with V the coefficients' block of (X' W_g X)^-1 from the feature's own",
      "weights; p from F(2, df.prior + df.residual) as limma's; BH over all tested features"
    ),
    max_rel_diff = max(rel),
    feature_at_max = rownames(b)[which.max(rel)],
    median_rel_diff = median(rel),
    stdev_unscaled_max_rel_diff = su_diff,
    n_significant_exact = sum(adj <= alpha),
    n_changed_calls = length(changed),
    changed_calls = lapply(changed, function(g) list(
      feature = rownames(b)[g], limma_adj_p = limma_f$adj[g], exact_adj_p = adj[g],
      limma_F = limma_f$F[g], exact_F = fx[g]
    ))
  )
}

fit_model <- function(dir) {
  t0 <- proc.time()[["elapsed"]]
  settings <- read_settings(dir)
  where <- paste0("stage_de_limma.R: ", settings$omic, " ", settings$model)
  defaults <- check_defaults(settings$method, where)
  counts <- read_counts(file.path(dir, "counts.csv"), where)
  export_tmm <- read.csv(file.path(dir, "export_tmm.csv"), colClasses = c("character", "numeric"))
  design <- read_design(dir, where)
  check_ids(counts, export_tmm, design$X, design$table, where)
  tmm <- tmm_check(counts, export_tmm, settings, where)
  dcheck <- design_check(design, settings, where)
  X <- design$X
  w <- settings$method$weights
  e <- settings$method$moderation
  t1 <- proc.time()[["elapsed"]]
  v <- voomWithQualityWeights(
    tmm$y, design = X, normalize.method = w$normalize.method, span = w$span, adaptive.span = w$adaptive.span,
    var.design = NULL, var.group = NULL, method = w$method, maxiter = w$maxiter, tol = w$tol
  )
  fit <- lmFit(v, X, method = settings$method$fit$method)
  efit <- eBayes(fit, proportion = e$proportion, stdev.coef.lim = unlist(e$stdev.coef.lim), trend = e$trend,
                 robust = e$robust, legacy = NULL)
  t2 <- proc.time()[["elapsed"]]
  features <- rownames(counts)
  alpha <- settings$alpha
  tests <- list()
  for (tst in settings$tested) {
    tt <- topTable(efit, coef = tst$column, number = Inf, sort.by = "none", adjust.method = "BH")
    if (!identical(rownames(tt), features)) fail(where, ": topTable reordered the features of ", tst$label)
    se <- sqrt(efit$s2.post) * efit$stdev.unscaled[, tst$column]
    write_table(list(feature = features, AveExpr = tt$AveExpr, logFC = tt$logFC, SE = se, t = tt$t,
                     P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, B = tt$B),
                file.path(dir, paste0("coef__", tst$label, ".csv")))
    tests[[tst$label]] <- list(column = tst$column, n_significant = sum(tt$adj.P.Val <= alpha),
                               n_adj_missing = sum(is.na(tt$adj.P.Val)), min_adj_p = min(tt$adj.P.Val))
  }
  joint <- list()
  for (jt in settings$joint) {
    columns <- unlist(jt$columns)
    tt <- topTable(efit, coef = columns, number = Inf, sort.by = "none", adjust.method = "BH")
    if (!identical(rownames(tt), features)) fail(where, ": topTable reordered the features of ", jt$label)
    # limma's F as Context describes it: the moderated t's combined through the
    # correlation of the unweighted design's coefficients
    tm <- efit$t[, columns]
    Q <- solve(cov2cor(efit$cov.coefficients[columns, columns]))
    f_re <- rowSums((tm %*% Q) * tm) / length(columns)
    out <- list(feature = features, AveExpr = tt$AveExpr)
    for (k in seq_along(columns)) out[[paste0("logFC_", jt$coefficients[[k]])]] <- efit$coefficients[, columns[k]]
    out <- c(out, list(F = tt$F, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val))
    write_table(out, file.path(dir, paste0("joint__", jt$label, ".csv")))
    df2 <- efit$df.prior + efit$df.residual
    joint[[jt$label]] <- list(
      columns = columns, coefficients = unlist(jt$coefficients), df1 = length(columns), df2 = unique(df2),
      n_significant = sum(tt$adj.P.Val <= alpha), n_adj_missing = sum(is.na(tt$adj.P.Val)),
      min_adj_p = min(tt$adj.P.Val),
      limma_f_recomputed_max_rel_diff = max(abs(f_re - tt$F) / tt$F),
      exact = exact_f(efit, v, X, columns, list(F = tt$F, adj = tt$adj.P.Val), alpha)
    )
  }
  t3 <- proc.time()[["elapsed"]]
  sw <- v$targets$sample.weights
  write_table(list(donor_id = colnames(counts), weight = sw), file.path(dir, "sample_weights.csv"))
  write_table(tmm$table, file.path(dir, "tmm.csv"))
  dfr <- efit$df.residual
  write_json(list(
    versions = versions(),
    model = settings$model,
    omic = settings$omic,
    formula = settings$formula,
    n_features = nrow(counts),
    n_donors = ncol(counts),
    installed_defaults_equal_method = defaults,
    tmm = tmm$result,
    design_check = dcheck,
    voom = list(span = v$span, adaptive_span = w$adaptive.span, log_cpm = "log2((count + 0.5) / (lib.size * norm.factor + 1) * 1e6)"),
    sample_weights = list(min = min(sw), median = median(sw), max = max(sw)),
    ebayes = list(
      df_prior = efit$df.prior, s2_prior = efit$s2.prior, df_residual = unique(dfr),
      df_total = unique(efit$df.total), legacy_requested = NULL,
      legacy_resolved = identical(min(dfr[dfr > 0]), max(dfr[dfr > 0])), var_prior = efit$var.prior
    ),
    tests = tests,
    joint = joint,
    runtime_seconds = list(read_and_check = t1 - t0, fit = t2 - t1, tests_and_exact_f = t3 - t2, total = t3 - t0)
  ), file.path(dir, "fit.json"))
}

check_omic <- function(omic_dir, model_dirs) {
  first <- read_settings(model_dirs[1])
  where <- paste0("stage_de_limma.R: ", first$omic)
  counts <- read_counts(file.path(omic_dir, "counts.csv"), where)
  export_tmm <- read.csv(file.path(omic_dir, "export_tmm.csv"), colClasses = c("character", "numeric"))
  check_ids(counts, export_tmm, NULL, NULL, where)
  defaults <- check_defaults(first$method, where)
  tmm <- tmm_check(counts, export_tmm, first, where)
  write_table(tmm$table, file.path(omic_dir, "tmm.csv"))
  designs <- list()
  for (d in model_dirs) {
    s <- read_settings(d)
    w <- paste0("stage_de_limma.R: ", s$omic, " ", s$model)
    des <- read_design(d, w)
    check_ids(counts, export_tmm, des$X, des$table, w)
    designs[[s$model]] <- design_check(des, s, w)
  }
  write_json(list(versions = versions(), omic = first$omic, n_features = nrow(counts), n_donors = ncol(counts),
                  installed_defaults_equal_method = defaults, tmm = tmm$result, design_checks = designs),
             file.path(omic_dir, "check.json"))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && args[1] == "--versions") {
  cat(jsonlite::toJSON(versions(), auto_unbox = TRUE), "\n")
} else if (length(args) >= 3 && args[1] == "--check") {
  check_omic(args[2], args[-(1:2)])
} else if (length(args) == 1) {
  fit_model(args[1])
} else {
  stop("usage: stage_de_limma.R <model_dir> | --check <omic_dir> <model_dir>... | --versions")
}
