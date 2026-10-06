# Shared helpers for the replication scripts. Run scripts from the project root:
#   Rscript scripts/02_application_bfi.R
# Set ORDPP_QUICK=1 for a fast smoke test (tiny budgets; numbers not for the paper).
suppressPackageStartupMessages(library(ordpp))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

QUICK <- identical(Sys.getenv("ORDPP_QUICK"), "1")
budget <- function(full, quick) if (QUICK) quick else full
if (QUICK) message("ORDPP_QUICK=1: smoke-test budgets; do not report these numbers")

# ---- output helpers --------------------------------------------------------
tex_escape <- function(s) gsub("([_%&#$])", "\\\\\\1", s)

fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, digits = d, format = "f"))
fmt_p <- function(p) ifelse(is.na(p), "--",
  ifelse(p < 0.001, "<0.001", formatC(p, digits = 3, format = "f")))

# Writes LaTeX macros (\newcommand) read by the manuscript before placeholders.
write_macros <- function(values, file) {
  stopifnot(all(grepl("^[A-Za-z]+$", names(values))))
  vals <- vapply(values, function(v) as.character(v)[1], character(1))
  writeLines(sprintf("\\newcommand{\\%s}{%s}", names(values), vals), file)
}

# Writes a booktabs tabular (no float); the caption lives in the manuscript.
write_tex_table <- function(df, file, digits = 3, align = NULL) {
  num <- vapply(df, is.numeric, logical(1))
  body <- df
  for (j in which(num)) {
    v <- df[[j]]
    whole <- all(is.na(v) | abs(v - round(v)) < 1e-9)
    body[[j]] <- fmt(v, if (whole) 0 else digits)
  }
  for (j in which(!num)) body[[j]] <- tex_escape(as.character(df[[j]]))
  body[] <- lapply(body, function(v) ifelse(grepl("^<", v), paste0("$", v, "$"), v))
  if (is.null(align)) align <- paste0(ifelse(num, "r", "l"), collapse = "")
  head <- paste(tex_escape(names(df)), collapse = " & ")
  rows <- apply(as.matrix(body), 1, paste, collapse = " & ")
  writeLines(c(sprintf("\\begin{tabular}{%s}", align), "\\toprule",
               paste0(head, " \\\\"), "\\midrule", paste0(rows, " \\\\"),
               "\\bottomrule", "\\end{tabular}"), file)
}

capture_to <- function(file, expr) {
  out <- capture.output(expr)
  writeLines(out, file)
  invisible(out)
}

# ---- comparator maps -------------------------------------------------------
rank_codes <- function(x) as.data.frame(lapply(x, function(v)
  as.integer(if (is.ordered(v)) v else factor(v, levels = sort(unique(v))))))

equal_scores <- function(x) lapply(x, function(v)
  if (is.ordered(v)) seq_along(levels(v)) else seq_along(sort(unique(v))))

scores_from_transformed <- function(codes, Tm)
  lapply(seq_len(ncol(codes)), function(j) as.numeric(tapply(Tm[, j], codes[[j]], mean)))

pick_transformed <- function(obj, n, candidates = c("transform", "Z", "transformed")) {
  for (nm in candidates) {
    v <- obj[[nm]]
    if (!is.null(v)) { v <- as.matrix(v); if (nrow(v) == n) return(v) }
  }
  stop("transformed data not found; inspect names() of the fitted object")
}

pc1_map <- function(x, scores, label) {
  m0 <- as_ordpp_map(x, scores, rep(1, ncol(x)))
  a <- prcomp(m0$transformed, center = FALSE)$rotation[, 1]
  as_ordpp_map(x, scores, a, label = label)
}

map_pca <- function(x) pc1_map(x, equal_scores(x), "PCA, equal scores")

# Polychoric loadings applied to equal scores: an explicit, stated scoring rule.
map_polychoric <- function(x) {
  if (!requireNamespace("psych", quietly = TRUE)) return(NULL)
  R <- suppressWarnings(psych::polychoric(rank_codes(x))$rho)
  a <- eigen(R, symmetric = TRUE)$vectors[, 1]
  as_ordpp_map(x, equal_scores(x), a, label = "Polychoric PCA")
}

# Ordinal optimal-scaling PCA (Gifi). Check field names for your Gifi version.
map_gifi <- function(x) {
  if (!requireNamespace("Gifi", quietly = TRUE)) return(NULL)
  codes <- rank_codes(x)
  g <- Gifi::princals(codes, ndim = 1, levels = "ordinal")
  Tm <- pick_transformed(g, nrow(codes))
  pc1_map(x, scores_from_transformed(codes, Tm), "Gifi princals")
}

# Penalized ordinal PCA (ordPens; Hoshiyar et al.). Check ?ordPens::ordPCA.
map_ordpens <- function(x, lambda = 1) {
  if (!requireNamespace("ordPens", quietly = TRUE)) return(NULL)
  H <- as.matrix(rank_codes(x))
  o <- ordPens::ordPCA(H, p = 1, lambda = lambda, Ks = apply(H, 2, max),
                       constr = rep(TRUE, ncol(H)))
  sc <- if (!is.null(o$qs)) lapply(o$qs, function(q) as.numeric(as.matrix(q)[, 1]))
        else scores_from_transformed(rank_codes(x), pick_transformed(o, nrow(H)))
  pc1_map(x, sc, "ordPens ordPCA")
}

safe_map <- function(f, ..., label = deparse(substitute(f))) tryCatch(f(...), error = function(e) {
  message("comparator skipped (", label, "): ", conditionMessage(e)); NULL })

# Evaluate a fitted or external map on held-out data.
evaluate_map <- function(map, xtest, B, latent = NULL, seed = 1) {
  ht <- heldout_test_ordpp(map, xtest, B = B, seed = seed, type = "full")
  dc <- decompose_ordpp(map, xtest)
  data.frame(method = map$scoring, train_index = map$index,
             test_index = dc$full, test_shape = dc$shape,
             variance_ratio = dc$variance_ratio, p_value = ht$p.value,
             abs_latent_cor = if (is.null(latent)) NA_real_ else
               abs(stats::cor(predict(map, xtest), latent)))
}
