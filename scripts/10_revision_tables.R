# Revision tables from saved results (no refitting; about a minute).
#  * bfi: test-set evaluation of the smoother gamma = 0.01 fit (sensitivity)
#  * sensitivity table reordered, frequency-grid rows marked
#  * benchmark tables with Monte Carlo standard errors
source("scripts/common.R")

## 1. bfi gamma = 0.01 as a labelled sensitivity analysis -------------------
b <- readRDS("results/bfi.rds")
data("bfi", package = "psychTools")
items <- paste0("N", 1:5)
x <- bfi[complete.cases(bfi[, items]), items]
x <- as.data.frame(lapply(x, ordered, levels = 1:6))
n <- nrow(x); set.seed(2026)   # identical split to 02_application_bfi.R
idx <- sample(rep(c("train", "valid", "test"), times = c(n %/% 2, n %/% 4, n - n %/% 2 - n %/% 4)))
xtr <- x[idx == "train", ]; xte <- x[idx == "test", ]
g <- b$grid
ialt <- which(g$scoring == g$scoring[b$best] & abs(g$gamma - 0.01) < 1e-12)
alt <- b$cands[[ialt]]; sel <- b$cands[[b$best]]
keyed <- as_ordpp_map(xtr, equal_scores(xtr), rep(1, length(items)), label = "Unit-weighted sum")
alt$scoring <- "ordpp item, gamma 0.01 (sensitivity)"
r <- evaluate_map(alt, xte, B = 999)
r$abs_latent_cor <- NULL
r$cor_keyed <- stats::cor(predict(alt, xte), predict(keyed, xte))
r$p_value <- fmt_p(r$p_value)
tt <- read.csv("results/bfi-test.csv", check.names = FALSE)
tt <- tt[tt$method != r$method, , drop = FALSE]
tt <- rbind(tt[1, ], r[, names(tt)], tt[-1, ])      # place right after the selected map
write.csv(tt, "results/bfi-test.csv", row.names = FALSE)
tt$p_value <- fmt_p(as.numeric(tt$p_value))
names(tt) <- c("Map", "Train index", "Test index", "Test shape", "Var ratio",
               "Held-out p", "cor(sum)")
write_tex_table(tt, "results/tab-bfi-test.tex", digits = 3)
write_macros(list(
  bfiAltTest = fmt(r$test_index, 4), bfiAltShape = fmt(r$test_shape, 4),
  bfiAltVarRatio = fmt(r$variance_ratio, 2),
  bfiAltCos = fmt(abs(sum(alt$direction * sel$direction)), 4),
  bfiAltScoreDiff = fmt(max(abs(unlist(alt$scores) - unlist(sel$scores))), 2)),
  "results/macros-bfi-alt.tex")

## 2. Sensitivity table: comparable columns first, grid rows marked -----------
se <- read.csv("results/sensitivity.csv", check.names = FALSE)
grid_row <- grepl("^grid", se$Variant)
lab <- ifelse(grid_row, paste0(tex_escape(se$Variant), "$^\\dagger$"), tex_escape(se$Variant))
lines <- c("\\begin{tabular}{lrrrr}", "\\toprule",
  "Variant & $|\\cos|$ with selected & Max score change & Var ratio & Valid.\\ index \\\\",
  "\\midrule")
for (i in seq_len(nrow(se))) {
  if (i > 1 && grid_row[i] != grid_row[i - 1]) lines <- c(lines, "\\addlinespace")
  lines <- c(lines, sprintf("%s & %s & %s & %s & %s \\\\", lab[i],
    fmt(se[["|cos| with selected"]][i], 3), fmt(se[["Max score change"]][i], 3),
    fmt(se[["Var ratio"]][i], 3), fmt(se[["Valid. index"]][i], 3)))
}
writeLines(c(lines, "\\bottomrule", "\\end{tabular}"), "results/tab-sensitivity.tex")

## 3. Benchmark tables with Monte Carlo standard errors ----------------------
bs <- read.csv("results/benchmark-summary.csv")
nmax <- max(bs$n); s <- bs[bs$n == nmax & bs$scenario != "independent", ]
scen <- c("mixture", "rare", "two_factor", "response_style")
cell <- function(m, se, d) paste0(fmt(m, d), " (", fmt(se, d), ")")
mk <- function(col, secol, file, d = 3) {
  meths <- sort(unique(s$method))
  rows <- vapply(meths, function(m) {
    v <- vapply(scen, function(sc) {
      r <- s[s$method == m & s$scenario == sc, ]
      if (!nrow(r)) "--" else cell(r[[col]], r[[secol]], d)
    }, character(1))
    paste(c(tex_escape(m), v), collapse = " & ")
  }, character(1))
  writeLines(c("\\begin{tabular}{lllll}", "\\toprule",
    "Method & mixture & rare & two factor & response style \\\\", "\\midrule",
    paste0(rows, " \\\\"), "\\bottomrule", "\\end{tabular}"), file)
}
mk("test_index", "test_index_mcse", "results/tab-bench-index.tex")
mk("abs_cor", "abs_cor_mcse", "results/tab-bench-cor.tex")
co <- s[s$method == "ordpp common", ]; it <- s[s$method == "ordpp item", ]
co <- co[match(scen, co$scenario), ]; it <- it[match(scen, it$scenario), ]
d <- co$test_index - it$test_index
write_macros(list(benchCIDiffMax = fmt(max(abs(d)), 3),
                  benchCIMcseMax = fmt(max(sqrt(co$test_index_mcse^2 + it$test_index_mcse^2)), 3)),
             "results/macros-bench-mcse.tex")
message("revision tables written")
