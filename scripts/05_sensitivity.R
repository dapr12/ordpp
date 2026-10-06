# Sensitivity of the bfi selected setting to numerical safeguards and grid.
source("scripts/common.R")
if (!file.exists("results/bfi.rds")) stop("run 02_application_bfi.R first")
b <- readRDS("results/bfi.rds")
data("bfi", package = "psychTools")
items <- paste0("N", 1:5)
x <- bfi[complete.cases(bfi[, items]), items]
x <- as.data.frame(lapply(x, ordered, levels = 1:6))
n <- nrow(x); set.seed(2026)   # identical split to 02_application_bfi.R
idx <- sample(rep(c("train", "valid", "test"), times = c(n %/% 2, n %/% 4, n - n %/% 2 - n %/% 4)))
xtr <- x[idx == "train", ]; xva <- x[idx == "valid", ]
ref <- b$cands[[b$best]]
base <- list(scoring = ref$scoring, gamma = ref$settings$gamma,
             nstart = ref$settings$nstart, maxit = ref$settings$maxit)
if (QUICK) { base$nstart <- 1; base$maxit <- 20 }
variants <- list(
  "default" = list(),
  "grid 0.1-2, L=15" = list(frequencies = seq(0.1, 2, length.out = 15)),
  "grid 0.5-5, L=15" = list(frequencies = seq(0.5, 5, length.out = 15)),
  "grid 0.2-3, L=30" = list(frequencies = seq(0.2, 3, length.out = 30)),
  "min gap 0.01" = list(min_gap = 0.01), "min gap 0.10" = list(min_gap = 0.10),
  "variance floor 0.01" = list(variance_floor = 0.01),
  "variance floor 0.20" = list(variance_floor = 0.20),
  "shape index" = list(index = "shape"))
rows <- lapply(names(variants), function(v) {
  message("variant: ", v)
  f <- suppressWarnings(do.call(ordpp, c(list(x = xtr, seed = 99), base, variants[[v]])))
  data.frame(Variant = v,
    "Valid. index" = validate_ordpp(f, xva, type = "full"),
    "Var ratio" = decompose_ordpp(f, xva)$variance_ratio,
    "|cos| with selected" = abs(sum(f$direction * ref$direction)),
    "Max score change" = max(abs(unlist(f$scores) - unlist(ref$scores))),
    check.names = FALSE)
})
tab <- do.call(rbind, rows)
write.csv(tab, "results/sensitivity.csv", row.names = FALSE)
write_tex_table(tab, "results/tab-sensitivity.tex", digits = 3)
write_macros(list(sensMinCos = fmt(min(tab[["|cos| with selected"]]), 3)),
             "results/macros-sensitivity.tex")
