# Real-data application: SAPA/IPIP items in psychTools::bfi.
source("scripts/common.R")
if (!requireNamespace("psychTools", quietly = TRUE))
  stop("install.packages('psychTools') first")
NSTART <- budget(5, 2); MAXIT <- budget(300, 40)
B_HELD <- budget(999, 49); B_STAB <- budget(50, 3)

data("bfi", package = "psychTools")
items <- paste0("N", 1:5)
x_all <- bfi[, items, drop = FALSE]
n_raw <- nrow(x_all)
x <- x_all[complete.cases(x_all), , drop = FALSE]
x <- as.data.frame(lapply(x, ordered, levels = 1:6))
n <- nrow(x)

# 50/25/25 split into training, validation and final test rows.
set.seed(2026)
idx <- sample(rep(c("train", "valid", "test"), times = c(n %/% 2, n %/% 4, n - n %/% 2 - n %/% 4)))
xtr <- x[idx == "train", ]; xva <- x[idx == "valid", ]; xte <- x[idx == "test", ]

# Candidate settings, all fitted on training rows only.
grid <- rbind(data.frame(scoring = "linear", gamma = 0),
              expand.grid(scoring = c("common", "item"), gamma = c(0.001, 0.01, 0.1),
                          stringsAsFactors = FALSE))
cands <- lapply(seq_len(nrow(grid)), function(i) {
  message(sprintf("candidate %d of %d: %s, gamma = %g", i, nrow(grid),
                  grid$scoring[i], grid$gamma[i]))
  ordpp(xtr, scoring = grid$scoring[i], gamma = grid$gamma[i],
        nstart = NSTART, maxit = MAXIT, seed = 300 + i)
})
grid$train_index <- vapply(cands, `[[`, numeric(1), "index")
grid$valid_index <- vapply(cands, validate_ordpp, numeric(1), newdata = xva)
grid$roughness <- vapply(cands, `[[`, numeric(1), "roughness")
grid$convergence <- vapply(cands, `[[`, integer(1), "convergence")
grid$at_bounds <- vapply(cands, `[[`, numeric(1), "at_bounds")
best <- which.max(grid$valid_index)
sel <- cands[[best]]
write.csv(grid, "results/bfi-grid.csv", row.names = FALSE)
gt <- grid; names(gt) <- c("Scoring", "gamma", "Train index", "Valid. index",
                           "Roughness", "Conv.", "On bounds")
gt$gamma <- ifelse(gt$Scoring == "linear", NA, gt$gamma)
write_tex_table(gt, "results/tab-bfi-grid.tex", digits = 4)
capture_to("results/bfi-selected.txt", print(sel))

# Final test: selected map versus simpler and external comparators.
sel$scoring <- paste0("ordpp ", sel$scoring, " (selected)")
lin <- cands[[1]]; lin$scoring <- "ordpp linear"
keyed <- as_ordpp_map(xtr, equal_scores(xtr), rep(1, length(items)), label = "Unit-weighted sum")
maps <- list(sel, lin, keyed, map_pca(xtr), safe_map(map_polychoric, xtr),
             safe_map(map_gifi, xtr), safe_map(map_ordpens, xtr))
maps <- maps[!vapply(maps, is.null, logical(1))]
message("held-out tests ...")
tt <- do.call(rbind, lapply(maps, evaluate_map, xtest = xte, B = B_HELD))
tt$abs_latent_cor <- NULL
tt$cor_keyed <- vapply(maps, function(m) stats::cor(predict(m, xte), predict(keyed, xte)), numeric(1))
tt$p_value <- fmt_p(tt$p_value)
write.csv(tt, "results/bfi-test.csv", row.names = FALSE)
names(tt) <- c("Map", "Train index", "Test index", "Test shape", "Var ratio",
               "Held-out p", "cor(sum)")
write_tex_table(tt, "results/tab-bfi-test.tex", digits = 3)

# Stability of the selected map on training rows.
message("stability bootstrap ...")
stab <- stability_ordpp(cands[[best]], xtr, B = B_STAB, seed = 7, trace = TRUE)
capture_to("results/bfi-stability.txt", print(stab))

# Figures: selected score curves with bootstrap bands; held-out projections.
pdf("figures/fig-bfi.pdf", width = 8, height = 3.4)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
cols <- seq_along(items)
plot(NA, xlim = c(1, 6), ylim = range(unlist(sel$scores), unlist(stab$scores), na.rm = TRUE),
     xlab = "Response category", ylab = "Standardized score", main = "(a) Selected scores")
for (j in cols) {
  sj <- stab$scores[[j]]
  band <- apply(sj, 2, stats::quantile, c(.05, .95), na.rm = TRUE)
  if (all(is.finite(band))) polygon(c(1:6, 6:1), c(band[1, ], rev(band[2, ])),
          col = adjustcolor(j, .15), border = NA)
  lines(1:6, sel$scores[[j]], type = "b", col = j, pch = 19)
}
legend("topleft", items, col = cols, lty = 1, pch = 19, bty = "n", cex = .8)
zs <- predict(sel, xte); zl <- predict(lin, xte)
plot(stats::density(zs), main = "(b) Held-out projections", xlab = "Projected score",
     ylim = c(0, max(stats::density(zs)$y, stats::density(zl)$y)))
lines(stats::density(zl), lty = 2)
legend("topright", c("Selected", "Linear"), lty = 1:2, bty = "n", cex = .8)
dev.off()

sel_row <- tt[1, ]
write_macros(list(
  bfiNraw = n_raw, bfiN = n, bfiDropped = n_raw - n,
  bfiNtrain = nrow(xtr), bfiNvalid = nrow(xva), bfiNtest = nrow(xte),
  bfiSelScoring = grid$scoring[best],
  bfiSelGamma = if (grid$scoring[best] == "linear") "--" else format(grid$gamma[best]),
  bfiSelTest = fmt(as.numeric(sel_row[["Test index"]]), 4),
  bfiSelP = sel_row[["Held-out p"]],
  bfiSelVarRatio = fmt(as.numeric(sel_row[["Var ratio"]]), 2),
  bfiCorKeyed = fmt(as.numeric(sel_row[["cor(sum)"]]), 3),
  bfiStabCong = fmt(stats::median(stab$congruence), 3),
  bfiStabFail = stab$failures, bfiStabB = B_STAB, bfiHeldB = B_HELD),
  "results/macros-bfi.tex")
saveRDS(list(grid = grid, cands = cands, best = best, stab = stab, test = tt),
        "results/bfi.rds")
message("bfi application done")
