# Worked synthetic example (Section "A worked synthetic example").
# The code blocks shown in the article are reproduced verbatim below.
source("scripts/common.R")
NSTART <- budget(5, 2); MAXIT <- budget(300, 40)
B_HELD <- budget(999, 49); B_PERM <- budget(99, 9)

## ---- article block: fit ----
sim <- simulate_ordpp(n = 600, p = 6, seed = 103)
set.seed(104)
train <- sample.int(nrow(sim$x), 300)
test <- setdiff(seq_len(nrow(sim$x)), train)
fit <- ordpp(sim$x[train, ], scoring = "item", gamma = 0.01,
             nstart = NSTART, maxit = MAXIT, seed = 105)
capture_to("results/syn-fit.txt", { print(fit); print(fit$starts) })

## ---- article block: held-out ----
ht <- heldout_test_ordpp(fit, sim$x[test, ], B = B_HELD, seed = 106)
dc <- decompose_ordpp(fit, sim$x[test, ])
capture_to("results/syn-heldout.txt", { print(ht); print(dc, digits = 3) })

## ---- article block: refitted permutation test ----
fit_lin <- ordpp(sim$x[train, ], scoring = "linear", nstart = 3,
                 maxit = budget(150, 30), seed = 107)
perm <- permutation_ordpp(sim$x[train, ], B = B_PERM, seed = 108, fit = fit_lin,
                          trace = TRUE)
capture_to("results/syn-perm.txt", print(perm))

## ---- ablation table (same split, same budget) ----
fits <- list(
  linear = fit_lin,
  common = ordpp(sim$x[train, ], scoring = "common", nstart = NSTART, maxit = MAXIT, seed = 105),
  item = fit,
  item_shape = ordpp(sim$x[train, ], scoring = "item", index = "shape",
                     nstart = NSTART, maxit = MAXIT, seed = 105))
fits$item_shape$scoring <- "item (shape index)"
maps <- c(fits, list(pca = map_pca(sim$x[train, ]),
                     poly = safe_map(map_polychoric, sim$x[train, ]),
                     gifi = safe_map(map_gifi, sim$x[train, ]),
                     ordpens = safe_map(map_ordpens, sim$x[train, ])))
maps <- maps[!vapply(maps, is.null, logical(1))]
ab <- do.call(rbind, lapply(maps, evaluate_map, xtest = sim$x[test, ],
                            B = B_HELD, latent = sim$latent[test]))
ab$p_value <- fmt_p(ab$p_value)
write.csv(ab, "results/syn-ablation.csv", row.names = FALSE)
names(ab) <- c("Method", "Train index", "Test index", "Test shape",
               "Var ratio", "Held-out p", "|cor| latent")
write_tex_table(ab, "results/tab-syn-ablation.tex", digits = 3)

## ---- figures ----
pdf("figures/fig-syn.pdf", width = 8, height = 3.4)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(fit, type = "scores", main = "(a) Learned category scores")
z_test <- predict(fit, sim$x[test, ])
plot(stats::density(z_test), main = "(b) Held-out projection", xlab = "Projected score")
graphics::rug(z_test)
dev.off()

write_macros(list(
  synTestIndex = fmt(dc$full, 4), synShape = fmt(dc$shape, 4),
  synVarRatio = fmt(dc$variance_ratio, 2), synHeldP = fmt_p(ht$p.value),
  synHeldB = B_HELD, synPermP = fmt_p(perm$p.value), synPermB = B_PERM,
  synLatentCor = fmt(abs(stats::cor(z_test, sim$latent[test])), 2),
  synConv = fit$convergence, synBounds = fit$at_bounds),
  "results/macros-synthetic.tex")
saveRDS(list(fits = fits, maps = maps), "results/synthetic.rds")
message("synthetic example done")
