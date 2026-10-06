# Figures 1-2: explanatory (not fitted) illustrations, produced in R.
source("scripts/common.R")

# Figure 1: equal versus unequal spacing under uniform category probabilities.
std <- function(s) { s <- s - mean(s); s / sqrt(mean(s^2)) }
gaps_u <- c(2.15, 0.23, 0.23, 0.23, 2.15)
eq <- std(1:6); uq <- std(c(0, cumsum(gaps_u)))
pdf("figures/fig-scoring.pdf", width = 8, height = 3.2)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
matplot(1:6, cbind(eq, uq), type = "b", pch = 19, lty = 1,
        col = c("#2c5f8a", "#c0622b"), xlab = "Ordered response category",
        ylab = "Standardized score", main = "(a) Same order, different spacing")
legend("topleft", c("Equal spacing", "Unequal spacing"), bty = "n",
       col = c("#2c5f8a", "#c0622b"), lty = 1, pch = 19, cex = .8)
g <- rbind(rep(1, 5), gaps_u / mean(gaps_u))
barplot(g, beside = TRUE, col = c("#2c5f8a", "#c0622b"),
        names.arg = c("1-2", "2-3", "3-4", "4-5", "5-6"),
        xlab = "Adjacent categories", ylab = "Gap / mean gap",
        main = "(b) Gaps penalized by R(s)")
abline(h = 1, lty = 3)
dev.off()

# Figure 2: an exact dependent table, its product-marginal reference and the
# induced laws of the raw sum.
tab_dep <- matrix(c(.4, .1, .1, .4), 2)
tab_ind <- matrix(.25, 2, 2)
draw_tab <- function(tb, main) {
  plot(NA, xlim = c(0, 2), ylim = c(0, 2), axes = FALSE, asp = 1,
       xlab = "Item 2", ylab = "Item 1", main = main)
  for (i in 1:2) for (j in 1:2) {
    v <- tb[i, j]
    rect(j - 1, 2 - i, j, 3 - i, col = gray(1 - v * 1.6), border = "white")
    text(j - .5, 2.5 - i, sprintf("%.2f", v), col = if (v > .3) "white" else "black")
  }
  axis(1, at = c(.5, 1.5), labels = 0:1, tick = FALSE)
  axis(2, at = c(1.5, .5), labels = 0:1, tick = FALSE, las = 1)
}
pdf("figures/fig-reference.pdf", width = 9, height = 3)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
draw_tab(tab_dep, "(a) Dependent joint law")
draw_tab(tab_ind, "(b) Product of marginals")
barplot(rbind(c(.4, .2, .4), c(.25, .5, .25)), beside = TRUE, names.arg = 0:2,
        col = c("#2c5f8a", "#c0622b"), ylim = c(0, .6), xlab = "Raw sum X1 + X2",
        ylab = "Probability", main = "(c) Projected laws")
legend("topright", c("Joint", "Independence"), fill = c("#2c5f8a", "#c0622b"),
       bty = "n", cex = .9)
dev.off()
message("figures/fig-scoring.pdf and figures/fig-reference.pdf written")
