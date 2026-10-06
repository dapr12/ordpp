# Redraws the runtime figure from results/runtime.csv (no refitting).
grid <- read.csv("results/runtime.csv")
ns <- sort(unique(grid$n)); ps <- sort(unique(grid$p))
pdf("figures/fig-runtime.pdf", width = 5.5, height = 3.8)
par(mar = c(4, 4.5, 1, 7.5), xpd = FALSE)
plot(NA, xlim = range(ns), ylim = range(grid$seconds), log = "xy",
     xlab = "Sample size n", ylab = "Seconds (one start)", xaxt = "n", yaxt = "n")
yt <- c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000)
yt <- yt[yt >= min(grid$seconds) / 1.5 & yt <= max(grid$seconds) * 1.5]
axis(2, at = yt, labels = format(yt, scientific = FALSE, drop0trailing = TRUE, trim = TRUE), las = 1)
axis(1, at = ns, labels = format(ns, big.mark = ",", trim = TRUE))
for (sc in c("linear", "item")) for (p in ps) {
  g <- grid[grid$scoring == sc & grid$p == p, ]
  g <- g[order(g$n), ]
  lines(g$n, g$seconds, type = "b", pch = if (sc == "item") 19 else 1,
        lty = if (sc == "item") 1 else 2, col = match(p, ps))
}
par(xpd = TRUE)
legend("topleft", inset = c(1.02, 0),
       c(paste("p =", ps), "item", "linear"),
       col = c(seq_along(ps), 1, 1), lty = c(rep(1, length(ps)), 1, 2),
       pch = c(rep(NA, length(ps)), 19, 1), bty = "n", cex = .8)
dev.off()
