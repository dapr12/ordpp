# Runtime scaling of a single fit (one start, fixed iteration budget).
source("scripts/common.R")
ns <- budget(c(250, 500, 1000, 2000, 4000), c(100, 200))
ps <- budget(c(5, 10, 20), c(4, 6))
MAXIT <- budget(100, 10)
grid <- expand.grid(n = ns, p = ps, scoring = c("linear", "item"),
                    stringsAsFactors = FALSE)
grid$seconds <- NA_real_; grid$iterations_cap <- MAXIT
for (i in seq_len(nrow(grid))) {
  sim <- simulate_ordpp(grid$n[i], grid$p[i], seed = i)
  grid$seconds[i] <- system.time(suppressWarnings(
    ordpp(sim$x, scoring = grid$scoring[i], nstart = 1, maxit = MAXIT, seed = i)))["elapsed"]
  message(sprintf("n=%d p=%d %s: %.1fs", grid$n[i], grid$p[i], grid$scoring[i], grid$seconds[i]))
}
grid$seconds <- pmax(grid$seconds, 0.01)  # log axis needs positive times
write.csv(grid, "results/runtime.csv", row.names = FALSE)
source("scripts/04b_runtime_plot.R")
big <- grid[grid$n == max(ns) & grid$p == max(ps) & grid$scoring == "item", "seconds"]
write_macros(list(rtMaxN = max(ns), rtMaxP = max(ps), rtMaxSec = fmt(big, 0),
                  rtMaxit = MAXIT, rtCPU = tex_escape(R.version$platform)),
             "results/macros-runtime.tex")
