# Simulation benchmark. Usage: Rscript scripts/03_benchmark.R [reps] [cores] [maxit]
# Default 50 replications per cell. Parallel over replications on Unix via
# parallel::mclapply (cores = 1 on Windows).
source("scripts/common.R")
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args) >= 1) as.integer(args[1]) else budget(50, 2)
CORES <- if (length(args) >= 2) as.integer(args[2]) else 1L
NSTART <- budget(3, 1)
MAXIT <- if (length(args) >= 3) as.integer(args[3]) else budget(500, 20)
B_HELD <- budget(199, 19)
scenarios <- c("independent", "mixture", "rare", "two_factor", "response_style")
ns <- budget(c(300, 800), 200)
cells <- expand.grid(rep = seq_len(REPS), n = ns, scenario = scenarios,
                     stringsAsFactors = FALSE)

one_cell <- function(i) {
  cl <- cells[i, ]
  if (i %% 10 == 0) message(sprintf("cell %d of %d (%s)", i, nrow(cells), format(Sys.time(), "%H:%M")))
  seed <- 1e5 + i
  sim <- simulate_ordpp(cl$n, p = 8, separation = 1, seed = seed, scenario = cl$scenario)
  set.seed(seed); tr <- sample.int(cl$n, cl$n %/% 2); te <- setdiff(seq_len(cl$n), tr)
  xtr <- sim$x[tr, ]; xte <- sim$x[te, ]
  lat <- if (cl$scenario == "independent") NULL else sim$latent[te]
  rows <- list()
  add <- function(map, secs) {
    # Permutation seed kept separate from the data-generating seed.
    r <- tryCatch(evaluate_map(map, xte, B = B_HELD, latent = lat, seed = seed + 500000L),
                  error = function(e) NULL)   # e.g. unseen category in test rows
    if (is.null(r)) return(invisible())
    r$seconds <- secs
    r$convergence <- if (is.na(map$convergence)) NA_integer_ else map$convergence
    rows[[length(rows) + 1]] <<- r
  }
  for (sc in c("linear", "common", "item")) {
    tm <- system.time(f <- tryCatch(suppressWarnings(
      ordpp(xtr, scoring = sc, nstart = NSTART, maxit = MAXIT, seed = seed)),
      error = function(e) NULL))["elapsed"]
    if (!is.null(f)) { f$scoring <- paste("ordpp", sc); add(f, tm) }
  }
  comps <- list(pca = map_pca, polychoric = map_polychoric, gifi = map_gifi,
                ordpens = map_ordpens)
  for (nm in names(comps)) {
    tm <- system.time(m <- safe_map(comps[[nm]], xtr, label = paste(nm, cl$scenario)))["elapsed"]
    if (!is.null(m)) add(m, tm)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  cbind(scenario = cl$scenario, n = cl$n, rep = cl$rep, out)
}

res <- if (CORES > 1) {
  parallel::mclapply(seq_len(nrow(cells)), one_cell, mc.cores = CORES)
} else {
  lapply(seq_len(nrow(cells)), one_cell)
}
lost <- sum(!vapply(res, is.data.frame, logical(1)))  # e.g. unseen test category
raw <- do.call(rbind, res[vapply(res, is.data.frame, logical(1))])
write.csv(raw, "results/benchmark-raw.csv", row.names = FALSE)

# Summaries with Monte Carlo standard errors.
agg <- function(v) c(mean = mean(v, na.rm = TRUE),
                     mcse = stats::sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))
key <- interaction(raw$scenario, raw$n, raw$method, drop = TRUE, sep = "|")
summ <- do.call(rbind, lapply(split(raw, key), function(d) data.frame(
  scenario = d$scenario[1], n = d$n[1], method = d$method[1], runs = nrow(d),
  reject = mean(d$p_value <= 0.05),
  reject_mcse = sqrt(mean(d$p_value <= 0.05) * (1 - mean(d$p_value <= 0.05)) / nrow(d)),
  test_index = agg(d$test_index)[1], test_index_mcse = agg(d$test_index)[2],
  var_ratio = mean(d$variance_ratio), abs_cor = agg(d$abs_latent_cor)[1],
  abs_cor_mcse = agg(d$abs_latent_cor)[2], seconds = stats::median(d$seconds),
  nonconverged = sum(d$convergence != 0, na.rm = TRUE))))
write.csv(summ, "results/benchmark-summary.csv", row.names = FALSE)

# Wide tables for the largest n: held-out index and latent correlation.
nmax <- max(summ$n); s <- summ[summ$n == nmax, ]
wide <- function(col) {
  w <- reshape(s[, c("scenario", "method", col)], idvar = "method",
               timevar = "scenario", direction = "wide")
  names(w) <- sub(paste0("^", col, "\\."), "", names(w))
  w[, c("method", intersect(scenarios, names(w)))]
}
ix <- wide("test_index"); ix$independent <- NULL; names(ix)[1] <- "Method"
names(ix) <- gsub("_", " ", names(ix))
write_tex_table(ix, "results/tab-bench-index.tex", digits = 3)
co <- wide("abs_cor"); co$independent <- NULL; names(co)[1] <- "Method"
names(co) <- gsub("_", " ", names(co))
write_tex_table(co, "results/tab-bench-cor.tex", digits = 3)

dep <- s[s$scenario != "independent", ]
mean_ix <- function(m) mean(dep$test_index[dep$method == m])
it <- summ[summ$method == "ordpp item", ]
write_macros(list(benchReps = REPS, benchNmax = nmax, benchB = B_HELD,
                  benchMaxit = MAXIT, benchLost = lost, benchCells = nrow(cells),
                  benchItemNonconv = sum(it$nonconverged), benchItemRuns = sum(it$runs),
                  benchGainCommon = fmt(100 * (mean_ix("ordpp common") / mean_ix("ordpp linear") - 1), 0),
                  benchGainItem = fmt(100 * (mean_ix("ordpp item") / mean_ix("ordpp linear") - 1), 0)),
             "results/macros-benchmark.tex")
message("benchmark done: ", nrow(raw), " method evaluations; ", lost, " replications lost")
