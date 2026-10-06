# Size and power of the held-out permutation test.
# Usage: Rscript scripts/08_size_power.R [cores]
source("scripts/common.R")
args <- commandArgs(trailingOnly = TRUE)
CORES <- if (length(args) >= 1) as.integer(args[1]) else 1L
REPS_SIZE <- budget(1000, 20); REPS_POW <- budget(200, 5)
B <- budget(199, 19); MAXIT <- budget(500, 20)
strengths <- c(0.1, 0.2, 0.3)

fitters <- list(
  "ordpp linear" = function(x, s) suppressWarnings(ordpp(x, scoring = "linear", maxit = MAXIT, seed = s)),
  "ordpp common" = function(x, s) suppressWarnings(ordpp(x, scoring = "common", maxit = MAXIT, seed = s)),
  "ordpp item"   = function(x, s) suppressWarnings(ordpp(x, scoring = "item", maxit = MAXIT, seed = s)),
  "PCA, equal scores" = function(x, s) map_pca(x),
  "Polychoric PCA" = function(x, s) map_polychoric(x),
  "Gifi princals" = function(x, s) map_gifi(x))
cheap <- c("ordpp linear", "PCA, equal scores", "Polychoric PCA")

one_rep <- function(job) {
  sim <- simulate_ordpp(job$n, p = 8, separation = 1, strength = job$strength,
                        seed = job$seed, scenario = "mixture")
  set.seed(job$seed + 1L)
  tr <- sample.int(job$n, job$n %/% 2); te <- setdiff(seq_len(job$n), tr)
  out <- lapply(job$methods, function(m) {
    map <- tryCatch(fitters[[m]](sim$x[tr, ], job$seed), error = function(e) NULL)
    if (is.null(map)) return(NULL)
    ht <- tryCatch(heldout_test_ordpp(map, sim$x[te, ], B = B,
                                      seed = job$seed + 500000L),   # separate stream
                   error = function(e) NULL)
    if (is.null(ht)) return(NULL)
    data.frame(method = m, p_value = ht$p.value)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) return(NULL)
  cbind(study = job$study, n = job$n, strength = job$strength, out)
}

jobs <- list(); k <- 0
for (n in budget(c(300, 800), 200)) for (r in seq_len(REPS_SIZE)) {
  k <- k + 1
  jobs[[k]] <- list(study = "size", n = n, strength = 0, seed = 2e6 + k,
                    methods = c(cheap, if (r <= REPS_POW) c("ordpp item", "Gifi princals")))
}
for (st in strengths) for (r in seq_len(REPS_POW)) {
  k <- k + 1
  jobs[[k]] <- list(study = "power", n = 300, strength = st, seed = 2e6 + k,
                    methods = names(fitters))
}
message(length(jobs), " jobs")
run <- function(i) {
  if (i %% 100 == 0) message(sprintf("job %d of %d (%s)", i, length(jobs), format(Sys.time(), "%H:%M")))
  one_rep(jobs[[i]])
}
res <- if (CORES > 1) {
  parallel::mclapply(seq_along(jobs), run, mc.cores = CORES)
} else {
  lapply(seq_along(jobs), run)
}
raw <- do.call(rbind, res[vapply(res, is.data.frame, logical(1))])
write.csv(raw, "results/size-power-raw.csv", row.names = FALSE)

rate <- function(d) data.frame(rate = mean(d$p_value <= 0.05), runs = nrow(d),
                               mcse = sqrt(mean(d$p_value <= 0.05) * (1 - mean(d$p_value <= 0.05)) / nrow(d)))
summ <- do.call(rbind, lapply(split(raw, interaction(raw$study, raw$n, raw$strength, raw$method, drop = TRUE)),
  function(d) cbind(d[1, c("study", "n", "strength", "method")], rate(d))))
write.csv(summ, "results/size-power-summary.csv", row.names = FALSE)

sz <- summ[summ$study == "size", ]
szw <- reshape(sz[, c("method", "n", "rate")], idvar = "method", timevar = "n", direction = "wide")
names(szw) <- c("Method", paste0("n = ", sub("rate\\.", "", names(szw)[-1])))
runs_per <- tapply(sz$runs, sz$method, max)
szw$Replications <- as.numeric(runs_per[szw$Method])
write_tex_table(szw, "results/tab-size.tex", digits = 3)

pw <- summ[summ$study == "power", ]
pww <- reshape(pw[, c("method", "strength", "rate")], idvar = "method", timevar = "strength", direction = "wide")
names(pww) <- c("Method", paste0("c = ", sub("rate\\.", "", names(pww)[-1])))
write_tex_table(pww, "results/tab-power.tex", digits = 2)

cheap_sz <- sz[sz$method %in% cheap, ]
write_macros(list(sizeReps = REPS_SIZE, sizeItemReps = REPS_POW, powReps = REPS_POW, spB = B,
                  sizeMin = fmt(min(sz$rate), 3), sizeMax = fmt(max(sz$rate), 3),
                  sizeMcse = fmt(max(cheap_sz$mcse), 3)),
             "results/macros-sizepower.tex")
message("size/power done")
