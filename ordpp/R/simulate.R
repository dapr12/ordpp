simulate_ordpp <- function(n = 500, p = 8, separation = 1.5, seed = NULL,
                           strength = 1,
                           scenario = c("mixture", "independent", "rare",
                                        "two_factor", "response_style")) {
  scenario <- match.arg(scenario)
  rng <- .op_set_seed(seed)
  on.exit(.op_restore_seed(rng), add = TRUE)
  if (n < 10 || p < 4 || separation < 0 || strength < 0)
    stop("invalid simulation arguments")
  group <- rbinom(n, 1, .5)
  latent <- rnorm(n, (2 * group - 1) * separation, 1)
  loading <- strength * c(rep(1, 3), -.8, rep(0, p - 4))
  if (scenario == "independent") loading <- rep(0, p)
  y <- outer(latent, loading) + matrix(rnorm(n * p), n, p)
  if (scenario == "two_factor") {
    load2 <- c(rep(0, 4), rep(1, min(3, p - 4)), rep(0, max(0, p - 7)))
    y <- y + outer(rnorm(n), strength * load2)
  }
  if (scenario == "response_style") {
    extreme <- rbinom(n, 1, 0.3) == 1   # extreme responders stretch all items
    y[extreme, ] <- 1.8 * y[extreme, ]
  }
  cuts <- c(-Inf, -1.7, -.4, .15, 1.8, Inf)
  rare_cuts <- c(-Inf, -1.7, -.4, .15, 3.2, Inf)
  x <- vapply(seq_len(p), function(j) {
    cj <- if (scenario == "rare" && loading[j] != 0) rare_cuts else cuts
    as.integer(cut(y[, j], cj, labels = FALSE))
  }, integer(n))
  x <- matrix(x, n, p)
  colnames(x) <- paste0("item", seq_len(p))
  list(x = as.data.frame(x), latent = latent, group = group,
       loading = loading, thresholds = cuts, scenario = scenario)
}
