library(ordpp)
s <- simulate_ordpp(100, 4, seed = 2)
f <- suppressWarnings(ordpp(s$x, scoring = "common", nstart = 1, maxit = 15, seed = 3))
stopifnot(abs(sum(f$direction^2) - 1) < 1e-9,
          all(vapply(f$scores, function(x) all(diff(x) > 0), logical(1))),
          max(abs(predict(f, s$x) - f$projection)) < 1e-10)

# Independence gives an exactly zero index on a balanced Cartesian table,
# for both index types; perfect dependence gives a positive index.
y <- as.matrix(expand.grid(a = c(-1, 1), b = c(-1, 1)))
stopifnot(ordpp_index(y, c(1, 1)) < 1e-20,
          ordpp_index(y, c(1, 1), type = "shape") < 1e-20,
          ordpp_index(cbind(c(-1, 1), c(-1, 1)), c(1, 1)) > .001)

# Single-item direction: joint and reference laws coincide.
stopifnot(ordpp_index(cbind(c(-1, 0, 2), c(1, 0, -1)), c(1, 0)) < 1e-20)

# Weight rescaling leaves the criterion unchanged.
w <- seq_len(nrow(s$x))
f1 <- suppressWarnings(ordpp(s$x, weights = w, scoring = "linear", nstart = 1, maxit = 10, seed = 3))
f2 <- suppressWarnings(ordpp(s$x, weights = 100 * w, scoring = "linear", nstart = 1, maxit = 10, seed = 3))
stopifnot(abs(f1$index - f2$index) < 1e-8)

# Unseen categories are rejected; column order is irrelevant at prediction.
bad <- s$x; bad[1, 1] <- 99
stopifnot(inherits(try(predict(f, bad), silent = TRUE), "try-error"))
stopifnot(max(abs(predict(f, s$x[, 4:1]) - f$projection)) < 1e-10)
stopifnot(abs(predict(f, s$x[1, , drop = FALSE]) - f$projection[1]) < 1e-10)

# Validation on the training rows reproduces the training index.
stopifnot(abs(validate_ordpp(f, s$x) - f$index) < 1e-10)

# A frozen map built from equal scores reproduces a linear fit exactly.
fl <- suppressWarnings(ordpp(s$x, scoring = "linear", nstart = 1, maxit = 10, seed = 3))
m <- as_ordpp_map(s$x, lapply(fl$levels, seq_along), fl$direction)
stopifnot(abs(m$index - fl$index) < 1e-10,
          max(abs(predict(m, s$x) - fl$projection)) < 1e-10)

# Decomposition: variance ratio equals stored projection variance in training.
dc <- decompose_ordpp(fl)
stopifnot(abs(dc$variance_ratio - fl$projection_variance) < 1e-8,
          abs(dc$full - fl$index) < 1e-10, dc$shape >= 0)

# Held-out permutation test returns a valid Monte Carlo p-value.
ht <- heldout_test_ordpp(fl, s$x, B = 19, seed = 1)
stopifnot(ht$p.value > 0, ht$p.value <= 1, length(ht$null) == 19)

# Refitted test reuses the settings of a supplied fit.
pt <- permutation_ordpp(s$x, B = 3, seed = 1, fit = fl)
stopifnot(pt$observed_fit$scoring == "linear", length(pt$null) == 3)

# Stability: one row of summaries per item; congruence in [0, 1].
st <- stability_ordpp(fl, s$x, B = 3, seed = 1)
stopifnot(nrow(st$summary) == 4, all(st$congruence <= 1 + 1e-12))

# Shape-index fit runs and is recorded.
fs <- suppressWarnings(ordpp(s$x, scoring = "linear", index = "shape",
                             nstart = 1, maxit = 10, seed = 3))
stopifnot(fs$settings$index == "shape",
          abs(validate_ordpp(fs, s$x) - fs$index) < 1e-10)

# Simulation scenarios produce complete integer data.
for (sc in c("mixture", "independent", "rare", "two_factor", "response_style")) {
  d <- simulate_ordpp(60, 7, seed = 1, scenario = sc)$x
  stopifnot(nrow(d) == 60, ncol(d) == 7, !anyNA(d))
}

# Tie-collapsed computation equals the direct n-row computation.
naive <- function(y, a, w, f) {
  w <- w / sum(w); yc <- sweep(y, 2, colSums(w * y), "-")
  v0 <- sum(a^2 * colSums(w * yc^2)); t0 <- f / sqrt(v0); z <- as.vector(yc %*% a)
  joint <- colSums(w * exp(1i * outer(z, t0)))
  ref <- rep(1 + 0i, length(f))
  for (j in seq_along(a)) ref <- ref * colSums(w * exp(1i * outer(yc[, j] * a[j], t0)))
  u <- exp(-f^2 / 2); u <- u / sum(u); sum(u * Mod(joint - ref)^2)
}
set.seed(9)
yy <- matrix(sample(1:4, 600, TRUE), 200, 3); ww <- runif(200); aa <- c(.6, -.3, .5)
ff <- seq(0.2, 3, length.out = 15)
stopifnot(abs(ordpp_index(yy, aa, ww, ff) - naive(yy, aa, ww, ff)) < 1e-12)

# Seeded calls leave the caller's random number state untouched.
set.seed(42); before <- .Random.seed
invisible(simulate_ordpp(50, 4, seed = 1))
invisible(heldout_test_ordpp(fl, s$x, B = 9, seed = 2))
stopifnot(identical(before, .Random.seed))
