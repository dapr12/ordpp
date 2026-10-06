ordpp <- function(x, weights = NULL, scoring = c("item", "common", "linear"),
                  lambda = 0, gamma = 0.01, nstart = 3L, maxit = 500L,
                  frequencies = seq(0.2, 3, length.out = 15),
                  min_gap = 0.03, variance_floor = 0.05,
                  index = c("full", "shape"), seed = NULL) {
  cl <- match.call()
  scoring <- match.arg(scoring)
  index <- match.arg(index)
  rng <- .op_set_seed(seed)
  on.exit(.op_restore_seed(rng), add = TRUE)
  if (lambda < 0 || gamma < 0 || min_gap <= 0 || min_gap > 1 ||
      variance_floor <= 0 || nstart < 1 || maxit < 1)
    stop("invalid tuning parameters")
  if (!length(frequencies) || any(!is.finite(frequencies)) || any(frequencies <= 0))
    stop("frequencies must be positive and finite")
  en <- .op_encode(x); ix <- en$x; n <- nrow(ix); d <- ncol(ix)
  if (n < 2) stop("fitting requires at least two rows")
  if (d < 2) stop("fitting requires at least two items")
  w <- .op_weights(weights, n)
  K <- lengths(en$levels)
  if (any(K < 2)) stop("all items must have at least two categories")
  q <- .op_marginals(ix, K, w)
  if (any(unlist(q) <= 0)) stop("every declared category must occur in training data")
  if (scoring == "common" && length(unique(K)) != 1)
    stop("common scoring requires identical category counts and meanings")
  ng <- if (scoring == "linear") integer(0) else
    if (scoring == "common") K[1] - 2L else K - 2L
  decode <- function(par) {
    a <- par[seq_len(d)]
    if (sum(a^2) < 1e-10) return(NULL)
    a <- a / sqrt(sum(a^2)); pos <- d
    raw <- list()
    if (scoring != "linear") for (r in seq_along(ng)) {
      h <- c(0, if (ng[r] > 0) par[pos + seq_len(ng[r])] else numeric(0))
      pos <- pos + ng[r]
      g <- exp(h - max(h)); g <- min_gap + (1 - min_gap) * g / mean(g)
      raw[[r]] <- c(0, cumsum(g))
    }
    if (scoring == "linear") raw <- lapply(K, function(k) seq_len(k) - 1)
    if (scoring == "common") raw <- rep(raw, d)
    scores <- lapply(seq_len(d), function(j) .op_standardize(raw[[j]], q[[j]]))
    y <- .op_map(ix, scores)
    rough <- mean(vapply(raw, function(s) {
      g <- diff(s); mean((g / mean(g) - 1)^2)
    }, numeric(1)))
    list(a = a, scores = scores, y = y, rough = rough,
         variance = sum(w * as.vector(y %*% a)^2))
  }
  objective <- function(par) {
    z <- decode(par)
    if (is.null(z)) return(1e6)
    discrepancy <- ordpp_index(z$y, z$a, w, frequencies, index)
    -discrepancy + lambda * sum(abs(z$a)) + gamma * z$rough +
      100 * max(0, variance_floor - z$variance)^2
  }
  equal <- lapply(seq_len(d), function(j) .op_standardize(seq_len(K[j]), q[[j]]))
  yy <- .op_map(ix, equal)
  pc <- prcomp(sweep(yy, 1, sqrt(w), "*"), center = FALSE)$rotation[, 1]
  fits <- vector("list", nstart)
  for (s in seq_len(nstart)) {
    a0 <- if (s == 1) pc else rnorm(d)
    par <- c(a0 / sqrt(sum(a0^2)), rep(0, sum(ng)))
    fits[[s]] <- optim(par, objective, method = "L-BFGS-B",
      lower = c(rep(-1, d), rep(-3, sum(ng))),
      upper = c(rep(1, d), rep(3, sum(ng))), control = list(maxit = maxit))
  }
  best <- which.min(vapply(fits, function(f) f$value, numeric(1)))
  z <- decode(fits[[best]]$par)
  if (z$a[which.max(abs(z$a))] < 0) z$a <- -z$a
  names(z$a) <- en$names; names(z$scores) <- en$names
  for (j in seq_len(d)) names(z$scores[[j]]) <- as.character(en$levels[[j]])
  gap_par <- fits[[best]]$par[-seq_len(d)]
  out <- list(direction = z$a, scores = z$scores, levels = en$levels,
    variable_names = en$names, projection = as.vector(z$y %*% z$a),
    transformed = z$y, weights = w, weights_supplied = !is.null(weights),
    index = ordpp_index(z$y, z$a, w, frequencies, index),
    objective = -fits[[best]]$value, roughness = z$rough,
    projection_variance = z$variance, n = n, scoring = scoring,
    convergence = fits[[best]]$convergence,
    at_bounds = sum(abs(gap_par) >= 3 - 1e-6),
    starts = data.frame(
      objective = -vapply(fits, function(f) f$value, numeric(1)),
      convergence = vapply(fits, function(f) f$convergence, integer(1))),
    settings = list(lambda = lambda, gamma = gamma, frequencies = frequencies,
                    min_gap = min_gap, variance_floor = variance_floor,
                    index = index, nstart = nstart, maxit = maxit),
    call = cl)
  class(out) <- "ordpp"
  if (out$convergence != 0)
    warning("best optimization did not converge; increase maxit/nstart")
  if (z$variance < variance_floor)
    warning("projection variance below requested floor")
  if (out$at_bounds > 0)
    warning(out$at_bounds, " gap parameter(s) on their bounds; report and check sensitivity")
  out
}

# Wrap scores and loadings produced elsewhere (e.g. optimal-scaling PCA) as a
# frozen map so they can be evaluated with the same held-out tools.
as_ordpp_map <- function(x, scores, direction, weights = NULL,
                         frequencies = seq(0.2, 3, length.out = 15),
                         index = c("full", "shape"), label = "external") {
  index <- match.arg(index)
  en <- .op_encode(x); d <- ncol(en$x)
  w <- .op_weights(weights, nrow(en$x))
  if (length(scores) != d || length(direction) != d)
    stop("scores and direction need one entry per item")
  K <- lengths(en$levels)
  if (any(lengths(scores) != K))
    stop("each score vector needs one value per observed category, in order")
  q <- .op_marginals(en$x, K, w)
  sc <- lapply(seq_len(d), function(j) .op_standardize(as.numeric(scores[[j]]), q[[j]]))
  a <- as.numeric(direction)
  if (any(!is.finite(a)) || sum(a^2) < 1e-12) stop("invalid direction")
  a <- a / sqrt(sum(a^2))
  if (a[which.max(abs(a))] < 0) a <- -a
  names(a) <- en$names; names(sc) <- en$names
  for (j in seq_len(d)) names(sc[[j]]) <- as.character(en$levels[[j]])
  y <- .op_map(en$x, sc)
  out <- list(direction = a, scores = sc, levels = en$levels,
    variable_names = en$names, projection = as.vector(y %*% a),
    transformed = y, weights = w, weights_supplied = !is.null(weights),
    index = ordpp_index(y, a, w, frequencies, index),
    objective = NA_real_, roughness = NA_real_,
    projection_variance = sum(w * as.vector(y %*% a)^2), n = nrow(y),
    scoring = label, convergence = NA_integer_, at_bounds = NULL,
    starts = NULL, settings = list(frequencies = frequencies, index = index),
    call = match.call())
  class(out) <- "ordpp"
  out
}
