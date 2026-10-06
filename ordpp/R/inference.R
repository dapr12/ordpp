# Exact Monte Carlo permutation test on held-out data along a frozen map.
# Under iid rows and item independence, independently permuting each column
# leaves the joint law of the held-out sample unchanged; because the map was
# estimated on other observations, no refitting is needed.
heldout_test_ordpp <- function(object, newdata, B = 999L, seed = NULL, type = NULL) {
  if (B < 1 || B != as.integer(B)) stop("B must be a positive integer")
  y <- .op_newy(object, newdata)
  if (nrow(y) < 2) stop("at least two held-out rows are required")
  if (is.null(type)) type <- .op_type(object)
  f <- object$settings$frequencies; a <- object$direction
  stat <- ordpp_index(y, a, NULL, f, type)
  rng <- .op_set_seed(seed)
  on.exit(.op_restore_seed(rng), add = TRUE)
  n <- nrow(y)
  null <- vapply(seq_len(B), function(b) {
    yb <- y
    for (j in seq_len(ncol(y))) yb[, j] <- y[sample.int(n), j]
    ordpp_index(yb, a, NULL, f, type)
  }, numeric(1))
  structure(list(statistic = stat, null = null,
                 p.value = (1 + sum(null >= stat)) / (B + 1), B = B, n = n,
                 type = type,
                 method = "Frozen-map held-out permutation test",
                 hypothesis = "iid held-out rows; mutual item independence"),
            class = "ordpp_test")
}

# Refitted permutation test: repeats the complete fixed-setting search.
permutation_ordpp <- function(x, B = 99L, seed = NULL, fit = NULL, trace = FALSE, ...) {
  if (B < 1 || B != as.integer(B)) stop("B must be a positive integer")
  dots <- list(...)
  if (any(c("weights", "seed") %in% names(dots)))
    stop("iid unweighted test only; use design-specific resampling for surveys")
  x <- as.data.frame(x)
  if (!is.null(fit)) {
    if (!inherits(fit, "ordpp")) stop("fit must be an ordpp object")
    if (length(dots)) stop("supply either fit or fitting arguments, not both")
    if (isTRUE(fit$weights_supplied))
      stop("fit used weights; this iid permutation test does not apply")
    dots <- .op_refit_args(fit)
    x <- x[, fit$variable_names, drop = FALSE]
  }
  rng <- .op_set_seed(seed)
  on.exit(.op_restore_seed(rng), add = TRUE)
  perms <- lapply(seq_len(B), function(b) {
    xb <- x
    xb[] <- lapply(x, function(v) v[sample.int(length(v))])
    xb
  })
  seeds <- sample.int(.Machine$integer.max, B + 1)
  obs <- suppressWarnings(do.call(ordpp, c(list(x = x, seed = seeds[1]), dots)))
  null <- numeric(B); conv <- integer(B)
  for (b in seq_len(B)) {
    fb <- suppressWarnings(do.call(ordpp,
            c(list(x = perms[[b]], seed = seeds[b + 1]), dots)))
    null[b] <- fb$objective; conv[b] <- fb$convergence
    if (trace && (b %% 10 == 0 || b == B)) message("permutation ", b, " of ", B)
  }
  structure(list(statistic = obs$objective, null = null,
                 p.value = (1 + sum(null >= obs$objective)) / (B + 1), B = B,
                 n = nrow(x), type = .op_type(obs),
                 null_convergence = conv, observed_fit = obs,
                 method = "Refitted permutation test (entire search repeated)",
                 hypothesis = "iid rows; mutual item independence"),
            class = "ordpp_test")
}

# Bootstrap stability of loadings and category scores.
stability_ordpp <- function(fit, x, B = 50L, seed = NULL, weights = NULL,
                            cluster = NULL, trace = FALSE) {
  if (!inherits(fit, "ordpp")) stop("fit must be an ordpp object")
  if (B < 1 || B != as.integer(B)) stop("B must be a positive integer")
  args <- .op_refit_args(fit)
  x <- as.data.frame(x)
  miss <- setdiff(fit$variable_names, names(x))
  if (length(miss)) stop("x lacks item(s): ", paste(miss, collapse = ", "))
  x <- x[, fit$variable_names, drop = FALSE]
  # Freeze the fitted category dictionary so every refit is aligned.
  for (j in seq_along(x)) x[[j]] <- factor(as.character(x[[j]]),
    levels = as.character(fit$levels[[j]]), ordered = TRUE)
  if (anyNA(x)) stop("x has missing values or categories absent from the fit")
  n <- nrow(x)
  if (!is.null(weights) && length(weights) != n) stop("weights must match nrow(x)")
  if (!is.null(cluster) && length(cluster) != n) stop("cluster must match nrow(x)")
  rng <- .op_set_seed(seed)
  on.exit(.op_restore_seed(rng), add = TRUE)
  resamples <- lapply(seq_len(B), function(b) {
    if (is.null(cluster)) return(sample.int(n, n, replace = TRUE))
    ids <- unique(cluster)
    pick <- ids[sample.int(length(ids), length(ids), replace = TRUE)]
    unlist(lapply(pick, function(g) which(cluster == g)), use.names = FALSE)
  })
  seeds <- sample.int(.Machine$integer.max, B)
  p <- length(fit$direction)
  dirs <- matrix(NA_real_, B, p, dimnames = list(NULL, fit$variable_names))
  scores <- lapply(fit$scores, function(s)
    matrix(NA_real_, B, length(s), dimnames = list(NULL, names(s))))
  index <- rep(NA_real_, B); conv <- rep(NA_integer_, B)
  errors <- rep(NA_character_, B)
  for (b in seq_len(B)) {
    rows <- resamples[[b]]
    wb <- if (is.null(weights)) NULL else weights[rows]
    fb <- tryCatch(suppressWarnings(do.call(ordpp,
            c(list(x = x[rows, , drop = FALSE], weights = wb, seed = seeds[b]),
              args))), error = function(e) e)
    if (trace && (b %% 5 == 0 || b == B)) message("bootstrap refit ", b, " of ", B)
    if (inherits(fb, "error")) { errors[b] <- conditionMessage(fb); next }
    a <- fb$direction
    if (sum(a * fit$direction) < 0) a <- -a   # sign alignment
    dirs[b, ] <- a
    for (j in seq_len(p)) scores[[j]][b, ] <- fb$scores[[j]]
    index[b] <- fb$index; conv[b] <- fb$convergence
  }
  ok <- !is.na(index)
  dk <- dirs[ok, , drop = FALSE]
  summary <- data.frame(item = fit$variable_names, estimate = fit$direction,
    mean = colMeans(dk), sd = apply(dk, 2, stats::sd),
    lower = apply(dk, 2, stats::quantile, 0.025, na.rm = TRUE, names = FALSE),
    upper = apply(dk, 2, stats::quantile, 0.975, na.rm = TRUE, names = FALSE),
    row.names = NULL)
  structure(list(summary = summary, directions = dirs, scores = scores,
                 index = index, convergence = conv,
                 congruence = abs(as.vector(dk %*% fit$direction)),
                 failures = sum(!ok), errors = errors[!ok], B = B,
                 clustered = !is.null(cluster)),
            class = "ordpp_stability")
}
