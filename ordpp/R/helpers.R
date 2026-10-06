# Internal helpers. None of these are exported.

.op_weights <- function(w, n) {
  if (is.null(w)) w <- rep(1, n)
  if (length(w) != n || any(!is.finite(w)) || any(w <= 0))
    stop("weights must be finite, positive and match nrow(x)")
  w / sum(w)
}

.op_encode <- function(x, levels = NULL) {
  x <- as.data.frame(x)
  if (!ncol(x) || nrow(x) < 1 || anyNA(x))
    stop("x must have at least one complete row and one column; recode missing values first")
  if (is.null(names(x))) names(x) <- paste0("V", seq_len(ncol(x)))
  if (is.null(levels)) levels <- lapply(x, function(v) {
    if (is.ordered(v)) levels(v) else {
      if (!is.numeric(v)) stop("use numeric categories or ordered factors")
      if (any(!is.finite(v))) stop("non-finite categories")
      sort(unique(v))
    }
  })
  if (length(levels) != ncol(x)) stop("incompatible category specification")
  ix <- mapply(function(v, lv) match(as.character(v), as.character(lv)),
               x, levels, SIMPLIFY = FALSE)
  if (any(vapply(ix, anyNA, logical(1))))
    stop("new or unrecognised category; prediction uses the training dictionary")
  list(x = matrix(unlist(ix), nrow(x), ncol(x)), levels = levels,
       names = names(x))
}

.op_map <- function(ix, scores) {
  matrix(unlist(lapply(seq_along(scores), function(j) scores[[j]][ix[, j]])),
         nrow(ix), ncol(ix))
}

.op_standardize <- function(s, q) {
  s <- s - sum(q * s)
  v <- sum(q * s^2)
  if (v < 1e-12) stop("constant item in positive-weight training data")
  s / sqrt(v)
}

.op_marginals <- function(ix, K, w) {
  lapply(seq_along(K), function(j)
    vapply(seq_len(K[j]), function(k) sum(w[ix[, j] == k]), numeric(1)))
}

# Transformed held-out matrix using the frozen dictionary of a fitted map.
.op_newy <- function(object, newdata) {
  x <- as.data.frame(newdata)
  miss <- setdiff(object$variable_names, names(x))
  if (length(miss))
    stop("newdata lacks item(s): ", paste(miss, collapse = ", "))
  x <- x[, object$variable_names, drop = FALSE]
  en <- .op_encode(x, object$levels)
  .op_map(en$x, object$scores)
}

.op_type <- function(object) {
  if (is.null(object$settings$index)) "full" else object$settings$index
}

# Var(Z) / Var(Z0) for weighted data; Var(Z0) is the independence variance.
.op_variance_ratio <- function(y, a, w) {
  yc <- sweep(y, 2, colSums(w * y), "-")
  v0 <- sum(a^2 * colSums(w * yc^2))
  vz <- sum(w * as.vector(yc %*% a)^2)
  vz / v0
}

# Settings needed to repeat a fit exactly (used by resampling helpers).
.op_refit_args <- function(fit) {
  s <- fit$settings
  if (is.null(s$nstart))
    stop("this object cannot be refitted (it was not produced by ordpp())")
  args <- list(scoring = fit$scoring, lambda = s$lambda, gamma = s$gamma,
               nstart = s$nstart, maxit = s$maxit,
               frequencies = s$frequencies, min_gap = s$min_gap,
               variance_floor = s$variance_floor, index = s$index)
  args[!vapply(args, is.null, logical(1))]
}

# Distinct values of v with summed weights (exact ties only).
.op_collapse <- function(v, w) {
  u <- unique(v)
  if (length(u) == length(v)) return(list(v = v, w = w))
  g <- factor(match(v, u), levels = seq_along(u))
  list(v = u, w = vapply(split(w, g), sum, numeric(1), USE.NAMES = FALSE))
}

# Seed the random number generator and return the caller's previous state, so
# that functions can restore it on exit (CRAN: do not alter global RNG state).
.op_set_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- list(had = had, state = if (had)
    get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL)
  set.seed(seed)
  old
}

.op_restore_seed <- function(old) {
  if (is.null(old)) return(invisible(NULL))
  if (old$had) assign(".Random.seed", old$state, envir = globalenv())
  else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())
  invisible(NULL)
}
