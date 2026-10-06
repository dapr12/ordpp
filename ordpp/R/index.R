ordpp_index <- function(y, direction, weights = NULL,
                        frequencies = seq(0.2, 3, length.out = 15),
                        type = c("full", "shape")) {
  type <- match.arg(type)
  y <- as.matrix(y); a <- as.numeric(direction)
  if (length(a) != ncol(y) || any(!is.finite(y)) || any(!is.finite(a)) ||
      sum(a^2) < 1e-12 || !length(frequencies) ||
      any(!is.finite(frequencies)) || any(frequencies <= 0))
    stop("invalid transformed data, direction or frequencies")
  w <- .op_weights(weights, nrow(y))
  mu <- colSums(w * y)
  yc <- sweep(y, 2, mu, "-")
  # Frequencies are expressed on the scale of the independence SD, so the
  # reference projection always has unit variance.
  v0 <- sum(a^2 * colSums(w * yc^2))
  if (v0 < 1e-12) return(0)
  z <- as.vector(yc %*% a)
  t0 <- frequencies / sqrt(v0)
  tz <- t0
  if (type == "shape") {
    # Compare Z / sd(Z) with Z0 / sd(Z0): removes the variance (pairwise
    # covariance) component, leaving differences in distributional shape.
    vz <- sum(w * z^2)
    if (vz < 1e-12) return(0)
    tz <- frequencies / sqrt(vz)
  }
  # Ordinal data take few distinct values: collapse ties before evaluating
  # characteristic functions (exactly the same value, much less work).
  zc <- .op_collapse(z, w)
  joint <- colSums(zc$w * exp(1i * outer(zc$v, tz)))
  ref <- rep(1 + 0i, length(frequencies))
  for (j in seq_along(a)) {
    if (a[j] == 0) next   # factor is exactly one
    mj <- .op_collapse(yc[, j], w)
    ref <- ref * colSums(mj$w * exp(1i * outer(mj$v * a[j], t0)))
  }
  u <- exp(-frequencies^2 / 2); u <- u / sum(u)
  sum(u * Mod(joint - ref)^2)
}

validate_ordpp <- function(object, newdata, weights = NULL, type = NULL) {
  y <- .op_newy(object, newdata)
  if (is.null(type)) type <- .op_type(object)
  ordpp_index(y, object$direction, weights, object$settings$frequencies, type)
}

decompose_ordpp <- function(object, newdata = NULL, weights = NULL) {
  if (is.null(newdata)) {
    y <- object$transformed; w <- object$weights
    if (is.null(y)) stop("object does not store its transformed training data")
  } else {
    y <- .op_newy(object, newdata); w <- .op_weights(weights, nrow(y))
  }
  f <- object$settings$frequencies; a <- object$direction
  data.frame(variance_ratio = .op_variance_ratio(y, a, w),
             full = ordpp_index(y, a, w, f, "full"),
             shape = ordpp_index(y, a, w, f, "shape"))
}
