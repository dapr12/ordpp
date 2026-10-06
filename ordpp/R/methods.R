predict.ordpp <- function(object, newdata, ...) {
  as.vector(.op_newy(object, newdata) %*% object$direction)
}

print.ordpp <- function(x, ...) {
  cat("Ordinal projection pursuit\n")
  cat("  n =", x$n, "; scoring =", x$scoring, "; index type =", .op_type(x), "\n")
  cat("  index =", format(x$index, digits = 4))
  if (!is.na(x$objective)) cat("; penalized objective =", format(x$objective, digits = 4))
  cat("; Var(Z)/Var(Z0) =", format(x$projection_variance, digits = 3), "\n")
  if (!is.na(x$convergence))
    cat("  convergence code =", x$convergence, "; starts =", NROW(x$starts),
        "; gap parameters on bounds =", x$at_bounds, "\n")
  cat("Loadings:\n")
  print(round(x$direction, 3))
  invisible(x)
}

plot.ordpp <- function(x, type = c("projection", "scores"), ...) {
  type <- match.arg(type)
  if (type == "projection") {
    w <- if (isTRUE(x$weights_supplied)) x$weights else NULL
    dens <- suppressWarnings(stats::density(x$projection, weights = w))
    plot(dens, main = if (is.null(w)) "Training projection" else
           "Training projection (weighted density)",
         xlab = "Projected score", ...)
    graphics::rug(x$projection)
  } else {
    k <- max(lengths(x$scores))
    z <- sapply(x$scores, function(s) c(s, rep(NA_real_, k - length(s))))
    cols <- seq_len(ncol(z))
    matplot(seq_len(k), z, type = "b", lty = 1, pch = 1, col = cols,
            xlab = "Ordered category index", ylab = "Standardized score", ...)
    legend("topleft", legend = names(x$scores), col = cols, lty = 1,
           cex = .7, bty = "n")
  }
  invisible(x)
}

print.ordpp_test <- function(x, ...) {
  cat(x$method, "\n")
  cat("  H0:", x$hypothesis, "\n")
  cat("  n =", x$n, "; B =", x$B, "; index type =", x$type, "\n")
  cat("  statistic =", format(x$statistic, digits = 4),
      "; p-value =", format(x$p.value, digits = 3), "\n")
  if (!is.null(x$null_convergence) && any(x$null_convergence != 0))
    cat("  null fits with nonzero convergence code:",
        sum(x$null_convergence != 0), "\n")
  invisible(x)
}

print.ordpp_stability <- function(x, ...) {
  cat("Bootstrap stability of an ordinal projection\n")
  cat("  B =", x$B, "; failed refits =", x$failures,
      if (x$clustered) "; cluster resampling" else "", "\n")
  if (length(x$congruence))
    cat("  median |cosine| with original direction =",
        format(stats::median(x$congruence), digits = 3), "\n")
  print(x$summary, digits = 3, row.names = FALSE)
  invisible(x)
}
