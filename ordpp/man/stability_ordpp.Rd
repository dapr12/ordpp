\name{stability_ordpp}
\alias{stability_ordpp}
\alias{print.ordpp_stability}
\title{Bootstrap Stability of Loadings and Category Scores}
\description{Refits the projection on bootstrap resamples of rows (or of
clusters), aligns signs with the original direction and summarizes variability.}
\usage{
stability_ordpp(fit, x, B = 50L, seed = NULL, weights = NULL, cluster = NULL,
                trace = FALSE)
\method{print}{ordpp_stability}(x, ...)
}
\arguments{
\item{fit}{A fitted \code{ordpp} object; its settings are reused.}
\item{x}{For \code{stability_ordpp}, the data used for \code{fit}; for the
print method, a stability result.}
\item{B}{Number of bootstrap refits.}
\item{seed}{Optional seed for resampling and optimizer starts; the caller's random number state is restored on exit.}
\item{weights}{Optional weights, resampled with their rows.}
\item{cluster}{Optional cluster identifiers; whole clusters are resampled.}
\item{trace}{Print progress every 5 refits.}
\item{...}{Unused.}
}
\value{An object of class \code{"ordpp_stability"} with a loading
\code{summary}, the matrix of aligned \code{directions}, per-item score
matrices \code{scores}, \code{index}, \code{convergence}, absolute cosine
\code{congruence} with the original direction and the number of
\code{failures} (for example, resamples missing a rare category).}
\details{The summaries describe algorithmic stability. Percentile intervals
after data-driven selection need not have nominal coverage, and resampling
rows ignores stratification and other design features.}
\examples{
\donttest{
s <- simulate_ordpp(150, 4, seed = 1)
f <- suppressWarnings(ordpp(s$x, scoring = "linear", nstart = 1,
                            maxit = 20, seed = 2))
stability_ordpp(f, s$x, B = 5)
}
}
