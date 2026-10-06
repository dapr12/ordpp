\name{permutation_ordpp}
\alias{permutation_ordpp}
\title{Refitted Permutation Test of Item Independence}
\description{Independently permutes each item and repeats the entire
fixed-setting optimization, so that the null distribution includes the
search over scores and directions.}
\usage{permutation_ordpp(x, B = 99L, seed = NULL, fit = NULL, trace = FALSE, ...)}
\arguments{
\item{x}{Complete ordinal data frame.}
\item{B}{Number of permutations.}
\item{seed}{Optional seed for the permutations and optimizer starts; the caller's random number state is restored on exit.}
\item{fit}{Optional \code{ordpp} object whose settings (scoring, penalties,
starts, iterations, frequencies, safeguards, index type) are reused; this
guarantees that observed and null fits use the same procedure.}
\item{trace}{Print progress every 10 permutations.}
\item{...}{Fitting arguments passed to \code{\link{ordpp}} when \code{fit}
is not supplied; \code{weights} and \code{seed} are not allowed.}
}
\value{An object of class \code{"ordpp_test"}; the observed statistic is the
maximized penalized objective of a refit to \code{x}, returned as
\code{observed_fit}. \code{null_convergence} records optimizer codes of the
null fits.}
\details{Valid for iid unweighted rows under mutual item independence when the
settings are fixed independently of \code{x}. If tuning parameters were chosen
using \code{x}, the selection must be repeated inside each permutation. The
test is computationally expensive; \code{\link{heldout_test_ordpp}} is a fast
alternative when a held-out sample is available.}
\examples{
\donttest{
s <- simulate_ordpp(120, 4, seed = 1)
f <- suppressWarnings(ordpp(s$x, scoring = "linear", nstart = 1,
                            maxit = 20, seed = 2))
permutation_ordpp(s$x, B = 9, seed = 3, fit = f)
}
}
