\name{heldout_test_ordpp}
\alias{heldout_test_ordpp}
\alias{print.ordpp_test}
\title{Frozen-Map Held-out Permutation Test}
\description{Tests mutual independence of the items in held-out data along a
map (scores and loadings) estimated on separate observations.}
\usage{
heldout_test_ordpp(object, newdata, B = 999L, seed = NULL, type = NULL)
\method{print}{ordpp_test}(x, ...)
}
\arguments{
\item{object}{A fitted \code{ordpp} object or frozen map.}
\item{newdata}{Held-out data not used to estimate \code{object}.}
\item{B}{Number of permutations.}
\item{seed}{Optional seed for reproducible permutations; the caller's random number state is restored on exit.}
\item{type}{Index type; defaults to the fitted type.}
\item{x}{A test result.}
\item{...}{Unused.}
}
\value{An object of class \code{"ordpp_test"} with \code{statistic},
\code{null}, \code{p.value}, \code{B}, \code{n}, \code{type}, \code{method}
and \code{hypothesis}.}
\details{Each item column of the held-out data is permuted independently and
the index is recomputed with the map held fixed. Under iid rows and mutual
item independence the permuted samples are exchangeable with the observed one,
so the Monte Carlo p-value \eqn{(1+\sum_b 1\{T_b \ge T\})/(B+1)} is valid in finite
samples. Validity requires that \code{newdata} played no part in fitting or
tuning. No refitting is performed, so the test is fast. It is not valid for
clustered or complex survey samples.}
\examples{
s <- simulate_ordpp(300, 5, seed = 1)
f <- suppressWarnings(ordpp(s$x[1:150, ], scoring = "linear", nstart = 1,
                            maxit = 30, seed = 2))
heldout_test_ordpp(f, s$x[151:300, ], B = 99, seed = 3)
}
