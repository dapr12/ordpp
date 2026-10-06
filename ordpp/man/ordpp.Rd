\name{ordpp}
\alias{ordpp}
\alias{predict.ordpp}
\alias{print.ordpp}
\alias{plot.ordpp}
\title{Fit an Exploratory Ordinal Projection}
\description{Jointly estimates monotone category scores and a unit-length
projection direction by maximizing a characteristic-function discrepancy
between the projected joint distribution and its counterpart under item
independence with the observed marginals.}
\usage{
ordpp(x, weights = NULL, scoring = c("item", "common", "linear"),
      lambda = 0, gamma = 0.01, nstart = 3L, maxit = 500L,
      frequencies = seq(0.2, 3, length.out = 15),
      min_gap = 0.03, variance_floor = 0.05,
      index = c("full", "shape"), seed = NULL)
\method{predict}{ordpp}(object, newdata, ...)
\method{print}{ordpp}(x, ...)
\method{plot}{ordpp}(x, type = c("projection", "scores"), ...)
}
\arguments{
\item{x}{For \code{ordpp}, a complete data frame of numeric category codes or
ordered factors (numeric labels are used as ranks only). For methods, a
fitted object.}
\item{weights}{Optional positive observation weights; their scale is irrelevant.}
\item{scoring}{\code{"item"}: item-specific monotone curves;
\code{"common"}: one raw curve shared by items with identical response scales;
\code{"linear"}: fixed equally spaced scores.}
\item{lambda}{Nonnegative L1 penalty on loadings. Encourages concentration;
it does not perform exact variable selection.}
\item{gamma}{Nonnegative penalty on unequal adjacent score gaps.}
\item{nstart}{Number of optimizer starts (the first is the leading principal
component of equally scored items).}
\item{maxit}{Iteration budget per start for L-BFGS-B.}
\item{frequencies}{Positive frequencies, on the scale of the independence
standard deviation, at which characteristic functions are compared.}
\item{min_gap}{Lower bound on each raw gap as a fraction of the mean gap.}
\item{variance_floor}{Soft lower bound on the projection variance (the
reference variance is one).}
\item{index}{\code{"full"} compares the projected law with the independence
reference; \code{"shape"} first rescales the projection to the reference
variance, removing the pairwise-covariance component. See
\code{\link{decompose_ordpp}}.}
\item{seed}{Optional seed for reproducible starts; the caller's random number state is restored on exit.}
\item{object}{A fitted \code{ordpp} object.}
\item{newdata}{Complete data with the same named items and no unseen categories.}
\item{type}{Kernel density of the training projection, or category-score curves.}
\item{...}{Further graphical arguments, or unused.}
}
\value{An object of class \code{"ordpp"}: a list with \code{direction},
\code{scores}, \code{levels}, \code{variable_names}, \code{projection},
\code{transformed}, \code{weights}, \code{index} (unpenalized),
\code{objective} (penalized), \code{roughness}, \code{projection_variance},
\code{n}, \code{scoring}, \code{convergence}, \code{at_bounds},
\code{starts}, \code{settings} and \code{call}. \code{predict} returns the
projected scores of \code{newdata} using the frozen training dictionary.}
\details{All declared categories must occur in the training data; missing
values and unseen categories are rejected. Scores are strictly increasing,
centred and of unit variance under the (weighted) training marginals. Weights
change the target distribution but do not provide design-based standard
errors. Inspect \code{convergence}, \code{starts} and \code{at_bounds}
before interpretation.}
\seealso{\code{\link{validate_ordpp}}, \code{\link{heldout_test_ordpp}},
\code{\link{permutation_ordpp}}, \code{\link{stability_ordpp}}}
\examples{
s <- simulate_ordpp(150, 5, seed = 1)
f <- ordpp(s$x, scoring = "linear", nstart = 2, seed = 2)
f
predict(f, s$x[1:3, ])
\donttest{
# Item-specific monotone scores (slower)
fi <- ordpp(s$x, scoring = "item", nstart = 2, seed = 2)
fi$scores
}
}
