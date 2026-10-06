\name{as_ordpp_map}
\alias{as_ordpp_map}
\title{Frozen Map from External Scores and Loadings}
\description{Wraps category scores and a direction obtained by another method
(for example optimal-scaling or polychoric PCA) as an \code{ordpp} object, so
that it can be evaluated with \code{\link{validate_ordpp}},
\code{\link{heldout_test_ordpp}} and \code{\link{decompose_ordpp}}.}
\usage{
as_ordpp_map(x, scores, direction, weights = NULL,
             frequencies = seq(0.2, 3, length.out = 15),
             index = c("full", "shape"), label = "external")
}
\arguments{
\item{x}{Training data used to standardize the scores.}
\item{scores}{List with one numeric vector per item, one value per observed
category in increasing category order.}
\item{direction}{Loading vector.}
\item{weights}{Optional positive training weights.}
\item{frequencies, index}{As in \code{\link{ordpp}}.}
\item{label}{Name stored in the \code{scoring} field.}
}
\value{An object of class \code{"ordpp"} that can be evaluated but not refitted.}
\examples{
s <- simulate_ordpp(150, 4, seed = 1)
m <- as_ordpp_map(s$x, lapply(s$x, function(v) sort(unique(v))),
                  direction = c(1, 1, 1, -1), label = "sum score")
validate_ordpp(m, s$x)
}
