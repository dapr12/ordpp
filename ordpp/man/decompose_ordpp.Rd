\name{decompose_ordpp}
\alias{decompose_ordpp}
\title{Variance and Shape Components of the Projection Index}
\description{Reports the ratio of the projection variance to its independence
variance together with the full and shape indices.}
\usage{decompose_ordpp(object, newdata = NULL, weights = NULL)}
\arguments{
\item{object}{A fitted \code{ordpp} object or frozen map.}
\item{newdata}{Optional held-out data; if omitted the training data are used.}
\item{weights}{Optional positive weights for \code{newdata}.}
}
\value{A one-row data frame with \code{variance_ratio}, \code{full} and
\code{shape}.}
\details{At low frequencies the full index is dominated by
\eqn{(\mathrm{Var}(Z)/\mathrm{Var}(Z^0)-1)^2}, which is driven by pairwise
correlations between transformed items and rewards both inflated and deflated
projection variance. The shape index removes this component. A variance ratio
well below one indicates a contrast among positively correlated items.}
\examples{
s <- simulate_ordpp(150, 5, seed = 1)
f <- suppressWarnings(ordpp(s$x, scoring = "linear", nstart = 1,
                            maxit = 30, seed = 2))
decompose_ordpp(f)
}
