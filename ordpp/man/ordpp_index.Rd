\name{ordpp_index}
\alias{ordpp_index}
\alias{validate_ordpp}
\title{Characteristic-Function Discrepancy and Held-out Evaluation}
\description{Compares the characteristic function of a projection with the
product of the characteristic functions of its item components, i.e. with the
projection's law under item independence with the same marginals.}
\usage{
ordpp_index(y, direction, weights = NULL,
            frequencies = seq(0.2, 3, length.out = 15),
            type = c("full", "shape"))
validate_ordpp(object, newdata, weights = NULL, type = NULL)
}
\arguments{
\item{y}{Finite numeric matrix of transformed items.}
\item{direction}{Nonzero loading vector, one entry per column of \code{y}.}
\item{weights}{Optional positive observation weights.}
\item{frequencies}{Positive frequencies on the independence-SD scale.}
\item{type}{\code{"full"} or \code{"shape"}; for \code{validate_ordpp}
the default is the type used in fitting.}
\item{object}{A fitted \code{ordpp} object.}
\item{newdata}{Complete held-out data with known categories and named items.}
}
\value{A nonnegative scalar.}
\details{\code{validate_ordpp} freezes the fitted scores and loadings; the
reference marginals and the independence variance are re-estimated from the
held-out rows. The index measures departure from item independence along the
projection; it is not a test of polarization or of a particular factor
structure. A finite frequency grid need not detect every alternative. The
empirical index is positive in finite samples even under independence; use
\code{\link{heldout_test_ordpp}} for calibration.}
\examples{
y <- as.matrix(expand.grid(a = c(-1, 1), b = c(-1, 1)))
ordpp_index(y, c(1, 1))                         # exact independence: 0
ordpp_index(cbind(c(-1, 1), c(-1, 1)), c(1, 1))  # perfect dependence
}
