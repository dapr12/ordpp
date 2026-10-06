\name{simulate_ordpp}
\alias{simulate_ordpp}
\title{Synthetic Ordinal Data}
\description{Five-category items obtained by thresholding latent continuous
responses, with informative and noise items.}
\usage{
simulate_ordpp(n = 500, p = 8, separation = 1.5, seed = NULL, strength = 1,
               scenario = c("mixture", "independent", "rare",
                            "two_factor", "response_style"))
}
\arguments{
\item{n}{Number of observations (at least 10).}
\item{p}{Number of items (at least 4).}
\item{separation}{Half-distance between the means of a two-component normal
mixture for the latent factor; 0 gives a single normal factor.}
\item{seed}{Optional seed; the caller's random number state is restored on exit.}
\item{strength}{Multiplier of all factor loadings; small values give weak
dependence for power studies, and 0 gives independent items.}
\item{scenario}{\code{"mixture"}: items 1--3 load positively and item 4
negatively on the latent factor; remaining items are noise.
\code{"independent"}: all items independent with skewed marginals.
\code{"rare"}: as \code{"mixture"} with a rarely used top category on the
informative items. \code{"two_factor"}: adds a second normal factor on items
5--7. \code{"response_style"}: 30 percent of respondents stretch all latent
responses, inducing dependence even among noise items.}
}
\value{A list with \code{x}, \code{latent}, \code{group}, \code{loading},
\code{thresholds} and \code{scenario}.}
\details{Zero separation still yields dependent items through a shared factor
and is not an independence null; use \code{scenario = "independent"} for that.}
\examples{
s <- simulate_ordpp(100, 6, seed = 1, scenario = "independent")
sapply(s$x, table)
}
