# ordpp

Projection pursuit with monotone scoring for ordered categorical (ordinal) data in R,
together with the replication materials for the article
*ordpp: Projection Pursuit with Monotone Scoring for Ordinal Data in R*
(D. A. Perez Ruiz, submitted to The R Journal).

## Install

From CRAN (once available):

```r
install.packages("ordpp")
```

Development version from this repository (the package lives in the `ordpp/` folder):

```r
# install.packages("remotes")
remotes::install_github("GITHUB-USERNAME/ordpp", subdir = "ordpp")
```

## Quick start

```r
library(ordpp)
sim  <- simulate_ordpp(n = 600, p = 6, seed = 103)
fit  <- ordpp(sim$x[1:300, ], scoring = "item", seed = 105)
fit
heldout_test_ordpp(fit, sim$x[301:600, ], B = 999, seed = 106)
```

## Repository layout

    ordpp/               the R package (GPL-3)
    scripts/             replication scripts, run from the repository root
    results/, figures/   outputs used by the article
    perez-ruiz.tex       article; RJwrapper.tex is the R Journal wrapper
    ordpp.bib, packages.bib

## Reproducing the article

```bash
R CMD INSTALL ordpp
install.packages(c("psychTools", "psych", "Gifi"))   # in R; ordPens from the CRAN archive
Rscript scripts/reproduce_article.R                  # every number, table and figure; a few minutes
FULL=1 Rscript scripts/reproduce_article.R           # also re-runs the long simulations (hours)
pdflatex RJwrapper; bibtex RJwrapper; pdflatex RJwrapper; pdflatex RJwrapper; pdflatex RJwrapper
```

The SAPA (`bfi`) data are loaded from the psychTools package and are not
redistributed here.
