# User-downloaded ESS CSV only. No data are bundled or retrieved by this script.
library(ordpp)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript scripts/application_ess.R file.csv country_code")
d <- read.csv(args[1]); d <- d[d$cntry == args[2], ]
items <- c("trstprl", "trstlgl", "trstplc", "trstplt", "trstprt")
if (!all(c(items, "anweight") %in% names(d))) stop("verify round-specific variables/weights")
x <- d[, items]
for (j in seq_along(x)) {
  x[[j]][!x[[j]] %in% 0:10] <- NA_real_ # Includes refusal/don't know/no answer.
  x[[j]] <- ordered(x[[j]], levels = 0:10)
}
ok <- complete.cases(x) & is.finite(d$anweight) & d$anweight > 0
x <- x[ok, ]; w <- d$anweight[ok]
# Household/PSU-grouped split should replace row split for final analysis.
set.seed(401); tr <- sample.int(nrow(x), floor(.7*nrow(x)))
te <- setdiff(seq_len(nrow(x)), tr)
f <- ordpp(x[tr, ], weights = w[tr], scoring = "item", nstart = 5, maxit = 200, seed = 402)
print(f); print(validate_ordpp(f, x[te, ], weights = w[te]))
dir.create("results", showWarnings = FALSE)
saveRDS(f, "results/ess_fit.rds")
write.csv(data.frame(projection = predict(f, x[te, ]), weight = w[te]),
          "results/ess_holdout.csv", row.names = FALSE)
# Fitting weights are not design-based standard errors. No iid permutation test here.
