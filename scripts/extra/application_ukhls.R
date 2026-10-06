# GHQ item names and coding must be checked in the chosen wave's documentation.
# No licensed observations are distributed with the package.
library(ordpp)
if (!requireNamespace("haven", quietly = TRUE)) stop("install haven to read Stata")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript scripts/application_ukhls.R indresp.dta prefix weight_name")
d <- haven::read_dta(args[1]); prefix <- args[2]
items <- paste0(prefix, "_scghq", letters[1:12])
if (!all(c(items, args[3]) %in% names(d))) stop("verify wave/item/weight names")
x <- as.data.frame(lapply(d[, items], as.numeric))
for (j in seq_along(x)) {
  x[[j]][!x[[j]] %in% 1:4] <- NA_real_
  x[[j]] <- ordered(x[[j]], levels = 1:4)
}
# Verify each response order against questionnaire. Never reverse merely by item wording.
w <- as.numeric(d[[args[3]]]); ok <- complete.cases(x) & is.finite(w) & w > 0
fit <- ordpp(x[ok, ], weights = w[ok], scoring = "item", nstart = 5, maxit = 200, seed = 401)
print(fit)
dir.create("results", showWarnings = FALSE); saveRDS(fit, "results/ukhls_fit.rds")
# Choose an eligible self-completion weight, not automatically adult-interview weight.
# Final longitudinal validation must freeze scores/loadings and account for households.
