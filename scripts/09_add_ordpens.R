# Adds the ordPens comparator to the synthetic and bfi test tables, using the
# same fixed splits as scripts 01 and 02 (no refitting of ordpp models).
source("scripts/common.R")
if (!requireNamespace("ordPens", quietly = TRUE)) stop("ordPens is not installed")
B_HELD <- budget(999, 49)

append_row <- function(file, row) {
  d <- read.csv(file, check.names = FALSE)
  d <- d[d$method != row$method, , drop = FALSE]   # avoid duplicates on rerun
  d <- rbind(d, row[, names(d)])
  write.csv(d, file, row.names = FALSE)
  d
}

## Synthetic example (split as in 01_synthetic_example.R)
sim <- simulate_ordpp(n = 600, p = 6, seed = 103)
set.seed(104)
train <- sample.int(nrow(sim$x), 300)
test <- setdiff(seq_len(nrow(sim$x)), train)
m <- map_ordpens(sim$x[train, ])
r <- evaluate_map(m, sim$x[test, ], B = B_HELD, latent = sim$latent[test])
r$p_value <- fmt_p(r$p_value)
ab <- append_row("results/syn-ablation.csv", r)
ab$p_value <- fmt_p(as.numeric(ab$p_value))
names(ab) <- c("Method", "Train index", "Test index", "Test shape",
               "Var ratio", "Held-out p", "|cor| latent")
write_tex_table(ab, "results/tab-syn-ablation.tex", digits = 3)

## bfi (split as in 02_application_bfi.R)
data("bfi", package = "psychTools")
items <- paste0("N", 1:5)
x <- bfi[complete.cases(bfi[, items]), items]
x <- as.data.frame(lapply(x, ordered, levels = 1:6))
n <- nrow(x); set.seed(2026)
idx <- sample(rep(c("train", "valid", "test"), times = c(n %/% 2, n %/% 4, n - n %/% 2 - n %/% 4)))
xtr <- x[idx == "train", ]; xte <- x[idx == "test", ]
keyed <- as_ordpp_map(xtr, equal_scores(xtr), rep(1, length(items)), label = "Unit-weighted sum")
m <- map_ordpens(xtr)
r <- evaluate_map(m, xte, B = B_HELD)
r$abs_latent_cor <- NULL
r$cor_keyed <- stats::cor(predict(m, xte), predict(keyed, xte))
r$p_value <- fmt_p(r$p_value)
tt <- append_row("results/bfi-test.csv", r)
tt$p_value <- fmt_p(as.numeric(tt$p_value))
names(tt) <- c("Map", "Train index", "Test index", "Test shape", "Var ratio",
               "Held-out p", "cor(sum)")
write_tex_table(tt, "results/tab-bfi-test.tex", digits = 3)
print(tt[, 1:3])
message("ordPens rows added")
