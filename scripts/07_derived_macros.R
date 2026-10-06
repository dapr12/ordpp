# Derived comparisons quoted in the text, computed from saved results (fast).
# Run after 01 and 02.
source("scripts/common.R")
g  <- read.csv("results/bfi-grid.csv")
tt <- read.csv("results/bfi-test.csv")
sy <- read.csv("results/syn-ablation.csv")
pick <- function(d, m) d[d$method == m, , drop = FALSE]
pct <- function(a, b) fmt(100 * (a / b - 1), 0)

best <- which.max(g$valid_index)
alt <- g[g$scoring == g$scoring[best] & g$gamma > g$gamma[best], , drop = FALSE]
alt <- alt[which.min(alt$gamma), , drop = FALSE]
ngap <- if (g$scoring[best] == "item") 5 * (6 - 2) else if (g$scoring[best] == "common") 6 - 2 else 0
sel <- tt[grepl("selected", tt$method), , drop = FALSE]
lin <- pick(tt, "ordpp linear"); gifi <- pick(tt, "Gifi princals")
ordp <- pick(tt, "ordPens ordPCA")

s_lin <- pick(sy, "linear"); s_item <- pick(sy, "item")
s_shape <- pick(sy, "item (shape index)")
syn_train <- readRDS("results/synthetic.rds")$fits$item_shape   # training shape index

write_macros(list(
  bfiSelValid = fmt(g$valid_index[best], 4), bfiSelBounds = g$at_bounds[best],
  bfiNgap = ngap,
  bfiAltGamma = if (nrow(alt)) format(alt$gamma) else "--",
  bfiAltValid = if (nrow(alt)) fmt(alt$valid_index, 4) else "--",
  bfiAltBounds = if (nrow(alt)) alt$at_bounds else "--",
  bfiLinTest = fmt(lin$test_index, 4), bfiGainPct = pct(sel$test_index, lin$test_index),
  bfiSelShape = fmt(sel$test_shape, 4), bfiLinShape = fmt(lin$test_shape, 4),
  bfiLinVarRatio = fmt(lin$variance_ratio, 2),
  bfiGifiTest = if (nrow(gifi)) fmt(gifi$test_index, 4) else "--",
  bfiOrdpensTest = if (nrow(ordp)) fmt(ordp$test_index, 4) else "--",
  synLinTest = fmt(s_lin$test_index, 4), synItemTest = fmt(s_item$test_index, 4),
  synGainPct = pct(s_item$test_index, s_lin$test_index),
  synLinCor = fmt(s_lin$abs_latent_cor, 3), synItemCor = fmt(s_item$abs_latent_cor, 3),
  synShapeFitVar = fmt(s_shape$variance_ratio, 2),
  synShapeFitCor = fmt(s_shape$abs_latent_cor, 2),
  synShapeFitTrain = fmt(syn_train$index, 4),
  synShapeFitTest = fmt(s_shape$test_shape, 4),
  synItemTestShape = fmt(s_item$test_shape, 4)),
  "results/macros-derived.tex")

# Sensitivity (written separately so 07 still runs if 05 has not).
if (file.exists("results/sensitivity.csv")) {
  se <- read.csv("results/sensitivity.csv", check.names = FALSE)
  num <- se[!se$Variant %in% c("default", "shape index"), ]
  sh <- se[se$Variant == "shape index", ]
  dflt <- se[se$Variant == "default", ]
  write_macros(list(
    sensMinCosNum = fmt(min(num[["|cos| with selected"]]), 5),
    sensMaxScore = fmt(max(num[["Max score change"]]), 2),
    sensMinGapValid = fmt(se[se$Variant == "min gap 0.10", "Valid. index"], 4),
    sensDefaultScore = fmt(dflt[["Max score change"]], 3),
    sensShapeCos = fmt(sh[["|cos| with selected"]], 2),
    sensShapeVar = fmt(sh[["Var ratio"]], 2)),
    "results/macros-derived-sens.tex")
}

# Size and power (from 08_size_power.R).
if (file.exists("results/size-power-summary.csv")) {
  sp <- read.csv("results/size-power-summary.csv")
  pw <- sp[sp$study == "power", ]
  pow <- function(m, c = 0.3) fmt(pw$rate[pw$method == m & pw$strength == c], 2)
  weak <- pw[pw$strength < 0.3, ]
  sz <- sp[sp$study == "size", ]
  write_macros(list(
    powPCA = pow("PCA, equal scores"), powPoly = pow("Polychoric PCA"),
    powLin = pow("ordpp linear"), powCommon = pow("ordpp common"),
    powItem = pow("ordpp item"), powGifi = pow("Gifi princals"),
    powWeakMin = fmt(min(weak$rate), 2), powWeakMax = fmt(max(weak$rate), 2),
    sizeLost = max(sz$runs[sz$method == "ordpp linear"]) - min(sz$runs[sz$method == "ordpp linear"])),
    "results/macros-derived-power.tex")
}

# Package versions for "Computational details".
ver <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "--"
write_macros(list(verPsych = ver("psych"), verGifi = ver("Gifi"), verOrdPens = ver("ordPens")),
             "results/macros-versions.tex")
message("derived macros written")
