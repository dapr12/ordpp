# Regenerates packages.bib from citation() for the versions actually used.
# Only software manuals here; packages described in journal articles (tourr,
# mirt, lavaan, FactoMineR, ca) are cited through those articles in ordpp.bib.
pk <- c("ordpp", "Gifi", "ordPens", "psych", "psychTools", "energy")
keys <- c(ordpp = "ordpp", Gifi = "Gifi", ordPens = "ordPens", psych = "psych",
          psychTools = "psychTools", energy = "energy")
out <- character()
for (p in pk) {
  if (!requireNamespace(p, quietly = TRUE)) { message("not installed: ", p); next }
  ci <- utils::citation(p)[[1]]          # first entry of each package citation
  bib <- format(utils::toBibtex(ci))
  bib[1] <- sub("\\{[^,]*,", paste0("{", keys[[p]], ","), bib[1])
  out <- c(out, bib, "")
}
writeLines(out, "packages.bib")
message("packages.bib written; check entries before submission")
