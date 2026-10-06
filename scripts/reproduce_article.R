# Reproduces every number, table and figure in the article in a few minutes.
#
#   Rscript scripts/reproduce_article.R          # fast path (< 10 minutes)
#   FULL=1 Rscript scripts/reproduce_article.R   # also re-runs the long analyses (hours)
#
# Fast path: re-runs the explanatory figures and the worked synthetic example,
# re-evaluates the frozen bfi maps on the test rows, and rebuilds all tables,
# figures and in-text numbers from the saved outputs of the long-running
# analyses (bfi tuning and bootstrap, sensitivity, benchmark, size/power,
# runtime), which are stored in results/. At the end it compares the
# regenerated CSV files with the saved ones and reports any difference.
t0 <- Sys.time()
full <- identical(Sys.getenv("FULL"), "1")
run <- function(s) {
  message("==== ", s)
  local(source(file.path("scripts", s), local = TRUE))
}

# Keep a copy of the saved results to compare against.
dir.create("results_saved", showWarnings = FALSE)
csvs <- list.files("results", pattern = "\\.csv$")
file.copy(file.path("results", csvs), "results_saved", overwrite = TRUE)

if (full) {
  for (s in c("02_application_bfi.R", "05_sensitivity.R", "03_benchmark.R",
              "08_size_power.R", "04_runtime.R")) run(s)
}
for (s in c("00_explanatory_figures.R", "01_synthetic_example.R",
            "09_add_ordpens.R", "04b_runtime_plot.R", "07_derived_macros.R",
            "10_revision_tables.R")) run(s)
writeLines(sprintf("\\newcommand{\\Rversion}{%s}",
                   paste(R.version$major, R.version$minor, sep = ".")),
           "results/macros-session.tex")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

# Compare regenerated tables with the saved ones (numeric tolerance 1e-6).
message("==== comparison with saved results")
same <- TRUE
for (f in csvs) {
  a <- read.csv(file.path("results_saved", f)); b <- read.csv(file.path("results", f))
  num <- intersect(names(a)[vapply(a, is.numeric, TRUE)], names(b))
  ok <- nrow(a) == nrow(b) && all(vapply(num, function(v)
    isTRUE(all.equal(a[[v]], b[[v]], tolerance = 1e-6)), TRUE))
  if (!ok) { same <- FALSE; message("DIFFERS: ", f) }
}
if (same) message("all regenerated tables match the saved results")
message("finished in ", format(round(Sys.time() - t0, 1)))
