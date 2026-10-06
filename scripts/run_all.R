# Runs the full replication. From the project root:
#   R CMD INSTALL ordpp
#   Rscript scripts/run_all.R            # full budgets (hours)
#   ORDPP_QUICK=1 Rscript scripts/run_all.R   # smoke test (minutes)
for (s in c("00_explanatory_figures.R", "01_synthetic_example.R",
            "02_application_bfi.R", "03_benchmark.R", "04_runtime.R",
            "05_sensitivity.R", "06_package_bib.R", "08_size_power.R",
            "09_add_ordpens.R", "07_derived_macros.R", "10_revision_tables.R")) {
  message("==== ", s, " ====")
  t0 <- Sys.time()
  local(source(file.path("scripts", s), local = TRUE))
  message("finished in ", format(round(Sys.time() - t0, 1)))
}
dir.create("results", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
writeLines(sprintf("\\newcommand{\\Rversion}{%s}",
                   paste(R.version$major, R.version$minor, sep = ".")),
           "results/macros-session.tex")
