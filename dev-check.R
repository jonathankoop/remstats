#!/usr/bin/env Rscript
# Fast correctness check after touching src/*.cpp.
#
#   R CMD INSTALL . && Rscript dev-check.R
#
# Runs the highest-signal subset of inst/tinytest rather than the whole suite:
# the sampled-vs-full comparisons (which cross-check tomstats.cpp against
# tomstats_sampled.cpp on active/manual risksets) plus the endogenous and
# memory tests. ~47k assertions in a few seconds. Use tinytest::test_all()
# before committing.
suppressMessages({library(tinytest); library(remify); library(remstats)})

files <- c(
  "test_compare_sampled_stats4.R",          # decay + typed + riskset="active"
  "test_compare_sampled_stats6.R",          # riskset="manual", undirected
  "test-compare-sampled-typed-ext-FALSE.R", # extend_riskset_by_type = FALSE
  "test-compare-sampled-typed-ext-TRUE.R",  # extend_riskset_by_type = TRUE
  "test-endogenous-stats3.R",
  "test-endogenous-stats4.R",               # undirected
  "test-memory.R",                          # full / interval / decay
  "test-memory-custom.R",                   # custom (tabulated) decay kernel
  "test-remstats-typed-events.R",
  "test-weights.R"
)

dir <- "inst/tinytest"
total <- 0L; failed <- 0L
for (f in files) {
  path <- file.path(dir, f)
  if (!file.exists(path)) { cat(sprintf("%-42s MISSING\n", f)); next }
  res <- suppressWarnings(run_test_file(path, verbose = 0))
  n  <- length(res)
  nf <- sum(!vapply(res, isTRUE, logical(1)))
  total <- total + n; failed <- failed + nf
  cat(sprintf("%-42s %6d tests %4d fails%s\n", f, n, nf,
              if (nf > 0) "   <<<" else ""))
}
cat(sprintf("\n%s: %d assertions, %d failures\n",
            if (failed == 0) "PASS" else "FAIL", total, failed))
if (failed > 0) quit(status = 1)
