# Run from the directory containing run_audit.R.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, file.exists("run_audit.R"))
suppressPackageStartupMessages(library(data.table))
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
all_runs <- lapply(1:20, function(seed) {
  tmp <- tempfile()
  r <- file.path(R.home("bin"), "Rscript")
  status <- system2(r, c("--vanilla", "run_audit.R", shQuote(args[1]),
                         shQuote(args[2]), shQuote(tmp), seed), stdout = FALSE)
  stopifnot(status == 0L)
  res <- fread(file.path(tmp, "candidate_summary.csv"))
  res[, seed := seed]
  unlink(tmp, recursive = TRUE)
  res
})
all <- rbindlist(all_runs)
fwrite(all, file.path(args[3], "seed_candidates.csv"))
fwrite(all[, .(matched = sum(positive_matched == 3L), mast = sum(positive_mast == 3L),
               combined_mean_margin = sum(retains_slope & specific),
               combined_each_margin = sum(retains_slope & specific_all_above_margin)),
           by = .(seed, disease)], file.path(args[3], "seed_summary.csv"))
