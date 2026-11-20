# Reproduce the six-section exploratory attribution audit.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 3:4)
suppressPackageStartupMessages({library(data.table); library(Matrix); library(RANN)})
inp <- readRDS(args[1])
cases <- fread(args[2])
output <- args[3]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
stopifnot(identical(inp$cells$cell, colnames(inp$expression)),
          length(inp$sections) == 6L, nrow(unique(cases)) == 20L)

nearest_other <- function(xy, refs) {
  stopifnot(length(refs) >= 2L)
  near <- nn2(xy[refs, , drop = FALSE], xy, k = 2)
  self <- refs[near$nn.idx[, 1L]] == seq_len(nrow(xy))
  ifelse(self, near$nn.dists[, 2L], near$nn.dists[, 1L])
}

fit_distance <- function(design, y, valid) {
  use <- valid & complete.cases(design)
  keep <- apply(design[use, , drop = FALSE], 2, sd) > 0
  keep[1:2] <- TRUE
  fit <- lm.fit(design[use, keep, drop = FALSE], y[use])
  c(slope = unname(fit$coefficients[2]), n = sum(use))
}

estimates <- list()
quality <- list()
pairs <- list()
seed <- if (length(args) == 4L) as.integer(args[4]) else 42L
set.seed(seed)
for (sect in inp$sections) {
  a <- inp$cells[tissue == sect]
  group <- if (sect %in% c("18_72_TB", "Lung_TB_18-78", "Lung_TB_18-79")) "MTB" else "NTM"
  xy <- as.matrix(a[, .(x, y)])
  neural <- which(a$broad == "Neural")
  mast <- which(a$fine == "Mast_cells")
  nb <- nn2(xy, xy[neural, , drop = FALSE], k = 25, searchtype = "radius", radius = 100)
  selected <- apply(nb$nn.idx, 1L, function(row) {
    ix <- row[row > 0L]
    eligible <- ix[a$broad[ix] != "Neural"]
    if (!length(eligible)) return(NA_integer_)
    eligible[sample.int(length(eligible), 1L)]
  })
  hit <- !is.na(selected)
  offsets <- sqrt(rowSums((xy[neural[hit], , drop = FALSE] - xy[selected[hit], , drop = FALSE])^2))
  stopifnot(all(offsets <= 100), all(a$broad[selected[hit]] != "Neural"))
  refs <- unique(selected[hit])
  # Preserve the November sequence of 20 abundance-matched control draws per section.
  pool <- which(a$broad != "Neural")
  for (draw in seq_len(20L)) invisible(sample(pool, length(neural)))
  pairs[[sect]] <- data.table(tissue = sect, anchor = a$cell[neural[hit]],
                              matched = a$cell[selected[hit]], offset_px = offsets)
  d <- list(neural = nearest_other(xy, neural), matched = nearest_other(xy, refs),
            mast = nearest_other(xy, mast))
  eligible_types <- cases[disease == group, cell_type]
  target <- a$fine %in% eligible_types
  quality[[sect]] <- data.table(tissue = sect, cells = nrow(a), neural = length(neural),
                                matched_anchors = sum(hit), unique_matches = length(refs),
                                mast = length(mast), max_match_offset_px = max(offsets),
                                max_target_distance_cache_diff_px =
                                  max(abs(a$cached_distance[target] - d$neural[target])))
  for (case in which(cases$disease == group)) {
    gene <- cases$gene[case]
    ct <- cases$cell_type[case]
    idx <- which(a$fine == ct)
    s <- a[idx]
    y <- as.numeric(inp$expression[gene, s$cell])
    rho <- vapply(d, function(v) cor(y, v[idx], method = "spearman"), numeric(1))
    valid <- is.finite(s$area) & s$area > 0 & is.finite(s$ncount) & s$ncount > 0 &
             is.finite(s$cached_distance) & s$cached_distance >= 0
    X0 <- cbind(intercept = 1, d100 = s$cached_distance / 100)
    X3 <- cbind(X0, log1p(s$ncount), log1p(s$area), s$neg_per_probe,
                as.matrix(s[, .(Epithelial, Immune, Stromal, Unclassified)]), log1p(s$n_neighbours))
    m0 <- fit_distance(X0, y, valid)
    m3 <- fit_distance(X3, y, valid)
    estimates[[length(estimates) + 1L]] <- data.table(tissue = sect, disease = group, gene, cell_type = ct,
      rho_neural = rho["neural"], rho_matched = rho["matched"], rho_mast = rho["mast"],
      M0 = m0["slope"], M3 = m3["slope"], n_M0 = m0["n"], n_M3 = m3["n"])
  }
}
e <- rbindlist(estimates)
stopifnot(nrow(e) == 60L, all(complete.cases(e)))
res <- e[, .(sections = .N, positive_M3 = sum(M3 > 0), retention = mean(M3)/mean(M0),
             positive_neural = sum(rho_neural > 0), positive_matched = sum(rho_matched > 0),
             positive_mast = sum(rho_mast > 0), mean_delta = mean(rho_neural - rho_matched),
             positive_delta = sum(rho_neural > rho_matched),
             delta_above_margin = sum(rho_neural - rho_matched > 0.01)), by = .(disease, gene, cell_type)]
res[, retains_slope := positive_M3 == 3L & retention >= 0.4]
res[, specific := mean_delta > 0.01 & positive_delta == 3L]
res[, specific_all_above_margin := delta_above_margin == 3L]
ref <- res[, .(candidates = .N, neural = sum(positive_neural == 3L), matched = sum(positive_matched == 3L),
               mast = sum(positive_mast == 3L), retains_slope = sum(retains_slope), specific = sum(specific)),
           by = disease]
fwrite(e, file.path(output, "section_estimates.csv"))
fwrite(res, file.path(output, "candidate_summary.csv"))
fwrite(ref, file.path(output, "reference_summary.csv"))
fwrite(rbindlist(quality), file.path(output, "reference_qc.csv"))
fwrite(rbindlist(pairs), file.path(output, "matched_pairs.csv"))
print(ref)

