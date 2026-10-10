# Session 7: does radial truncation steepen the RF like-to-like slope within each cell type?
#
# Run from the repo root after scripts/truncate_like_to_like.py:
#   .r-env/bin/Rscript scripts/celltype_truncation_check.R 200
# V1->V1, Nsyn Poisson GLMM (main-variant covariates, random intercept per presynaptic neuron), fitted
# separately to 23P and 4P presynaptic axons at radial f = 1, 0.7, 0.5, 0.3. CIs: paired bootstrap over
# presynaptic neurons of that cell type (seeds 30000 + b). Writes data/truncation/celltype_boot_rf_nsyn.csv.
# First run 2026-10-09 in the cloud test environment (glmmTMB 1.1.8); see notebooks/07_results.md.
suppressPackageStartupMessages(source("scripts/like_to_like.R"))
SET_DIR <- "data/truncation/sets_v1822"
ids <- c("radial_1.00_0", "radial_0.70_0", "radial_0.50_0", "radial_0.30_0")
P <- lapply(ids, function(id) { t <- as.data.frame(read_parquet(file.path(SET_DIR, paste0(id, ".parquet")))); p <- pairs; p$Ld_um[t$row + 1] <- t$Ld_um; p$Nsyn[t$row + 1] <- t$Nsyn; p$sim <- real_sim$rf_dist_deg / 10; p[p$proj == "V1->V1", ] })
names(P) <- ids
slope <- function(d) { f <- suppressWarnings(glmmTMB(Nsyn ~ sim + post_dendrite_level + (1 | pre_proj), family = poisson, data = d)); fixef(f)$cond[["sim"]] }
one <- function(b, ct) {
  set.seed(30000 + b)
  pres <- unique(P[[1]]$pre_root_id[P[[1]]$pre_cell_type == ct])
  draw <- sample(pres, length(pres), replace = TRUE)
  sapply(ids, function(id) {
    p <- P[[id]]
    d <- do.call(rbind, lapply(seq_along(draw), function(k) { x <- p[p$pre_root_id == draw[k], ]; x$pre_proj <- paste(k, x$pre_root_id); x }))
    d$pre_proj <- factor(d$pre_proj); slope(d)
  })
}
B <- as.integer(commandArgs(TRUE)[1])
res <- list()
for (ct in c("23P", "4P")) {
  est <- sapply(ids, function(id) { d <- P[[id]][P[[id]]$pre_cell_type == ct, ]; d$pre_proj <- droplevels(d$pre_proj); slope(d) })
  bs <- do.call(rbind, parallel::mclapply(seq_len(B), one, ct = ct, mc.cores = as.integer(Sys.getenv("WORKERS", "4"))))
  sh <- bs - bs[, 1]
  res[[ct]] <- data.frame(cell_type = ct, level = c(1, .7, .5, .3), slope = est, shift = est - est[1],
                          lo = apply(sh, 2, quantile, .025), hi = apply(sh, 2, quantile, .975), n_boot = sum(complete.cases(bs)))
}
r <- do.call(rbind, res); print(r, row.names = FALSE, digits = 3)
write.csv(r, "data/truncation/celltype_boot_rf_nsyn.csv", row.names = FALSE)
