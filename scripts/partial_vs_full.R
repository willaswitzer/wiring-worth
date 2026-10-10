# Session 7: like-to-like on real partially vs. fully extended axons, within cell type (third curve).
#
# Run from the repo root after scripts/prepare_partial.py:
#   .r-env/bin/Rscript scripts/partial_vs_full.R [B] [n_perm] [workers]      (defaults 200 100 2)
#
# For 23P and 4P presynaptic axons, V1->V1, each similarity (RF distance, orientation) and model (Nsyn;
# Nsyn with offset log Ld on Ld > 0), fits
#     y ~ sim * axon_level + post_dendrite_level [+ soma_dist_100um] + (1 | pre)       (Poisson)
# so the slope for full axons is `sim` and the partial - full difference is `sim:axon_levelpartial`.
# Uncertainty: bootstrap over presynaptic neurons, resampled within cell type x axon level (B samples);
# shuffled-wiring control: postsynaptic functional properties permuted within cell type x area (as
# like_to_like.R, seeds 10000 + i), giving a null for each group's slope.
# Variants: main (dendrite level) and soma (+ soma-to-soma distance per 100 um; Session 5 found partial
# and full 4P somata placed differently relative to the co-registered population).
# Writes data/truncation/partial_vs_full.csv, partial_vs_full_boot.csv, partial_vs_full_perm.csv.

suppressPackageStartupMessages(source("scripts/like_to_like.R"))
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 200
N_PERM <- if (length(args) >= 2) as.integer(args[2]) else 100
W <- if (length(args) >= 3) as.integer(args[3]) else 2
OUTD <- "data/truncation"

fp <- as.data.frame(read_parquet(sprintf("data/l2l_pairs_fullpartial_v%d.parquet", VERSION)))
fp$pre_root_id <- as.character(fp$pre_root_id); fp$post_root_id <- as.character(fp$post_root_id)
fp <- fp[fp$proj == "V1->V1" & fp$pre_cell_type %in% c("23P", "4P"), ]
fp$axon_level <- factor(fp$pre_axon_level, levels = c("full", "partial"))
fp$post_dendrite_level <- factor(fp$post_dendrite_level)
fp$soma_dist_100um <- fp$soma_dist_um / 100
pi_ <- match(fp$pre_root_id, neurons$root_id); po_ <- match(fp$post_root_id, neurons$root_id)
sims_for <- function(props) {
  a <- rf_unit(neurons$readout_loc_x[pi_], neurons$readout_loc_y[pi_])
  b <- rf_unit(props$readout_loc_x[po_], props$readout_loc_y[po_])
  rf <- acos(pmin(pmax(rowSums(a * b), -1), 1)) * 180 / pi
  d <- abs(neurons$pref_ori[pi_] - props$pref_ori[po_]) %% 180
  ori <- pmin(d, 180 - d)
  ori[!(neurons$gOSI[pi_] > GOSI_MIN & props$gOSI[po_] > GOSI_MIN)] <- NA
  list(rf_dist_deg = rf, ori_diff_deg = ori)
}
real <- sims_for(neurons)
stopifnot(max(abs(real$rf_dist_deg - fp$rf_dist_deg)) < 1e-4)

COV <- list(main = "post_dendrite_level", soma = c("post_dendrite_level", "soma_dist_100um"))

fit_ct <- function(d, model, covars) {
  if (model == "density") d <- d[d$Ld_um > 0, ]
  rhs <- paste(c("sim * axon_level", covars, if (model == "density") "offset(log(Ld_um))", "(1 | pre_root_id)"), collapse = " + ")
  f <- suppressWarnings(glmmTMB(as.formula(paste("Nsyn ~", rhs)), family = poisson, data = d))
  b <- fixef(f)$cond; V <- vcov(f)$cond; k <- "sim:axon_levelpartial"
  data.frame(full = b[["sim"]], partial = b[["sim"]] + b[[k]], diff = b[[k]], SE_diff = sqrt(V[k, k]),
             converged = f$fit$convergence == 0 && isTRUE(f$sdr$pdHess))
}

run_all <- function(d0, sim_vals, variants = names(COV)) {
  out <- list()
  for (ct in c("23P", "4P")) for (s in names(sim_vals)) for (model in c("Nsyn", "density")) for (v in variants) {
    d <- d0; d$sim <- sim_vals[[s]] / 10
    d <- d[d$pre_cell_type == ct & !is.na(d$sim), ]
    r <- tryCatch(fit_ct(d, model, COV[[v]]), error = function(e) data.frame(full = NA, partial = NA, diff = NA, SE_diff = NA, converged = FALSE))
    out[[length(out) + 1]] <- cbind(cell_type = ct, sim = s, model = model, variant = v, r)
  }
  do.call(rbind, out)
}

t0 <- Sys.time()
est <- run_all(fp, real)
est$n_full <- sapply(est$cell_type, function(ct) length(unique(fp$pre_root_id[fp$pre_cell_type == ct & fp$axon_level == "full"])))
est$n_partial <- sapply(est$cell_type, function(ct) length(unique(fp$pre_root_id[fp$pre_cell_type == ct & fp$axon_level == "partial"])))
write.csv(est, file.path(OUTD, "partial_vs_full.csv"), row.names = FALSE)
print(est, digits = 3, row.names = FALSE)

# bootstrap over presynaptic neurons within cell type x axon level (main variant)
grp <- unique(fp[, c("pre_root_id", "pre_cell_type", "axon_level")])
rows_by_pre <- split(seq_len(nrow(fp)), fp$pre_root_id)
boot_one <- function(b) {
  set.seed(40000 + b)
  draw <- unlist(lapply(split(grp$pre_root_id, paste(grp$pre_cell_type, grp$axon_level)), function(x) sample(x, length(x), replace = TRUE)))
  idx <- unlist(rows_by_pre[draw]); cp <- rep(seq_along(draw), lengths(rows_by_pre[draw]))
  d <- fp[idx, ]; d$pre_root_id <- paste(cp, d$pre_root_id)
  cbind(boot = b, run_all(d, lapply(real, function(x) x[idx]), "main"))
}
bt <- do.call(rbind, mclapply(seq_len(B), boot_one, mc.cores = W))
write.csv(bt, file.path(OUTD, "partial_vs_full_boot.csv"), row.names = FALSE)

# shuffled-wiring control (main variant)
pm <- do.call(rbind, mclapply(seq_len(N_PERM), function(i) cbind(perm = i, run_all(fp, sims_for(permute_props(10000 + i)), "main")),
                              mc.cores = W))
write.csv(pm, file.path(OUTD, "partial_vs_full_perm.csv"), row.names = FALSE)
cat(sprintf("done: %d boot, %d perms, %.1f min\n", B, N_PERM, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
