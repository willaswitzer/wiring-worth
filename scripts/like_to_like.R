# Session 6: like-to-like GLMMs (Ding et al. 2025) with a shuffled-wiring control.
#
# Run from the repo root after scripts/prepare_like_to_like.py:
#   .r-env/bin/Rscript scripts/like_to_like.R real                     # real fits (all variants)
#   .r-env/bin/Rscript scripts/like_to_like.R perm <variant> <track> <n> <min> [workers]
#     permutations, resumable; track = poisson (Nsyn + density models) or tweedie (Ld models).
#     The Tweedie Ld fits cost ~10x more time and ~1.7 GB each, so the tracks run separately
#     with different worker counts and end up with different permutation counts.
#
# Models (per similarity measure; fixed effects sim x projection type; random intercept per
# presynaptic neuron within projection type, as Ding's u_k(i,j),i; covariate postsynaptic
# dendrite level):
#   Ld       : Ld_um ~ sim * proj + covars + (1 | pre:proj), Tweedie (log link), all pairs
#   Nsyn     : Nsyn  ~ sim * proj + covars + (1 | pre:proj), Poisson, all pairs
#   density  : Nsyn  ~ sim * proj + covars + offset(log(Ld_um)) + (1 | pre:proj), Poisson, Ld > 0
# Similarity is in units of 10 degrees; slopes are on the log (link) scale per 10 degrees.
# Projection types are included only with > 30 synapses, > 5 presynaptic neurons with synapses,
# and no presynaptic neuron contributing > 50% of synapses (Ding et al. 2025), evaluated on the
# real data for each similarity x variant.
# Shuffle: postsynaptic functional properties (readout location, pref_ori, gOSI together) are
# permuted among all usable neurons within cell type x area; presynaptic properties stay real.

suppressPackageStartupMessages({
  library(glmmTMB)
  library(emmeans)
  library(arrow)
  library(parallel)
})

VERSION <- 1822
OUT_DIR <- "data/like_to_like"
dir.create(OUT_DIR, showWarnings = FALSE)
MONITOR_W_CM <- 56.5; MONITOR_H_CM <- 31.8; EYE_CM <- 15.0  # MICrONS Consortium 2025, Methods
GOSI_MIN <- 0.25
SIMS <- c("rf_dist_deg", "ori_diff_deg")
MODELS <- c("Ld", "Nsyn", "density")
WORKERS <- 3  # each permutation peaks at ~1.7 GB (Tweedie Ld fit); sized for an 8 GB machine

VARIANTS <- list(
  main    = list(covars = "post_dendrite_level", ext_only = FALSE),
  soma    = list(covars = c("post_dendrite_level", "soma_dist_100um"), ext_only = FALSE),
  ext_dend = list(covars = character(0), ext_only = TRUE)
)

# ---------------------------------------------------------------- data
pairs <- as.data.frame(read_parquet(sprintf("data/l2l_pairs_v%d.parquet", VERSION)))
neurons <- as.data.frame(read_parquet(sprintf("data/l2l_neurons_v%d.parquet", VERSION)))
pairs$pre_root_id <- as.character(pairs$pre_root_id)
pairs$post_root_id <- as.character(pairs$post_root_id)
neurons$root_id <- as.character(neurons$root_id)
pairs$pre_proj <- factor(paste(pairs$pre_root_id, pairs$proj))
pairs$post_dendrite_level <- factor(pairs$post_dendrite_level)
pairs$soma_dist_100um <- pairs$soma_dist_um / 100
pre_idx <- match(pairs$pre_root_id, neurons$root_id)
post_idx <- match(pairs$post_root_id, neurons$root_id)

rf_unit <- function(x, y) {
  v <- cbind(x * MONITOR_W_CM / 2, y * MONITOR_H_CM / 2, EYE_CM)
  v / sqrt(rowSums(v^2))
}
similarity <- function(post_props) {
  # post_props: neuron table (possibly permuted), rows aligned with `neurons`
  a <- rf_unit(neurons$readout_loc_x[pre_idx], neurons$readout_loc_y[pre_idx])
  b <- rf_unit(post_props$readout_loc_x[post_idx], post_props$readout_loc_y[post_idx])
  rf <- acos(pmin(pmax(rowSums(a * b), -1), 1)) * 180 / pi
  d <- abs(neurons$pref_ori[pre_idx] - post_props$pref_ori[post_idx]) %% 180
  ori <- pmin(d, 180 - d)
  tuned <- neurons$gOSI[pre_idx] > GOSI_MIN & post_props$gOSI[post_idx] > GOSI_MIN
  ori[!tuned] <- NA
  list(rf_dist_deg = rf, ori_diff_deg = ori)
}
real_sim <- similarity(neurons)
stopifnot(max(abs(real_sim$rf_dist_deg - pairs$rf_dist_deg)) < 1e-4,  # float32 inputs
          identical(is.na(real_sim$ori_diff_deg), is.na(pairs$ori_diff_deg)),
          max(abs(real_sim$ori_diff_deg - pairs$ori_diff_deg), na.rm = TRUE) < 1e-3)

permute_props <- function(seed) {
  set.seed(seed)
  grp <- paste(neurons$cell_type, neurons$area)
  idx <- seq_len(nrow(neurons))
  for (g in unique(grp)) {
    w <- which(grp == g)
    idx[w] <- w[sample.int(length(w))]
  }
  p <- neurons
  cols <- c("readout_loc_x", "readout_loc_y", "pref_ori", "gOSI")
  p[, cols] <- neurons[idx, cols]
  p
}

# ---------------------------------------------------------------- inclusion + fitting
base_data <- function(variant, sim_name, sim_values) {
  v <- VARIANTS[[variant]]
  d <- pairs
  d$sim <- sim_values[[sim_name]] / 10
  keep <- !is.na(d$sim)
  if (v$ext_only) keep <- keep & d$post_dendrite_level == "dendrite_extended"
  d[keep, ]
}

included_projs <- function(d) {
  s <- aggregate(Nsyn ~ proj + pre_root_id, d, sum)
  out <- c()
  for (k in unique(s$proj)) {
    x <- s[s$proj == k, ]
    tot <- sum(x$Nsyn)
    if (tot > 30 && sum(x$Nsyn > 0) > 5 && max(x$Nsyn) / tot <= 0.5) out <- c(out, k)
  }
  sort(out)
}

formula_for <- function(model, covars, n_proj) {
  # with a single included projection type the model has no proj term
  rhs <- paste(c(if (n_proj > 1) "sim * proj" else "sim", covars, if (model == "density") "offset(log(Ld_um))", "(1 | pre_proj)"),
               collapse = " + ")
  as.formula(paste(if (model == "Ld") "Ld_um" else "Nsyn", "~", rhs))
}

fit_one <- function(model, d, covars, start = NULL) {
  if (model == "density") d <- d[d$Ld_um > 0, ]
  fam <- if (model == "Ld") tweedie(link = "log") else poisson
  args <- list(formula_for(model, covars, nlevels(d$proj)), family = fam, data = d,
               control = glmmTMBControl(parallel = 1))
  if (!is.null(start)) args$start <- start
  fit <- do.call(glmmTMB, args)
  attr(fit, "proj_levels") <- d$proj
  fit
}

slopes_from_fixef <- function(fit) {
  b <- fixef(fit)$cond
  lv <- levels(fit$frame$proj)
  if (is.null(lv)) lv <- levels(attr(fit, "proj_levels"))
  s <- sapply(lv, function(k) b[["sim"]] + if (k == lv[1]) 0 else b[[paste0("sim:proj", k)]])
  data.frame(proj = lv, slope = unname(s))
}

prep_variant <- function(variant, sim_name) {
  d <- base_data(variant, sim_name, real_sim)
  projs <- included_projs(d)
  d <- d[d$proj %in% projs, ]
  d$proj <- factor(d$proj, levels = projs)
  d$pre_proj <- droplevels(d$pre_proj)
  list(d = d, projs = projs)
}

# ---------------------------------------------------------------- real fits
run_real <- function() {
  rows <- list(); incl <- list(); starts <- list()
  for (variant in names(VARIANTS)) for (sim_name in SIMS) {
    pv <- prep_variant(variant, sim_name)
    d <- pv$d
    s_all <- aggregate(Nsyn ~ proj, base_data(variant, sim_name, real_sim), sum)
    incl[[length(incl) + 1]] <- data.frame(variant, sim = sim_name, proj = s_all$proj, Nsyn = s_all$Nsyn,
                                           included = s_all$proj %in% pv$projs)
    for (model in MODELS) {
      t0 <- Sys.time()
      fit <- fit_one(model, d, VARIANTS[[variant]]$covars)
      et <- if (nlevels(d$proj) > 1) {
        as.data.frame(summary(emtrends(fit, ~ proj, var = "sim"), infer = c(TRUE, TRUE)))
      } else {
        cbind(proj = levels(d$proj), as.data.frame(summary(emtrends(fit, ~ 1, var = "sim"), infer = c(TRUE, TRUE))))
      }
      chk <- slopes_from_fixef(fit)
      stopifnot(max(abs(chk$slope - et$sim.trend)) < 1e-6)
      r2 <- if (variant == "main") tryCatch(unname(unlist(performance::r2(fit))), error = function(e) c(NA, NA)) else c(NA, NA)
      rows[[length(rows) + 1]] <- data.frame(
        variant, sim = sim_name, model, proj = as.character(et$proj), slope = et$sim.trend, SE = et$SE,
        CI_low = et$asymp.LCL, CI_high = et$asymp.UCL, p_wald = et$p.value,
        n_pairs = nrow(if (model == "density") d[d$Ld_um > 0, ] else d),
        converged = fit$fit$convergence == 0, pdHess = fit$sdr$pdHess,
        R2_conditional = r2[1], R2_marginal = r2[2],
        fit_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
      starts[[paste(variant, sim_name, model)]] <- list(beta = fixef(fit)$cond, betadisp = fit$fit$par[names(fit$fit$par) == "betadisp"],
                                                       psi = fit$fit$par[names(fit$fit$par) == "psi"],
                                                       theta = fit$fit$par[names(fit$fit$par) == "theta"])
      cat(sprintf("%-8s %-13s %-8s %5.1fs  %s\n", variant, sim_name, model, rows[[length(rows)]]$fit_seconds[1],
                  paste(sprintf("%s=%.3f", et$proj, et$sim.trend), collapse = "  ")))
    }
  }
  write.csv(do.call(rbind, rows), file.path(OUT_DIR, "real_slopes.csv"), row.names = FALSE)
  write.csv(do.call(rbind, incl), file.path(OUT_DIR, "inclusion.csv"), row.names = FALSE)
  saveRDS(starts, file.path(OUT_DIR, "real_starts.rds"))
}

# ---------------------------------------------------------------- permutations
TRACKS <- list(poisson = c("Nsyn", "density"), tweedie = "Ld")

perm_one <- function(i, variant, prepped, starts, models = MODELS) {
  props <- permute_props(10000 + i)
  sim_p <- similarity(props)
  out <- list()
  for (sim_name in SIMS) {
    projs <- prepped[[sim_name]]$projs
    d <- base_data(variant, sim_name, sim_p)
    d <- d[d$proj %in% projs, ]
    d$proj <- factor(d$proj, levels = projs)
    d$pre_proj <- droplevels(d$pre_proj)
    for (model in models) {
      st <- starts[[paste(variant, sim_name, model)]]
      st <- st[lengths(st) > 0]
      res <- tryCatch({
        fit <- suppressWarnings(fit_one(model, d, VARIANTS[[variant]]$covars, start = st))
        s <- slopes_from_fixef(fit)
        s$converged <- fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess)
        rm(fit); gc(verbose = FALSE)
        s
      }, error = function(e) data.frame(proj = projs, slope = NA_real_, converged = FALSE))
      out[[length(out) + 1]] <- cbind(perm = i, sim = sim_name, model = model, res)
    }
  }
  do.call(rbind, out)
}

run_perm <- function(variant, track, n_target, budget_min) {
  f <- file.path(OUT_DIR, sprintf("perm_slopes_%s_%s.csv", variant, track))
  done <- if (file.exists(f)) unique(read.csv(f)$perm) else integer(0)
  todo <- setdiff(seq_len(n_target), done)
  starts <- readRDS(file.path(OUT_DIR, "real_starts.rds"))
  prepped <- setNames(lapply(SIMS, function(s) prep_variant(variant, s)), SIMS)
  t0 <- Sys.time()
  while (length(todo) > 0 && as.numeric(difftime(Sys.time(), t0, units = "mins")) < budget_min) {
    batch <- head(todo, WORKERS * 2)
    res <- mclapply(batch, perm_one, variant = variant, prepped = prepped, starts = starts, models = TRACKS[[track]],
                    mc.cores = WORKERS, mc.preschedule = FALSE)
    ok <- vapply(res, is.data.frame, TRUE)
    if (!all(ok)) cat("worker errors:", sum(!ok), "\n")
    res <- do.call(rbind, res[ok])
    write.table(res, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
    todo <- setdiff(todo, unique(res$perm))
    n_done <- length(unique(read.csv(f)$perm))
    cat(sprintf("%s/%s: %d perms done, %.1f min elapsed\n", variant, track, n_done,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || args[1] == "real") {
  run_real()
} else if (args[1] == "perm") {
  if (length(args) >= 6) WORKERS <- as.integer(args[6])
  run_perm(args[2], args[3], as.integer(args[4]), as.numeric(args[5]))
}
