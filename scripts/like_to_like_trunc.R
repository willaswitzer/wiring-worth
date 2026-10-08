# Session 7: like-to-like slopes as a function of simulated axon completeness.
#
# Run from the repo root after scripts/truncate_like_to_like.py:
#   .r-env/bin/Rscript scripts/like_to_like_trunc.R fit  <track> [workers] [max_sets]
#       refit the main-variant models on every (cond, level, rep) truncation set
#   .r-env/bin/Rscript scripts/like_to_like_trunc.R perm <track> <n> [workers] [levels]
#       shuffled-wiring control on the radial sets (default levels 1,0.5,0.3), same permutation
#       seeds as like_to_like.R, so level 1 reproduces the Session 6 null
#   .r-env/bin/Rscript scripts/like_to_like_trunc.R boot <track> <B> [workers]
#       paired bootstrap over presynaptic neurons: every bootstrap sample is refitted at all radial
#       levels and at rep 0 of both nulls, so CIs are for the shift slope(f) - slope(1)
# track = poisson (Nsyn + density models) or tweedie (Ld model). The Tweedie fits cost ~8x more,
# so the tweedie track uses only levels in LD_LEVELS (and rep 0 of the nulls).
# All modes resume from their output CSV in data/truncation/.
#
# Models, similarity measures, covariates (post dendrite level) and random effects are those of the
# Session 6 main variant (scripts/like_to_like.R, which this file sources). Projection types are
# those that pass the inclusion rule at every radial level (see FULL_PROJS), fixed across levels.

suppressPackageStartupMessages(source("scripts/like_to_like.R"))

TRUNC_DIR <- "data/truncation"
LD_LEVELS <- c(1, 0.7, 0.5, 0.3)
VARIANT <- "main"
COVARS <- VARIANTS[[VARIANT]]$covars

SET_DIR <- sprintf("%s/sets_v%d", TRUNC_DIR, VERSION)
sets <- read.csv(file.path(SET_DIR, "index.csv"), stringsAsFactors = FALSE)
sets$level <- round(sets$level, 2)
sets <- sets[order(sets$cond != "radial", -sets$level, sets$cond, sets$rep), ]

# Projection types: Ding's inclusion rule (> 30 synapses, > 5 presynaptic neurons with synapses,
# none > 50%) must hold on the real data at EVERY radial level, so no slope is fitted on a
# projection that truncation has emptied. V1->HVA (feedforward) fails below f = 1: its synapses are
# distal (3 of 119 left at f = 0.5), so it is excluded; see notebooks/07_results.md.
radial_ids <- sets$set[sets$cond == "radial"]
set_pairs_raw <- function(set_id) {
  t <- as.data.frame(read_parquet(file.path(SET_DIR, paste0(set_id, ".parquet"))))
  p <- pairs; p$Ld_um[t$row + 1] <- t$Ld_um; p$Nsyn[t$row + 1] <- t$Nsyn; p
}
radial_pairs <- lapply(radial_ids, set_pairs_raw)
FULL_PROJS <- setNames(lapply(SIMS, function(s) {
  ok <- lapply(radial_pairs, function(p) {
    d <- p; d$sim <- real_sim[[s]]; included_projs(d[!is.na(d$sim), ])
  })
  sort(Reduce(intersect, ok))
}), SIMS)
rm(radial_pairs)
cat("projection types fitted:", paste(names(FULL_PROJS), sapply(FULL_PROJS, paste, collapse = ", "), sep = ": ", collapse = "; "), "\n")
STARTS <- readRDS(file.path(OUT_DIR, "real_starts.rds"))
TRACK_MODELS <- list(poisson = c("Nsyn", "density"), tweedie = "Ld")

set_pairs <- function(set_id) {
  # pairs with Ld_um and Nsyn replaced by one truncation set (read from its own file)
  t <- as.data.frame(read_parquet(file.path(SET_DIR, paste0(set_id, ".parquet"))))
  stopifnot(nrow(t) == nrow(pairs))
  p <- pairs
  r <- t$row + 1
  p$Ld_um[r] <- t$Ld_um
  p$Nsyn[r] <- t$Nsyn
  p
}

model_data <- function(p, sim_name, sim_values) {
  d <- p
  d$sim <- sim_values[[sim_name]] / 10
  d <- d[!is.na(d$sim) & d$proj %in% FULL_PROJS[[sim_name]], ]
  d$proj <- factor(d$proj, levels = FULL_PROJS[[sim_name]])
  d$pre_proj <- droplevels(d$pre_proj)
  d
}

slopes_se <- function(fit) {
  # per-projection slope and Wald SE from fixef + vcov (equal to emtrends for this design)
  b <- fixef(fit)$cond; V <- vcov(fit)$cond
  lv <- levels(fit$frame$proj)
  if (is.null(lv)) lv <- levels(attr(fit, "proj_levels"))
  do.call(rbind, lapply(lv, function(k) {
    w <- setNames(rep(0, length(b)), names(b)); w["sim"] <- 1
    if (k != lv[1]) w[paste0("sim:proj", k)] <- 1
    data.frame(proj = k, slope = sum(w * b), SE = sqrt(drop(t(w) %*% V %*% w)))
  }))
}

fit_set <- function(p, models, sim_values = real_sim) {
  out <- list()
  for (sim_name in SIMS) {
    d_all <- model_data(p, sim_name, sim_values)
    for (model in models) {
      # a bootstrap sample can leave a projection with no synapses (or no Ld > 0 pairs); its slope
      # is then not estimable, so that projection is dropped from this one fit
      d <- if (model == "density") d_all[d_all$Ld_um > 0, ] else d_all
      has <- tapply(d$Nsyn, d$proj, sum) > 0
      d <- d[d$proj %in% names(has)[!is.na(has) & has], ]
      d$proj <- droplevels(d$proj); d$pre_proj <- droplevels(d$pre_proj)
      st <- STARTS[[paste(VARIANT, sim_name, model)]]
      st <- st[lengths(st) > 0]
      t0 <- Sys.time()
      res <- tryCatch({
        # start values from the Session 6 fit; retry without them if they don't fit this glmmTMB build
        fit <- tryCatch(suppressWarnings(fit_one(model, d, COVARS, start = st)),
                        error = function(e) suppressWarnings(fit_one(model, d, COVARS)))
        s <- slopes_se(fit)
        s$converged <- fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess)
        rm(fit); gc(verbose = FALSE)
        s
      }, error = function(e) {
        message(sprintf("fit failed (%s, %s): %s", sim_name, model, conditionMessage(e)))
        data.frame(proj = FULL_PROJS[[sim_name]], slope = NA_real_, SE = NA_real_, converged = FALSE)
      })
      res$fit_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      out[[length(out) + 1]] <- cbind(sim = sim_name, model = model, res)
    }
  }
  do.call(rbind, out)
}

track_sets <- function(track) {
  if (track == "tweedie") sets[sets$level %in% LD_LEVELS & sets$rep == 0, ] else sets
}

append_csv <- function(res, f) {
  write.table(res, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
}

# ---------------------------------------------------------------- fit: real similarity, every set
run_fit <- function(track, workers, max_sets = Inf) {
  f <- file.path(TRUNC_DIR, sprintf("trunc_slopes_%s.csv", track))
  done <- if (file.exists(f)) unique(read.csv(f)$set) else character(0)
  todo <- track_sets(track)
  todo <- head(todo[!todo$set %in% done, ], max_sets)
  t0 <- Sys.time()
  for (b in split(seq_len(nrow(todo)), ceiling(seq_len(nrow(todo)) / (workers * 2)))) {
    res <- mclapply(b, function(i) {
      s <- todo[i, ]
      cbind(s[, c("set", "cond", "level", "rep")], fit_set(set_pairs(s$set), TRACK_MODELS[[track]]),
            row.names = NULL)
    }, mc.cores = workers, mc.preschedule = FALSE)
    ok <- vapply(res, is.data.frame, TRUE)
    if (!all(ok)) cat("worker errors:", sum(!ok), "\n", paste(unlist(res[!ok]), collapse = "\n"), "\n")
    append_csv(do.call(rbind, res[ok]), f)
    cat(sprintf("fit/%s: %d of %d sets, %.1f min\n", track, length(unique(read.csv(f)$set)),
                nrow(track_sets(track)), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

# ---------------------------------------------------------------- perm: shuffled control on radial sets
run_perm_trunc <- function(track, n, workers, levels) {
  f <- file.path(TRUNC_DIR, sprintf("trunc_perm_%s.csv", track))
  done <- if (file.exists(f)) unique(read.csv(f)[, c("level", "perm")]) else data.frame(level = numeric(0), perm = integer(0))
  jobs <- expand.grid(perm = seq_len(n), level = levels)
  jobs <- jobs[!paste(jobs$level, jobs$perm) %in% paste(round(done$level, 2), done$perm), ]
  t0 <- Sys.time()
  for (b in split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / (workers * 2)))) {
    res <- mclapply(b, function(i) {
      j <- jobs[i, ]
      p <- set_pairs(sprintf("radial_%.2f_0", j$level))
      cbind(level = j$level, perm = j$perm,
            fit_set(p, TRACK_MODELS[[track]], sim_values = similarity(permute_props(10000 + j$perm))))
    }, mc.cores = workers, mc.preschedule = FALSE)
    ok <- vapply(res, is.data.frame, TRUE)
    if (!all(ok)) cat("worker errors:", sum(!ok), "\n")
    append_csv(do.call(rbind, res[ok]), f)
    cat(sprintf("perm/%s: %d of %d jobs, %.1f min\n", track, nrow(unique(read.csv(f)[, c("level", "perm")])),
                n * length(levels), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

# ---------------------------------------------------------------- boot: paired bootstrap over pre neurons
boot_sets <- function(track) {
  s <- sets[sets$rep == 0, ]
  if (track == "tweedie") s <- s[s$level %in% LD_LEVELS, ]
  s
}

boot_one <- function(b, track) {
  set.seed(20000 + b)
  pres <- unique(pairs$pre_root_id)
  draw <- sample(pres, length(pres), replace = TRUE)
  rows <- unlist(lapply(seq_along(draw), function(k) which(pairs$pre_root_id == draw[k])))
  copy <- rep(seq_along(draw), vapply(draw, function(x) sum(pairs$pre_root_id == x), 1L))
  out <- list()
  for (i in seq_len(nrow(boot_sets(track)))) {
    s <- boot_sets(track)[i, ]
    p <- set_pairs(s$set)[rows, ]
    p$pre_proj <- factor(paste(copy, p$pre_proj))  # each drawn copy is its own random-effect level
    sim_b <- lapply(real_sim, function(x) x[rows])
    out[[i]] <- cbind(boot = b, s[, c("set", "cond", "level", "rep")], fit_set(p, TRACK_MODELS[[track]], sim_b),
                      row.names = NULL)
  }
  do.call(rbind, out)
}

run_boot <- function(track, B, workers) {
  f <- file.path(TRUNC_DIR, sprintf("trunc_boot_%s.csv", track))
  done <- if (file.exists(f)) unique(read.csv(f)$boot) else integer(0)
  todo <- setdiff(seq_len(B), done)
  t0 <- Sys.time()
  for (b in split(todo, ceiling(seq_along(todo) / workers))) {
    res <- mclapply(b, boot_one, track = track, mc.cores = workers, mc.preschedule = FALSE)
    ok <- vapply(res, is.data.frame, TRUE)
    if (!all(ok)) cat("worker errors:", sum(!ok), "\n")
    append_csv(do.call(rbind, res[ok]), f)
    cat(sprintf("boot/%s: %d of %d, %.1f min\n", track, length(unique(read.csv(f)$boot)), B,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]; track <- args[2]
if (mode == "fit") {
  run_fit(track, if (length(args) >= 3) as.integer(args[3]) else 2,
          if (length(args) >= 4) as.integer(args[4]) else Inf)
} else if (mode == "perm") {
  lv <- if (length(args) >= 5) as.numeric(strsplit(args[5], ",")[[1]]) else c(1, 0.5, 0.3)
  run_perm_trunc(track, as.integer(args[3]), if (length(args) >= 4) as.integer(args[4]) else 2, lv)
} else if (mode == "boot") {
  run_boot(track, as.integer(args[3]), if (length(args) >= 4) as.integer(args[4]) else 2)
}
