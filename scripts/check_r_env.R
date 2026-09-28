# Session 5: confirm the R environment for the like-to-like models (no model fitting).
# Run from the repo root:  .r-env/bin/Rscript scripts/check_r_env.R
# Environment: r-environment.yml (spec) / r-environment.lock.txt (exact versions).

pkgs <- c("glmmTMB", "emmeans", "performance", "arrow")
for (p in pkgs) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  cat(sprintf("%-12s %s\n", p, as.character(packageVersion(p))))
}
cat(R.version.string, "\n")
stopifnot(arrow::arrow_info()$capabilities[["parquet"]])

files <- c("usable_set_v1822_s2", "synapses_usable_v1822", "pairs_usable_v1822")
for (f in files) {
  d <- arrow::read_parquet(file.path("data", paste0(f, ".parquet")))
  cat(sprintf("\n%s: %d rows x %d cols\n", f, nrow(d), ncol(d)))
  cls <- vapply(d, function(x) class(x)[1], "")
  print(table(cls))
}

# 64-bit root ids must survive the round trip (they exceed 2^53, so doubles would corrupt them).
pairs <- arrow::read_parquet("data/pairs_usable_v1822.parquet", col_select = c("pre_root_id", "post_root_id"))
cat("\nroot id class:", class(pairs$pre_root_id), "\n")
stopifnot(inherits(pairs$pre_root_id, "integer64"))
cat("example root id:", format(pairs$pre_root_id[1]), "\n")
# Smoke test of the compiled stack on SIMULATED data only (not a like-to-like model):
# Poisson counts with a log-length offset and a random intercept.
set.seed(1)
sim <- data.frame(g = factor(rep(1:30, each = 20)), x = factor(rep(c("a", "b"), 300)),
                  len = runif(600, 5, 100))
sim$y <- rpois(600, sim$len / 100 * exp(0.3 * (sim$x == "b") + rnorm(30, 0, 0.3)[sim$g]))
fit <- glmmTMB(y ~ x + offset(log(len)) + (1 | g), family = poisson, data = sim)
stopifnot(fit$fit$convergence == 0)
cat("\nsmoke test (simulated): fixed effect b =", round(fixef(fit)$cond[["xb"]], 3), "(true 0.3)\n")
print(emmeans(fit, ~ x, type = "response", offset = 0))
cat("R2 (Nakagawa):", round(unlist(performance::r2(fit)), 3), "\n")

cat("\nR environment OK\n")
