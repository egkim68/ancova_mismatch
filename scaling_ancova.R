# SCALING EXPERIMENT (H8): multiply every group size by 1, 2, 4, 8, 16 for designs chosen by rule.
# Hypothesis mismatch: error tends to 1 as n grows. Variance mismatch: error settles at a fixed level.
# Needs explore_ancova.R in the same folder (its functions are reused). Base R only. Writes CSV only.
Sys.setenv(SIM_RUN = "0")
source("explore_ancova.R")                                           # defines draw_design(), run_design(), METHODS ...
SC_REPS <- as.integer(Sys.getenv("SC_REPS", unset = "300")); SC_DESIGNS <- as.integer(Sys.getenv("SC_DESIGNS", unset = "12"))
SC_SEED <- as.integer(Sys.getenv("SC_SEED", unset = "9100000"))       # pool of design seeds, separate from the main runs
FACTORS <- c(1, 2, 4, 8, 16)
SC_DIR <- file.path(OUTPUT_ROOT, "ancova_scaling"); dir.create(SC_DIR, recursive = TRUE, showWarnings = FALSE)
N_REPS <- SC_REPS                                                    # run_design() reads this global

rules <- list(
  set1_raw_mismatch  = function(d) d$a_effect == 0 & d$slope_range > 0 & d$d0 >= 2 & d$sd_ratio == 1,
  set2_cov_imbalance = function(d) d$a_effect == 0 & d$slope_range == 0 & d$cov_imbalance > 0 & d$sd_ratio == 1,
  set3_variance      = function(d) d$a_effect == 0 & d$sd_ratio >= 2 & d$pairing == -1 & d$imbalance >= 0.4)
chosen <- list(); i <- 0L
while (any(sapply(names(rules), function(s) length(chosen[[s]]) < SC_DESIGNS))) {
  i <- i + 1L; d <- draw_design(i, SC_SEED + i)
  for (s in names(rules)) if (length(chosen[[s]]) < SC_DESIGNS && isTRUE(rules[[s]](d))) chosen[[s]][[length(chosen[[s]]) + 1L]] <- d
}
out <- list(); t0 <- Sys.time()
for (s in names(chosen)) for (k in seq_along(chosen[[s]])) for (f in FACTORS) {
  d <- chosen[[s]][[k]]; d$N_target <- d$N_target * f
  out[[length(out) + 1L]] <- cbind(set = s, factor = f, run_design(d))
}
res <- do.call(rbind, out); write.csv(res, file.path(SC_DIR, "scaling_results.csv"), row.names = FALSE)
cat(sprintf("Written %d rows in %.1f min: %s\n", nrow(res), as.numeric(difftime(Sys.time(), t0, units = "mins")),
            normalizePath(file.path(SC_DIR, "scaling_results.csv"))))
