# Writes 100 random data sets (one replication each, from the study's generator) and the R p values to CSV,
# so validate_python_ancova.py can recompute them with statsmodels. Paste into R or run with Rscript.
Sys.setenv(SIM_RUN = "0", SIM_MODE = "pilot", SIM_OUTPUT_ROOT = tempdir())
source("explore_ancova.R")
dir.create("check", showWarnings = FALSE); set.seed(11); data <- list(); pv <- list()
for (i in 1:100) {
  d <- designs[sample(nrow(designs), 1), ]; gr <- build_groups(d); A <- rep(gr$g, gr$n); N <- length(A)
  xc0 <- rep(gr$mx, gr$n) + rnorm(N); xr <- gr$L[1] + xc0
  y <- rep(gr$tau, gr$n) + rep(gr$beta, gr$n) * xc0 + rep(gr$sd, gr$n) * draw_err(N, d$shape)
  p <- test_all(y, A, xr, d$a)
  data[[i]] <- data.frame(case = i, A = A, x = xr, y = y); pv[[i]] <- data.frame(case = i, t(p))
}
write.csv(do.call(rbind, data), "check/datasets.csv", row.names = FALSE)
write.csv(do.call(rbind, pv), "check/r_pvalues.csv", row.names = FALSE)
cat("Written check/datasets.csv and check/r_pvalues.csv\n")
