# Checks the hand-built tests in explore_ancova.R against base R nested-model comparisons (anova of lm fits)
# on random designs from the same generator. Reports the largest absolute difference in p values and F-scale
# relative differences. Paste into R or run with Rscript from the folder that holds explore_ancova.R.
Sys.setenv(SIM_RUN = "0", SIM_MODE = "pilot", SIM_OUTPUT_ROOT = tempdir())
source("explore_ancova.R")
set.seed(1); worst <- c(); nchk <- 300
for (i in seq_len(nchk)) {
  d <- designs[sample(nrow(designs), 1), ]; gr <- build_groups(d); a <- d$a
  A <- rep(gr$g, gr$n); N <- length(A); xc0 <- rep(gr$mx, gr$n) + rnorm(N); xr <- gr$L[1] + xc0
  y <- rep(gr$tau, gr$n) + rep(gr$beta, gr$n) * xc0 + rep(gr$sd, gr$n) * draw_err(N, d$shape)
  mine <- test_all(y, A, xr, a); f <- factor(A); xc <- xr - mean(xr); full <- lm(y ~ f * xc)
  ref <- c(Oneway = anova(lm(y ~ f))[1, "Pr(>F)"], R_default = anova(full)[1, "Pr(>F)"],
           TypeII = anova(lm(y ~ xc + f * xc))[2, "Pr(>F)"], ANCOVA = anova(lm(y ~ xc + f))[2, "Pr(>F)"],
           TypeIII_raw = anova(lm(y ~ xr + f:xr), full)[2, "Pr(>F)"], TypeIII_ctr = anova(lm(y ~ xc + f:xc), full)[2, "Pr(>F)"])
  worst <- rbind(worst, abs(mine[names(ref)] - ref))
}
cat(sprintf("%d random designs; largest absolute difference in p over all analyses: %.2e\n", nchk, max(worst)))
print(signif(apply(worst, 2, max), 3))

# ---- packages: car::Anova (Types II and III, dummy and sum coding) and HC3 (white.adjust = "hc3") ----------
if (requireNamespace("car", quietly = TRUE)) {
  set.seed(2); w2 <- c()
  for (i in seq_len(nchk)) {
    d <- designs[sample(nrow(designs), 1), ]; gr <- build_groups(d); a <- d$a
    A <- rep(gr$g, gr$n); N <- length(A); xc0 <- rep(gr$mx, gr$n) + rnorm(N); xr <- gr$L[1] + xc0
    y <- rep(gr$tau, gr$n) + rep(gr$beta, gr$n) * xc0 + rep(gr$sd, gr$n) * draw_err(N, d$shape)
    mine <- test_all(y, A, xr, a); f <- factor(A); xc <- xr - mean(xr)
    fs <- f; contrasts(fs) <- contr.sum(a)
    pr <- function(tab, row) tab[row, "Pr(>F)"]
    ref <- c(TypeII = pr(car::Anova(lm(y ~ f * xc), type = 2, test.statistic = "F"), "f"),
             TypeIII_raw_trt = pr(car::Anova(lm(y ~ f * xr), type = 3, test.statistic = "F"), "f"),
             TypeIII_raw_sum = pr(car::Anova(lm(y ~ fs * xr), type = 3, test.statistic = "F"), "fs"),
             TypeIII_ctr = pr(car::Anova(lm(y ~ fs * xc), type = 3, test.statistic = "F"), "fs"),
             TypeIII_ctr_HC3 = pr(car::Anova(lm(y ~ fs * xc), type = 3, white.adjust = "hc3"), "fs"))
    m <- c(TypeII = mine[["TypeII"]], TypeIII_raw_trt = mine[["TypeIII_raw"]], TypeIII_raw_sum = mine[["TypeIII_raw"]],
           TypeIII_ctr = mine[["TypeIII_ctr"]], TypeIII_ctr_HC3 = mine[["TypeIII_ctr_HC3"]])
    w2 <- rbind(w2, abs(m - ref))
  }
  cat(sprintf("\ncar::Anova, %d random designs; largest absolute difference in p: %.2e\n", nchk, max(w2)))
  print(signif(apply(w2, 2, max), 3))
}
