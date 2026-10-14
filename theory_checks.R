# Checks requested in review: (a) exact error of the sample-centered Type III test vs the same test centered at the
# true population grand mean (oracle); (b) adequacy of 30 covariate draws. Paste into R or source().
Sys.setenv(TH_RUN = "0"); source("theory_ancova.R")
NDES <- as.integer(Sys.getenv("TC_NDES", "24")); NX <- as.integer(Sys.getenv("TC_NX", "150"))
dir.create("output", showWarnings = FALSE)
res <- read.csv("output/ancova_full/ancova_results.csv")
sel <- subset(res, shape == "normal" & a_effect == 0)
set.seed(1); sel <- sel[sample(nrow(sel), min(40, nrow(sel))), ][seq_len(NDES), ]
one_design <- function(i) {
  d <- draw_design(sel$design_id[i], sel$seed[i]); gr <- build_groups(d); a <- d$a; A <- rep(gr$g, gr$n); N <- length(A)
  mxv <- rep(gr$mx, gr$n); sdv <- rep(gr$sd, gr$n); dA <- a - 1; sA <- contr.sum(a)[A, , drop = FALSE]; one <- rep(1, N)
  dfF <- N - 2 * a; set.seed(d$seed + 99L)
  out <- matrix(NA_real_, NX, 8, dimnames = list(NULL, c("Oneway","R_default","TypeII","ANCOVA","TypeIII_raw","TypeIII_ctr","TypeIII_pop","d_shift")))
  for (k in seq_len(NX)) {
    xc0 <- mxv + rnorm(N); xr <- gr$L[1] + xc0; xc <- xr - mean(xr); mu <- rep(gr$beta, gr$n) * xc0
    Pf <- proj(cbind(one, sA, xc, sA * xc)); P1 <- proj(cbind(one)); PA <- proj(cbind(one, sA)); PX <- proj(cbind(one, xc))
    Padd <- proj(cbind(one, xc, sA)); Praw <- proj(cbind(one, xr, sA * xr)); Pctr <- proj(cbind(one, xc, sA * xc))
    Ppop <- proj(cbind(one, xc0, sA * xc0)); Pfp <- proj(cbind(one, sA, xc0, sA * xc0))
    out[k, ] <- c(exact_p(PA - P1, PA, dA, N - a, mu, sdv), exact_p(PA - P1, Pf, dA, dfF, mu, sdv),
                  exact_p(Padd - PX, Pf, dA, dfF, mu, sdv), exact_p(Padd - PX, Padd, dA, N - a - 1, mu, sdv),
                  exact_p(Pf - Praw, Pf, dA, dfF, mu, sdv), exact_p(Pf - Pctr, Pf, dA, dfF, mu, sdv),
                  exact_p(Pfp - Ppop, Pfp, dA, N - 2 * a, mu, sdv), mean(xr) - gr$L[1])
  }
  out
}
cores <- max(1L, min(parallel::detectCores() - 1L, 8L))
L <- parallel::mclapply(seq_len(nrow(sel)), one_design, mc.cores = cores)
nm <- colnames(L[[1]])[1:7]
tab <- do.call(rbind, lapply(seq_along(L), function(i) { o <- L[[i]]
  data.frame(design_id = sel$design_id[i], N = sel$N[i], a = sel$a[i], sd_ratio = sel$sd_ratio[i],
    t(setNames(colMeans(o[1:30, nm]), paste0("m30_", nm))), t(setNames(colMeans(o[, nm]), paste0("mAll_", nm))),
    t(setNames(apply(o[, nm], 2, sd) / sqrt(30), paste0("se30_", nm)))) }))
write.csv(tab, "output/theory_checks.csv", row.names = FALSE)
cat("designs:", nrow(tab), " draws for the long run:", NX, "\n")
for (a in nm) cat(sprintf("%-12s |m30 - mAll|: mean %.4f max %.4f ; mean SE of the 30-draw mean %.4f\n", a,
   mean(abs(tab[[paste0("m30_", a)]] - tab[[paste0("mAll_", a)]])), max(abs(tab[[paste0("m30_", a)]] - tab[[paste0("mAll_", a)]])), mean(tab[[paste0("se30_", a)]])))
cat(sprintf("sample-centered vs population-centered (long run): mean diff %.4f  max abs diff %.4f  mean(sample) %.4f mean(pop) %.4f\n",
   mean(tab$mAll_TypeIII_ctr - tab$mAll_TypeIII_pop), max(abs(tab$mAll_TypeIII_ctr - tab$mAll_TypeIII_pop)), mean(tab$mAll_TypeIII_ctr), mean(tab$mAll_TypeIII_pop)))

# Table 6 of the manuscript: analytical error (long run) against the saved simulation, and 30 versus long-run draws.
# Use TC_NDES=40 TC_NX=300 to match the manuscript (about 10 to 20 minutes). Output: output/table6_theory.csv
nm6 <- nm[1:6]
simr <- res[match(tab$design_id, res$design_id), paste0("rej_", nm6)]
t6 <- do.call(rbind, lapply(nm6, function(a) {
  ex <- tab[[paste0("mAll_", a)]]; sm <- simr[[paste0("rej_", a)]]; d30 <- abs(tab[[paste0("m30_", a)]] - ex)
  data.frame(analysis = a, mean_abs_diff_sim = round(mean(abs(ex - sm)), 3), max_abs_diff_sim = round(max(abs(ex - sm)), 3),
             corr_sim = round(cor(ex, sm), 3), mean_diff_30v300 = round(mean(d30), 4), max_diff_30v300 = round(max(d30), 4)) }))
print(t6); write.csv(t6, "output/table6_theory.csv", row.names = FALSE)
