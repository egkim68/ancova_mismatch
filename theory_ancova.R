# Exact Type I error (Imhof 1961) conditional on covariate values, normal errors, then averaged over covariate draws.
# Sources the frozen design generator without running it.
Sys.setenv(SIM_RUN = "0", SIM_MODE = "full", SIM_SEED = "8419", SIM_OUTPUT_ROOT = tempdir())
src <- readLines("explore_ancova.R"); src <- src[!grepl("^dir.create", src)]
eval(parse(text = src[seq_len(grep("^# ---- run", src) - 1)]))

proj <- function(X) { q <- qr.Q(qr(X)); tcrossprod(q) }
imhof_pos <- function(S, delta_vec_fun = NULL, w) {
  ev <- eigen(S, symmetric = TRUE); keep <- abs(ev$values) > 1e-9 * max(abs(ev$values))
  lam <- ev$values[keep]; del <- drop(crossprod(ev$vectors[, keep, drop = FALSE], w))
  f <- function(u) sapply(u, function(uu) {
    th <- 0.5 * sum(atan(lam * uu) + del^2 * lam * uu / (1 + lam^2 * uu^2))
    lr <- 0.25 * sum(log1p(lam^2 * uu^2)) + 0.5 * sum(del^2 * lam^2 * uu^2 / (1 + lam^2 * uu^2))
    sin(th) / (uu * exp(lr)) })
  0.5 + integrate(f, 0, Inf, subdivisions = 2000L, rel.tol = 1e-8, abs.tol = 1e-10, stop.on.error = FALSE)$value / pi
}
# exact P(F > crit) for F = (y' (Pn) y / df1) / (y' (I - Pd) y / df2), y = mu + diag(sd) z
exact_p <- function(Pn, Pd, df1, df2, mu, sdv) {
  N <- length(mu); crit <- qf(1 - ALPHA, df1, df2)
  B <- Pn - (crit * df1 / df2) * (diag(N) - Pd)
  Sh <- sdv; S <- (Sh * B) * rep(Sh, each = N)
  imhof_pos(S, w = mu / sdv)
}
exact_all <- function(d, nx = 30) {
  gr <- build_groups(d); a <- d$a; A <- rep(gr$g, gr$n); N <- length(A)
  mxv <- rep(gr$mx, gr$n); sdv <- rep(gr$sd, gr$n); dA <- a - 1
  sA <- contr.sum(a)[A, , drop = FALSE]; one <- rep(1, N)
  out <- matrix(NA_real_, nx, 6, dimnames = list(NULL, c("Oneway","R_default","TypeII","ANCOVA","TypeIII_raw","TypeIII_ctr")))
  set.seed(d$seed + 99L)
  for (k in seq_len(nx)) {
    xc0 <- mxv + rnorm(N); xr <- gr$L[1] + xc0; xc <- xr - mean(xr)
    mu <- rep(gr$beta, gr$n) * xc0
    Pf <- proj(cbind(one, sA, xc, sA * xc)); dfF <- N - 2 * a
    P1 <- proj(cbind(one)); PA <- proj(cbind(one, sA)); PX <- proj(cbind(one, xc)); Padd <- proj(cbind(one, xc, sA))
    Praw <- proj(cbind(one, xr, sA * xr)); Pctr <- proj(cbind(one, xc, sA * xc))
    out[k, ] <- c(exact_p(PA - P1, PA, dA, N - a, mu, sdv),
                  exact_p(PA - P1, Pf, dA, dfF, mu, sdv),
                  exact_p(Padd - PX, Pf, dA, dfF, mu, sdv),
                  exact_p(Padd - PX, Padd, dA, N - a - 1, mu, sdv),
                  exact_p(Pf - Praw, Pf, dA, dfF, mu, sdv),
                  exact_p(Pf - Pctr, Pf, dA, dfF, mu, sdv))
  }
  colMeans(out)
}
if (identical(Sys.getenv("TH_RUN", "1"), "1")) {
  res <- read.csv("output/ancova_full/ancova_results.csv")
  sel <- subset(res, shape == "normal" & a_effect == 0)
  set.seed(1); sel <- sel[sample(nrow(sel), min(40, nrow(sel))), ]
  designs <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) draw_design(sel$design_id[i], sel$seed[i])))
  nm <- c("Oneway","R_default","TypeII","ANCOVA","TypeIII_raw","TypeIII_ctr")
  ex <- t(sapply(seq_len(nrow(designs)), function(i) exact_all(designs[i, ])))
  sim <- as.matrix(sel[, paste0("rej_", nm)])
  cmp <- data.frame(design_id = sel$design_id, a = sel$a, N = sel$N, exact = round(ex, 4), sim = round(sim, 4))
  write.csv(cmp, "output/theory_vs_sim.csv", row.names = FALSE)
  cat("n designs:", nrow(sel), "\n")
  for (j in seq_along(nm)) cat(sprintf("%-12s mean |exact-sim| %.4f  max %.4f  corr %.4f\n", nm[j],
      mean(abs(ex[, j] - sim[, j])), max(abs(ex[, j] - sim[, j])), cor(ex[, j], sim[, j])))
}
