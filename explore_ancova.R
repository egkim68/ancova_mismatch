# EXPLORER: one categorical factor A and one numeric covariate x, with a possible A-by-x interaction.
# Base R only. Writes CSV files only (figures come from a separate script). Paste into R or source().
# Output folders are created automatically. Finished designs are checkpointed, so an interrupted run
# resumes when you run it again. See docs/PROTOCOL_ancova_v0.md for the question and the analyses.
#
# Each replication draws one data set and tests the effect of A in eight ways:
#   Oneway          one-way ANOVA on A, ignoring x
#   R_default       anova(lm(y ~ A*x)): sequential, A first (numerator ignores x)
#   TypeII          A adjusted for x, interaction ignored, error from the full model
#   ANCOVA          y ~ x + A, error from the additive model
#   TypeIII_raw     Type III with x on its raw scale: equality of the group lines at x = 0
#   TypeIII_ctr     Type III with x centered at the sample grand mean (targets the intended null)
#   TypeIII_ctr_HC3 the previous test with HC3 standard errors
#   Drop_int        test the interaction; if p > .05 use ANCOVA, otherwise TypeIII_raw
# The null is "equal adjusted means at the population grand mean of x" (a_effect = 0).

MODE <- Sys.getenv("SIM_MODE", unset = "pilot")             # "pilot", "medium", "full"
N_DESIGNS <- switch(MODE, pilot = 160L, medium = 400L, full = 2000L, 20L)
N_REPS    <- switch(MODE, pilot = 200L, medium = 200L, full = 1000L, 100L)
BASE_SEED <- as.integer(Sys.getenv("SIM_SEED", unset = "7301"))   # exploratory 7301; confirmatory 8419
ALPHA <- 0.05
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-ancova")
OUT_DIR <- file.path(OUTPUT_ROOT, paste0("ancova_", MODE)); CK_DIR <- file.path(OUT_DIR, "checkpoints")
dir.create(CK_DIR, recursive = TRUE, showWarnings = FALSE)
METHODS <- c("Oneway", "R_default", "TypeII", "ANCOVA", "TypeIII_raw", "TypeIII_ctr", "TypeIII_ctr_HC3", "Drop_int")

# ---- design generator -------------------------------------------------------------------------------
draw_design <- function(id, seed) {
  set.seed(seed)
  a <- sample(2:4, 1, prob = c(0.40, 0.35, 0.25))
  N <- round(exp(runif(1, log(40), log(300))))
  imb <- if (runif(1) < 0.2) 0 else runif(1, 0.2, 1.2)
  sdr <- if (runif(1) < 0.4) 1 else exp(runif(1, log(1.5), log(4)))
  pairing <- sample(c(-1, 0, 1), 1)
  shape <- sample(c("normal", "skewed", "heavy"), 1)
  cov_imb <- if (runif(1) < 0.3) 0 else runif(1, 0.2, 1.0)           # range of group means of x (SD units)
  d0 <- if (runif(1) < 0.25) 0 else exp(runif(1, log(0.5), log(6)))   # grand mean of x, distance from 0 (SD units)
  slope <- runif(1, 0.2, 0.8)
  s_range <- if (runif(1) < 0.4) 0 else runif(1, 0.1, 0.6)           # slope range across groups = interaction
  a_eff <- if (runif(1) < 0.5) 0 else runif(1, 0.2, 0.8)             # range of adjusted means
  data.frame(design_id = id, seed = seed, a = a, N_target = N, imbalance = imb, sd_ratio = sdr, pairing = pairing,
             shape = shape, cov_imbalance = cov_imb, d0 = d0, slope = slope, slope_range = s_range,
             a_effect = a_eff, stringsAsFactors = FALSE)
}

spread <- function(k) { v <- runif(k); if (k > 1 && diff(range(v)) > 0) (v - min(v)) / diff(range(v)) else rep(0, k) }
build_groups <- function(d) {
  set.seed(d$seed + 1L); a <- d$a
  w <- exp(d$imbalance * rnorm(a)); n <- pmax(3L, round(d$N_target * w / sum(w))); pi <- n / sum(n)
  u <- runif(a); lsd <- if (d$sd_ratio == 1) rep(0, a) else (u - min(u)) / (max(u) - min(u)) * log(d$sd_ratio)
  if (d$pairing != 0 && d$sd_ratio > 1) { ord <- order(n); lsd <- sort(lsd, decreasing = d$pairing < 0)[order(ord)] }
  sdv <- exp(lsd); sdv <- sdv / exp(mean(log(sdv)))
  cen <- function(v) v - sum(pi * v)                                   # zero weighted mean: grand mean of x stays at L
  mx <- cen(d$cov_imbalance * spread(a))
  tau <- if (d$a_effect == 0) rep(0, a) else d$a_effect * spread(a)    # adjusted means at the grand mean
  sb <- spread(a); beta <- d$slope + d$slope_range * (sb - mean(sb))
  L <- d$d0 * sample(c(-1, 1), 1)
  data.frame(g = seq_len(a), n = n, sd = sdv, mx = mx, tau = tau, beta = beta, L = L)
}

draw_err <- function(n, shape) switch(shape,
  normal = rnorm(n),
  skewed = { s <- 0.8; (exp(s * rnorm(n)) - exp(s^2 / 2)) / sqrt((exp(s^2) - 1) * exp(s^2)) },
  heavy  = rt(n, 3) / sqrt(3))

# ---- tests ------------------------------------------------------------------------------------------
rss <- function(X, y) sum(lm.fit(X, y)$residuals^2)
pF <- function(ss, df1, rssE, df2) pf((ss / df1) / (rssE / df2), df1, df2, lower.tail = FALSE)

test_all <- function(y, A, xr, a) {
  N <- length(y); dA <- a - 1
  sA <- contr.sum(a)[A, , drop = FALSE]; one <- rep(1, N); xc <- xr - mean(xr)
  sAxc <- sA * xc; sAxr <- sA * xr                                     # (N x (a-1)) interaction columns
  Xf <- cbind(one, sA, xc, sAxc); fit <- lm.fit(Xf, y); e <- fit$residuals
  rssF <- sum(e^2); dfF <- N - 2 * a
  rss1 <- rss(cbind(one), y); rssA <- rss(cbind(one, sA), y); rssX <- rss(cbind(one, xc), y); rssAdd <- rss(cbind(one, xc, sA), y)
  p <- numeric(length(METHODS)); names(p) <- METHODS
  p["Oneway"] <- pF(rss1 - rssA, dA, rssA, N - a)
  p["R_default"] <- pF(rss1 - rssA, dA, rssF, dfF)
  p["TypeII"] <- pF(rssX - rssAdd, dA, rssF, dfF)
  p["ANCOVA"] <- pF(rssX - rssAdd, dA, rssAdd, N - a - 1)
  p["TypeIII_raw"] <- pF(rss(cbind(one, xr, sAxr), y) - rssF, dA, rssF, dfF)
  p["TypeIII_ctr"] <- pF(rss(cbind(one, xc, sAxc), y) - rssF, dA, rssF, dfF)
  Q <- qr.Q(fit$qr); h <- rowSums(Q^2); XtXi <- chol2inv(qr.R(fit$qr))
  piv <- fit$qr$pivot; XtXi <- XtXi[order(piv), order(piv)]
  meat <- crossprod(Xf * (e / (1 - pmin(h, 0.9999))))
  V <- XtXi %*% meat %*% XtXi; idx <- 1 + seq_len(dA); bA <- fit$coefficients[idx]
  W <- tryCatch(drop(t(bA) %*% solve(V[idx, idx, drop = FALSE], bA)), error = function(err) NA_real_)
  p["TypeIII_ctr_HC3"] <- pf(W / dA, dA, dfF, lower.tail = FALSE)
  p_int <- pF(rssAdd - rssF, dA, rssF, dfF)
  p["Drop_int"] <- if (p_int > ALPHA) p["ANCOVA"] else p["TypeIII_raw"]
  c(p, p_int = p_int)
}

run_design <- function(d) {
  gr <- build_groups(d); a <- d$a
  A <- rep(gr$g, gr$n); mxv <- rep(gr$mx, gr$n); sdv <- rep(gr$sd, gr$n); N <- length(A)
  set.seed(d$seed + 2L)
  P <- matrix(NA_real_, N_REPS, length(METHODS) + 1)
  for (r in seq_len(N_REPS)) {
    xc0 <- mxv + rnorm(N)                                              # within-group SD 1; population grand mean 0
    xr <- gr$L[1] + xc0                                                   # raw scale: grand mean L (d0 SD from 0)
    y <- rep(gr$tau, gr$n) + rep(gr$beta, gr$n) * xc0 + sdv * draw_err(N, d$shape)
    P[r, ] <- test_all(y, A, xr, a)
  }
  sig <- P[, seq_along(METHODS)] < ALPHA
  rej <- colMeans(sig, na.rm = TRUE); names(rej) <- paste0("rej_", METHODS)
  flip_any <- mean(apply(sig[, 1:7], 1, function(s) length(unique(s)) > 1), na.rm = TRUE)
  data.frame(d, N = N, min_group = min(gr$n), group_ratio = max(gr$n) / min(gr$n), as.list(rej),
             rej_interaction = mean(P[, ncol(P)] < ALPHA), flip_any = flip_any, check.names = FALSE)
}

# ---- run ----------------------------------------------------------------------------------------------
designs <- do.call(rbind, lapply(seq_len(N_DESIGNS), function(i) draw_design(i, BASE_SEED * 1000L + i)))
if (identical(Sys.getenv("SIM_RUN", "1"), "1")) {
  t0 <- Sys.time()
  for (i in seq_len(N_DESIGNS)) {
    f <- file.path(CK_DIR, sprintf("design_%05d.rds", i))
    if (file.exists(f)) next
    saveRDS(run_design(designs[i, ]), f)
    if (i %% 20 == 0) cat(sprintf("%d / %d designs (%.1f min)\n", i, N_DESIGNS, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  res <- do.call(rbind, lapply(sprintf("design_%05d.rds", seq_len(N_DESIGNS)), function(f) readRDS(file.path(CK_DIR, f))))
  write.csv(res, file.path(OUT_DIR, "ancova_results.csv"), row.names = FALSE)
  cat("Written:", normalizePath(file.path(OUT_DIR, "ancova_results.csv")), "\n")
}
