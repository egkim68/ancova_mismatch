# Independent checks (seed 11): (1) the Type II contrast equals the group-line difference at t; (2) Imhof errors vs direct Monte Carlo from nested least-squares fits (3 levels, unequal variances and slopes).
ALPHA <- .05
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
set.seed(11)
# ---- Part 1: Type II contrast (two groups) vs weighted covariate point t
n <- c(14, 22); N <- sum(n); g <- rep(1:2, n)
x <- c(rnorm(n[1], 0.5), rnorm(n[2], -0.3)) ; al <- c(1, 0.4); be <- c(1.2, 0.2)
Ey <- al[g] + be[g] * x
fit <- lm(Ey ~ x + factor(g))                      # additive model fitted to noise-free means
d_add <- -coef(fit)[3]                              # group1 - group2 adjusted difference (sign)
Sxx <- tapply(x, g, function(v) sum((v - mean(v))^2)); w <- Sxx / sum(Sxx)
xb <- tapply(x, g, mean); tpt <- w[2] * xb[1] + w[1] * xb[2]
d_form <- (al[1] - al[2]) + (be[1] - be[2]) * tpt
cat("Type II additive difference:", round(unname(coef(fit)[3]), 6), " formula (g2-g1):", round(-d_form, 6), "\n")
# ---- Part 2: independent Monte Carlo vs Imhof, a=3, unequal sd and slopes, N=54
n <- c(10, 16, 28); N <- sum(n); a <- 3; g <- rep(1:3, n); sd <- c(2.2, 1, 0.6)[g]
x0 <- rnorm(N) + c(.6, 0, -.5)[g]; xr <- x0 + 2; xc <- xr - mean(xr); be <- c(1.0, 0.2, -0.6)
mu <- be[g] * x0
sA <- contr.sum(a)[g, ]; one <- rep(1, N)
Pf <- proj(cbind(one, sA, xc, sA * xc)); P1 <- proj(cbind(one)); PA <- proj(cbind(one, sA))
PX <- proj(cbind(one, xc)); Padd <- proj(cbind(one, xc, sA)); Praw <- proj(cbind(one, xr, sA * xr)); Pctr <- proj(cbind(one, xc, sA * xc))
ex <- c(TypeII = exact_p(Padd - PX, Pf, a - 1, N - 2 * a, mu, sd),
        raw = exact_p(Pf - Praw, Pf, a - 1, N - 2 * a, mu, sd),
        ctr = exact_p(Pf - Pctr, Pf, a - 1, N - 2 * a, mu, sd))
R <- 40000; rej <- c(0, 0, 0)
Xf <- cbind(1, sA, xc, sA * xc); Xadd <- cbind(1, xc, sA); Xc <- cbind(1, xc); Xr <- cbind(1, xr, sA * xr)
rss <- function(X, Y) colSums(qr.resid(qr(X), Y)^2)
Y <- mu + sd * matrix(rnorm(N * R), N, R)
Ff <- rss(Xf, Y); dff <- N - 2 * a
fII <- ((rss(Xc, Y) - rss(Xadd, Y)) / (a - 1)) / (Ff / dff)
# raw type III: compare full vs model without A at x=0 : drop A main effect with raw x => columns (1, xr, sA*xr)
fr <- ((rss(Xr, Y) - Ff) / (a - 1)) / (Ff / dff)
Xct <- cbind(1, xc, sA * xc)
fc <- ((rss(Xct, Y) - Ff) / (a - 1)) / (Ff / dff)
cr <- qf(.95, a - 1, dff)
sim <- c(TypeII = mean(fII > cr), raw = mean(fr > cr), ctr = mean(fc > cr))
print(rbind(imhof = round(ex, 4), mc = round(sim, 4), mc_se = round(sqrt(sim * (1 - sim) / R), 4)))
