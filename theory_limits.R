Sys.setenv(TH_RUN="0"); source("theory_ancova.R")
NX <- as.integer(Sys.getenv("TL_NX", "60"))   # covariate draws per design and multiplier
res <- read.csv("output/ancova_full/ancova_results.csv")
pick <- function(cond) { r <- res[cond, ][1, ]; draw_design(r$design_id, r$seed) }
ok <- res$shape=="normal" & res$a_effect==0 & res$N_target<70
sets <- list(
 raw_mismatch = pick(ok & res$slope_range>0 & res$d0>=2 & res$sd_ratio==1 & res$a==2),
 cov_imbalance = pick(ok & res$slope_range==0 & res$cov_imbalance>0.5 & res$sd_ratio==1 & res$a==2),
 variance = pick(ok & res$sd_ratio>=2 & res$pairing==-1 & res$imbalance>=.4 & res$a==2))
run <- function(d, m) {
  gr <- build_groups(d); gr$n <- gr$n * m; a <- d$a; A <- rep(gr$g, gr$n); N <- length(A)
  mxv <- rep(gr$mx, gr$n); sdv <- rep(gr$sd, gr$n); one <- rep(1,N); sA <- contr.sum(a)[A,,drop=FALSE]; dA <- a-1
  set.seed(d$seed+7); nx <- NX; out <- matrix(NA,nx,6)
  for (k in 1:nx) { xc0 <- mxv+rnorm(N); xr <- gr$L[1]+xc0; xc <- xr-mean(xr); mu <- rep(gr$beta,gr$n)*xc0; dfF <- N-2*a
    Pf <- proj(cbind(one,sA,xc,sA*xc)); P1 <- proj(one); PA <- proj(cbind(one,sA)); PX <- proj(cbind(one,xc)); Padd <- proj(cbind(one,xc,sA))
    out[k,] <- c(exact_p(PA-P1,PA,dA,N-a,mu,sdv), exact_p(PA-P1,Pf,dA,dfF,mu,sdv), exact_p(Padd-PX,Pf,dA,dfF,mu,sdv),
                 exact_p(Padd-PX,Padd,dA,N-a-1,mu,sdv), exact_p(Pf-proj(cbind(one,xr,sA*xr)),Pf,dA,dfF,mu,sdv), exact_p(Pf-proj(cbind(one,xc,sA*xc)),Pf,dA,dfF,mu,sdv)) }
  c(colMeans(out), apply(out, 2, sd) / sqrt(nx)) }
nm <- c("Oneway","R_default","TypeII","ANCOVA","raw","ctr"); allres <- list()
for (s in names(sets)) { cat("\n", s, " (N at 1x =", sum(build_groups(sets[[s]])$n), ")\n")
  tab <- t(sapply(c(1,2,4,8), function(m) run(sets[[s]], m))); rownames(tab) <- paste0(c(1,2,4,8),"x"); colnames(tab) <- c(nm, paste0("se_", nm))
  print(round(tab, 3)); allres[[s]] <- data.frame(set = s, mult = rownames(tab), round(tab, 4)) }
dir.create("output", showWarnings = FALSE); write.csv(do.call(rbind, allres), "output/theory_limits.csv", row.names = FALSE)
