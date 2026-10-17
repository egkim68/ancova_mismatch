# TABLES 3-6 of the ANCOVA manuscript, rebuilt from the saved result files. Base R only; writes CSV files to <OUTPUT_ROOT>/tables.
# Paste into R or source(). Needs the three result files written by the simulation scripts (or the copies in output/).
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "output")
TAB_DIR <- file.path(OUTPUT_ROOT, "tables"); dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
M <- c("Oneway", "R_default", "TypeII", "ANCOVA", "TypeIII_raw", "TypeIII_ctr", "TypeIII_ctr_HC3", "Drop_int")
LAB <- c("One-way ANOVA", "R default (Type I, A first)", "Type II", "ANCOVA", "Type III, raw",
         "Type III, centered", "Type III, centered + HC3", "Drop interaction if not significant")
SHORT <- setNames(c("One-way ANOVA", "R default", "Type II", "ANCOVA", "Type III, raw", "Type III, centered",
                    "Type III, centered + HC3", "Drop interaction"), M)
SHARE <- c(`2` = .40, `3` = .35, `4` = .25)            # confirmatory shares of 2-, 3-, 4-level designs
f1 <- function(v) sprintf("%.1f", v)
share <- function(d, m) 100 * mean(d[[paste0("rej_", m)]] > 0.10)
harm <- function(d) sapply(M, function(m) share(d, m))
wmean <- function(d, f) { v <- sapply(names(SHARE), function(a) f(d[d$a == as.integer(a), , drop = FALSE])); sum(v * SHARE) / sum(SHARE) }
wshare <- function(d, m) wmean(d, function(g) share(g, m))
wharm <- function(d) sapply(M, function(m) wshare(d, m))
wdis <- function(d) wmean(d, function(g) 100 * mean(g$flip_any))
full <- read.csv(file.path(OUTPUT_ROOT, "ancova_full", "ancova_results.csv")); fn <- subset(full, a_effect == 0)
ext <- read.csv(file.path(OUTPUT_ROOT, "ancova_ext", "ancova_results.csv")); en <- subset(ext, a_effect == 0)
e24 <- subset(en, a <= 4); eh <- subset(e24, design_id %% 4 == 0); fh <- subset(fn, design_id %% 4 == 0)

# ---- Table 3: share of null designs with Type I error above .10 -------------------------------------------------
t3 <- data.frame(Analysis = c(LAB, "Disagreement (% of replications)"),
  Confirmatory = round(c(harm(fn), 100 * mean(fn$flip_any)), 1), Extension_weighted = round(c(wharm(e24), wdis(e24)), 1),
  Extension_heldout = round(c(wharm(eh), wdis(eh)), 1))
names(t3)[2:4] <- sprintf("%s (n = %s)", c("Confirmatory", "Extension, weighted", "Extension, held-out"), prettyNum(c(nrow(fn), nrow(e24), nrow(eh)), big.mark = ","))
write.csv(t3, file.path(TAB_DIR, "Table_3_harmful_share.csv"), row.names = FALSE)

# ---- Table 4: frozen criteria (protocol H1-H7) and observed values --------------------------------------------
# sh(d, m) is the unweighted share for the confirmatory columns and the weighted share for the extension columns.
crit <- function(d, sh, dis) {
  s <- setNames(sapply(M, function(m) sh(d, m)), M)
  imb <- subset(d, slope_range == 0 & cov_imbalance > 0); far <- subset(d, slope_range > 0 & d0 >= 2)
  eq <- subset(d, sd_ratio == 1); ne <- subset(d, sd_ratio > 1)
  r <- round(c(s["R_default"], sh(imb, "Oneway"), sh(imb, "R_default"), sh(imb, "TypeII"), sh(imb, "ANCOVA"),
               sh(far, "TypeIII_raw"), sh(far, "TypeIII_ctr"), sh(eq, "TypeIII_ctr"), sh(ne, "TypeIII_ctr"),
               s["TypeIII_ctr_HC3"], min(s[names(s) != "TypeIII_ctr_HC3"]), dis(d), s["Drop_int"]), 1)
  lead1 <- r[1] - round(max(s[names(s) != "R_default"]), 1)   # lead computed from the displayed (rounded) shares
  c(H1 = sprintf("%s; %s", f1(r[1]), f1(lead1)), H2 = sprintf("%s, %s; %s, %s", f1(r[2]), f1(r[3]), f1(r[4]), f1(r[5])),
    H3 = sprintf("%s; %s", f1(r[6]), f1(r[6] - r[7])), H4 = sprintf("%s; %s", f1(r[8]), f1(r[9])),
    H5 = sprintf("%s; %s", f1(r[10]), f1(r[11])), H6 = f1(r[12]), H7 = f1(r[13]))
}
dis_u <- function(d) 100 * mean(d$flip_any)
t4 <- data.frame(Hypothesis = paste0("H", 1:7),
  Criterion = c("R default >= 50%; lead >= 10 points", "One-way, R default >= 50%; Type II, ANCOVA <= 20%",
                "Raw Type III >= 60%; lead >= 30 points", "Centered <= 5% (equal variances); >= 15% (unequal)",
                "HC3 <= 10%; below every other analysis", "Disagreement >= 25% of replications", "Drop interaction >= 30%"),
  Confirmatory = crit(fn, share, dis_u), Heldout = crit(fh, share, dis_u),
  Extension = crit(e24, wshare, wdis), Extension_heldout = crit(eh, wshare, wdis))
names(t4)[3:6] <- sprintf("%s (n = %s)", c("Confirmatory", "Held-out", "Extension", "Extension held-out"), prettyNum(c(nrow(fn), nrow(fh), nrow(e24), nrow(eh)), big.mark = ","))
write.csv(t4, file.path(TAB_DIR, "Table_4_criteria.csv"), row.names = FALSE)

# ---- Table 5: mean Type I error as group sizes are multiplied --------------------------------------------------
sc <- read.csv(file.path(OUTPUT_ROOT, "ancova_scaling", "scaling_results.csv"))
pick <- list(c("set1_raw_mismatch", "TypeIII_raw"), c("set1_raw_mismatch", "TypeII"), c("set1_raw_mismatch", "ANCOVA"),
             c("set1_raw_mismatch", "TypeIII_ctr"), c("set2_cov_imbalance", "Oneway"), c("set2_cov_imbalance", "R_default"),
             c("set2_cov_imbalance", "TypeII"), c("set2_cov_imbalance", "ANCOVA"), c("set2_cov_imbalance", "TypeIII_ctr"),
             c("set3_variance", "TypeIII_ctr"), c("set3_variance", "TypeIII_ctr_HC3"))
t5 <- do.call(rbind, lapply(pick, function(p) { d <- subset(sc, set == p[1]); y <- tapply(d[[paste0("rej_", p[2])]], d$factor, mean)
  data.frame(Set_and_analysis = paste(sub("set(\\d).*", "Set \\1", p[1]), SHORT[p[2]], sep = ", "), t(round(y + 1e-10, 3))) }))
names(t5)[-1] <- paste0("x", c(1, 2, 4, 8, 16)); write.csv(t5, file.path(TAB_DIR, "Table_5_scaling.csv"), row.names = FALSE)

# ---- Table 6: centered Type III with and without HC3 by number of levels (extension, unweighted) ------------------
t6 <- do.call(rbind, lapply(2:5, function(a) { g <- en[en$a == a, ]
  data.frame(Levels = paste(a, "levels"), Null_designs = nrow(g), TypeIII_centered = round(share(g, "TypeIII_ctr"), 1),
             TypeIII_centered_HC3 = round(share(g, "TypeIII_ctr_HC3"), 1)) }))
write.csv(t6, file.path(TAB_DIR, "Table_6_HC3_by_levels.csv"), row.names = FALSE)

print(t3); print(t4); print(t5); print(t6); cat("Tables written to", normalizePath(TAB_DIR), "\n")
