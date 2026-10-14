# TABLES 3-5 of the ANCOVA manuscript, rebuilt from the saved result files. Base R only; writes CSV files to <OUTPUT_ROOT>/tables.
# Paste into R or source(). Needs the three result files written by the simulation scripts (or the copies in results/).
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-ancova")
TAB_DIR <- file.path(OUTPUT_ROOT, "tables"); dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
M <- c("Oneway", "R_default", "TypeII", "ANCOVA", "TypeIII_raw", "TypeIII_ctr", "TypeIII_ctr_HC3", "Drop_int")
LAB <- c("One-way ANOVA", "R default (Type I, A first)", "Type II", "ANCOVA", "Type III, raw covariate",
         "Type III, centered", "Centered + HC3", "Drop interaction if n.s.")
SHARE <- c(`2` = .40, `3` = .35, `4` = .25)
harm <- function(d) sapply(M, function(m) 100 * mean(d[[paste0("rej_", m)]] > 0.10))
wmean <- function(d, f) { v <- sapply(names(SHARE), function(a) f(d[d$a == as.integer(a), , drop = FALSE])); sum(v * SHARE) / sum(SHARE) }
wharm <- function(d) sapply(M, function(m) wmean(d, function(g) 100 * mean(g[[paste0("rej_", m)]] > 0.10)))
wdis <- function(d) wmean(d, function(g) 100 * mean(g$flip_any))
full <- read.csv(file.path(OUTPUT_ROOT, "ancova_full", "ancova_results.csv")); fn <- subset(full, a_effect == 0)
ext <- read.csv(file.path(OUTPUT_ROOT, "ancova_ext", "ancova_results.csv")); en <- subset(ext, a_effect == 0)
eh <- subset(en, design_id %% 4 == 0)
t3 <- data.frame(Analysis = c(LAB, "Disagreement (% of replications)"),
  Confirmatory = round(c(harm(fn), 100 * mean(fn$flip_any)), 1), Extension_weighted = round(c(wharm(en), wdis(en)), 1),
  Extension_heldout = round(c(wharm(eh), wdis(eh)), 1))
names(t3)[2:4] <- sprintf("%s (n = %s)", c("Confirmatory", "Extension, weighted", "Extension, held-out"), format(c(nrow(fn), sum(en$a <= 4), sum(eh$a <= 4)), big.mark = ","))
write.csv(t3, file.path(TAB_DIR, "Table_3_harmful_share.csv"), row.names = FALSE)
sc <- read.csv(file.path(OUTPUT_ROOT, "ancova_scaling", "scaling_results.csv"))
pick <- list(c("set1_raw_mismatch", "TypeIII_raw"), c("set1_raw_mismatch", "TypeII"), c("set1_raw_mismatch", "TypeIII_ctr"),
             c("set2_cov_imbalance", "Oneway"), c("set2_cov_imbalance", "R_default"), c("set2_cov_imbalance", "TypeII"),
             c("set3_variance", "TypeIII_ctr"), c("set3_variance", "TypeIII_ctr_HC3"))
t4 <- do.call(rbind, lapply(pick, function(p) { d <- subset(sc, set == p[1]); y <- tapply(d[[paste0("rej_", p[2])]], d$factor, mean)
  data.frame(Set_and_analysis = paste(sub("set(\\d).*", "Set \\1", p[1]), LAB[match(p[2], M)], sep = ", "), t(round(y + 1e-10, 3))) }))
names(t4)[-1] <- paste0("x", c(1, 2, 4, 8, 16)); write.csv(t4, file.path(TAB_DIR, "Table_4_scaling.csv"), row.names = FALSE)
bins <- cut(en$min_group, c(0, 3, 5, 9, Inf), labels = c("Smallest group 3", "Smallest group 4 to 5", "Smallest group 6 to 9", "Smallest group 10 or more"))
row <- function(g, lab) data.frame(Group = lab, Null_designs = nrow(g), Centered_classical = round(100 * mean(g$rej_TypeIII_ctr > 0.10), 1),
                                   Centered_HC3 = round(100 * mean(g$rej_TypeIII_ctr_HC3 > 0.10), 1))
t5 <- rbind(do.call(rbind, lapply(2:5, function(a) row(en[en$a == a, ], paste(a, "levels")))),
            do.call(rbind, lapply(levels(bins), function(b) row(en[bins == b, ], b))))
write.csv(t5, file.path(TAB_DIR, "Table_5_HC3_by_layout_and_group.csv"), row.names = FALSE)
print(t3); print(t4); print(t5); cat("Tables written to", normalizePath(TAB_DIR), "\n")
