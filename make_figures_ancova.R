# FIGURES for the ANCOVA study. Base R only. Reads the CSV files written by explore_ancova.R and
# scaling_ancova.R and writes image files (not PDF) named Fig1 to Fig5 to <OUTPUT_ROOT>/figures. Independent of the
# simulation: rerun any time to change sizes or the harmful line. Paste into R or source().
# Journal layout: TIFF at 600 dpi (FIG_FORMAT = "png" writes PNG instead), sans-serif lettering (Arial on Windows),
# no titles inside the artwork, panels lettered a, b, c.
# Colors: eight categorical slots in a fixed order, one per analysis (validated: lightness band, chroma,
# adjacent colorblind separation, normal-vision separation pass; low-contrast slots are directly labeled).

MODE <- Sys.getenv("SIM_MODE", unset = "full")
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "N:/paper-simulation-ancova")
FIG_DIR <- file.path(OUTPUT_ROOT, "figures"); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
HARM <- 0.10; DPI <- as.integer(Sys.getenv("FIG_DPI", unset = "600"))
FMT <- tolower(Sys.getenv("FIG_FORMAT", unset = "tiff")); FONT <- Sys.getenv("FIG_FONT", unset = "sans")
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"
ANALYSES <- c(Oneway = "One-way ANOVA", R_default = "R default (Type I, A first)", TypeII = "Type II",
              ANCOVA = "ANCOVA (no interaction)", TypeIII_raw = "Type III, raw covariate",
              TypeIII_ctr = "Type III, centered covariate", TypeIII_ctr_HC3 = "Type III, centered + HC3",
              Drop_int = "Drop interaction if not significant")
COL <- c(Oneway = "#6A3D9A", R_default = "#D55E00", TypeII = "#56B4E9", ANCOVA = "#E69F00", TypeIII_raw = "#CC79A7",
         TypeIII_ctr = "#0072B2", TypeIII_ctr_HC3 = "#009E73", Drop_int = "#000000")
PCH <- c(Oneway = 16, R_default = 15, TypeII = 17, ANCOVA = 18, TypeIII_raw = 16, TypeIII_ctr = 15, TypeIII_ctr_HC3 = 17, Drop_int = 18)

r <- read.csv(file.path(OUTPUT_ROOT, paste0("ancova_", MODE), "ancova_results.csv"), stringsAsFactors = FALSE)
null <- subset(r, a_effect == 0)
rej <- function(d, m) d[[paste0("rej_", m)]]
png_open <- function(name, w, h) {
  f <- file.path(FIG_DIR, paste0(name, if (FMT == "png") ".png" else ".tif")); ty <- if (capabilities("cairo")) "cairo" else getOption("bitmapType")
  if (FMT == "png") png(f, width = w, height = h, units = "in", res = DPI, type = ty, family = FONT)
  else tiff(f, width = w, height = h, units = "in", res = DPI, compression = "lzw", type = ty, family = FONT) }
dodge <- function(y, gap = 0.075) { o <- order(y); z <- y[o]; for (i in seq_along(z)[-1]) if (z[i] - z[i - 1] < gap) z[i] <- z[i - 1] + gap; z[order(o)] }
base_par <- function(...) par(las = 1, col.axis = INK2, col.lab = INK, col.main = INK, fg = INK2, cex.axis = 1.05, cex.lab = 1.2, ...)

# Theme: soft gray panel, white grid lines, no frame, no tick marks (applied to the line, scatter, and bar graphs)
th_plot <- function(xlim, ylim, xlab = "", ylab = "", xat = NULL, xlabs = xat, yat = seq(0, 1, 0.2), ylabs = yat, xaxs = "r", yline = 3.6, xline = 3.2) {
  plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "", xaxs = xaxs)
  if (!is.null(xat)) axis(1, at = xat, labels = xlabs, col = INK, col.axis = INK, lwd = 1.6)
  axis(2, at = yat, labels = ylabs, col = INK, col.axis = INK, lwd = 1.6)
  mtext(ylab, side = 2, line = yline, las = 0, col = INK, cex = 0.92); mtext(xlab, side = 1, line = xline, col = INK, cex = 0.92)
}
# ---- Figure 1: share of harmful designs by situation ----------------------------------------------------
groups <- list("All null\ndesigns" = rep(TRUE, nrow(null)),
               "No interaction,\ncovariate imbalance" = with(null, slope_range == 0 & cov_imbalance > 0),
               "Interaction,\ngrand mean far from 0" = with(null, slope_range > 0 & d0 >= 2),
               "Equal\nvariances" = with(null, sd_ratio == 1), "Unequal\nvariances" = with(null, sd_ratio > 1))
m <- sapply(groups, function(g) sapply(names(ANALYSES), function(a) mean(rej(null[g, ], a) > HARM)))
m <- m[order(-m[, 1]), , drop = FALSE]  # rows ordered from most to least often harmful (all null designs)
png_open("Fig2", 7.6, 5.7)
# heatmap: one hue (blue), light = few harmful designs, dark = many; every cell carries its number
SEQ <- c("#FFF3DB", "#FFD08A", "#FFA040", "#F26B0F", "#C23E00", "#8A2800", "#4F1500")
ramp <- colorRampPalette(SEQ, space = "Lab")(101)
nr <- nrow(m); nc <- ncol(m)
# text color with the higher contrast against the cell (white or near-black); contrast stays at 4.5:1 or above
rel_lum <- function(h) { v <- col2rgb(h)[, 1] / 255; v <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4); sum(v * c(0.2126, 0.7152, 0.0722)) }
txtcol <- function(h) { L <- rel_lum(h); if ((1.05) / (L + 0.05) > (L + 0.05) / (rel_lum(INK) + 0.05)) "white" else INK }
par(mar = c(0.4, 9.8, 8.6, 4.4), las = 1, fg = INK, col.axis = INK, ps = 11)
plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(nr + 0.5, 0.5), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
gap <- 0.04
for (i in seq_len(nr)) for (j in seq_len(nc)) {
  v <- m[i, j]
  rect(j - 0.5 + gap, i - 0.5 + gap, j + 0.5 - gap, i + 0.5 - gap, col = ramp[round(v * 100) + 1], border = NA)
  text(j, i, sprintf("%.0f%%", 100 * v), cex = 1.15, col = txtcol(ramp[round(v * 100) + 1]), font = 2)
}
for (i in seq_len(nr)) { a <- rownames(m)[i]
  if (a == "R_default") text(0.4, i, expression(atop("R default (Type I,", italic(A) * " first)")), adj = c(1, 0.5), cex = 0.88, col = INK, xpd = NA)
  else text(0.4, i, paste(strwrap(ANALYSES[a], 22), collapse = "\n"), adj = c(1, 0.5), cex = 0.88, col = INK, xpd = NA) }
HL <- list("All null\ndesigns", "No\ninteraction,\ncovariate\nimbalance", "Interaction,\ngrand mean\nfar from 0", "Equal\nvariances", "Unequal\nvariances")
UN <- par("pin")[2] / nr; LS <- 0.2 / UN; NY <- 0.5 - 0.2 / UN  # header spacing set in inches so it holds for any row height
for (j in seq_len(nc)) text(j, NY, bquote(italic(n) == .(sum(groups[[j]]))), adj = c(0.5, 0.5), cex = 0.85, col = INK, xpd = NA)
for (j in seq_len(nc)) { ln <- strsplit(HL[[j]], "\n")[[1]]; for (k in seq_along(ln)) text(j, NY - LS * (length(ln) - k + 1), ln[k], adj = c(0.5, 0.5), cex = 0.9, col = INK, xpd = NA) }
text(nc + 0.7, NY - LS * 0.5, "Harmful\ndesigns", adj = c(0, 0.5), cex = 0.9, col = INK, xpd = NA)
# color key
kx <- nc + 0.7; ky <- seq(1, nr, length.out = 101)
for (k in 1:100) rect(kx, ky[k] - 0.5 * (nr - 1) / 100 - 0.001, kx + 0.12, ky[k + 1] - 0.5 * (nr - 1) / 100 + 0.001, col = ramp[k], border = NA, xpd = NA)
text(kx + 0.17, c(1, nr), c("0%", "100%"), adj = 0, cex = 0.9, col = INK, xpd = NA)
dev.off()

# ---- Figure 3: Type I error against the distance of the grand mean of x from zero (large fonts, black text, no grey) ----
s <- subset(null, slope_range > 0); show <- c("TypeIII_raw", "TypeIII_ctr", "TypeIII_ctr_HC3")
png_open("Fig3", 6.5, 5.8); par(las = 1, ps = 12, fg = INK, col.lab = INK, col.axis = INK, mar = c(4.6, 5.0, 0.4, 0.5), cex.lab = 1, cex.axis = 1)
th_plot(c(0, 6), c(0, 1), expression("Distance of the grand mean of " * italic(x) * " from 0 (standard deviation units)"), expression("Type I error for the effect of " * italic(A)), xat = 0:6)
abline(h = c(0.05, HARM), lty = c(2, 3), col = INK, lwd = 1.4)
lo <- lapply(show, function(a) lowess(s$d0, rej(s, a), f = 0.5)); names(lo) <- show
for (a in show) {
  points(s$d0, rej(s, a), pch = PCH[a], cex = 0.8, col = adjustcolor(COL[a], 0.6)); lines(lo[[a]], col = COL[a], lwd = 5.5)
}
legend("topleft", legend = c("Type III, raw covariate", "Type III, centered", "Type III, centered + HC3"), col = COL[c("TypeIII_raw", "TypeIII_ctr", "TypeIII_ctr_HC3")],
       lwd = 5.5, bty = "n", cex = 1, text.col = INK, seg.len = 1.6, y.intersp = 1.35)
dev.off()

# ---- Figure 3: scaling experiment -------------------------------------------------------------------------
sc <- read.csv(file.path(OUTPUT_ROOT, "ancova_scaling", "scaling_results.csv"), stringsAsFactors = FALSE)
SHORT <- c(Oneway = "One-way", R_default = "R default", TypeII = "Type II", ANCOVA = "ANCOVA", TypeIII_raw = "Type III,\nraw",
           TypeIII_ctr = "Type III,\ncentered", TypeIII_ctr_HC3 = "Centered\n+ HC3", Drop_int = "Drop\ninteraction")
panels <- list(set1_raw_mismatch = list(title = "A. Interaction, grand mean far from 0,\nequal variances", show = c("TypeIII_raw", "TypeII", "TypeIII_ctr")),
               set2_cov_imbalance = list(title = "B. Covariate imbalance, no interaction,\nequal variances", show = c("R_default", "Oneway", "TypeII")),
               set3_variance = list(title = "C. Unequal variances, small groups\nwith large spread", show = c("TypeIII_ctr", "TypeIII_ctr_HC3")))
png_open("Fig4", 6.5, 9.4); par(las = 1, ps = 19, fg = INK, col.lab = INK, col.axis = INK, cex.lab = 1, cex.axis = 1); layout(matrix(1:3, ncol = 1), heights = c(1, 1, 1.2))
ptitle <- c(set1_raw_mismatch = "a", set2_cov_imbalance = "b", set3_variance = "c")  # panel letters only; the panel descriptions are in the figure note
SHORT <- c(Oneway = "One-way", R_default = "R default", TypeII = "Type II", ANCOVA = "ANCOVA", TypeIII_raw = "Type III, raw",
           TypeIII_ctr = "Type III, centered", TypeIII_ctr_HC3 = "Centered + HC3", Drop_int = "Drop interaction")
for (nm in names(panels)) {
  p <- panels[[nm]]; d <- subset(sc, set == nm); lastp <- nm == names(panels)[3]
  par(mar = c(if (lastp) 4.8 else 3.4, 6.0, 1.6, 0.5))
  YL <- if (nm == "set3_variance") 0.4 else 1
  th_plot(c(0, 4), c(0, YL), if (lastp) "Group sizes multiplied by" else "", "Mean Type I error", xat = 0:4, yline = 4.0, xline = 3.4,
          xlabs = if (lastp) c(1, 2, 4, 8, 16) else rep("", 5), yat = if (YL == 1) seq(0, 1, 0.2) else seq(0, 0.4, 0.1))
  mtext(ptitle[nm], side = 3, line = 0.2, adj = 0, font = 2, cex = 1.1, col = INK)
  segments(0, 0.05, 4, 0.05, lty = 2, col = INK, lwd = 1.4)
  for (k in seq_along(p$show)) { a <- p$show[k]; y <- tapply(rej(d, a), d$factor, mean); lines(0:4, y, col = COL[a], lwd = 4.5); points(0:4, y, pch = PCH[a], col = COL[a], cex = 1.7) }
  legend("topleft", inset = c(0.12, 0.02), legend = SHORT[p$show], col = COL[p$show], pch = PCH[p$show], lwd = 4.5, bty = "n", cex = 1, text.col = INK, seg.len = 1.6, y.intersp = 1.3)
}
dev.off()

# ---- Figure 4: centered Type III with and without HC3, by the smallest group ---------------------------------
bins <- cut(null$min_group, c(0, 3, 5, 9, Inf), labels = c("3", "4 to 5", "6 to 9", "10 or more"))
h4 <- sapply(c("TypeIII_ctr", "TypeIII_ctr_HC3"), function(a) tapply(rej(null, a) > HARM, bins, mean))
png_open("Fig5", 6.5, 5.2); par(las = 1, ps = 12, fg = INK, col.lab = INK, col.axis = INK, cex.lab = 1, cex.axis = 1, mar = c(3.9, 4.8, 0.5, 0.5))
bp <- barplot(t(h4), beside = TRUE, plot = FALSE)
th_plot(range(bp) + c(-1, 1), c(0, 0.72), "Smallest group size", "Share of harmful designs", yline = 3.4, xline = 2.6, xaxs = "i", yat = seq(0, 0.6, 0.2), ylabs = paste0(seq(0, 60, 20), "%"))
barplot(t(h4), beside = TRUE, col = COL[c("TypeIII_ctr", "TypeIII_ctr_HC3")], border = NA, axes = FALSE, add = TRUE, names.arg = rep("", 4))
axis(1, at = colMeans(bp), labels = c("3", "4 to 5", "6 to 9", "10 or more"), tick = FALSE, line = -0.3)
legend("topright", legend = c("Type III, centered", "Type III, centered + HC3"), fill = COL[c("TypeIII_ctr", "TypeIII_ctr_HC3")], border = NA, bty = "n", cex = 0.95, text.col = INK)
dev.off()


# ---- Figure 1: workflow of how the simulation is conducted (colored, three panels, few words) ------------------------
FILL <- "#DCEBF7"; EDGE <- "#0072B2"; OR <- "#E69F00"; ORF <- "#FCEFD2"
png_open("Fig1", 6.5, 3.9); par(mar = c(0, 0, 0, 0), ps = 10, fg = INK)
plot(NA, xlim = c(0, 100), ylim = c(0, 69.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
bx <- function(x0, y0, x1, y1, txt, fill = FILL, edge = EDGE, cex = 1.05, font = 1) {
  rect(x0, y0, x1, y1, col = fill, border = edge, lwd = 2); text((x0 + x1) / 2, (y0 + y1) / 2, txt, cex = cex, col = INK, lheight = 1.0, font = font) }
ar <- function(x0, y0, x1, y1) arrows(x0, y0, x1, y1, length = 0.1, lwd = 2.2, col = EDGE)
PH <- 22; BH <- 11
pn <- function(y0, lab, ttl) { y1 <- y0 + PH; rect(0.5, y0, 99.5, y1, col = NA, border = EDGE, lwd = 1.6)
  text(2.5, y1 - 4, lab, font = 2, cex = 1.3, adj = c(0, 0.5), col = EDGE); invisible(ttl) }  # panel titles are described in the figure note
row4 <- function(y0, labs) { w <- 20.5; g <- 4.2; x <- 3 + (w + g) * 0:3; yb <- y0 + 3
  for (i in 1:4) bx(x[i], yb, x[i] + w, yb + BH, labs[[i]]); for (i in 1:3) ar(x[i] + w + 0.3, yb + BH / 2, x[i + 1] - 0.3, yb + BH / 2) }
pn(47.5, "a", ""); row4(47.5, list("Draw\ndesign", "Generate\ndata", "Run 8\nanalyses", expression(atop("Test " * italic(A) * ",", "record " * italic(p)))))
pn(24, "b", ""); row4(24, c("Exploratory", "Frozen\nprotocol", "Confirmatory", "Extension,\nscaling"))
pn(0.5, "c", ""); w <- 27; g <- 6.5; x <- 3.5 + (w + g) * 0:2; yb <- 3.5
bx(x[1], yb, x[1] + w, yb + BH, expression(atop("Condition", "on " * italic(x)))); bx(x[2], yb, x[2] + w, yb + BH, "Imhof error,\naverage"); bx(x[3], yb, x[3] + w, yb + BH, "Compare with\nsimulation")
for (i in 1:2) ar(x[i] + w + 0.3, yb + BH / 2, x[i + 1] - 0.3, yb + BH / 2)
dev.off()
cat("Figures written to", normalizePath(FIG_DIR), "\n")
