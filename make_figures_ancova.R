# FIGURES for the ANCOVA study. Base R only. Reads the CSV files written by explore_ancova.R and
# scaling_ancova.R and writes image files (not PDF) named Fig1 to Fig5 to <OUTPUT_ROOT>/figures. Independent of the
# simulation: rerun any time to change sizes or the harmful line. Paste into R or source().
# Journal layout: TIFF at 600 dpi (FIG_FORMAT = "png" writes PNG instead), sans-serif lettering (Arial on Windows),
# no titles inside the artwork, panels lettered a, b, c.
# Colors: eight categorical slots in a fixed order, one per analysis; every analysis also has its own plotting symbol,
# and the legends show both, so the figures do not depend on color alone.

MODE <- Sys.getenv("SIM_MODE", unset = "full")
OUTPUT_ROOT <- Sys.getenv("SIM_OUTPUT_ROOT", unset = "output")
FIG_DIR <- file.path(OUTPUT_ROOT, "figures"); dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
HARM <- 0.10; DPI <- as.integer(Sys.getenv("FIG_DPI", unset = "600"))
FMT <- tolower(Sys.getenv("FIG_FORMAT", unset = "tiff")); FONT <- Sys.getenv("FIG_FONT", unset = "sans")
INK <- "#0b0b0b"; AXT <- 1.15   # axis title size (25% larger than the earlier 0.92)
# Labels used in every figure; they match the row labels of the tables in the manuscript
ANALYSES <- c(Oneway = "One-way ANOVA", R_default = "R default (Type I, A first)", TypeII = "Type II",
              ANCOVA = "ANCOVA", TypeIII_raw = "Type III, raw",
              TypeIII_ctr = "Type III, centered", TypeIII_ctr_HC3 = "Type III, centered + HC3",
              Drop_int = "Drop interaction if not significant")
SHORT <- ANALYSES; SHORT["R_default"] <- "R default"  # shorter label where italics are not drawn
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

# Theme: no frame, axes drawn separately (applied to the line, scatter, and bar graphs)
th_plot <- function(xlim, ylim, xlab = "", ylab = "", xat = NULL, xlabs = xat, yat = seq(0, 1, 0.2), ylabs = yat, xaxs = "r", yline = 3.6, xline = 3.2, axt = AXT) {
  plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "", xaxs = xaxs)
  if (!is.null(xat)) axis(1, at = xat, labels = xlabs, col = INK, col.axis = INK, lwd = 1.6)
  axis(2, at = yat, labels = ylabs, col = INK, col.axis = INK, lwd = 1.6)
  mtext(ylab, side = 2, line = yline, las = 0, col = INK, cex = axt); mtext(xlab, side = 1, line = xline, col = INK, cex = axt)
}

# ---- Figure 2: share of harmful designs by situation (heatmap) ---------------------------------------------
groups <- list("All null\ndesigns" = rep(TRUE, nrow(null)),
               "No interaction,\ncovariate imbalance" = with(null, slope_range == 0 & cov_imbalance > 0),
               "Interaction,\ngrand mean far from 0" = with(null, slope_range > 0 & d0 >= 2),
               "Equal\nvariances" = with(null, sd_ratio == 1), "Unequal\nvariances" = with(null, sd_ratio > 1))
m <- sapply(groups, function(g) sapply(names(ANALYSES), function(a) mean(rej(null[g, ], a) > HARM)))
m <- m[order(-m[, 1]), , drop = FALSE]  # rows ordered from most to least often harmful (all null designs)
png_open("Fig2", 7.6, 6.5)
LB <- 1.25    # Fig. 2 label size: 25% larger than the original row, column, and key labels
# heatmap: one hue (orange), light = few harmful designs, dark = many; every cell carries its number
SEQ <- c("#FFF3DB", "#FFD08A", "#FFA040", "#F26B0F", "#C23E00", "#8A2800", "#4F1500")
ramp <- colorRampPalette(SEQ, space = "Lab")(101)
nr <- nrow(m); nc <- ncol(m)
# text color with the higher contrast against the cell (white or near-black)
rel_lum <- function(h) { v <- col2rgb(h)[, 1] / 255; v <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4); sum(v * c(0.2126, 0.7152, 0.0722)) }
txtcol <- function(h) { L <- rel_lum(h); if ((1.05) / (L + 0.05) > (L + 0.05) / (rel_lum(INK) + 0.05)) "white" else INK }
par(mar = c(4.2, 9.6, 9.0, 0.6), las = 1, fg = INK, col.axis = INK, ps = 11)
plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(nr + 0.5, 0.5), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
UN <- par("pin")[2] / nr                  # inches per row, so spacing below is set in inches
gap <- 0.04
for (i in seq_len(nr)) for (j in seq_len(nc)) {
  v <- m[i, j]
  rect(j - 0.5 + gap, i - 0.5 + gap, j + 0.5 - gap, i + 0.5 - gap, col = ramp[round(v * 100) + 1], border = NA)
  text(j, i, sprintf("%.0f%%", 100 * v), cex = 1.15, col = txtcol(ramp[round(v * 100) + 1]), font = 2)
}
# row labels, at most two lines
g <- 0.10 / UN                            # half the gap between two label lines (inches converted to rows)
for (i in seq_len(nr)) { a <- rownames(m)[i]
  if (a == "R_default") {
    text(0.4, i - g, "R default (Type I,", adj = c(1, 0.5), cex = 0.88 * LB, col = INK, xpd = NA)
    text(0.4, i + g, expression(italic(A) * " first)"), adj = c(1, 0.5), cex = 0.88 * LB, col = INK, xpd = NA)
  } else {
    ln <- strwrap(ANALYSES[a], 19)
    if (length(ln) == 1) text(0.4, i, ln, adj = c(1, 0.5), cex = 0.88 * LB, col = INK, xpd = NA)
    else { text(0.4, i - g, ln[1], adj = c(1, 0.5), cex = 0.88 * LB, col = INK, xpd = NA)
           text(0.4, i + g, ln[2], adj = c(1, 0.5), cex = 0.88 * LB, col = INK, xpd = NA) }
  } }
# column headers with n
HL <- list("All null\ndesigns", "No\ninteraction,\ncovariate\nimbalance", "Interaction,\ngrand mean\nfar from 0", "Equal\nvariances", "Unequal\nvariances")
LS <- 0.2 * LB / UN; NY <- 0.5 - 0.2 * LB / UN
for (j in seq_len(nc)) text(j, NY, bquote(italic(n) == .(sum(groups[[j]]))), adj = c(0.5, 0.5), cex = 0.85 * LB, col = INK, xpd = NA)
for (j in seq_len(nc)) { ln <- strsplit(HL[[j]], "\n")[[1]]; for (k in seq_along(ln)) text(j, NY - LS * (length(ln) - k + 1), ln[k], adj = c(0.5, 0.5), cex = 0.9 * LB, col = INK, xpd = NA) }
# color key: horizontal bar under the grid, 0% on the left and 100% on the right
ky <- nr + 0.5 + 0.42 / UN; kh <- 0.13 / UN; kx0 <- 1.6; kx1 <- 4.4; xs <- seq(kx0, kx1, length.out = 102)
for (k in 1:101) rect(xs[k], ky - kh / 2, xs[k + 1] + 0.002, ky + kh / 2, col = ramp[k], border = NA, xpd = NA)
text(kx0 - 0.08, ky, "0%", adj = c(1, 0.5), cex = 0.9 * LB, col = INK, xpd = NA)
text(kx1 + 0.08, ky, "100%", adj = c(0, 0.5), cex = 0.9 * LB, col = INK, xpd = NA)
text(0.4, ky, "Harmful designs", adj = c(1, 0.5), cex = 0.9 * LB, col = INK, xpd = NA)
dev.off()

# ---- Figure 3: Type I error against the distance of the grand mean of x from zero ---------------------------
s <- subset(null, slope_range > 0); show <- c("TypeIII_raw", "TypeIII_ctr", "TypeIII_ctr_HC3")
png_open("Fig3", 6.5, 5.8); par(las = 1, ps = 12, fg = INK, col.lab = INK, col.axis = INK, mar = c(6.4, 5.8, 0.4, 0.5), cex.lab = 1, cex.axis = 1)
AXT3 <- 1.3    # Fig. 3 axis titles
th_plot(c(0, 6), c(0, 1), "", "", xat = 0:6)
mtext(expression("Type I error for the effect of " * italic(A)), side = 2, line = 3.7, las = 0, col = INK, cex = AXT3)
# x title in two lines, centered on the whole figure
mtext(expression("Distance of the grand mean of " * italic(x) * " from 0"), side = 1, line = 3.0, col = INK, cex = AXT3, at = grconvertX(0.5, "ndc", "user"))
mtext("(standard deviation units)", side = 1, line = 4.7, col = INK, cex = AXT3, at = grconvertX(0.5, "ndc", "user"))
abline(h = c(0.05, HARM), lty = c(2, 3), col = INK, lwd = 1.4)
lo <- lapply(show, function(a) lowess(s$d0, rej(s, a), f = 0.5)); names(lo) <- show
for (a in show) {
  points(s$d0, rej(s, a), pch = PCH[a], cex = 0.8, col = adjustcolor(COL[a], 0.6)); lines(lo[[a]], col = COL[a], lwd = 5.5)
}
legend("topleft", inset = c(0.03, 0), legend = ANALYSES[show], col = COL[show], pch = PCH[show], pt.cex = 1.4,
       lwd = 5.5, bty = "n", cex = 1, text.col = INK, seg.len = 2.2, y.intersp = 1.35)
dev.off()

# ---- Figure 4: scaling experiment ---------------------------------------------------------------------------
sc <- read.csv(file.path(OUTPUT_ROOT, "ancova_scaling", "scaling_results.csv"), stringsAsFactors = FALSE)
# panels: letter plus a short subtitle naming the set (equal variances in sets 1 and 2 are stated in the figure note)
panels <- list(set1_raw_mismatch = c("TypeIII_raw", "TypeII", "TypeIII_ctr"),
               set2_cov_imbalance = c("R_default", "Oneway", "TypeII"),
               set3_variance = c("TypeIII_ctr", "TypeIII_ctr_HC3"))
ptitle <- c(set1_raw_mismatch = "a", set2_cov_imbalance = "b", set3_variance = "c")
psub <- c(set1_raw_mismatch = "Set 1: interaction, distant grand mean", set2_cov_imbalance = "Set 2: covariate imbalance, no interaction",
          set3_variance = "Set 3: unequal variances")
AXT4 <- 0.8    # Fig. 4 axis titles (layout() shrinks the tick labels, so these are set separately)
png_open("Fig4", 6.5, 10.6); par(las = 1, ps = 19, fg = INK, col.lab = INK, col.axis = INK, cex.lab = 1, cex.axis = 1); layout(matrix(1:3, ncol = 1), heights = c(1, 1, 1.14))
for (nm in names(panels)) {
  show <- panels[[nm]]; d <- subset(sc, set == nm); lastp <- nm == names(panels)[3]
  par(mar = c(if (lastp) 4.2 else 3.4, 5.2, 3.4, 0.5))
  YL <- if (nm == "set3_variance") 0.4 else 1
  th_plot(c(0, 4), c(0, YL), if (lastp) "Group sizes multiplied by" else "", "Mean Type I error", xat = 0:4, yline = 3.6, xline = 2.9, axt = AXT4,
          xlabs = if (lastp) c(1, 2, 4, 8, 16) else rep("", 5), yat = if (YL == 1) seq(0, 1, 0.2) else seq(0, 0.4, 0.1))
  # bold panel title, letter and description in one size: "a. Set 1: ..."
  mtext(paste0(ptitle[nm], ". ", psub[nm]), side = 3, line = 1.1, adj = 0, at = grconvertX(0, "nfc", "user") + 0.02 * diff(par("usr")[1:2]),
        font = 2, cex = 0.85, col = INK)
  segments(0, 0.05, 4, 0.05, lty = 2, col = INK, lwd = 1.4)
  for (a in show) { y <- tapply(rej(d, a), d$factor, mean); lines(0:4, y, col = COL[a], lwd = 4.5); points(0:4, y, pch = PCH[a], col = COL[a], cex = 1.7) }
  legend("topleft", inset = c(0.06, 0.02), legend = SHORT[show], col = COL[show], pch = PCH[show], pt.cex = 1.7, lwd = 4.5,
         bty = "n", cex = 1, text.col = INK, seg.len = 2.4, y.intersp = 1.3)
}
dev.off()

# ---- Figure 5: centered Type III with and without HC3, by the smallest group ----------------------------------
bins <- cut(null$min_group, c(0, 3, 5, 9, Inf), labels = c("3", "4 to 5", "6 to 9", "10 or more"))
h4 <- sapply(c("TypeIII_ctr", "TypeIII_ctr_HC3"), function(a) tapply(rej(null, a) > HARM, bins, mean))
png_open("Fig5", 6.5, 5.2); par(las = 1, ps = 12, fg = INK, col.lab = INK, col.axis = INK, cex.lab = 1, cex.axis = 1, mar = c(3.9, 4.8, 0.5, 0.5))
bp <- barplot(t(h4), beside = TRUE, plot = FALSE)
th_plot(range(bp) + c(-1, 1), c(0, 0.72), "Smallest group size", "Share of harmful designs", yline = 3.4, xline = 2.6, xaxs = "i", yat = seq(0, 0.6, 0.2), ylabs = paste0(seq(0, 60, 20), "%"))
barplot(t(h4), beside = TRUE, col = COL[c("TypeIII_ctr", "TypeIII_ctr_HC3")], border = NA, axes = FALSE, add = TRUE, names.arg = rep("", 4))
axis(1, at = colMeans(bp), labels = c("3", "4 to 5", "6 to 9", "10 or more"), tick = FALSE, line = -0.3)
legend("topright", legend = ANALYSES[c("TypeIII_ctr", "TypeIII_ctr_HC3")], fill = COL[c("TypeIII_ctr", "TypeIII_ctr_HC3")], border = NA, bty = "n", cex = 0.95, text.col = INK)
dev.off()

# ---- Figure 1: workflow of the study (three panels, few words) -------------------------------------------------
FILL <- "#DCEBF7"; EDGE <- "#0072B2"
png_open("Fig1", 6.5, 3.9); par(mar = c(0, 0, 0, 0), ps = 10, fg = INK)
plot(NA, xlim = c(0, 100), ylim = c(0, 69.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
# one or two lines of text, each line a list of (string, font) pieces, centered at (xc, yc) with a fixed line gap,
# so that labels with italic symbols have the same spacing as plain labels
mixed_line <- function(xc, y, parts, cex) {
  w <- sapply(parts, function(p) strwidth(p[[1]], cex = cex, font = p[[2]])); x <- xc - sum(w) / 2
  for (k in seq_along(parts)) { text(x, y, parts[[k]][[1]], adj = c(0, 0.5), cex = cex, font = parts[[k]][[2]], col = INK); x <- x + w[k] }
}
bx <- function(x0, y0, x1, y1, lines, cex = 1.05) {
  rect(x0, y0, x1, y1, col = FILL, border = EDGE, lwd = 2)
  xc <- (x0 + x1) / 2; yc <- (y0 + y1) / 2; gap <- 2.1 * strheight("Mg", cex = cex)
  ys <- if (length(lines) == 1) yc else yc + c(gap / 2, -gap / 2)
  for (k in seq_along(lines)) { L <- lines[[k]]; if (is.character(L)) L <- list(list(L, 1)); mixed_line(xc, ys[k], L, cex) }
}
ar <- function(x0, y0, x1, y1) arrows(x0, y0, x1, y1, length = 0.1, lwd = 2.2, col = EDGE)
PH <- 22; BH <- 11
pn <- function(y0, lab) { y1 <- y0 + PH; rect(0.5, y0, 99.5, y1, col = NA, border = EDGE, lwd = 1.6)
  text(2.5, y1 - 4, lab, font = 2, cex = 1.3, adj = c(0, 0.5), col = EDGE) }  # panel contents are described in the figure note
row4 <- function(y0, labs) { w <- 20.5; g <- 4.2; x <- 3 + (w + g) * 0:3; yb <- y0 + 3
  for (i in 1:4) bx(x[i], yb, x[i] + w, yb + BH, labs[[i]]); for (i in 1:3) ar(x[i] + w + 0.3, yb + BH / 2, x[i + 1] - 0.3, yb + BH / 2) }
pn(47.5, "a"); row4(47.5, list(c("Draw", "design"), c("Generate", "data"), c("Run 8", "analyses"),
                               list(list(list("Test ", 1), list("A", 3), list(",", 1)), list(list("record ", 1), list("p", 3)))))
pn(24, "b"); row4(24, list("Exploratory", c("Frozen", "protocol"), "Confirmatory", c("Extension,", "scaling")))
pn(0.5, "c"); w <- 27; g <- 6.5; x <- 3.5 + (w + g) * 0:2; yb <- 3.5
bx(x[1], yb, x[1] + w, yb + BH, list("Condition", list(list("on ", 1), list("x", 3))))
bx(x[2], yb, x[2] + w, yb + BH, c("Imhof error,", "average")); bx(x[3], yb, x[3] + w, yb + BH, c("Compare with", "simulation"))
for (i in 1:2) ar(x[i] + w + 0.3, yb + BH / 2, x[i + 1] - 0.3, yb + BH / 2)
dev.off()
cat("Figures written to", normalizePath(FIG_DIR), "\n")
