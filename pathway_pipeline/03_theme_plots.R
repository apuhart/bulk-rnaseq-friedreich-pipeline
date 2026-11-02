#!/usr/bin/env Rscript
# Two plots coloured by pathway THEME (needs the CSVs from pathway_ABC.R):
#  1) "Rescued" pathways : NES(A: FA vs ctrl) vs NES(C: FA vs edited)   [diagonal y = x]
#  2) "Edited-vs-ctrl shift": NES(B: edited vs ctrl) vs NES(C: FA vs edited) [anti-diagonal y = -x, since A ~ 0]
# Packages: ggplot2, ggrepel, patchwork
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(patchwork)})

# ---- shared configuration (00_config.R lives next to this script) ----
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()      # in RStudio: setwd() to the script folder first
source(file.path(HERE, "00_config.R"))

## ---- parameters specific to this script ------------------------------
LAB_PER_THEME <- 2      # labelled pathways per theme (largest |NES_C|)
## ---------------------------------------------------------------------

panel <- function(cell, klass, xcol, ycol, xlab, ylab, diag_slope, note) {
  m <- read.csv(file.path(OUT_DIR, sprintf("Pathways_ABC_%s_all.csv", cell)), stringsAsFactors = FALSE)
  m <- m[m$class == klass, ]
  if (DROP_EMT) m <- m[!m$EMT_flag, ]
  if (nrow(m) == 0) return(NULL)
  m$theme <- factor(theme_of(m$pathway), levels = names(THEME_COL))
  m$x <- m[[xcol]]; m$y <- m[[ycol]]; m$lab <- clean(m$pathway)
  m <- m[order(m$theme == "Other", decreasing = TRUE), ]
  lab <- do.call(rbind, lapply(split(m, m$theme), function(g) head(g[order(-abs(g$NES_C)), ], LAB_PER_THEME)))
  # zoom on the bulk of the points (2-98 % quantiles) so that a lone outlier does not squash the plot
  qx <- quantile(m$x, c(0.02, 0.98)); qy <- quantile(m$y, c(0.02, 0.98))
  padx <- 0.25 * diff(qx) + 0.3; pady <- 0.25 * diff(qy) + 0.3
  xl <- range(c(qx, 0)) + c(-padx, padx); yl <- range(c(qy, 0)) + c(-pady, pady)
  nout <- sum(m$x < xl[1] | m$x > xl[2] | m$y < yl[1] | m$y > yl[2])
  if (nout > 0) note <- paste0(note, sprintf("  [%d point(s) outside the zoom]", nout))
  ggplot(m, aes(x, y)) +
    geom_abline(slope = diag_slope, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_point(aes(colour = theme, shape = EMT_flag), size = 2.4, alpha = 0.85) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    geom_text_repel(data = lab, aes(label = substr(lab, 1, 42), colour = theme), size = 2.9, fontface = "bold",
                    min.segment.length = 0, segment.size = 0.25, box.padding = 0.5, max.overlaps = Inf,
                    force = 6, seed = 1, show.legend = FALSE) +
    scale_colour_manual(values = THEME_COL, drop = FALSE, name = "Theme") +
    coord_cartesian(xlim = xl, ylim = yl) +
    labs(title = sprintf("%s  -  %s (n = %d)", cell, klass, nrow(m)), x = xlab, y = ylab, caption = note) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13), plot.caption = element_text(size = 7, colour = "grey40"),
          legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 8))
}

p1 <- wrap_plots(Filter(Negate(is.null), lapply(CELLS, panel, klass = "Rescued", xcol = "NES_A", ycol = "NES_C",
                        xlab = "NES  A: FA vs ctrl", ylab = "NES  C: FA vs edited", diag_slope = 1,
                        note = "diagonal = effect fully reversed by editing")), nrow = 1, guides = "collect")
p2 <- wrap_plots(Filter(Negate(is.null), lapply(CELLS, panel, klass = "Edited-vs-ctrl shift", xcol = "NES_B", ycol = "NES_C",
                        xlab = "NES  B: edited vs ctrl", ylab = "NES  C: FA vs edited", diag_slope = -1,
                        note = "A not significant; anti-diagonal = edited differs from both FA and ctrl")), nrow = 1, guides = "collect")
ggsave(file.path(OUT_DIR, "theme_plot_rescued.png"), p1, width = 16, height = 7.5, dpi = 170)
ggsave(file.path(OUT_DIR, "theme_plot_rescued.pdf"), p1, width = 16, height = 7.5)
ggsave(file.path(OUT_DIR, "theme_plot_edited_shift.png"), p2, width = 10, height = 7.5, dpi = 170)
ggsave(file.path(OUT_DIR, "theme_plot_edited_shift.pdf"), p2, width = 10, height = 7.5)
# theme table, useful to check the keyword mapping
for (cell in CELLS) {
  m <- read.csv(file.path(OUT_DIR, sprintf("Pathways_ABC_%s_all.csv", cell)), stringsAsFactors = FALSE)
  m <- m[m$class %in% c("Rescued", "Edited-vs-ctrl shift"), ]; m$theme <- theme_of(m$pathway)
  write.csv(m[, c("class", "theme", "collection", "pathway", "NES_A", "NES_B", "NES_C", "EMT_flag")],
            file.path(OUT_DIR, sprintf("theme_assignment_%s.csv", cell)), row.names = FALSE)
  print(table(m$class, m$theme))
}
cat("Done ->", OUT_DIR, "\n")
