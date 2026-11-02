#!/usr/bin/env Rscript
# Overview heatmap: pathway THEME x contrast (A, B, C) per cell type. Needs Pathways_ABC_<cell>_all.csv from pathway_ABC.R
# Tile colour = mean NES of the theme's pathways (those significant in >=1 contrast); number = pathways significant in that contrast / pathways in theme.
suppressPackageStartupMessages({library(ggplot2)})

# ---- shared configuration (00_config.R lives next to this script) ----
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()      # in RStudio: setwd() to the script folder first
source(file.path(HERE, "00_config.R"))

## ---- parameters specific to this script ------------------------------
MIN_N    <- 5            # hide themes with fewer pathways than this (per cell)
## ---------------------------------------------------------------------

rows <- list()
for (cell in CELLS) {
  m <- read.csv(file.path(OUT_DIR, sprintf("Pathways_ABC_%s_all.csv", cell)), stringsAsFactors = FALSE)
  if (DROP_EMT) m <- m[!m$EMT_flag, ]
  m$theme <- theme_of(m$pathway)
  for (th in unique(m$theme)) {
    g <- m[m$theme == th, ]; if (nrow(g) < MIN_N) next
    for (ct in c("A", "B", "C")) {
      nes <- g[[paste0("NES_", ct)]]; pad <- g[[paste0("padj_", ct)]]
      rows[[length(rows) + 1]] <- data.frame(cell = cell, theme = th, contrast = ct, n_theme = nrow(g),
        meanNES = mean(nes, na.rm = TRUE), n_sig = sum(pad < PADJ, na.rm = TRUE))
    }
  }
}
d <- do.call(rbind, rows)
d$contrast <- factor(d$contrast, levels = c("A", "B", "C"),
                     labels = c("A\nFA vs ctrl", "B\nedited vs ctrl", "C\nFA vs edited"))
d$cell <- factor(d$cell, levels = CELLS)
ord <- aggregate(abs(meanNES) ~ theme, d, max); d$theme <- factor(d$theme, levels = rev(ord$theme[order(-ord[[2]])]))
d$label <- sprintf("%.1f\n(%d/%d)", d$meanNES, d$n_sig, d$n_theme)
lim <- ceiling(max(abs(d$meanNES)) * 2) / 2

p <- ggplot(d, aes(contrast, theme, fill = meanNES)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.7, lineheight = 0.85,
            colour = ifelse(abs(d$meanNES) > 0.6 * lim, "white", "grey15")) +
  facet_wrap(~cell, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim),
                       name = "mean NES\n(red = higher in\nfirst-named group)") +
  labs(x = NULL, y = NULL,
       caption = "number = significant pathways (padj<0.05) in this contrast / pathways in the theme (those significant in >=1 contrast); EMT/matrix-driven pathways removed") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 9), plot.caption = element_text(size = 7, colour = "grey40"))
ggsave(file.path(OUT_DIR, "theme_heatmap.png"), p, width = 12, height = 5.6, dpi = 170)
ggsave(file.path(OUT_DIR, "theme_heatmap.pdf"), p, width = 12, height = 5.6)
write.csv(d[, c("cell", "theme", "contrast", "n_theme", "n_sig", "meanNES")], file.path(OUT_DIR, "theme_heatmap_data.csv"), row.names = FALSE)
cat("Done ->", OUT_DIR, "\n")
