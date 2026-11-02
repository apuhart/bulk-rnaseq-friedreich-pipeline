#!/usr/bin/env Rscript
# Rescue map: NES(A: FA vs ctrl) vs NES(C: FA vs edited), one panel per cell type.
# Needs the CSVs written by pathway_ABC.R. Packages: ggplot2, ggrepel, patchwork.
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(patchwork)})

# ---- shared configuration (00_config.R lives next to this script) ----
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()      # in RStudio: setwd() to the script folder first
source(file.path(HERE, "00_config.R"))

## ---- parameters specific to this script ------------------------------
N_LABEL <- c(iPSC = 6, iNeurons = 8)      # number of cluster representatives to label per panel
OUT     <- file.path(OUT_DIR, "rescue_map_NES_A_vs_C")
## ---------------------------------------------------------------------

COL <- c("Rescued" = "#0E7A5A", "Partly rescued" = "#E08A00", "Over-corrected" = "#1F4E5F",
         "Persistent" = "#8B1F14", "Intermediate" = "#7A8FA6", "Edited-vs-ctrl shift" = "#8E5BB5",
         "Other" = "#CCCCCC", "Discordant (A vs C)" = "#CCCCCC")
ORDER <- c("Other", "Discordant (A vs C)", "Intermediate", "Edited-vs-ctrl shift",
           "Persistent", "Partly rescued", "Over-corrected", "Rescued")   # drawing order (rescued on top)

short <- function(x) tolower(gsub("_", " ", sub("^(GOBP|GOCC|GOMF)_", "", sub("^HALLMARK_", "H: ", x))))

panel <- function(cell) {
  m  <- read.csv(file.path(OUT_DIR, sprintf("Pathways_ABC_%s_all.csv", cell)), stringsAsFactors = FALSE)
  cl <- read.csv(file.path(OUT_DIR, sprintf("Pathways_ABC_%s_clusters.csv", cell)), stringsAsFactors = FALSE)
  m <- m[m$padj_A < 0.05 | m$padj_C < 0.05, ]
  m$class <- factor(m$class, levels = ORDER)
  m <- m[order(m$class), ]                                   # rescued drawn last
  n <- table(m$class); lev <- paste0(names(n), " (", n, ")")
  m$class_lab <- factor(paste0(m$class, " (", n[as.character(m$class)], ")"), levels = lev)
  cols <- setNames(COL[names(n)], lev)
  lab <- head(cl, N_LABEL[[cell]]); lab$lab <- short(lab$cluster_rep)
  lim <- ceiling(max(abs(c(m$NES_A, m$NES_C)))) + 0.2
  ggplot(m, aes(NES_A, NES_C)) +
    geom_abline(linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_point(data = subset(m, !EMT_flag), aes(colour = class_lab), size = 1.8, alpha = 0.85) +
    geom_point(data = subset(m, EMT_flag),  aes(colour = class_lab), size = 1.8, shape = 1, alpha = 0.8) +
    geom_text_repel(data = lab, aes(NES_A, NES_C, label = lab), size = 2.8, min.segment.length = 0,
                    segment.size = 0.25, box.padding = 0.6, max.overlaps = Inf, force = 6, seed = 1) +
    scale_colour_manual(values = cols, name = NULL) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(title = cell, x = "NES  A: FA vs ctrl", y = "NES  C: FA vs edited",
         caption = "diagonal = effect fully reversed by editing; hollow = EMT/matrix-driven") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 14), legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1), legend.key.size = unit(0.35, "cm"),
          legend.text = element_text(size = 7), legend.background = element_blank(),
          plot.caption = element_text(size = 7, colour = "grey40"))
}

p <- wrap_plots(lapply(CELLS, panel), nrow = 1)
ggsave(paste0(OUT, ".png"), p, width = 16, height = 7.5, dpi = 170)
ggsave(paste0(OUT, ".pdf"), p, width = 16, height = 7.5)
cat("Saved", OUT, "\n")
