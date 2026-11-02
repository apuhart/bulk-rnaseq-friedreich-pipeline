#!/usr/bin/env Rscript
# Gene-level check: heatmap of selected genes per sample (control / FA / FA_edited), row z-scored VST.
# Answers e.g. "is the edited line more proliferative?", "is the proteasome / mito import signal carried by many genes?".
# Needs DESeq2, ggplot2.
suppressPackageStartupMessages({library(ggplot2)})

# ---- shared configuration (00_config.R lives next to this script) ----
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()      # in RStudio: setwd() to the script folder first
source(file.path(HERE, "00_config.R"))

## ---- parameters specific to this script ------------------------------
CELL <- CELL_DDS            # cell type (set in 00_config.R)
GENES <- list(
  "Proliferation / cell cycle" = c("MKI67", "TOP2A", "CCNB1", "CCNA2", "CDK1", "PCNA", "MCM2", "MCM7", "E2F1", "CENPF", "AURKB"),
  "Proteasome"                 = c("PSMA1", "PSMA2", "PSMA3", "PSMB1", "PSMB2", "PSMB5", "PSMC1", "PSMC3", "PSMD1", "PSMD2", "PSMD11"),
  "Mitochondrial import"       = c("TOMM20", "TOMM22", "TOMM40", "TIMM10", "TIMM13", "TIMM17A", "TIMM23", "TIMM44", "PAM16", "HSPD1"),
  "Fe-S / frataxin"            = c("FXN", "NFS1", "ISCU", "LYRM4", "HSCB", "GLRX5", "NFU1", "ISCA1", "FDX2"),
  "OXPHOS (examples)"          = c("NDUFA4", "NDUFB8", "NDUFS1", "SDHB", "UQCRC1", "UQCRFS1", "COX4I1", "COX7A2", "ATP5F1A", "ATP5F1B"),
  "Neuronal maturation"        = c("MAP2", "TUBB3", "RBFOX3", "SYP", "SNAP25", "SYT1", "DLG4", "GRIN1", "GRIA1", "STMN2"),
  "Progenitor / stem"          = c("SOX2", "NES", "PAX6", "VIM", "HES1", "POU5F1", "LIN28A"))
## ---------------------------------------------------------------------

em <- load_expr_meta(CELL, fake_genes = unlist(GENES)); expr <- em$expr; meta <- em$meta
meta <- meta[order(meta$group, meta$line), ]; expr <- expr[, meta$sample]

gl <- do.call(rbind, lapply(names(GENES), function(k) data.frame(set = k, gene = GENES[[k]], stringsAsFactors = FALSE)))
miss <- setdiff(gl$gene, rownames(expr)); if (length(miss)) message("Not found: ", paste(miss, collapse = ", "))
gl <- gl[gl$gene %in% rownames(expr), ]
z <- t(scale(t(expr[gl$gene, ])))
d <- data.frame(gene = rep(gl$gene, ncol(z)), set = rep(gl$set, ncol(z)), sample = rep(colnames(z), each = nrow(z)), z = as.vector(z))
d$set <- factor(d$set, levels = names(GENES)); d$gene <- factor(d$gene, levels = rev(gl$gene))
d <- merge(d, meta[, c("sample", "group", "line")], by = "sample")
d$col <- factor(paste(d$group, d$sample, sep = "|"), levels = paste(meta$group, meta$sample, sep = "|"))
lab <- setNames(paste0(meta$sample), paste(meta$group, meta$sample, sep = "|"))

p <- ggplot(d, aes(col, gene, fill = z)) + geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(set ~ group, scales = "free", space = "free", switch = "y") +
  scale_x_discrete(labels = function(x) sub(".*\\|", "", x)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-2, 2), oob = scales::squish, name = "z-score\n(per gene)") +
  labs(x = NULL, y = NULL, title = paste0(CELL, ": selected genes per sample"),
       caption = "VST expression, z-scored per gene across the 9 samples; FA and FA_edited samples are ordered by line (F4259, F4676, F86)") +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 7.5), strip.text.x = element_text(face = "bold"),
    strip.placement = "outside", plot.title = element_text(face = "bold"), plot.caption = element_text(size = 7, colour = "grey40"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT_DIR, sprintf("gene_heatmap_%s.png", CELL)), p, width = 7.5, height = 0.17 * nrow(gl) + 2.6, dpi = 170)
ggsave(file.path(OUT_DIR, sprintf("gene_heatmap_%s.pdf", CELL)), p, width = 7.5, height = 0.17 * nrow(gl) + 2.6)
write.csv(d[, c("set", "gene", "sample", "group", "line", "z")], file.path(OUT_DIR, sprintf("gene_heatmap_%s.csv", CELL)), row.names = FALSE)
cat("Done ->", OUT_DIR, "\n")
