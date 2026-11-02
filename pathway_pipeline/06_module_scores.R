#!/usr/bin/env Rscript
# Module score per sample (control / FA / FA_edited), FA joined to its edited line.
# Module genes = leading-edge genes of chosen fgsea pathways (union over A, B, C of the cell).
# Score = mean over module genes of per-gene z-scored VST expression.  Needs: DESeq2, ggplot2, patchwork.
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

# ---- shared configuration (00_config.R lives next to this script) ----
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()      # in RStudio: setwd() to the script folder first
source(file.path(HERE, "00_config.R"))

## ---- parameters specific to this script ------------------------------
CELL <- CELL_DDS            # cell type (set in 00_config.R)
# module name = regex on the fgsea pathway name (all matches are pooled)
MODULES <- c("OXPHOS / respiratory chain" = "OXIDATIVE_PHOSPHORYLATION$|^GOCC_RESPIRATORY_CHAIN_COMPLEX$",
             "Synapse (synaptic membrane)" = "^GOCC_SYNAPTIC_MEMBRANE$",
             "Cytosolic ribosome"          = "^GOCC_CYTOSOLIC_RIBOSOME$",
             "Fe-S cluster assembly"       = "^GOCC_IRON_SULFUR_CLUSTER_ASSEMBLY_COMPLEX$")
## ---------------------------------------------------------------------

read_fgsea <- function(con) do.call(rbind, lapply(COLLS, function(cl) {
  f <- file.path(RES_DIR, CELL, con, sprintf("%s_%s_fgsea.csv", con, cl))
  d <- read.csv(f, stringsAsFactors = FALSE); d[, c("pathway", "leadingEdge")] }))
fg <- do.call(rbind, lapply(c(CON_A, CON_B, CON_C), read_fgsea))
mod_genes <- lapply(MODULES, function(rx) {
  hit <- unique(fg$pathway[grepl(rx, fg$pathway)])
  if (!length(hit)) stop("No pathway matches: ", rx)
  unique(unlist(strsplit(fg$leadingEdge[fg$pathway %in% hit], ";", fixed = TRUE))) })
cat("Module sizes (leading-edge genes):\n"); print(vapply(mod_genes, length, 1L))

# ---- expression: VST of all samples of the cell (control + FA + edited) ----
em <- load_expr_meta(CELL, fake_genes = unlist(mod_genes)); expr <- em$expr; meta <- em$meta
cat("Genes in expression matrix:", nrow(expr), "| example rownames:", paste(head(rownames(expr), 3), collapse = ", "), "\n")
for (nm in names(mod_genes)) cat(sprintf("  %-30s %d / %d module genes found\n", nm, sum(mod_genes[[nm]] %in% rownames(expr)), length(mod_genes[[nm]])))
if (any(vapply(mod_genes, function(g) sum(g %in% rownames(expr)) < 5, TRUE))) stop("A module has <5 genes in the expression matrix: gene identifiers do not match (see rownames above).")


z <- t(scale(t(expr[rowSums(is.finite(expr)) == ncol(expr), ])))
sc <- do.call(rbind, lapply(names(mod_genes), function(nm) {
  g <- intersect(mod_genes[[nm]], rownames(z))
  data.frame(module = sprintf("%s (%d genes)", nm, length(g)), sample = colnames(z), score = colMeans(z[g, , drop = FALSE])) }))
sc <- merge(sc, meta[, c("sample", "group", "line")], by = "sample")
sc$group  <- factor(as.character(sc$group), levels = c("control", "FA", "FA_edited"))
sc$module <- factor(sc$module, levels = unique(sc$module[order(match(sub(" \\(.*", "", sc$module), names(MODULES)))]))
write.csv(sc, file.path(OUT_DIR, sprintf("module_scores_%s.csv", CELL)), row.names = FALSE)

COLG <- c(control = "#555555", FA = "#8B1F14", FA_edited = "#0E7A5A")
p <- ggplot(sc, aes(group, score)) +
  geom_line(data = subset(sc, group != "control"), aes(group = line), colour = "grey55", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "grey30") +
  geom_point(aes(colour = group), size = 3.2) +
  geom_text(aes(label = line), size = 2.4, nudge_x = 0.13, hjust = 0, colour = "grey30") +
  scale_colour_manual(values = COLG, guide = "none") +
  scale_x_discrete(limits = c("control", "FA", "FA_edited"), expand = expansion(add = c(0.5, 0.9))) +
  facet_wrap(~module, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "module score (mean z-scored VST)", title = paste0(CELL, ": module scores per sample"),
       caption = "grey lines join each FA line to its edited line; bar = group mean; n = 3 per group, no test shown") +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9), plot.title = element_text(face = "bold"),
        plot.caption = element_text(size = 7, colour = "grey40"), axis.text.x = element_text(size = 9))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT_DIR, sprintf("module_scores_%s.png", CELL)), p, width = 15, height = 4.8, dpi = 170)
ggsave(file.path(OUT_DIR, sprintf("module_scores_%s.pdf", CELL)), p, width = 15, height = 4.8)
cat("Done ->", OUT_DIR, "\n")
