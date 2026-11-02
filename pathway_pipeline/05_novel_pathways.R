#!/usr/bin/env Rscript
# New pathways that follow the A/B/C trend of the reference pathways you highlight.
# Needs Pathways_ABC_<cell>_all.csv from pathway_ABC.R (with leadingEdge_C).  Packages: ggplot2, ggrepel, igraph
# Outputs (OUT_DIR/novel/): New_pathways_following_reference_trends.xlsx (+ .csv) ; heatmap_<cell>_<group>.png ; profiles.png ; network_<cell>_<group>.png
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(igraph)})

# ---- shared configuration (00_config.R lives next to this script) ----
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()      # in RStudio: setwd() to the script folder first
source(file.path(HERE, "00_config.R"))

## ---- parameters specific to this script ------------------------------
COS_MIN     <- 0.9          # cosine similarity of the (NES_A, NES_B, NES_C) vectors
OVERLAP_MAX <- 0.3          # "novel" = < 30 % of leading-edge genes shared with the reference pathways
N_SHOW      <- 15           # novel pathways drawn per heatmap / network
## ---------------------------------------------------------------------

OUT <- file.path(OUT_DIR, "novel"); dir.create(OUT, showWarnings = FALSE)
le_of <- function(s) strsplit(s, ";", fixed = TRUE)
cosine <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))

res <- list()
for (cell in CELLS) {
  m <- read.csv(file.path(OUT_DIR, sprintf("Pathways_ABC_%s_all.csv", cell)), stringsAsFactors = FALSE)
  if (!"leadingEdge_C" %in% names(m)) stop("Re-run the updated pathway_ABC.R first (leadingEdge_C missing).")
  m$LE <- ifelse(abs(m$NES_C) > abs(m$NES_A) & !is.na(m$leadingEdge_C), m$leadingEdge_C, m$leadingEdge_A)
  for (g in names(REFS)) {
    r <- m[m$pathway %in% REFS[[g]], ]
    if (nrow(r) < 2) { message(cell, " | ", g, ": <2 reference pathways significant, skipped"); next }
    v <- colMeans(r[, c("NES_A", "NES_B", "NES_C")]); genes <- unique(unlist(le_of(paste(r$LE, collapse = ";"))))
    cand <- m[!(m$pathway %in% r$pathway) & !m$EMT_flag & m$size >= SIZE_MIN & m$size <= SIZE_MAX, ]
    P <- as.matrix(cand[, c("NES_A", "NES_B", "NES_C")]); pad <- as.matrix(cand[, c("padj_A", "padj_B", "padj_C")])
    ok <- rep(TRUE, nrow(cand))
    for (i in 1:3) ok <- ok & if (abs(v[i]) > 1) (sign(P[, i]) == sign(v[i]) & pad[, i] < 0.05) else (pad[, i] >= 0.05 | abs(P[, i]) < 1)
    cand <- cand[ok, ]; P <- P[ok, , drop = FALSE]
    if (!nrow(cand)) next
    cand$cos  <- apply(P, 1, cosine, b = v); cand$dist <- sqrt(rowSums((P - matrix(v, nrow(P), 3, byrow = TRUE))^2))
    cand$overlap_with_ref_genes <- vapply(le_of(cand$LE), function(l) mean(l %in% genes), 0)
    cand <- cand[cand$cos >= COS_MIN, ]
    if (!nrow(cand)) next
    cand$novel <- cand$overlap_with_ref_genes < OVERLAP_MAX
    cand$cell <- cell; cand$ref_group <- g; cand$ref_profile <- paste(round(v, 1), collapse = " / ")
    res[[length(res) + 1]] <- list(cell = cell, group = g, ref = r, cand = cand[order(cand$dist), ])
  }
}
if (!length(res)) stop("No candidate found.")
allc <- do.call(rbind, lapply(res, function(x) x$cand[, c("cell", "ref_group", "ref_profile", "collection", "pathway", "size", "NES_A", "padj_A", "NES_B", "padj_B",
                                                       "NES_C", "padj_C", "class", "cos", "dist", "overlap_with_ref_genes", "novel")]))
allc$pathway_short <- clean(allc$pathway)
allc <- allc[order(allc$cell, allc$ref_group, !allc$novel, allc$dist), ]
allc[sapply(allc, is.numeric)] <- lapply(allc[sapply(allc, is.numeric)], round, 3)
write.csv(allc, file.path(OUT, "new_pathways_candidates.csv"), row.names = FALSE)
summ <- aggregate(novel ~ cell + ref_group, allc, function(x) paste0(sum(x), " novel / ", length(x)))
names(summ)[3] <- "novel / candidates"; print(summ)
method <- data.frame(Note = c(
  paste0("Reference groups = the pathways you highlighted (", paste(names(REFS), collapse = "; "), "). A group is used only if >= 2 of its pathways are significant in the cell type."),
  "ref_profile = mean NES of the reference pathways in contrasts A, B, C (A = FA vs ctrl, B = edited vs ctrl, C = FA vs edited).",
  sprintf("Candidates: same sign and significance pattern as the reference group, cosine similarity >= %.2f of the (A,B,C) NES vector, %d-%d genes, EMT/matrix-driven pathways and reference pathways excluded.", COS_MIN, SIZE_MIN, SIZE_MAX),
  sprintf("dist = Euclidean distance between candidate and reference profile (smaller = more similar). overlap_with_ref_genes = fraction of the candidate's leading-edge genes also in the reference pathways' leading edges; novel = overlap < %.0f%%.", 100 * OVERLAP_MAX),
  "Caution: related GO terms share genes, so these are not independent evidence; n = 3 per group."))
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "summary");    openxlsx::writeData(wb, "summary", summ)
  openxlsx::addWorksheet(wb, "candidates"); openxlsx::writeData(wb, "candidates", allc)
  openxlsx::addWorksheet(wb, "Method");     openxlsx::writeData(wb, "Method", method)
  openxlsx::freezePane(wb, "candidates", firstRow = TRUE); openxlsx::setColWidths(wb, "candidates", cols = 1:ncol(allc), widths = "auto")
  openxlsx::setColWidths(wb, "Method", cols = 1, widths = 140)
  openxlsx::saveWorkbook(wb, file.path(OUT, "New_pathways_following_reference_trends.xlsx"), overwrite = TRUE)
} else message("openxlsx not installed -> CSV only")

long <- function(d, type) {
  do.call(rbind, lapply(c("A", "B", "C"), function(k)
    data.frame(path = clean(d$pathway), type = type, contrast = k, NES = d[[paste0("NES_", k)]], padj = d[[paste0("padj_", k)]],
               overlap = if (type == "New") d$overlap_with_ref_genes else NA, stringsAsFactors = FALSE))) }
CLAB <- c(A = "A\nFA vs ctrl", B = "B\nedited vs ctrl", C = "C\nFA vs edited")

# ---------- 1) heatmaps: references on top, new pathways below ----------
for (x in res) {
  nv <- head(x$cand[x$cand$novel, ], N_SHOW); if (!nrow(nv)) next
  d <- rbind(long(x$ref, "Reference"), long(nv, "New"))
  d$type <- factor(d$type, levels = c("Reference", "New"))
  ord <- unique(d$path[order(d$type, d$contrast == "A", decreasing = FALSE)]); d$path <- factor(d$path, levels = rev(unique(d$path)))
  d$contrast <- factor(CLAB[d$contrast], levels = CLAB); d$star <- ifelse(d$padj < 0.05, "*", "")
  lim <- ceiling(max(abs(d$NES))); ov <- unique(d[d$type == "New", c("path", "type", "overlap")])
  p <- ggplot(d, aes(contrast, path, fill = NES)) + geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.1f%s", NES, star)), size = 2.9, colour = ifelse(abs(d$NES) > 0.6 * lim, "white", "grey15")) +
    geom_text(data = ov, aes(x = 3.6, y = path, label = sprintf("%.0f%%", 100 * overlap)), inherit.aes = FALSE, size = 2.7, colour = "grey35", hjust = 0) +
    geom_text(data = data.frame(type = factor("New", levels = c("Reference", "New"))), aes(x = 3.6, y = Inf, label = "gene overlap"), inherit.aes = FALSE, size = 2.4, colour = "grey35", hjust = 0, vjust = -0.6) +
    facet_grid(type ~ ., scales = "free_y", space = "free_y") + coord_cartesian(clip = "off", xlim = c(0.5, 3.5)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-lim, lim), name = "NES") +
    labs(x = NULL, y = NULL, title = sprintf("%s - trend of '%s'", x$cell, x$group),
         caption = "* padj < 0.05 ; gene overlap = share of the pathway's leading-edge genes also in the reference pathways; EMT/matrix-driven pathways removed") +
    theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), plot.margin = margin(5, 45, 5, 5),
        plot.title = element_text(face = "bold"), plot.caption = element_text(size = 7, colour = "grey40"), strip.text.y = element_text(face = "bold", angle = 0))
  ggsave(file.path(OUT, sprintf("heatmap_%s_%s.png", x$cell, gsub("[^A-Za-z]+", "_", x$group))), p, width = 8.5, height = 2 + 0.26 * (nrow(x$ref) + nrow(nv)), dpi = 170)
}

# ---------- 2) profiles A -> B -> C ----------
pr <- do.call(rbind, lapply(res, function(x) {
  nv <- head(x$cand[x$cand$novel, ], 25)
  rbind(cbind(long(x$ref, "Reference"), cell = x$cell, grp = x$group), if (nrow(nv)) cbind(long(nv, "New"), cell = x$cell, grp = x$group)) }))
pr$id <- paste(pr$cell, pr$grp, pr$path); pr$contrast <- factor(pr$contrast, levels = c("A", "B", "C"), labels = c("A\nFA vs ctrl", "B\nedited vs ctrl", "C\nFA vs edited"))
lab <- pr[pr$type == "Reference" & pr$contrast == levels(pr$contrast)[3], ]
p2 <- ggplot(pr, aes(contrast, NES, group = id)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(data = subset(pr, type == "New"), colour = "#4C8BC8", alpha = 0.55, linewidth = 0.5) +
  geom_line(data = subset(pr, type == "Reference"), colour = "#B2182B", linewidth = 1.1) +
  geom_text_repel(data = lab, aes(label = substr(path, 1, 30)), colour = "#B2182B", size = 2.5, hjust = 0, nudge_x = 0.1, direction = "y", segment.size = 0.2, max.overlaps = Inf, xlim = c(3.1, NA)) +
  scale_x_discrete(expand = expansion(add = c(0.2, 1.4))) + facet_grid(cell ~ grp, scales = "free_y") +
  labs(x = NULL, y = "NES", title = "Reference pathways (red) and new pathways following the same trend (blue)") +
  theme_classic(base_size = 10) + theme(strip.text = element_text(face = "bold", size = 8), plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "profiles.png"), p2, width = 4 * length(unique(pr$grp)) + 2, height = 7, dpi = 170, limitsize = FALSE)

# ---------- 3) network by shared leading-edge genes ----------
for (x in res) {
  nv <- head(x$cand[x$cand$novel, ], 30); if (!nrow(nv)) next
  nd <- rbind(data.frame(path = x$ref$pathway, type = "Reference", LE = x$ref$LE, NES_C = x$ref$NES_C),
              data.frame(path = nv$pathway, type = "New", LE = nv$LE, NES_C = nv$NES_C))
  L <- le_of(nd$LE); n <- nrow(nd); J <- matrix(0, n, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) J[i, j] <- J[j, i] <- length(intersect(L[[i]], L[[j]])) / length(union(L[[i]], L[[j]]))
  ed <- which(upper.tri(J) & J >= 0.10, arr.ind = TRUE)
  g <- graph_from_data_frame(data.frame(from = ed[, 1], to = ed[, 2], w = J[ed]), vertices = data.frame(name = 1:n), directed = FALSE)
  set.seed(1); xy <- layout_with_fr(g, weights = E(g)$w * 5); nd$x <- xy[, 1]; nd$y <- xy[, 2]
  seg <- data.frame(x = nd$x[ed[, 1]], y = nd$y[ed[, 1]], xend = nd$x[ed[, 2]], yend = nd$y[ed[, 2]], w = J[ed])
  nd$lab <- substr(clean(nd$path), 1, 32)
  p3 <- ggplot() + geom_segment(data = seg, aes(x, y, xend = xend, yend = yend, linewidth = w), colour = "grey70", alpha = 0.7) +
    geom_point(data = nd, aes(x, y, fill = type, size = type), shape = 21, colour = "white", stroke = 0.6) +
    geom_text_repel(data = nd, aes(x, y, label = lab, colour = type), size = 2.6, max.overlaps = Inf, box.padding = 0.35, segment.size = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = c(Reference = "#B2182B", New = "#4C8BC8"), name = NULL) + scale_colour_manual(values = c(Reference = "#B2182B", New = "#2B5C8A")) +
    scale_size_manual(values = c(Reference = 6, New = 3.5), guide = "none") + scale_linewidth(range = c(0.2, 2), name = "shared\ngenes (Jaccard)") +
    labs(title = sprintf("%s - pathways linked by shared leading-edge genes ('%s')", x$cell, x$group)) +
    theme_void(base_size = 10) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path(OUT, sprintf("network_%s_%s.png", x$cell, gsub("[^A-Za-z]+", "_", x$group))), p3, width = 9, height = 7.5, dpi = 170)
}
cat("Done ->", OUT, "\n")
