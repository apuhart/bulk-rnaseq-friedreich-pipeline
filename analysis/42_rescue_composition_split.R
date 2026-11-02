#!/usr/bin/env Rscript
# =============================================================================
# 42_rescue_composition_split.R  (iNeurons, bulk RNA v4)
# -----------------------------------------------------------------------------
# Question : parmi les gènes « rescued » des iNeurons, lesquels suivent la
# composition cellulaire (score mésenchyme) et lesquels sont un vrai rescue
# indépendant de la composition (ex. FXN) ?
#
# 1. Reconstruit la liste des gènes « Rescued » avec EXACTEMENT les règles du
#    pipeline bulk (run_rescue_scatter) :
#      FA-DEG     : padj(FA vs ctrl) < 0.05 et |log2FC| >= 1
#      Rescued    : FA-DEG, NON significatif en éditée vs ctrl, et |LFC éditée| < |LFC FA|
#      (+ Partial rescue : encore significatif mais |LFC éditée| < |LFC FA|)
# 2. Score mésenchyme par échantillon = moyenne log2(comptages normalisés + 1)
#    de COL1A1, COL1A2, FN1, THY1, ACTA2, TAGLN (même score que la slide).
# 3. Pour chaque gène rescued :
#      r_score  : corrélation (9 échantillons) entre son expression et le score,
#                 dans le sens de son effet FA (r_signed > 0 = suit la composition)
#      pairs_ok : dans combien des 3 paires FA -> éditée (même donneur) le gène
#                 bouge dans le sens du rescue (|delta log2| >= DELTA_MIN)
# 4. Classes :
#      composition    r_signed >= R_COMP
#      F4676 only     seule la paire F4676 va dans le bon sens
#      consistent     3/3 paires dans le bon sens et r_signed < R_COMP
#      partial        2/3 paires, r_signed < R_COMP
#      unclear        le reste
#
# Sorties : OUT_DIR (défaut results/rescue_composition_split)
#   42_rescued_genes_classified.csv
#   42_summary.csv / .pdf (barplot des classes + nuage r_signed vs paires)
#
# Options : RNA_V4_BASE (dossier iNeurons du bulk v4), INCLUDE_PARTIAL=1, SCORE_GENES,
#           R_COMP (0.7), DELTA_MIN (0.3), LFC_THR (1), PADJ (0.05)
# Usage   : Rscript scripts/42_rescue_composition_split.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
  library(DESeq2); library(AnnotationDbi); library(org.Hs.eg.db)
})

BASE <- Sys.getenv("RNA_V4_BASE",
  unset = file.path(getwd(), "results_v4", "iNeurons"))
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
OUT_DIR <- Sys.getenv("OUT_DIR", unset = file.path(PROJECT_ROOT, "results", "rescue_composition_split"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
PADJ <- as.numeric(Sys.getenv("PADJ", "0.05")); LFC_THR <- as.numeric(Sys.getenv("LFC_THR", "1"))
R_COMP <- as.numeric(Sys.getenv("R_COMP", "0.7")); DELTA_MIN <- as.numeric(Sys.getenv("DELTA_MIN", "0.3"))
INCLUDE_PARTIAL <- Sys.getenv("INCLUDE_PARTIAL", "0") == "1"
SCORE_GENES <- { e <- Sys.getenv("SCORE_GENES", ""); if (nzchar(e)) trimws(strsplit(e, ",")[[1]]) else c("COL1A1", "COL1A2", "FN1", "THY1", "ACTA2", "TAGLN") }
say <- function(...) cat(sprintf(...), "\n", sep = "")

# ---------------------------------------------------------------------------
# 1. DEG tables -> rescued genes (règles du pipeline)
# ---------------------------------------------------------------------------
read_deg <- function(ct) {
  f <- file.path(BASE, ct, paste0(ct, "_DEG_results.csv"))
  if (!file.exists(f)) stop("Introuvable : ", f, "  -> ajuste RNA_V4_BASE")
  d <- read_csv(f, show_col_types = FALSE)
  idc <- intersect(c("GENEID", "gene_id", "ensembl_id", "ENSEMBL", "...1", "X"), names(d))[1]
  if (is.na(idc)) stop("Colonne d'identifiant introuvable dans ", f, " : ", paste(names(d), collapse = ", "))
  d %>% transmute(GENEID = sub("\\..*$", "", .data[[idc]]), lfc = log2FoldChange, padj = padj)
}
fa <- read_deg("control_vs_FA"); ed <- read_deg("control_vs_FA_edited")
df <- inner_join(fa, ed, by = "GENEID", suffix = c("_FA", "_ED")) %>%
  filter(!is.na(lfc_FA), !is.na(lfc_ED)) %>%
  mutate(sig_FA = !is.na(padj_FA) & padj_FA < PADJ & abs(lfc_FA) >= LFC_THR,
         sig_ED = !is.na(padj_ED) & padj_ED < PADJ & abs(lfc_ED) >= LFC_THR,
         status = case_when(!sig_FA ~ "Not FA-DEG",
                            !sig_ED & abs(lfc_ED) < abs(lfc_FA) ~ "Rescued",
                            sig_ED & abs(lfc_ED) < abs(lfc_FA) ~ "Partial rescue",
                            sign(lfc_FA) != sign(lfc_ED) ~ "Reversed",
                            sig_ED ~ "Not rescued", TRUE ~ "Other"))
say("FA-DEGs : %d | Rescued : %d | Partial : %d", sum(df$sig_FA), sum(df$status == "Rescued"), sum(df$status == "Partial rescue"))
keep_status <- if (INCLUDE_PARTIAL) c("Rescued", "Partial rescue") else "Rescued"
resc <- df %>% filter(status %in% keep_status)

# ---------------------------------------------------------------------------
# 2. Normalized expression of the 9 samples + mesenchymal score
# ---------------------------------------------------------------------------
dds_files <- file.path(BASE, c("control_vs_FA", "FA_edited_vs_FA", "control_vs_FA_edited"),
                       paste0(c("control_vs_FA", "FA_edited_vs_FA", "control_vs_FA_edited"), "_dds.rds"))
dl <- lapply(dds_files, readRDS)
common <- Reduce(intersect, lapply(dl, rownames))
cnt <- do.call(cbind, lapply(dl, function(d) counts(d)[common, , drop = FALSE]))
cnt <- cnt[, !duplicated(colnames(cnt)), drop = FALSE]
cd <- do.call(rbind, lapply(dl, function(d) as.data.frame(colData(d))[, "condition", drop = FALSE]))
cond <- setNames(as.character(cd$condition), rownames(cd))[colnames(cnt)]
sf <- estimateSizeFactorsForMatrix(cnt)
lex <- log2(t(t(cnt) / sf) + 1)
rownames(lex) <- sub("\\..*$", "", rownames(lex))
lex <- lex[!duplicated(rownames(lex)), , drop = FALSE]
sym <- suppressMessages(mapIds(org.Hs.eg.db, rownames(lex), "SYMBOL", "ENSEMBL", multiVals = "first"))
sid <- names(sym)[!is.na(sym) & sym %in% SCORE_GENES]
if (length(sid) < 3) stop("Moins de 3 gènes du score trouvés.")
score <- colMeans(lex[sid, , drop = FALSE])
say("Score mésenchyme (%s) :", paste(sym[sid], collapse = ", "))
print(round(score, 2))

# paires FA -> éditée (même donneur)
donor <- sub("(_R|ed)$", "", colnames(lex))
pairs <- lapply(unique(donor[cond %in% c("FA", "FA_edited")]), function(d) {
  fa_s <- colnames(lex)[donor == d & cond == "FA"]; ed_s <- colnames(lex)[donor == d & cond == "FA_edited"]
  if (length(fa_s) == 1 && length(ed_s) == 1) c(donor = d, FA = fa_s, ED = ed_s) else NULL
})
pairs <- do.call(rbind, pairs)
say("Paires : %s", paste(sprintf("%s->%s", pairs[, "FA"], pairs[, "ED"]), collapse = ", "))

# ---------------------------------------------------------------------------
# 3. Per-gene metrics and classes
# ---------------------------------------------------------------------------
g <- intersect(resc$GENEID, rownames(lex))
res <- resc %>% filter(GENEID %in% g) %>% mutate(SYMBOL = unname(sym[GENEID]))
res$r_score <- sapply(res$GENEID, function(id) suppressWarnings(cor(lex[id, ], score)))
res$r_signed <- res$r_score * sign(res$lfc_FA)
for (k in seq_len(nrow(pairs))) {
  d <- pairs[k, ]
  delta <- lex[res$GENEID, d["ED"]] - lex[res$GENEID, d["FA"]]
  res[[paste0("delta_", d["donor"])]] <- unname(delta)
  res[[paste0("ok_", d["donor"])]] <- unname(sign(delta) == -sign(res$lfc_FA) & abs(delta) >= DELTA_MIN)
}
okc <- grep("^ok_", names(res), value = TRUE)
res$pairs_ok <- rowSums(res[, okc, drop = FALSE])
only4676 <- if ("ok_F4676" %in% okc) res$ok_F4676 & res$pairs_ok == 1 else FALSE
res <- res %>% mutate(class = case_when(
  r_signed >= R_COMP ~ "composition",
  only4676 ~ "F4676 only",
  pairs_ok == length(okc) ~ "consistent rescue",
  pairs_ok == length(okc) - 1 ~ "partial (2/3 pairs)",
  TRUE ~ "unclear"),
  class = factor(class, levels = c("consistent rescue", "partial (2/3 pairs)", "F4676 only", "composition", "unclear"))) %>%
  arrange(class, desc(abs(lfc_FA)))
write_csv(res, file.path(OUT_DIR, "42_rescued_genes_classified.csv"))

summ <- res %>% dplyr::count(class, .drop = FALSE) %>% mutate(pct = round(100 * n / sum(n), 1))
write_csv(summ, file.path(OUT_DIR, "42_summary.csv"))
say("\n== %d gènes %s ==", nrow(res), paste(keep_status, collapse = " + "))
print(as.data.frame(summ), row.names = FALSE)
say("\nRescue cohérent (indépendant de la composition) :")
print(as.data.frame(res %>% filter(class == "consistent rescue") %>%
        dplyr::select(SYMBOL, lfc_FA, lfc_ED, r_signed, pairs_ok) %>% mutate(across(where(is.numeric), ~ round(.x, 2)))), row.names = FALSE)
if (any(res$SYMBOL == "FXN", na.rm = TRUE)) say("\nFXN : classe = %s", as.character(res$class[which(res$SYMBOL == "FXN")]))

cls_col <- c("consistent rescue" = "#1B7837", "partial (2/3 pairs)" = "#7FBF7B", "F4676 only" = "#D08C2E",
             "composition" = "#B2182B", "unclear" = "grey70")
p1 <- ggplot(summ, aes(class, n, fill = class)) + geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", n, pct)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = cls_col, guide = "none") +
  labs(title = sprintf("iNeurons: %d rescued genes, split by cause", nrow(res)),
       subtitle = sprintf("composition = follows the mesenchymal score (r >= %.1f); consistent = rescued in all FA -> edited pairs", R_COMP),
       x = NULL, y = "genes") + theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
set.seed(42); res$y_jit <- res$pairs_ok + runif(nrow(res), -0.15, 0.15)
lab <- res %>% filter(class == "consistent rescue" | SYMBOL %in% c("FXN", SCORE_GENES)) %>% head(25)
p2 <- ggplot(res, aes(r_signed, y_jit, colour = class)) +
  geom_vline(xintercept = R_COMP, linetype = 2, colour = "grey50") + geom_point(size = 2) +
  scale_colour_manual(values = cls_col, name = NULL) +
  scale_y_continuous(breaks = 0:length(okc)) +
  labs(x = "correlation with mesenchymal score (in the direction of the FA effect)",
       y = "FA -> edited pairs moving towards control", title = "Each rescued gene") +
  theme_minimal(base_size = 11)
if (requireNamespace("ggrepel", quietly = TRUE))
  p2 <- p2 + ggrepel::geom_text_repel(data = lab, aes(label = SYMBOL), size = 3, show.legend = FALSE, max.overlaps = 30)
pdf(file.path(OUT_DIR, "42_summary.pdf"), width = 12, height = 5.5)
if (requireNamespace("patchwork", quietly = TRUE)) print(patchwork::wrap_plots(p1, p2, widths = c(1, 1.3))) else { print(p1); print(p2) }
invisible(dev.off())
say("\nSorties : %s", OUT_DIR)
