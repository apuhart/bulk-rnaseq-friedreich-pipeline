#!/usr/bin/env Rscript
# =============================================================================
# 52_paired_rescue_allgenes.R — Rescue apparié sur TOUS les gènes (iNeurons)
# -----------------------------------------------------------------------------
# Idée : au lieu de partir des DEG (|log2FC| >= 1, dominés par la composition
# et les clones), on cherche les gènes qui reviennent vers le contrôle après
# correction de FXN, de façon REPRODUCTIBLE chez les 3 donneurs (chaque lignée
# FA comparée à SA lignée éditée), sans seuil d'amplitude, et hors composition.
#
# Pour chaque gène exprimé :
#   lfc_dis / p_dis   FA vs contrôle (DESeq2, log2FC NON rétréci)
#   lfc_pair / p_pair édité vs FA, modèle apparié (dds FA_edited_vs_FA, ~ lignée + condition)
#   delta_<donneur>   édité - FA pour chaque paire (log2 comptages normalisés)
#   donors_back       nb de paires qui reviennent vers le contrôle (|delta| >= DELTA_MIN,
#                     sens opposé à l'effet FA)
#   rescue_frac       -lfc_pair / lfc_dis (1 = retour complet au niveau contrôle)
#   r_comp            corrélation (9 échantillons) avec le score mésenchyme, orientée
#                     dans le sens de l'effet FA
#   r_within          même corrélation APRÈS retrait de l'effet condition (intra-groupe) :
#                     composition = r_comp >= R_COMP ET r_within >= R_WITHIN. Sans ce second
#                     critère, tout vrai gène FA (FXN !) serait exclu, car le score est
#                     lui-même plus haut chez les FA.
#   fibro_ratio       (optionnel, FIBRO_DDS) expression fibroblastes / iNeurons contrôle :
#                     >= FIBRO_RATIO = gène « de fibroblaste », gonflable par contamination
# Niveaux :
#   Tier 1  p_dis < P_DIS, 3/3 donneurs, p_pair < P_PAIR, hors composition
#   Tier 2  p_dis < P_DIS, 3/3 donneurs, p_pair < 0.2  OU  2/3 donneurs et p_pair < P_PAIR
# Bruit de fond : les mêmes critères appliqués dans le sens INVERSE (l'édition
# éloigne encore plus du contrôle chez les 3 donneurs) donnent le nombre de
# gènes attendus par hasard / par effet de clone -> à comparer au Tier 1.
#
# Entrées : RNA_V4_ROOT/<DATASET>/<contraste>/<contraste>_dds.rds (3 contrastes)
# Sorties : OUT_DIR/52_paired_rescue_all_genes.csv, 52_candidates.xlsx,
#           52_paired_rescue.pdf, 52_verdict.txt
# Options : DATASET=iNeurons  RNA_V4_ROOT=...  OUT_DIR=...  P_DIS=0.05  P_PAIR=0.05
#           DELTA_MIN=0.2  R_COMP=0.6  R_WITHIN=0.5  MIN_COUNT=10  SCORE_GENES="COL1A1,..."
#           FIBRO_DDS=/chemin/fibro/control_vs_FA/control_vs_FA_dds.rds  FIBRO_RATIO=4
# Usage   : Rscript scripts/52_paired_rescue_allgenes.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
  library(DESeq2); library(AnnotationDbi); library(org.Hs.eg.db)
})
ROOT <- Sys.getenv("RNA_V4_ROOT", unset = file.path(getwd(), "results_v4"))
DS <- Sys.getenv("DATASET", "iNeurons")
OUT_DIR <- Sys.getenv("OUT_DIR", unset = file.path(Sys.getenv("PROJECT_ROOT", "."), "results", "paired_rescue", DS))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
num <- function(v, d) as.numeric(Sys.getenv(v, d))
P_DIS <- num("P_DIS", "0.05"); P_PAIR <- num("P_PAIR", "0.05"); DELTA_MIN <- num("DELTA_MIN", "0.2")
R_COMP <- num("R_COMP", "0.6"); MIN_COUNT <- num("MIN_COUNT", "10"); FIBRO_RATIO <- num("FIBRO_RATIO", "4")
SCORE_GENES <- { e <- Sys.getenv("SCORE_GENES", ""); if (nzchar(e)) trimws(strsplit(e, ",")[[1]]) else
  c("COL1A1", "COL1A2", "FN1", "THY1", "ACTA2", "TAGLN") }
FIBRO_DDS <- Sys.getenv("FIBRO_DDS", "")
CT <- c(DIS = "control_vs_FA", ED = "control_vs_FA_edited", PAIR = "FA_edited_vs_FA")
say <- function(...) cat(sprintf(...), "\n", sep = "")
strip <- function(x) sub("\\..*$", "", x)

# ---------------------------------------------------------------------------
# 1. Statistiques DESeq2 (log2FC non rétrécis, contrastes explicites)
# ---------------------------------------------------------------------------
load_dds <- function(ct) {
  f <- file.path(ROOT, DS, ct, paste0(ct, "_dds.rds"))
  if (!file.exists(f)) stop("Introuvable : ", f, "  -> ajuste RNA_V4_ROOT / DATASET")
  d <- readRDS(f)
  fitted <- tryCatch(length(resultsNames(d)) > 0, error = function(e) FALSE)
  if (!fitted) { say("  [%s] DESeq() (dds non ajusté)", ct); d <- DESeq(d, quiet = TRUE) }
  d
}
lvl <- function(d, pat) { l <- levels(factor(colData(d)$condition)); l[grepl(pat, l, ignore.case = TRUE)][1] }
get_res <- function(d, num_pat, den_pat) {
  num <- lvl(d, num_pat); den <- lvl(d, den_pat)
  if (is.na(num) || is.na(den)) stop("Niveaux de condition introuvables : ", paste(levels(factor(colData(d)$condition)), collapse = ", "))
  r <- results(d, contrast = c("condition", num, den))
  say("  contraste %s - %s : %d gènes testés", num, den, sum(!is.na(r$pvalue)))
  tibble(GENEID = strip(rownames(r)), baseMean = r$baseMean, lfc = r$log2FoldChange, p = r$pvalue, padj = r$padj)
}
say("== %s : statistiques ==", DS)
d_dis <- load_dds(CT[["DIS"]]); d_pair <- load_dds(CT[["PAIR"]]); d_ed <- readRDS(file.path(ROOT, DS, CT[["ED"]], paste0(CT[["ED"]], "_dds.rds")))
dis <- get_res(d_dis, "^FA$|^FRDA$", "control|ctrl")
pair <- get_res(d_pair, "edit", "^FA$|^FRDA$")
say("  modèle apparié : %s", paste(deparse(design(d_pair)), collapse = ""))

# ---------------------------------------------------------------------------
# 2. Expression normalisée des 9 échantillons, paires, score de composition
# ---------------------------------------------------------------------------
dl <- list(d_dis, d_ed, d_pair)
common <- Reduce(intersect, lapply(dl, rownames))
cnt <- do.call(cbind, lapply(dl, function(d) counts(d)[common, , drop = FALSE])); cnt <- cnt[, !duplicated(colnames(cnt)), drop = FALSE]
cd <- do.call(rbind, lapply(dl, function(d) as.data.frame(colData(d))[, "condition", drop = FALSE]))
cond <- setNames(as.character(cd$condition), rownames(cd))[colnames(cnt)]
nc <- t(t(cnt) / estimateSizeFactorsForMatrix(cnt)); rownames(nc) <- strip(rownames(nc)); nc <- nc[!duplicated(rownames(nc)), , drop = FALSE]
lex <- log2(nc + 1)
sym <- suppressMessages(mapIds(org.Hs.eg.db, rownames(lex), "SYMBOL", "ENSEMBL", multiVals = "first"))
say("Échantillons : %s", paste(sprintf("%s[%s]", colnames(lex), cond), collapse = ", "))

is_ed <- grepl("edit", cond, ignore.case = TRUE)
is_fa <- !is_ed & !grepl("control|ctrl", cond, ignore.case = TRUE)
key <- function(x) gsub("[^A-Z0-9]", "", gsub("[_.-]?(EDITED|EDIT|ED|E|R)[_.-]?[0-9]*$", "", toupper(x)))
digits <- function(x) sub("^\\D*(\\d+).*$", "\\1", x)
fa_s <- colnames(lex)[is_fa]; ed_s <- colnames(lex)[is_ed]
pairs <- do.call(rbind, lapply(fa_s, function(a) {
  b <- ed_s[key(ed_s) == key(a)]
  if (length(b) != 1) b <- ed_s[grepl("\\d", ed_s) & digits(ed_s) == digits(a)]
  if (length(b) == 1) c(donor = key(a), FA = a, ED = b) else NULL }))
if (is.null(pairs) || nrow(pairs) < 2) stop("Paires FA -> édité non reconnues (FA : ", paste(fa_s, collapse = ", "), " | édités : ", paste(ed_s, collapse = ", "), ")")
say("Paires : %s", paste(sprintf("%s -> %s", pairs[, "FA"], pairs[, "ED"]), collapse = ", "))

sid <- names(sym)[!is.na(sym) & sym %in% SCORE_GENES]
score <- if (length(sid) >= 3) colMeans(lex[sid, , drop = FALSE]) else NULL
if (is.null(score)) say("ATTENTION : moins de 3 gènes du score mésenchyme -> pas de filtre de composition")

# ---------------------------------------------------------------------------
# 3. Table par gène
# ---------------------------------------------------------------------------
tab <- dis %>% dplyr::rename(lfc_dis = lfc, p_dis = p, padj_dis = padj) %>%
  inner_join(pair %>% dplyr::select(GENEID, lfc_pair = lfc, p_pair = p, padj_pair = padj), by = "GENEID") %>%
  filter(GENEID %in% rownames(lex), !is.na(p_dis), !is.na(p_pair)) %>%
  filter(rowMeans(nc[GENEID, , drop = FALSE]) >= MIN_COUNT) %>%
  mutate(SYMBOL = coalesce(unname(sym[GENEID]), GENEID))
for (k in seq_len(nrow(pairs))) tab[[paste0("delta_", pairs[k, "donor"])]] <- unname(lex[tab$GENEID, pairs[k, "ED"]] - lex[tab$GENEID, pairs[k, "FA"]])
dcols <- grep("^delta_", names(tab), value = TRUE); D <- as.matrix(tab[, dcols])
sdir <- sign(tab$lfc_dis)
tab$donors_back <- rowSums(sign(D) == -sdir & abs(D) >= DELTA_MIN)
tab$donors_away <- rowSums(sign(D) == sdir & abs(D) >= DELTA_MIN)      # pour le bruit de fond
tab$min_delta_back <- apply(-sdir * D, 1, min)                           # le plus faible retour (> 0 = tous reviennent)
tab$rescue_frac <- -tab$lfc_pair / tab$lfc_dis
tab$r_comp <- if (!is.null(score)) as.numeric(cor(t(lex[tab$GENEID, names(score)]), score)) * sdir else NA_real_
# Corrélation INTRA-groupe : le score mésenchyme est lui-même plus haut chez les FA, donc tout vrai
# gène FA (ex. FXN) corrèle avec lui « par la condition ». On retire l'effet condition (résidus de
# gène ~ condition et de score ~ condition) : un gène de composition suit encore le score à
# l'intérieur des groupes (F4676 vs F86...), un vrai gène FA non.
if (!is.null(score)) {
  cf <- factor(cond[names(score)])
  res_score <- residuals(lm(score ~ cf))
  M <- lex[tab$GENEID, names(score), drop = FALSE]
  grp_mean <- sapply(levels(cf), function(l) rowMeans(M[, cf == l, drop = FALSE]))
  RM <- M - grp_mean[, as.character(cf), drop = FALSE]
  tab$r_within <- as.numeric(cor(t(RM), res_score)) * sdir
} else tab$r_within <- NA_real_

if (nzchar(FIBRO_DDS) && file.exists(FIBRO_DDS)) {
  fd <- readRDS(FIBRO_DDS); fc <- counts(fd); fn <- t(t(fc) / estimateSizeFactorsForMatrix(fc)); rownames(fn) <- strip(rownames(fn))
  fctrl <- grepl("control|ctrl", colData(fd)$condition, ignore.case = TRUE)
  fib <- rowMeans(fn[, fctrl, drop = FALSE]); ctrl_i <- rowMeans(nc[, grepl("control|ctrl", cond, ignore.case = TRUE), drop = FALSE])
  tab$fibro_ratio <- (fib[tab$GENEID] + 1) / (ctrl_i[tab$GENEID] + 1)
  say("Référence fibroblastes : %s", FIBRO_DDS)
} else tab$fibro_ratio <- NA_real_

np <- length(dcols)
R_WITHIN <- num("R_WITHIN", "0.5")
comp_flag <- (!is.na(tab$r_comp) & tab$r_comp >= R_COMP & !is.na(tab$r_within) & tab$r_within >= R_WITHIN) |
  (!is.na(tab$fibro_ratio) & tab$fibro_ratio >= FIBRO_RATIO)
tab <- tab %>% mutate(
  composition = comp_flag,
  direction_ok = sign(lfc_pair) == -sign(lfc_dis),
  tier = case_when(
    p_dis < P_DIS & direction_ok & donors_back == np & p_pair < P_PAIR & !composition ~ "Tier 1",
    p_dis < P_DIS & direction_ok & !composition & ((donors_back == np & p_pair < 0.2) | (donors_back == np - 1 & p_pair < P_PAIR)) ~ "Tier 2",
    p_dis < P_DIS & direction_ok & donors_back >= np - 1 & p_pair < 0.2 & composition ~ "composition",
    TRUE ~ ""),
  null_tier1 = p_dis < P_DIS & sign(lfc_pair) == sign(lfc_dis) & donors_away == np & p_pair < P_PAIR & !composition,
  score = ifelse(direction_ok, -log10(p_dis) - log10(p_pair), 0)) %>%
  arrange(factor(tier, levels = c("Tier 1", "Tier 2", "composition", "")), desc(score))
out <- tab %>% transmute(Gene = SYMBOL, Tier = tier, `LFC FA vs ctrl` = round(lfc_dis, 2), `p FA vs ctrl` = signif(p_dis, 2),
                         `LFC edited vs FA (paired)` = round(lfc_pair, 2), `p paired` = signif(p_pair, 2), `rescue fraction` = round(rescue_frac, 2),
                         `donors back (/3)` = donors_back, across(all_of(dcols), ~ round(.x, 2)), `r composition` = round(r_comp, 2), `r within groups` = round(r_within, 2),
                         `fibro / iNeuron ctrl` = round(fibro_ratio, 1), baseMean = round(baseMean), Ensembl = GENEID)
write_csv(out, file.path(OUT_DIR, "52_paired_rescue_all_genes.csv"))

n1 <- sum(tab$tier == "Tier 1"); n2 <- sum(tab$tier == "Tier 2"); nnull <- sum(tab$null_tier1); ncomp <- sum(tab$tier == "composition")
say("\n%d gènes exprimés testés | Tier 1 : %d | Tier 2 : %d | composition : %d | bruit de fond (sens inverse) : %d",
    nrow(tab), n1, n2, ncomp, nnull)
fxn <- tab %>% filter(SYMBOL == "FXN")
if (nrow(fxn)) say("FXN : %s, LFC FA %.2f (p %.1g), apparié %+.2f (p %.1g), %d/%d donneurs", ifelse(fxn$tier == "", "non classé", fxn$tier),
                   fxn$lfc_dis, fxn$p_dis, fxn$lfc_pair, fxn$p_pair, fxn$donors_back, np)
say("Tier 1 : %s", paste(head(tab$SYMBOL[tab$tier == "Tier 1"], 60), collapse = ", "))

# ---------------------------------------------------------------------------
# 4. Excel + figures + verdict
# ---------------------------------------------------------------------------
defs <- tibble(Term = c("Tier 1", "Tier 2", "composition", "rescue fraction", "donors back", "Background"),
               Definition = c(sprintf("FA vs ctrl p < %g; paired edited vs FA p < %g in the opposite direction; all %d donors move back (|delta log2| >= %g); not composition-related", P_DIS, P_PAIR, np, DELTA_MIN),
                              sprintf("same, but paired p < 0.2 with %d/%d donors, or %d/%d donors with paired p < %g", np, np, np - 1, np, P_PAIR),
                              sprintf("would be Tier 1/2 but follows the mesenchymal score across samples (r >= %g) AND within groups (r >= %g)%s", R_COMP, R_WITHIN, if (nzchar(FIBRO_DDS)) sprintf(" or is >= %gx higher in fibroblasts", FIBRO_RATIO) else ""),
                              "-paired LFC / FA LFC: 1 = back to control level, 0 = no change",
                              "number of FA -> edited pairs (same donor) moving back toward control",
                              sprintf("%d genes pass the Tier-1 criteria in the WRONG direction (editing moves them further from control in all donors): expected false positives / clone effects", nnull)))
sheets <- list(`Tier 1` = filter(out, Tier == "Tier 1"), `Tier 2` = filter(out, Tier == "Tier 2"),
               Composition = filter(out, Tier == "composition"), Definitions = defs)
xl <- file.path(OUT_DIR, "52_candidates.xlsx")
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook(); hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#1F4E5F", fontColour = "#FFFFFF")
  for (nm in names(sheets)) { openxlsx::addWorksheet(wb, nm); openxlsx::writeData(wb, nm, sheets[[nm]], headerStyle = hs)
    openxlsx::freezePane(wb, nm, firstRow = TRUE); openxlsx::setColWidths(wb, nm, seq_len(ncol(sheets[[nm]])), "auto") }
  openxlsx::saveWorkbook(wb, xl, overwrite = TRUE)
} else writexl::write_xlsx(sheets, xl)

pd <- tab %>% mutate(cls = ifelse(tier == "", "other", tier),
                     cls = factor(cls, levels = c("Tier 1", "Tier 2", "composition", "other")))
lab <- pd %>% filter(tier == "Tier 1" | SYMBOL == "FXN") %>% slice_max(score, n = 25, with_ties = FALSE) %>% bind_rows(filter(pd, SYMBOL == "FXN")) %>% distinct(GENEID, .keep_all = TRUE)
p1 <- ggplot(pd, aes(lfc_dis, lfc_pair)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(data = filter(pd, cls == "other"), colour = "grey85", size = 0.5) +
  geom_point(data = filter(pd, cls != "other"), aes(colour = cls), size = 1.6) +
  scale_colour_manual(values = c(`Tier 1` = "#0F6E56", `Tier 2` = "#74C476", composition = "#B2182B"), name = NULL) +
  coord_cartesian(xlim = quantile(pd$lfc_dis, c(0.001, 0.999)), ylim = quantile(pd$lfc_pair, c(0.001, 0.999))) +
  labs(title = sprintf("%s: genes moving back after FXN correction, all genes", DS),
       subtitle = sprintf("Dashed line = full rescue (edited back to control level). Tier 1: %d genes | background (wrong direction): %d", n1, nnull),
       x = "log2FC FA vs control", y = "log2FC edited vs FA (paired)") +
  theme_minimal(base_size = 11)
if (requireNamespace("ggrepel", quietly = TRUE)) p1 <- p1 + ggrepel::geom_text_repel(data = lab, aes(label = SYMBOL), size = 2.8, max.overlaps = 40, seed = 1)
top <- tab %>% filter(tier %in% c("Tier 1", "Tier 2")) %>% slice_max(score, n = 40, with_ties = FALSE)
p2 <- NULL
if (nrow(top)) {
  hm <- top %>% dplyr::select(SYMBOL, all_of(dcols), lfc_dis) %>% pivot_longer(all_of(dcols), names_to = "donor", values_to = "delta") %>%
    mutate(back = -sign(lfc_dis) * delta, donor = sub("^delta_", "", donor), SYMBOL = factor(SYMBOL, levels = rev(unique(top$SYMBOL))))
  p2 <- ggplot(hm, aes(donor, SYMBOL, fill = back)) + geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#0F6E56", midpoint = 0, name = "move back\n(log2)") +
    labs(title = "Top candidates: change toward control in each donor (edited - FA)", x = NULL, y = NULL) + theme_minimal(base_size = 10)
}
pdf(file.path(OUT_DIR, "52_paired_rescue.pdf"), width = 10, height = 8); print(p1); if (!is.null(p2)) print(p2); invisible(dev.off())

v <- c(sprintf("=== 52 : rescue apparié sur tous les gènes (%s) ===", DS), "",
       sprintf("Gènes exprimés testés : %d (moyenne >= %g comptes normalisés)", nrow(tab), MIN_COUNT),
       sprintf("Tier 1 (3/3 donneurs, p apparié < %g, hors composition) : %d", P_PAIR, n1),
       sprintf("Tier 2 : %d   |   composition (exclus) : %d", n2, ncomp),
       sprintf("Bruit de fond : %d gènes passent les mêmes critères dans le MAUVAIS sens", nnull),
       sprintf("  -> rapport Tier 1 / bruit de fond = %s (>> 1 : signal réel de la correction ; ~1 : surtout du bruit / des effets de clone)",
               if (nnull > 0) sprintf("%.1f", n1 / nnull) else "inf"), "",
       if (nrow(fxn)) sprintf("FXN (contrôle positif) : %s, apparié %+.2f (p %.1g), %d/%d donneurs, fraction de rescue %.2f",
                              ifelse(fxn$tier == "", "non classé", fxn$tier), fxn$lfc_pair, fxn$p_pair, fxn$donors_back, np, fxn$rescue_frac) else "FXN absent",
       "", "Tier 1 (classés par score) :", paste(" ", head(tab$SYMBOL[tab$tier == "Tier 1"], 80), collapse = ""))
writeLines(v, file.path(OUT_DIR, "52_verdict.txt")); cat("\n"); cat(v, sep = "\n")
say("\nSorties : %s", OUT_DIR)

