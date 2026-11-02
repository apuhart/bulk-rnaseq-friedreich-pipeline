#!/usr/bin/env Rscript
# =============================================================================
# 53_gsea_leave_one_donor_out.R — Les pathways GSEA tiennent-ils sans chaque lignée ? (v2, iNeurons + contribution par lignée)
# -----------------------------------------------------------------------------
# Pour un modèle (iPSC par défaut), refait la GSEA (Hallmark + GO BP + GO CC + GO MF) :
#   - FA vs contrôle      : toutes les lignées, puis en retirant UNE lignée FA à la fois
#                           (les 3 contrôles sont gardés) -> 4 GSEA
#   - édité vs FA apparié : toutes les paires, puis en retirant UNE paire à la fois -> 4 GSEA
#   - édité vs contrôle   : toutes les lignées éditées, puis en retirant UNE lignée à la fois (v2)
#   v2 : pour chaque pathway, CONTRIBUTION PAR LIGNÉE = NES(tout) - NES(sans la lignée)
#        (>0 : la lignée soutient le signal ; ~0 : neutre ; <0 : elle le tire dans l'autre sens)
#        + 'driver_line' (lignée dont le retrait réduit le plus |NES|) et 'single_line_driven'
#        (le pathway perd le sens ou la signif nominale quand on retire UNE lignée précise).
# Classement des gènes = statistique de Wald DESeq2 (même design que le dds d'origine).
# Un pathway est « robuste » s'il garde le même sens (signe du NES) dans les 3
# analyses sans un donneur, et reste au moins nominalement significatif (p < 0.05)
# dans chacune.
#
# Entrées : RNA_V4_ROOT/<DATASET>/<contraste>/<contraste>_dds.rds
# Sorties : OUT_DIR/53_loo_<contraste>.csv   (NES/padj complet + sans chaque donneur)
#           OUT_DIR/53_loo_focus.pdf          (pathways de la diapo, heatmap NES)
#           OUT_DIR/53_verdict.txt
# Options : DATASET=iPSC  RNA_V4_ROOT=...  OUT_DIR=...  MIN_SIZE=15  MAX_SIZE=500
#           GMT=fichier.gmt (si msigdbr absent)  FOCUS="GOBP_...,GOBP_..."
# Usage   : Rscript scripts/53_gsea_leave_one_donor_out.R
#           DATASET=iNeurons Rscript scripts/53_gsea_leave_one_donor_out.R
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(DESeq2); library(AnnotationDbi); library(org.Hs.eg.db) })
ROOT <- Sys.getenv("RNA_V4_ROOT", unset = file.path(getwd(), "results_v4"))
DS <- Sys.getenv("DATASET", "iNeurons")
OUT_DIR <- Sys.getenv("OUT_DIR", unset = file.path(Sys.getenv("PROJECT_ROOT", "."), "results", "gsea_loo", DS))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
MIN_SIZE <- as.numeric(Sys.getenv("MIN_SIZE", "15")); MAX_SIZE <- as.numeric(Sys.getenv("MAX_SIZE", "500"))
GMT <- Sys.getenv("GMT", "")
FOCUS_IPSC <- c(
  "GOBP_METALLO_SULFUR_CLUSTER_ASSEMBLY", "GOBP_IRON_SULFUR_CLUSTER_ASSEMBLY",
  "GOBP_MITOCHONDRIAL_TRANSLATION", "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
  "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_CHOLESTEROL_BIOSYNTHETIC_PROCESS", "GOBP_TRICARBOXYLIC_ACID_CYCLE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
# iNeurons : thèmes des slides (Fe-S / DNA-prolifération / OXPHOS / ribosome-traduction / synapse / inflammation)
FOCUS_INEU <- c(
  "GOBP_IRON_SULFUR_CLUSTER_ASSEMBLY", "GOBP_METALLO_SULFUR_CLUSTER_ASSEMBLY",
  "GOBP_DNA_REPLICATION", "GOBP_MITOTIC_CELL_CYCLE", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
  "GOBP_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "GOBP_MITOCHONDRIAL_TRANSLATION",
  "GOBP_CYTOPLASMIC_TRANSLATION", "GOBP_RIBOSOME_BIOGENESIS",
  "GOBP_SYNAPSE_ORGANIZATION", "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY", "GOBP_NEUROTRANSMITTER_TRANSPORT",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "GOBP_CHEMOKINE_MEDIATED_SIGNALING_PATHWAY", "GOBP_LYSOSOMAL_TRANSPORT")
FOCUS <- { e <- Sys.getenv("FOCUS", ""); if (nzchar(e)) trimws(strsplit(e, ",")[[1]]) else if (DS == "iPSC") FOCUS_IPSC else FOCUS_INEU }
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(53)

# --- gene sets ---------------------------------------------------------------------
sets <- list()
if (nzchar(GMT) && file.exists(GMT)) {
  for (l in readLines(GMT)) { x <- strsplit(l, "\t")[[1]]; sets[[x[1]]] <- x[-(1:2)] }
} else if (requireNamespace("msigdbr", quietly = TRUE)) {
  new_api <- utils::packageVersion("msigdbr") >= "10.0.0"
  get <- function(cat, sub) if (new_api) { if (is.null(sub)) msigdbr::msigdbr(species = "Homo sapiens", collection = cat)
    else msigdbr::msigdbr(species = "Homo sapiens", collection = cat, subcollection = sub) } else {
      if (is.null(sub)) msigdbr::msigdbr(species = "Homo sapiens", category = cat)
      else msigdbr::msigdbr(species = "Homo sapiens", category = cat, subcategory = sub) }
  m <- bind_rows(get("H", NULL), get("C5", if (new_api) "GO:BP" else "BP"),
                 get("C5", if (new_api) "GO:CC" else "CC"), get("C5", if (new_api) "GO:MF" else "MF"))   # v2 : + GO CC et MF
  sets <- split(m$gene_symbol, m$gs_name)
} else stop("Ni msigdbr ni GMT : installe msigdbr ou donne GMT=fichier.gmt")
use_fgsea <- requireNamespace("fgsea", quietly = TRUE)
say("%d gene sets ; méthode : %s", length(sets), if (use_fgsea) "fgseaMultilevel" else "limma::cameraPR (fgsea absent)")
run_gs <- function(stat) {
  stat <- stat[is.finite(stat)]
  if (use_fgsea) {
    r <- suppressWarnings(fgsea::fgseaMultilevel(sets, sort(stat, decreasing = TRUE), minSize = MIN_SIZE, maxSize = MAX_SIZE))
    tibble(pathway = r$pathway, NES = r$NES, pval = r$pval, padj = r$padj)
  } else {
    idx <- lapply(sets, function(s) which(names(stat) %in% s)); n <- lengths(idx); idx <- idx[n >= MIN_SIZE & n <= MAX_SIZE]
    r <- limma::cameraPR(stat, idx, sort = FALSE)
    tibble(pathway = rownames(r), NES = ifelse(r$Direction == "Up", 1, -1) * -log10(r$PValue), pval = r$PValue, padj = r$FDR)
  }
}

# --- DESeq2 sur un sous-ensemble d'échantillons -------------------------------------------
load_dds <- function(ct) { f <- file.path(ROOT, DS, ct, paste0(ct, "_dds.rds")); if (!file.exists(f)) stop("Introuvable : ", f); readRDS(f) }
lvl <- function(d, pat) { l <- levels(factor(colData(d)$condition)); l[grepl(pat, l, ignore.case = TRUE)][1] }
key <- function(x) { k <- sub("^\\D*(\\d+).*$", "\\1", x); ifelse(grepl("\\d", x), k, x) }   # F86, F86_R, F86ed, ED86 -> "86"
rank_stat <- function(d, keep, num_pat, den_pat) {
  d <- d[, keep]; colData(d) <- droplevels(colData(d))
  vars <- all.vars(design(d))
  ok <- vapply(vars, function(v) nlevels(factor(colData(d)[[v]])) > 1, logical(1))
  if (!all(ok)) design(d) <- as.formula(paste("~", paste(vars[ok], collapse = " + ")))
  d <- suppressMessages(DESeq(d, quiet = TRUE))
  r <- results(d, contrast = c("condition", lvl(d, num_pat), lvl(d, den_pat)))
  ens <- sub("\\..*$", "", rownames(r))
  sym <- suppressMessages(mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first"))
  tibble(SYMBOL = unname(sym), stat = r$stat) %>% filter(!is.na(SYMBOL), is.finite(stat)) %>%
    group_by(SYMBOL) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup() %>%
    { setNames(.$stat, .$SYMBOL) }
}

analyses <- list(
  list(name = "FA_vs_control", ct = "control_vs_FA", num = "^FA$|^FRDA$", den = "control|ctrl", drop_from = "^FA$|^FRDA$"),
  list(name = "edited_vs_control", ct = "control_vs_FA_edited", num = "edit", den = "control|ctrl", drop_from = "edit"),
  list(name = "edited_vs_FA", ct = "FA_edited_vs_FA", num = "edit", den = "^FA$|^FRDA$", drop_from = "pair"))
# noms de lignées (clé -> F4259, F4676, F86...) pris dans le dds FA vs control
{ d <- load_dds("control_vs_FA"); cn <- colnames(d); cd <- as.character(colData(d)$condition)
  fa_s <- cn[grepl("^FA$|^FRDA$", cd, ignore.case = TRUE)]
  line_name <- setNames(sub("[_.-]?R$", "", fa_s), key(fa_s)) }
summ <- list(); focus_all <- list(); donor_label <- c(all = "all donors")
for (a in analyses) {
  d0 <- load_dds(a$ct); cond <- as.character(colData(d0)$condition); s <- colnames(d0)
  say("\n== %s (%s) : %s", a$name, DS, paste(sprintf("%s[%s]", s, cond), collapse = ", "))
  donors <- if (a$drop_from == "pair") unique(key(s)) else unique(key(s[grepl(a$drop_from, cond, ignore.case = TRUE)]))
  runs <- list(all = s)
  for (dn in donors) {
    rm_s <- if (a$drop_from == "pair") s[key(s) == dn] else s[key(s) == dn & grepl(a$drop_from, cond, ignore.case = TRUE)]
    runs[[paste0("without_", dn)]] <- setdiff(s, rm_s)
    ln <- unname(line_name[dn]); donor_label[[paste0("without_", dn)]] <- paste("without", if (is.na(ln)) dn else ln)
  }
  res <- list()
  for (rn in names(runs)) {
    say("  %-14s %d échantillons", rn, length(runs[[rn]]))
    res[[rn]] <- run_gs(rank_stat(d0, runs[[rn]], a$num, a$den)) %>% mutate(run = rn)
  }
  R <- bind_rows(res)
  W <- R %>% dplyr::select(pathway, run, NES, pval, padj) %>%
    pivot_wider(names_from = run, values_from = c(NES, pval, padj), names_sep = "|")
  loo <- setdiff(names(runs), "all")
  sgn_all <- sign(W[["NES|all"]])
  same <- Reduce(`&`, lapply(loo, function(r) sign(W[[paste0("NES|", r)]]) == sgn_all))
  nom <- Reduce(`&`, lapply(loo, function(r) W[[paste0("pval|", r)]] < 0.05))
  W <- W %>% mutate(sig_all = `padj|all` < 0.05, same_sign_all_LOO = same, nominal_all_LOO = nom,
                    robust = sig_all & same & nom)
  # contribution par lignée : NES(tout) - NES(sans la lignée), signée dans le sens du NES complet
  for (r in loo) W[[paste0("contrib|", sub("without_", "", r))]] <- sgn_all * (W[["NES|all"]] - W[[paste0("NES|", r)]])
  cc <- as.matrix(W[paste0("contrib|", sub("without_", "", loo))])
  W$driver_line <- sub("contrib\\|", "", colnames(cc)[max.col(replace(cc, is.na(cc), -Inf), ties.method = "first")])
  W$driver_line <- unname(ifelse(is.na(line_name[W$driver_line]), W$driver_line, line_name[W$driver_line]))
  W$driver_share <- apply(cc, 1, function(x) { x <- pmax(x, 0); if (sum(x) == 0) NA else max(x) / sum(x) })   # ~1/3 = équilibré ; ~1 = une seule lignée
  lost <- Reduce(`|`, lapply(loo, function(r) sign(W[[paste0("NES|", r)]]) != sgn_all | W[[paste0("pval|", r)]] >= 0.05))
  W$single_line_driven <- W$sig_all & lost
  W <- W %>% arrange(desc(sig_all), `padj|all`)
  write_csv(W, file.path(OUT_DIR, sprintf("53_loo_%s.csv", a$name)))
  conc <- sapply(loo, function(r) suppressWarnings(cor(W[["NES|all"]], W[[paste0("NES|", r)]], method = "spearman", use = "complete.obs")))
  ns <- sum(W$sig_all, na.rm = TRUE)
  summ[[a$name]] <- tibble(analysis = a$name, sig_full = ns,
                           robust = sum(W$robust, na.rm = TRUE), n_single_line = sum(W$single_line_driven, na.rm = TRUE), same_sign = sum(W$sig_all & W$same_sign_all_LOO, na.rm = TRUE),
                           rho = paste(sprintf("%s %.2f", sub("without_", "-", loo), conc), collapse = " | "))
  focus_all[[a$name]] <- R %>% filter(pathway %in% FOCUS) %>% mutate(analysis = a$name)
}
summ <- bind_rows(summ); print(as.data.frame(summ), row.names = FALSE)

fa <- bind_rows(focus_all)
if (nrow(fa)) {
  run_lv <- unique(c("all", sort(unique(fa$run[fa$run != "all"]))))
  fa <- fa %>% mutate(
    label = ifelse(grepl("^HALLMARK_", pathway), paste0("Hallmark: ", tolower(gsub("_", " ", sub("^HALLMARK_", "", pathway)))),
                   tolower(gsub("_", " ", sub("^GO(BP|CC|MF)_", "", pathway)))),
    sig = ifelse(padj < 0.05, "**", ifelse(pval < 0.05, "*", "")), txt = sprintf("%+.2f%s", NES, sig),
    run = factor(unname(donor_label[run]), levels = unname(donor_label[run_lv])),
    analysis = factor(recode(analysis, FA_vs_control = "FA vs control", edited_vs_control = "Edited vs control", edited_vs_FA = "Edited vs FA (paired)"),
                      levels = c("FA vs control", "Edited vs control", "Edited vs FA (paired)")))
  ord <- fa %>% filter(analysis == "FA vs control", run == "all donors") %>% arrange(NES) %>% pull(label)
  fa$label <- factor(fa$label, levels = rev(unique(c(ord, as.character(fa$label)))))
  p <- ggplot(fa, aes(run, label, fill = NES)) + geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = txt, colour = abs(NES) > 1.8), size = 3.3, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-2.6, 2.6), oob = scales::squish) +
    facet_wrap(~ analysis) +
    labs(title = sprintf("%s: key pathways when leaving one donor out", DS),
         subtitle = "Numbers = NES (GSEA, GO BP / Hallmark); ** padj < 0.05, * p < 0.05 (nominal)", x = NULL, y = NULL) +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
                                          strip.text = element_text(face = "bold", size = 12))
  hh <- 2 + 0.35 * n_distinct(fa$label)
  ggsave(file.path(OUT_DIR, "53_loo_focus.pdf"), p, width = 16, height = hh)
  ggsave(file.path(OUT_DIR, "53_loo_focus.png"), p, width = 16, height = hh, dpi = 200, bg = "white")
}

v <- c(sprintf("=== 53 : GSEA sans chaque donneur (%s) ===", DS), "")
for (i in seq_len(nrow(summ))) with(summ[i, ], v <<- c(v,
                                                       sprintf("%s : %d pathways padj < 0.05 dans l'analyse complète", analysis, sig_full),
                                                       sprintf("   même sens sans chacun des donneurs : %d ; et p < 0.05 sans chacun : %d (robustes) ; perdus en retirant UNE lignée : %d", same_sign, robust, n_single_line),
                                                       sprintf("   corrélation des NES (complet vs sans donneur) : %s", rho)))
if (nrow(fa)) {
  v <- c(v, "", "Pathways de la diapo (NES ; ** padj<0.05, * p<0.05) :")
  for (an in unique(fa$analysis)) for (pw in unique(fa$label[fa$analysis == an])) {
    x <- fa %>% filter(analysis == an, label == pw) %>% arrange(run)
    v <- c(v, sprintf("  %-22s %-45s %s", an, substr(pw, 1, 45), paste(sprintf("%s %+.2f%s", sub("without_", "-", x$run), x$NES, x$sig), collapse = " | ")))
  }
}
v <- c(v, "", "Lecture : un pathway qui garde son signe et reste au moins nominal sans chacun des donneurs",
       "ne repose pas sur une seule lignée. Voir 53_loo_<contraste>.csv : contrib|<lignée>, driver_line, driver_share, single_line_driven.",
"NB : en retirant une lignée il ne reste que 2 lignées (n=2) : perte de puissance attendue, surtout pour edited_vs_FA.")
writeLines(v, file.path(OUT_DIR, "53_verdict.txt")); cat("\n"); cat(v, sep = "\n")
say("\nSorties : %s", OUT_DIR)