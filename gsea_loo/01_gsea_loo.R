#!/usr/bin/env Rscript
# =============================================================================
# 01_gsea_loo.R — UNE seule source pour toutes les analyses de pathways
# -----------------------------------------------------------------------------
# Pour chaque modèle (fibro, iPSC, iNeurons) et chaque contraste, refait la GSEA
# (fgsea ; Hallmark + GO BP + GO CC + GO MF) sur :
#   - tous les échantillons (« all »)
#   - puis en retirant UNE lignée à la fois (« without_<lignée> »)
# Classement des gènes = statistique de Wald DESeq2 (même design que le dds d'origine).
#
# Contrastes (toujours dans ce sens : NES > 0 = plus haut dans le PREMIER groupe) :
#   A = FA vs control               (retrait d'une lignée FA)
#   B = edited vs control           (retrait de la lignée éditée)         [pas en fibro]
#   C = FA vs edited                (retrait de la paire FA + éditée)     [pas en fibro]
#
# Entrées : ROOT/<DATASET>/<contraste>/<contraste>_dds.rds
#           dossiers de contraste : control_vs_FA, control_vs_FA_edited, FA_edited_vs_FA
# Sorties : OUT_DIR/<DATASET>_<A|B|C>.csv  (un fichier par modèle et contraste) :
#           pathway, collection, size, leadingEdge (analyse « all »),
#           NES|run, pval|run, padj|run pour chaque run
#           OUT_DIR/summary.txt
# Options (variables d'environnement) :
#   RNA_V4_ROOT=...  OUT_DIR=...  DATASETS=fibro,iPSC,iNeurons  MIN_SIZE=15  MAX_SIZE=500
#   GMT=fichier.gmt (si msigdbr absent)
# Usage : Rscript 01_gsea_loo.R
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(tibble)
  library(DESeq2); library(AnnotationDbi); library(org.Hs.eg.db); library(fgsea) })
ROOT <- Sys.getenv("RNA_V4_ROOT", unset = file.path(getwd(), "results_v4"))
OUT_DIR <- Sys.getenv("OUT_DIR", unset = "gsea_loo"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
DATASETS <- strsplit(Sys.getenv("DATASETS", "fibro,iPSC,iNeurons"), ",")[[1]]
MIN_SIZE <- as.numeric(Sys.getenv("MIN_SIZE", "15")); MAX_SIZE <- as.numeric(Sys.getenv("MAX_SIZE", "500"))
GMT <- Sys.getenv("GMT", "")
say <- function(...) cat(sprintf(...), "\n", sep = "")

# --- gene sets : Hallmark + GO BP + GO CC + GO MF -------------------------------------
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
                 get("C5", if (new_api) "GO:CC" else "CC"), get("C5", if (new_api) "GO:MF" else "MF"))
  sets <- split(m$gene_symbol, m$gs_name)
} else stop("Ni msigdbr ni GMT : installe msigdbr ou donne GMT=fichier.gmt")
collection_of <- function(p) ifelse(grepl("^HALLMARK_", p), "H", paste0("C5_", substr(p, 3, 4)))   # GOBP -> C5_BP
say("%d gene sets", length(sets))

run_gs <- function(stat) {
  stat <- sort(stat[is.finite(stat)], decreasing = TRUE)
  set.seed(53)                                   # même graine à chaque run -> résultats reproductibles
  r <- suppressWarnings(fgsea::fgseaMultilevel(sets, stat, minSize = MIN_SIZE, maxSize = MAX_SIZE))
  tibble(pathway = r$pathway, size = r$size, leadingEdge = vapply(r$leadingEdge, paste, "", collapse = ";"),
         NES = r$NES, pval = r$pval, padj = r$padj)
}

# --- DESeq2 sur un sous-ensemble d'échantillons -------------------------------------------
load_dds <- function(ds, ct) { f <- file.path(ROOT, ds, ct, paste0(ct, "_dds.rds")); if (!file.exists(f)) stop("Introuvable : ", f); readRDS(f) }
lvl <- function(d, pat) { l <- levels(factor(colData(d)$condition)); l[grepl(pat, l, ignore.case = TRUE)][1] }
key <- function(x) { k <- sub("^\\D*(\\d+).*$", "\\1", x); ifelse(grepl("\\d", x), k, x) }   # F86, F86_R, F86ed, ED86 -> "86"
rank_stat <- function(d, keep, num_pat, den_pat) {
  d <- d[, keep]; colData(d) <- droplevels(colData(d))
  vars <- all.vars(design(d))
  ok <- vapply(vars, function(v) nlevels(factor(colData(d)[[v]])) > 1, logical(1))
  if (!all(ok)) design(d) <- as.formula(paste("~", paste(vars[ok], collapse = " + ")))
  d <- tryCatch(suppressMessages(DESeq(d, quiet = TRUE)), error = function(e) {
    message("  DESeq2 impossible avec le design d'origine (", conditionMessage(e), ") -> ~ condition")
    design(d) <- ~ condition; suppressMessages(DESeq(d, quiet = TRUE)) })
  r <- results(d, contrast = c("condition", lvl(d, num_pat), lvl(d, den_pat)))
  ens <- sub("\\..*$", "", rownames(r))
  sym <- suppressMessages(mapIds(org.Hs.eg.db, ens, "SYMBOL", "ENSEMBL", multiVals = "first"))
  tibble(SYMBOL = unname(sym), stat = r$stat) %>% filter(!is.na(SYMBOL), is.finite(stat)) %>%
    group_by(SYMBOL) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup() %>%
    { setNames(.$stat, .$SYMBOL) }
}

FA <- "^FA$|^FRDA$"
ANALYSES <- list(
  A = list(ct = "control_vs_FA",        num = FA,     den = "control|ctrl", drop = FA),
  B = list(ct = "control_vs_FA_edited", num = "edit", den = "control|ctrl", drop = "edit"),
  C = list(ct = "FA_edited_vs_FA",      num = FA,     den = "edit",         drop = "pair"))
summ <- list()
for (ds in DATASETS) for (cn in names(ANALYSES)) {
  if (ds == "fibro" && cn != "A") next                       # pas de lignées éditées en fibroblastes
  a <- ANALYSES[[cn]]
  d0 <- load_dds(ds, a$ct); cond <- as.character(colData(d0)$condition); s <- colnames(d0)
  say("\n== %s contraste %s (%s) : %s", ds, cn, a$ct, paste(sprintf("%s[%s]", s, cond), collapse = ", "))
  donors <- if (a$drop == "pair") unique(key(s)) else unique(key(s[grepl(a$drop, cond, ignore.case = TRUE)]))
  runs <- list(all = s)
  for (dn in donors) {
    rm_s <- if (a$drop == "pair") s[key(s) == dn] else s[key(s) == dn & grepl(a$drop, cond, ignore.case = TRUE)]
    runs[[paste0("without_", dn)]] <- setdiff(s, rm_s)
  }
  res <- list()
  for (rn in names(runs)) { say("  %-16s %d échantillons", rn, length(runs[[rn]]))
    res[[rn]] <- run_gs(rank_stat(d0, runs[[rn]], a$num, a$den)) %>% mutate(run = rn) }
  R <- bind_rows(res)
  info <- R %>% filter(run == "all") %>% dplyr::select(pathway, size, leadingEdge) %>% mutate(collection = collection_of(pathway))
  W <- R %>% dplyr::select(pathway, run, NES, pval, padj) %>%
    pivot_wider(names_from = run, values_from = c(NES, pval, padj), names_sep = "|")
  W <- info %>% left_join(W, by = "pathway") %>% dplyr::select(pathway, collection, size, leadingEdge, everything())
  write_csv(W, file.path(OUT_DIR, sprintf("%s_%s.csv", ds, cn)))
  loo <- setdiff(names(runs), "all")
  rho <- sapply(loo, function(r) suppressWarnings(cor(W[["NES|all"]], W[[paste0("NES|", r)]], method = "spearman", use = "complete.obs")))
  summ[[paste(ds, cn)]] <- tibble(model = ds, contrast = cn, n_pathways = nrow(W), sig_all = sum(W[["padj|all"]] < 0.05, na.rm = TRUE),
    rho = paste(sprintf("%s %.2f", sub("without_", "-", loo), rho), collapse = " | "))
}
summ <- bind_rows(summ); print(as.data.frame(summ), row.names = FALSE)
writeLines(capture.output(print(as.data.frame(summ), row.names = FALSE)), file.path(OUT_DIR, "summary.txt"))
say("Terminé -> %s", OUT_DIR)
