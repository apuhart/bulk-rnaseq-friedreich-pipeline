#!/usr/bin/env Rscript
# =============================================================================
# 51_rescue_tables_ipsc_ineurons.R  — Tableau récapitulatif des gènes rescued /
#                                     non rescued, iPSC et iNeurons (bulk RNA v4)
# -----------------------------------------------------------------------------
# Mêmes règles que le pipeline bulk (run_rescue_scatter) :
#   FA-DEG          padj(FA vs ctrl) < 0.05 et |log2FC| >= 1
#   Rescued         FA-DEG, NON significatif en édité vs ctrl, |LFC édité| < |LFC FA|
#   Partial rescue  encore significatif en édité vs ctrl, mais |LFC édité| < |LFC FA|
#   Reversed        signe opposé et significatif en édité vs ctrl
#   Not rescued     encore significatif, |LFC édité| >= |LFC FA|
# Ajoute, pour chaque gène :
#   - LFC FA vs ctrl, édité vs ctrl, édité vs FA (apparié) et padj
#   - dans combien des 3 paires FA -> édité (même donneur) le gène va dans le
#     sens du rescue (|Δ log2| >= DELTA_MIN) ; quelle paire porte le signal
#   - corrélation avec le score mésenchyme (COL1A1, COL1A2, FN1, THY1, ACTA2,
#     TAGLN), dans le sens de l'effet FA : r >= R_COMP = « composition »
#   - une colonne « Confidence » lisible (robust / partial / one donor /
#     composition) pour les gènes rescued et partial
#
# Entrées : RNA_V4_ROOT/<dataset>/<contraste>/<contraste>_DEG_results.csv et _dds.rds
#           (contrastes control_vs_FA, control_vs_FA_edited, FA_edited_vs_FA)
#           Convention du pipeline : dossier « A_vs_B » -> log2FC = B - A.
# Sorties : OUT_DIR/51_rescue_tables.xlsx (Summary, Robust rescued, iPSC, iNeurons,
#           Definitions) + 51_rescue_<dataset>.csv
# Options : RNA_V4_ROOT, DATASETS="iPSC,iNeurons", OUT_DIR, PADJ=0.05, LFC_THR=1,
#           R_COMP=0.7, DELTA_MIN=0.3, COMP_DATASETS="iNeurons" (score mésenchyme seulement là)
# Usage   : Rscript scripts/51_rescue_tables_ipsc_ineurons.R
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(DESeq2); library(AnnotationDbi); library(org.Hs.eg.db)
})
ROOT <- Sys.getenv("RNA_V4_ROOT", unset = file.path(getwd(), "results_v4"))
DATASETS <- trimws(strsplit(Sys.getenv("DATASETS", "iPSC,iNeurons"), ",")[[1]])
OUT_DIR <- Sys.getenv("OUT_DIR", unset = file.path(Sys.getenv("PROJECT_ROOT", "."), "results", "rescue_tables"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
num <- function(v, d) as.numeric(Sys.getenv(v, d))
PADJ <- num("PADJ", "0.05"); LFC_THR <- num("LFC_THR", "1"); R_COMP <- num("R_COMP", "0.7"); DELTA_MIN <- num("DELTA_MIN", "0.3")
COMP_DATASETS <- trimws(strsplit(Sys.getenv("COMP_DATASETS", "iNeurons"), ",")[[1]])   # où la composition varie
SCORE_GENES <- c("COL1A1", "COL1A2", "FN1", "THY1", "ACTA2", "TAGLN")
CT <- c(FA = "control_vs_FA", ED = "control_vs_FA_edited", PAIR = "FA_edited_vs_FA")
STATUS_LEVELS <- c("Rescued", "Partial rescue", "Not rescued", "Reversed")
say <- function(...) cat(sprintf(...), "\n", sep = "")

read_deg <- function(ds, ct) {
  f <- file.path(ROOT, ds, ct, paste0(ct, "_DEG_results.csv"))
  if (!file.exists(f)) stop("Introuvable : ", f, "  -> ajuste RNA_V4_ROOT / DATASETS")
  d <- read_csv(f, show_col_types = FALSE)
  idc <- intersect(c("GENEID", "gene_id", "ensembl_id", "ENSEMBL", "...1", "X"), names(d))[1]
  if (is.na(idc)) stop("Colonne d'identifiant introuvable dans ", f)
  d %>% transmute(GENEID = sub("\\..*$", "", .data[[idc]]), lfc = as.numeric(log2FoldChange), padj = as.numeric(padj)) %>%
    filter(!is.na(lfc)) %>% distinct(GENEID, .keep_all = TRUE)
}

one_dataset <- function(ds) {
  say("\n== %s ==", ds)
  fa <- read_deg(ds, CT[["FA"]]); ed <- read_deg(ds, CT[["ED"]]); pr <- read_deg(ds, CT[["PAIR"]])
  df <- fa %>% inner_join(ed, by = "GENEID", suffix = c("_FA", "_ED")) %>%
    left_join(pr %>% transmute(GENEID, lfc_EDvsFA = -lfc, padj_EDvsFA = padj), by = "GENEID") %>%   # FA_edited_vs_FA = FA - édité
    mutate(sig_FA = !is.na(padj_FA) & padj_FA < PADJ & abs(lfc_FA) >= LFC_THR,
           sig_ED = !is.na(padj_ED) & padj_ED < PADJ & abs(lfc_ED) >= LFC_THR) %>%
    filter(sig_FA) %>%
    mutate(status = case_when(!sig_ED & abs(lfc_ED) < abs(lfc_FA) ~ "Rescued",
                              sig_ED & abs(lfc_ED) < abs(lfc_FA) ~ "Partial rescue",
                              sig_ED & sign(lfc_FA) != sign(lfc_ED) ~ "Reversed",
                              TRUE ~ "Not rescued"))
  # expression normalisée des 9 échantillons
  dl <- lapply(CT, function(ct) readRDS(file.path(ROOT, ds, ct, paste0(ct, "_dds.rds"))))
  common <- Reduce(intersect, lapply(dl, rownames))
  cnt <- do.call(cbind, lapply(unname(dl), function(d) counts(d)[common, , drop = FALSE])); cnt <- cnt[, !duplicated(colnames(cnt)), drop = FALSE]
  cd <- do.call(rbind, lapply(unname(dl), function(d) as.data.frame(colData(d))[, "condition", drop = FALSE]))
  cond <- setNames(as.character(cd$condition), rownames(cd))[colnames(cnt)]
  lex <- log2(t(t(cnt) / estimateSizeFactorsForMatrix(cnt)) + 1)
  rownames(lex) <- sub("\\..*$", "", rownames(lex)); lex <- lex[!duplicated(rownames(lex)), , drop = FALSE]
  sym_all <- suppressMessages(mapIds(org.Hs.eg.db, rownames(lex), "SYMBOL", "ENSEMBL", multiVals = "first"))
  sid <- names(sym_all)[!is.na(sym_all) & sym_all %in% SCORE_GENES]
  score <- if (ds %in% COMP_DATASETS && length(sid) >= 3) colMeans(lex[sid, , drop = FALSE]) else NULL
  # --- appariement FA -> édité (même donneur), tolérant aux formats de noms ---
  say("Échantillons : %s", paste(sprintf("%s[%s]", colnames(lex), cond), collapse = ", "))
  is_ed <- grepl("edit", cond, ignore.case = TRUE)
  is_fa <- !is_ed & grepl("^FA$|FRDA|^FA[^_]*$|patient", cond, ignore.case = TRUE) & !grepl("control|ctrl", cond, ignore.case = TRUE)
  key <- function(x) { x <- toupper(x)
  x <- gsub("[_.-]?(EDITED|EDIT|ED|E|R|RNA|IPSC|INEU(RONS)?|NEU|N)[_.-]?[0-9]*$", "", x)   # suffixes courants
  gsub("[^A-Z0-9]", "", x) }
  fa_s <- colnames(lex)[is_fa]; ed_s <- colnames(lex)[is_ed]
  pairs <- do.call(rbind, lapply(fa_s, function(a) {
    b <- ed_s[key(ed_s) == key(a)]
    if (length(b) != 1) b <- ed_s[startsWith(toupper(ed_s), toupper(sub("[_.-]?R$", "", a)))]   # repli : préfixe
    if (length(b) != 1) {                                   # repli : même numéro de lignée (F86 <-> ED86)
      digits <- function(x) sub("^\\D*(\\d+).*$", "\\1", x)
      b <- ed_s[grepl("\\d", ed_s) & digits(ed_s) == digits(a)]
    }
    if (length(b) == 1) c(donor = key(a), FA = a, ED = b) else NULL }))
  if (is.null(pairs)) {
    say("ATTENTION : aucune paire FA -> édité reconnue pour %s (FA : %s | édités : %s).",
        ds, paste(fa_s, collapse = ", "), paste(ed_s, collapse = ", "))
    say("  -> colonnes « donneurs » laissées vides ; envoie-moi la ligne « Échantillons » ci-dessus pour adapter.")
    pairs <- matrix(character(0), 0, 3, dimnames = list(NULL, c("donor", "FA", "ED")))
  } else say("Paires : %s", paste(sprintf("%s -> %s", pairs[, "FA"], pairs[, "ED"]), collapse = ", "))
  
  df <- df %>% filter(GENEID %in% rownames(lex)) %>%
    mutate(SYMBOL = coalesce(unname(sym_all[GENEID]), GENEID))
  okn <- character(0)
  for (k in seq_len(nrow(pairs))) {
    p <- pairs[k, ]; delta <- lex[df$GENEID, p[["ED"]]] - lex[df$GENEID, p[["FA"]]]
    df[[paste0("delta_", p[["donor"]])]] <- round(unname(delta), 2)
    ok <- paste0("ok_", p[["donor"]]); okn <- c(okn, ok)
    df[[ok]] <- unname(sign(delta) == -sign(df$lfc_FA) & abs(delta) >= DELTA_MIN)
  }
  if (length(okn)) {
    df$pairs_ok <- rowSums(df[, okn, drop = FALSE])
    df$pairs_carrying <- apply(df[, okn, drop = FALSE], 1, function(x) paste(sub("^ok_", "", okn)[x], collapse = ", "))
  } else { df$pairs_ok <- NA_integer_; df$pairs_carrying <- "" }
  df$r_composition <- if (!is.null(score)) round(sapply(df$GENEID, function(id) suppressWarnings(cor(lex[id, ], score))) * sign(df$lfc_FA), 2) else NA
  np <- length(okn)
  df <- df %>% mutate(
    Confidence = case_when(
      !status %in% c("Rescued", "Partial rescue") ~ "",
      !is.na(r_composition) & r_composition >= R_COMP ~ "composition-driven",
      np == 0 ~ "pairs not available",
      pairs_ok == np ~ sprintf("robust (%d/%d donors)", np, np),
      pairs_ok == np - 1 ~ sprintf("partial (%d/%d donors)", np - 1, np),
      pairs_ok == 1 ~ paste0("one donor only (", pairs_carrying, ")"),
      TRUE ~ "not reproduced in pairs"),
    Direction_FA = ifelse(lfc_FA > 0, "up in FA", "down in FA"),
    Model = ds, status = factor(status, levels = STATUS_LEVELS)) %>%
    arrange(status, desc(grepl("^robust", Confidence)), desc(abs(lfc_FA)))
  say("%d FA-DEGs : %s", nrow(df), paste(sprintf("%s %d", levels(df$status), table(df$status)), collapse = " | "))
  df
}

all <- bind_rows(lapply(DATASETS, one_dataset))
out_cols <- function(d) d %>% transmute(Gene = SYMBOL, Model, Status = as.character(status), Confidence, Direction_FA,
                                        `LFC FA vs ctrl` = round(lfc_FA, 2), `padj FA vs ctrl` = signif(padj_FA, 2),
                                        `LFC edited vs ctrl` = round(lfc_ED, 2), `padj edited vs ctrl` = signif(padj_ED, 2),
                                        `LFC edited vs FA (paired)` = round(lfc_EDvsFA, 2), `padj edited vs FA` = signif(padj_EDvsFA, 2),
                                        `Donors with rescue` = pairs_ok, `Donors` = pairs_carrying, `r mesenchymal score` = r_composition,
                                        across(starts_with("delta_")), Ensembl = GENEID)
for (ds in DATASETS) write_csv(out_cols(filter(all, Model == ds)), file.path(OUT_DIR, sprintf("51_rescue_%s.csv", ds)))

summ <- all %>% dplyr::count(Model, status, .drop = FALSE) %>% group_by(Model) %>% mutate(pct = round(100 * n / sum(n))) %>% ungroup() %>%
  mutate(cell = sprintf("%d (%d %%)", n, pct)) %>% dplyr::select(Model, status, cell) %>%
  pivot_wider(names_from = Model, values_from = cell) %>% dplyr::rename(Status = status)
conf <- all %>% filter(status %in% c("Rescued", "Partial rescue")) %>%
  mutate(Confidence = sub(" \\(.*", "", Confidence)) %>% dplyr::count(Model, Confidence) %>%
  pivot_wider(names_from = Model, values_from = n, values_fill = 0)
robust <- all %>% filter(status %in% c("Rescued", "Partial rescue"), grepl("^robust", Confidence)) %>% out_cols()
print(as.data.frame(summ)); print(as.data.frame(conf))
say("\nRescue robuste (tous les donneurs, hors composition) : %s", paste(robust$Gene, collapse = ", "))

defs <- tibble(Term = c("FA-DEG", "Rescued", "Partial rescue", "Not rescued", "Reversed",
                        "robust (3/3 donors)", "partial (2/3 donors)", "one donor only", "composition-driven",
                        "LFC edited vs FA (paired)", "delta_<donor>"),
               Definition = c(sprintf("padj < %g and |log2FC| >= %g in FA vs control", PADJ, LFC_THR),
                              "no longer significant in edited vs control, and closer to control than FA",
                              "still significant in edited vs control, but closer to control than FA",
                              "still significant in edited vs control, not closer to control",
                              "significant in edited vs control, opposite sign to FA",
                              sprintf("moves back towards control in all FA -> edited pairs (|delta log2| >= %g)", DELTA_MIN),
                              "moves back in 2 of the 3 pairs",
                              "only one pair moves back (named in the column)",
                              sprintf("expression follows the mesenchymal score (r >= %g, iNeurons): cell-mix effect", R_COMP),
                              "paired DESeq2 contrast (~ line + condition), positive = higher after editing",
                              "log2 change edited - FA for that donor (normalised counts)"))

xl <- file.path(OUT_DIR, "51_rescue_tables.xlsx")
sheets <- c(list(Summary = summ, `Confidence (rescued)` = conf, `Robust rescued` = robust),
            setNames(lapply(DATASETS, function(ds) out_cols(filter(all, Model == ds))), DATASETS),
            list(Definitions = defs))
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  wb <- createWorkbook(); hs <- createStyle(textDecoration = "bold", fgFill = "#1F4E5F", fontColour = "#FFFFFF", wrapText = TRUE)
  fills <- c(Rescued = "#D9EAD3", `Partial rescue` = "#FFF2CC", `Not rescued` = "#F4CCCC", Reversed = "#E6D0F0")
  for (nm in names(sheets)) {
    addWorksheet(wb, nm); writeData(wb, nm, sheets[[nm]], headerStyle = hs); freezePane(wb, nm, firstRow = TRUE)
    setColWidths(wb, nm, cols = seq_len(ncol(sheets[[nm]])), widths = "auto")
    if ("Status" %in% names(sheets[[nm]]) && nm != "Summary") for (s in names(fills)) {
      r <- which(sheets[[nm]]$Status == s) + 1
      if (length(r)) addStyle(wb, nm, createStyle(fgFill = fills[[s]]), rows = r, cols = seq_len(ncol(sheets[[nm]])), gridExpand = TRUE)
    }
  }
  saveWorkbook(wb, xl, overwrite = TRUE)
} else writexl::write_xlsx(sheets, xl)
say("\nSorties : %s", OUT_DIR)