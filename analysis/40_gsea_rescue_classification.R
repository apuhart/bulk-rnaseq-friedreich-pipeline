#!/usr/bin/env Rscript
# =============================================================================
# 40_gsea_rescue_classification.R  (tous types cellulaires du bulk RNA v4)
# -----------------------------------------------------------------------------
# Pour chaque voie GSEA, combine les 3 contrastes du pipeline bulk v4 pour dire
# si l'édition de FXN RESTAURE, CORRIGE PARTIELLEMENT ou DÉPASSE le niveau
# contrôle -- au lieu du seul « répond à la correction » (éditée vs FA).
#
#   FA vs contrôle       control_vs_FA         NES > 0 = plus haut dans FA
#   Éditée vs FA         FA_edited_vs_FA       fichier : NES > 0 = plus haut dans FA
#                                              -> inversé ici : NES > 0 = plus haut dans l'éditée
#   Éditée vs contrôle   control_vs_FA_edited  NES > 0 = plus haut dans l'éditée
#
# (convention du pipeline : "X_vs_Y" -> log2FC = Y - X ; le script vérifie le
#  sens du 3e contraste par cohérence : (E-C) ~ (F-C) + (E-F) sur toutes les voies)
#
# Classes (seuil PADJ, tendance = padj < PADJ_TREND) :
#   Déficit FA significatif (FA vs contrôle) :
#     restauré               éditée vs FA corrige ET éditée ≈ contrôle (n.s.)
#     partiellement corrigé  éditée vs FA corrige MAIS éditée encore du côté FA (sig.)
#     surcorrigé             éditée vs FA corrige ET dépasse le contrôle (sig., sens opposé)
#     tendance à la correction  éditée vs FA dans le bon sens, padj < PADJ_TREND
#     non corrigé            éditée encore différente du contrôle comme FA, sans correction
#     non corrigé (ambigu)   le reste
#   Pas de déficit FA significatif :
#     répond à l'édition, éditée ≈ contrôle   (correction plausible d'un déficit non détecté)
#     répond à l'édition, éditée ≠ contrôle   (effet propre à la lignée éditée / clonal)
#     non affecté
#
# Entrées : RNA_V4_ROOT/<dataset>/<contraste>/<contraste>_<COLL>_fgsea.csv
# Sorties : OUT_DIR (défaut results/gsea_rescue_classification)
#   40_<dataset>_<COLL>_classification.csv
#   40_classification_toutes.xlsx       (un onglet par type cellulaire, + résumé + voies suivies)
#   40_summary_counts.csv
#   40_focus_pathways.csv / .pdf        trajectoire Contrôle -> FA -> Éditée des voies suivies
#
# Options :
#   RNA_V4_ROOT   défaut ./results_v4 (dossier du projet)
#   DATASETS      liste (virgules) ; défaut : tous les dossiers qui ont les 3 contrastes
#   COLLECTIONS   défaut "H,C5_BP"
#   PADJ          défaut 0.05 ; PADJ_TREND défaut 0.25
#   SIGN_EC       +1 / -1 / auto (défaut auto)
#   OUT_DIR
# Usage :
#   Rscript scripts/40_gsea_rescue_classification.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(ggplot2)
})

ROOT <- Sys.getenv("RNA_V4_ROOT",
  unset = file.path(getwd(), "results_v4"))
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath("."))
OUT_DIR <- Sys.getenv("OUT_DIR", unset = file.path(PROJECT_ROOT, "results", "gsea_rescue_classification"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
COLLECTIONS <- trimws(strsplit(Sys.getenv("COLLECTIONS", unset = "H,C5_BP"), ",")[[1]])
PADJ       <- as.numeric(Sys.getenv("PADJ", unset = "0.05"))
PADJ_TREND <- as.numeric(Sys.getenv("PADJ_TREND", unset = "0.25"))
SIGN_EC    <- Sys.getenv("SIGN_EC", unset = "auto")

CT <- c(FC = "control_vs_FA", EF = "FA_edited_vs_FA", EC = "control_vs_FA_edited")

FOCUS_PATHWAYS <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION$",
  "IRON_SULFUR_CLUSTER_ASSEMBLY$",
  "RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY$",
  "(^|_)MITOCHONDRIAL_TRANSLATION$",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS$",
  "(^|_)STEROL_BIOSYNTHETIC_PROCESS$",
  "(^|_)CHOLESTEROL_BIOSYNTHETIC_PROCESS$")

say <- function(...) cat(sprintf(...), "\n", sep = "")

datasets <- {
  e <- Sys.getenv("DATASETS", unset = "")
  if (nzchar(e)) trimws(strsplit(e, ",")[[1]])
  else {
    d <- list.dirs(ROOT, recursive = FALSE, full.names = FALSE)
    d[vapply(d, function(x) all(dir.exists(file.path(ROOT, x, CT))), logical(1))]
  }
}
if (!length(datasets)) stop("Aucun dossier avec les 3 contrastes dans ", ROOT)
say("Types cellulaires : %s", paste(datasets, collapse = ", "))

read_fg <- function(ds, ct, coll) {
  f <- file.path(ROOT, ds, ct, sprintf("%s_%s_fgsea.csv", ct, coll))
  if (!file.exists(f)) return(NULL)
  read_csv(f, show_col_types = FALSE) %>%
    transmute(pathway, NES = as.numeric(NES), padj = as.numeric(padj), pval = as.numeric(pval),
              size = as.integer(size), leadingEdge = as.character(leadingEdge))
}

classify <- function(nFC, pFC, nEF, pEF, nEC, pEC) {
  # nEF est déjà orienté « éditée - FA »
  sFC <- !is.na(pFC) & pFC < PADJ
  d   <- sign(nFC)
  corr_sig   <- !is.na(pEF) & pEF < PADJ & sign(nEF) == -d
  corr_trend <- !is.na(pEF) & pEF < PADJ_TREND & sign(nEF) == -d
  ec_sig  <- !is.na(pEC) & pEC < PADJ
  resp    <- !is.na(pEF) & pEF < PADJ
  case_when(
    sFC & corr_sig & !ec_sig                     ~ "restauré",
    sFC & corr_sig & ec_sig & sign(nEC) == d     ~ "partiellement corrigé",
    sFC & corr_sig & ec_sig & sign(nEC) == -d    ~ "surcorrigé",
    sFC & corr_trend                             ~ "tendance à la correction",
    sFC & ec_sig & sign(nEC) == d                ~ "non corrigé",
    sFC                                          ~ "non corrigé (ambigu)",
    !sFC & resp & !ec_sig                        ~ "répond à l'édition, éditée ≈ contrôle",
    !sFC & resp & ec_sig                         ~ "répond à l'édition, éditée ≠ contrôle",
    TRUE                                         ~ "non affecté")
}
CLASS_LEVELS <- c("restauré", "partiellement corrigé", "surcorrigé", "tendance à la correction",
                  "non corrigé", "non corrigé (ambigu)",
                  "répond à l'édition, éditée ≈ contrôle", "répond à l'édition, éditée ≠ contrôle",
                  "non affecté")

all_tabs <- list(); counts <- list()
for (ds in datasets) for (coll in COLLECTIONS) {
  fc <- read_fg(ds, CT[["FC"]], coll); ef <- read_fg(ds, CT[["EF"]], coll); ec <- read_fg(ds, CT[["EC"]], coll)
  if (is.null(fc) || is.null(ef) || is.null(ec)) { say("  [%s %s] fichier(s) manquant(s) -- ignoré", ds, coll); next }
  tab <- fc %>% dplyr::select(pathway, size, NES_FA_vs_ctrl = NES, padj_FA_vs_ctrl = padj, LE_FA_vs_ctrl = leadingEdge) %>%
    full_join(ef %>% transmute(pathway, NES_edit_vs_FA = -NES, padj_edit_vs_FA = padj, LE_edit_vs_FA = leadingEdge), by = "pathway") %>%
    full_join(ec %>% transmute(pathway, NES_edit_vs_ctrl = NES, padj_edit_vs_ctrl = padj), by = "pathway")

  # Vérification du sens du contraste éditée vs contrôle
  ok <- with(tab, is.finite(NES_FA_vs_ctrl) & is.finite(NES_edit_vs_FA) & is.finite(NES_edit_vs_ctrl))
  r <- suppressWarnings(cor(tab$NES_edit_vs_ctrl[ok], (tab$NES_FA_vs_ctrl + tab$NES_edit_vs_FA)[ok]))
  sgn <- if (SIGN_EC == "auto") { if (is.finite(r) && r < -0.2) -1 else 1 } else as.numeric(SIGN_EC)
  say("  [%s %s] %d voies ; cohérence (E-C) vs (F-C)+(E-F) : r = %+.2f -> signe éditée vs contrôle %s%s",
      ds, coll, nrow(tab), r, ifelse(sgn > 0, "+1", "-1"),
      if (is.finite(r) && abs(r) < 0.2) "  (ATTENTION : cohérence faible, vérifier SIGN_EC)" else "")
  tab$NES_edit_vs_ctrl <- sgn * tab$NES_edit_vs_ctrl

  tab <- tab %>%
    mutate(classe = factor(classify(NES_FA_vs_ctrl, padj_FA_vs_ctrl, NES_edit_vs_FA, padj_edit_vs_FA,
                                    NES_edit_vs_ctrl, padj_edit_vs_ctrl), levels = CLASS_LEVELS),
           sens_FA = ifelse(NES_FA_vs_ctrl > 0, "↑ dans FA", "↓ dans FA"),
           dataset = ds, collection = coll) %>%
    arrange(classe, padj_FA_vs_ctrl) %>%
    dplyr::select(dataset, collection, pathway, classe, sens_FA, size,
                  NES_FA_vs_ctrl, padj_FA_vs_ctrl, NES_edit_vs_FA, padj_edit_vs_FA,
                  NES_edit_vs_ctrl, padj_edit_vs_ctrl, LE_FA_vs_ctrl, LE_edit_vs_FA)
  write_csv(tab, file.path(OUT_DIR, sprintf("40_%s_%s_classification.csv", ds, coll)))
  all_tabs[[paste(ds, coll)]] <- tab
  counts[[paste(ds, coll)]] <- tab %>% count(dataset, collection, classe, .drop = FALSE)
}
if (!length(all_tabs)) stop("Aucune table construite.")
all_df <- bind_rows(all_tabs)
counts <- bind_rows(counts) %>% pivot_wider(names_from = classe, values_from = n, values_fill = 0)
write_csv(counts, file.path(OUT_DIR, "40_summary_counts.csv"))
say("\nNombre de voies par classe :")
print(as.data.frame(counts))

# Voies suivies
hits <- unique(unlist(lapply(FOCUS_PATHWAYS, function(p) grep(p, unique(all_df$pathway), value = TRUE))))
focus <- all_df %>% filter(pathway %in% hits) %>%
  dplyr::select(dataset, pathway, classe, NES_FA_vs_ctrl, padj_FA_vs_ctrl, NES_edit_vs_FA, padj_edit_vs_FA,
                NES_edit_vs_ctrl, padj_edit_vs_ctrl) %>% arrange(pathway, dataset)
write_csv(focus, file.path(OUT_DIR, "40_focus_pathways.csv"))
say("\nVoies suivies :")
print(as.data.frame(focus %>% mutate(across(where(is.numeric), ~ signif(.x, 2)))), row.names = FALSE)

if (nrow(focus)) {
  short <- function(x) { x <- sub("^(HALLMARK|GOBP)_", "", x); tolower(gsub("_", " ", x)) }
  traj <- focus %>%
    transmute(dataset, voie = short(pathway), classe,
              `Contrôle` = 0, FA = NES_FA_vs_ctrl, `FA éditée` = NES_edit_vs_ctrl,
              sig_FA = padj_FA_vs_ctrl < PADJ, sig_E = padj_edit_vs_ctrl < PADJ) %>%
    pivot_longer(c(`Contrôle`, FA, `FA éditée`), names_to = "groupe", values_to = "NES") %>%
    mutate(groupe = factor(groupe, levels = c("Contrôle", "FA", "FA éditée")),
           sig = case_when(groupe == "Contrôle" ~ NA, groupe == "FA" ~ sig_FA, TRUE ~ sig_E))
  lab <- traj %>% filter(groupe == "FA éditée")
  p <- ggplot(traj, aes(groupe, NES, group = dataset, colour = dataset)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = ifelse(is.na(sig), "réf.", ifelse(sig, "padj < 0,05", "n.s."))), size = 2.6, fill = "white") +
    geom_text(data = lab, aes(label = classe), hjust = -0.08, size = 2.5, show.legend = FALSE) +
    scale_shape_manual(values = c("réf." = 15, "padj < 0,05" = 16, "n.s." = 21), name = NULL) +
    scale_x_discrete(expand = expansion(add = c(0.3, 1.6))) +
    facet_wrap(~ voie, scales = "free_y") +
    labs(title = "Trajectoire des voies suivies : contrôle -> FA -> FA éditée",
         subtitle = "NES de chaque groupe par rapport au contrôle (FA vs contrôle, éditée vs contrôle). Approximation visuelle : les NES de contrastes différents ne s'additionnent pas exactement.",
         x = NULL, y = "NES vs contrôle", colour = NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(OUT_DIR, "40_focus_pathways.pdf"), p, width = 13, height = 8, device = cairo_pdf)
}

# Excel
sheets <- c(list("Résumé" = counts, "Voies suivies" = focus),
            split(all_df, all_df$dataset))
sheets <- lapply(sheets, as.data.frame)
names(sheets) <- substr(names(sheets), 1, 31)
xl <- file.path(OUT_DIR, "40_classification_toutes.xlsx")
if (requireNamespace("writexl", quietly = TRUE)) writexl::write_xlsx(sheets, xl) else
  if (requireNamespace("openxlsx", quietly = TRUE)) openxlsx::write.xlsx(sheets, xl) else
    say("writexl / openxlsx absents : pas d'Excel (les CSV sont là).")
say("\nSorties : %s", OUT_DIR)
