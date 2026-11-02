# ============================================================
# R/extra_analyses.R — Analyses complémentaires
# ============================================================
#  - summarise_deg_counts()  : nombre de DEGs par contraste (LFC shrunk vs non shrunk)
#  - run_rescue_ratio()      : quantifie le rescue par le ratio des LFC (LFC non shrunk)
#  - classify_paired_degs()  : DEGs appariés = effet FA-spécifique ou effet de clone (éditée) ?
#  - run_stat_concordance()  : concordance entre datasets SANS seuil (stat de Wald)
#  - run_nes_concordance()   : concordance des NES fgsea entre datasets
#
# Conventions de signe (cf. run_deseq2_contrast : results(contrast = c("condition", B, A))) :
#   control_vs_FA        -> log2FC = FA - control
#   control_vs_FA_edited -> log2FC = FA_edited - control
#   FA_edited_vs_FA      -> log2FC = FA - FA_edited   (nom inversé par rapport au signe !)
# ============================================================


# ------------------------------------------------------------
# Résumé du nombre de DEGs (LFC shrunk vs non shrunk)
# ------------------------------------------------------------
# apeglm écrase certains effets réels vers 0 : on rapporte les deux versions.
summarise_deg_counts <- function(all_results_full, out_dir, dataset_name,
                                 padj_thr = 0.05, lfc_thr = 1) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  rows <- lapply(names(all_results_full), function(cn) {
    d <- all_results_full[[cn]]
    sig <- !is.na(d$padj) & d$padj < padj_thr
    data.frame(
      dataset = dataset_name, contrast = cn, n_genes_tested = nrow(d),
      n_padj = sum(sig),
      n_padj_lfc_shrunk = sum(sig & abs(d$log2FoldChange) > lfc_thr, na.rm = TRUE),
      n_padj_lfc_mle    = sum(sig & abs(d$lfc_mle) > lfc_thr, na.rm = TRUE),
      n_up_mle   = sum(sig & d$lfc_mle >  lfc_thr, na.rm = TRUE),
      n_down_mle = sum(sig & d$lfc_mle < -lfc_thr, na.rm = TRUE)
    )
  })
  tab <- do.call(rbind, rows)
  utils::write.csv(tab, file.path(out_dir, paste0(dataset_name, "_DEG_summary.csv")), row.names = FALSE)
  print(tab)
  invisible(tab)
}


# ------------------------------------------------------------
# Ratio de rescue
# ------------------------------------------------------------
# Pour chaque gène FA-DEG (control_vs_FA) :
#     r = LFC(FA vs edited) / LFC(FA vs control)
#   r ~ 1  : l'édition ramène le gène vers le niveau contrôle (rescue complet)
#   r ~ 0  : pas d'effet de l'édition
#   r < 0  : l'édition pousse dans le sens inverse du rescue
# Calculé sur le LFC NON shrunk (colonne lfc_mle) : le LFC shrunk par apeglm écrase vers 0
# les effets du contraste apparié et fausse le ratio.
# Biais résiduel : les gènes sont choisis parce que leur effet FA vs contrôle est grand
# (malédiction du gagnant) -> le ratio est tiré vers 0. À lire comme une distribution.
run_rescue_ratio <- function(
    all_results, out_dir, dataset_name,
    padj_thr = 0.05, lfc_thr = 1, lfc_col = "lfc_mle"
) {
  need <- c("control_vs_FA", "FA_edited_vs_FA")
  if (!all(need %in% names(all_results))) {
    message("  [rescue ratio] contrastes manquants pour ", dataset_name, " -> ignoré")
    return(invisible(NULL))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  fa  <- all_results[["control_vs_FA"]]
  ed  <- all_results[["FA_edited_vs_FA"]]
  if (!(lfc_col %in% colnames(fa)) || !(lfc_col %in% colnames(ed))) {
    warning("colonne ", lfc_col, " absente -> LFC shrunk utilisé")
    lfc_col <- "log2FoldChange"
  }
  fa$lfc <- fa[[lfc_col]]; ed$lfc <- ed[[lfc_col]]

  m <- merge(
    fa[, c("GENEID", "GENENAME", "lfc", "padj")],
    ed[, c("GENEID", "lfc", "padj")],
    by = "GENEID", suffixes = c("_FA_vs_ctrl", "_FA_vs_edited")
  )
  degs <- m[!is.na(m$padj_FA_vs_ctrl) & m$padj_FA_vs_ctrl < padj_thr &
              abs(m$lfc_FA_vs_ctrl) > lfc_thr, ]
  if (nrow(degs) == 0) {
    message("  [rescue ratio] aucun FA-DEG pour ", dataset_name)
    return(invisible(NULL))
  }

  degs$rescue_ratio <- degs$lfc_FA_vs_edited / degs$lfc_FA_vs_ctrl
  degs$rescue_class <- as.character(cut(
    degs$rescue_ratio,
    breaks = c(-Inf, 0, 0.25, 0.75, Inf),
    labels = c("Inverse (<0)", "Non rescue (0-0.25)", "Partiel (0.25-0.75)", "Rescue (>0.75)")
  ))
  degs <- degs[order(-degs$rescue_ratio), ]

  utils::write.csv(degs, file.path(out_dir, paste0(dataset_name, "_rescue_ratio.csv")),
                   row.names = FALSE)

  summ <- as.data.frame(table(rescue_class = degs$rescue_class))
  utils::write.csv(summ, file.path(out_dir, paste0(dataset_name, "_rescue_ratio_summary.csv")),
                   row.names = FALSE)
  message("  [rescue ratio] ", dataset_name, " : ", nrow(degs), " FA-DEGs (LFC : ", lfc_col,
          ") | ratio médian = ", round(stats::median(degs$rescue_ratio), 2))
  print(summ)

  with_pdf(file.path(out_dir, paste0(dataset_name, "_rescue_ratio_hist.pdf")), {
    r <- pmax(pmin(degs$rescue_ratio, 2), -1)     # écrêté pour la lisibilité
    graphics::hist(r, breaks = seq(-1, 2, by = 0.1), col = "steelblue", border = "white",
                   xlab = "ratio de rescue = LFC(FA vs edited) / LFC(FA vs control)  [ecrete a -1 / 2]",
                   main = paste0("Rescue des FA-DEGs — ", dataset_name,
                                 " (n = ", nrow(degs), ")"))
    graphics::abline(v = c(0, 1), lty = 2)
  })

  invisible(degs)
}


# ------------------------------------------------------------
# DEGs du contraste apparié : effet FA-spécifique ou effet de clone ?
# ------------------------------------------------------------
# Le contraste apparié FA vs édité mélange deux choses :
#   - un effet dépendant de FXN : FA s'écarte à la fois des contrôles et des éditées,
#     et l'éditée est proche des contrôles                        -> "FA_specifique"
#   - un effet de CLONE (sous-clonage/édition, état épigénétique) : FA est proche des
#     contrôles et c'est l'éditée qui s'en écarte (ex. MEG3, PCDH)  -> "clone_edite"
# Avec d = FA - edited, f_c = FA - control, e_c = edited - control (LFC non shrunk) :
#   FA_specifique : |e_c| < ratio * |d|  et  sign(f_c) == sign(d)
#   clone_edite   : |f_c| < ratio * |d|
#   ambigu        : ni l'un ni l'autre
#   non_testable  : gène absent d'un des contrastes contrôle (filtre, exclusion chrX/chrY/XIST)
# Les deux premières règles sont exclusives (d = f_c - e_c).
classify_paired_degs <- function(
    all_results, out_dir, dataset_name,
    padj_thr = 0.05,
    ratio    = get0("CLASSIFY_RATIO", ifnotfound = 0.5),
    lfc_col  = "lfc_mle"
) {
  need <- c("control_vs_FA", "control_vs_FA_edited", "FA_edited_vs_FA")
  if (!all(need %in% names(all_results))) {
    message("  [classification] contrastes manquants pour ", dataset_name, " -> ignoré")
    return(invisible(NULL))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  pick <- function(df) if (lfc_col %in% colnames(df)) df[[lfc_col]] else df$log2FoldChange
  pr <- all_results[["FA_edited_vs_FA"]]
  fc <- all_results[["control_vs_FA"]]
  ec <- all_results[["control_vs_FA_edited"]]

  sg <- pr[!is.na(pr$padj) & pr$padj < padj_thr, ]
  if (nrow(sg) == 0) {
    message("  [classification] aucun DEG apparié pour ", dataset_name)
    return(invisible(NULL))
  }
  d   <- pick(sg)
  f_c <- pick(fc)[match(sg$GENEID, fc$GENEID)]
  e_c <- pick(ec)[match(sg$GENEID, ec$GENEID)]

  ok  <- !is.na(f_c) & !is.na(e_c)
  cls <- rep("ambigu", nrow(sg))
  cls[!ok] <- "non_testable"
  cls[ok & abs(e_c) < ratio * abs(d) & sign(f_c) == sign(d)] <- "FA_specifique"
  cls[ok & abs(f_c) < ratio * abs(d)] <- "clone_edite"

  out <- data.frame(
    GENEID = sg$GENEID, GENENAME = sg$GENENAME, baseMean = sg$baseMean,
    lfc_FA_minus_edited = d, padj_paired = sg$padj,
    lfc_FA_minus_control = f_c, lfc_edited_minus_control = e_c,
    classe = cls, stringsAsFactors = FALSE
  )
  out <- out[order(out$padj_paired), ]

  utils::write.csv(out, file.path(out_dir, paste0(dataset_name, "_paired_deg_classification.csv")),
                   row.names = FALSE)
  summ <- as.data.frame(table(classe = out$classe))
  utils::write.csv(summ, file.path(out_dir, paste0(dataset_name, "_paired_deg_classification_summary.csv")),
                   row.names = FALSE)
  message("  [classification] ", dataset_name, " : ", nrow(out), " DEGs appariés (padj < ", padj_thr, ")")
  print(summ)

  # Graphique : effet FA (x) contre effet des éditées (y), par rapport aux contrôles
  pl <- out[!is.na(out$lfc_FA_minus_control) & !is.na(out$lfc_edited_minus_control), ]
  if (nrow(pl) > 0) {
    pl$label <- ifelse(is.na(pl$GENENAME) | pl$GENENAME == "", pl$GENEID, pl$GENENAME)
    lab <- pl[order(-abs(pl$lfc_FA_minus_edited)), ][seq_len(min(15, nrow(pl))), ]
    lab <- unique(rbind(lab, pl[pl$label == "FXN", ]))
    g <- ggplot2::ggplot(pl, ggplot2::aes(lfc_FA_minus_control, lfc_edited_minus_control, colour = classe)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey85") +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = label), size = 3,
                               max.overlaps = Inf, show.legend = FALSE) +
      ggplot2::labs(
        title    = paste0("DEGs apparies : FA-specifique ou clone ? — ", dataset_name),
        subtitle = "FA_specifique = sur l'axe y=0 (editee ~ controle) ; clone_edite = sur l'axe x=0 (FA ~ controle)",
        x = "LFC FA - controle (non shrunk)", y = "LFC editee - controle (non shrunk)"
      ) +
      ggplot2::theme_bw(base_size = 12)
    with_pdf(file.path(out_dir, paste0(dataset_name, "_paired_deg_classification.pdf")),
             print(g), width = 9, height = 7)
  }

  invisible(out)
}


# ------------------------------------------------------------
# Concordance sans seuil entre datasets (stat de Wald)
# ------------------------------------------------------------
# Corrélation de Spearman des stats de Wald sur tous les gènes communs.
# Plus robuste que l'intersection de listes de DEGs (3 gènes communs avec n = 3).
# Attention : les mêmes donneurs sont présents dans plusieurs datasets ; une partie
# de la concordance peut refléter le fond génétique des donneurs et pas la maladie.
# Le contraste apparié FA_edited_vs_FA est le seul à s'affranchir de cet effet donneur.
run_stat_concordance <- function(
    all_dataset_results, base_out,
    datasets  = names(all_dataset_results),
    contrasts = c("control_vs_FA", "control_vs_FA_edited", "FA_edited_vs_FA"),
    stat_thr  = 2
) {
  out_dir <- file.path(base_out, "CONCORDANCE")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (length(datasets) < 2) return(invisible(NULL))
  pairs <- utils::combn(datasets, 2, simplify = FALSE)
  rows  <- list()

  for (cn in contrasts) {
    for (pr in pairs) {
      a <- all_dataset_results[[pr[1]]][[cn]]
      b <- all_dataset_results[[pr[2]]][[cn]]
      if (is.null(a) || is.null(b)) next

      m <- merge(a[, c("GENEID", "stat")], b[, c("GENEID", "stat")],
                 by = "GENEID", suffixes = c("_a", "_b"))
      m <- m[is.finite(m$stat_a) & is.finite(m$stat_b), ]
      if (nrow(m) < 50) next

      ct     <- suppressWarnings(stats::cor.test(m$stat_a, m$stat_b, method = "spearman", exact = FALSE))
      either <- m[abs(m$stat_a) > stat_thr | abs(m$stat_b) > stat_thr, ]
      both   <- m[abs(m$stat_a) > stat_thr & abs(m$stat_b) > stat_thr, ]

      rows[[length(rows) + 1]] <- data.frame(
        contrast           = cn,
        dataset_a          = pr[1],
        dataset_b          = pr[2],
        n_genes            = nrow(m),
        spearman_rho       = unname(ct$estimate),
        spearman_p         = ct$p.value,
        n_either_sig       = nrow(either),
        rho_either_sig     = if (nrow(either) >= 10)
          suppressWarnings(stats::cor(either$stat_a, either$stat_b, method = "spearman")) else NA_real_,
        n_both_sig         = nrow(both),
        pct_same_sign_both = if (nrow(both) > 0)
          100 * mean(sign(both$stat_a) == sign(both$stat_b)) else NA_real_
      )

      g <- ggplot2::ggplot(m, ggplot2::aes(stat_a, stat_b)) +
        ggplot2::geom_point(alpha = 0.15, size = 0.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::labs(
          title    = paste0("Concordance des stats de Wald — ", cn),
          subtitle = sprintf("Spearman rho = %.3f (n = %d genes)", unname(ct$estimate), nrow(m)),
          x = paste0("stat de Wald — ", pr[1]),
          y = paste0("stat de Wald — ", pr[2])
        ) +
        ggplot2::theme_bw(base_size = 12)

      with_pdf(file.path(out_dir, paste0("stat_scatter_", cn, "_", pr[1], "_vs_", pr[2], ".pdf")),
               print(g))
    }
  }

  if (length(rows) > 0) {
    tab <- do.call(rbind, rows)
    utils::write.csv(tab, file.path(out_dir, "stat_concordance_summary.csv"), row.names = FALSE)
    message("\n  [concordance] résumé :")
    print(tab[, c("contrast", "dataset_a", "dataset_b", "n_genes", "spearman_rho")])
    invisible(tab)
  } else invisible(NULL)
}


# ------------------------------------------------------------
# Concordance des NES fgsea entre datasets
# ------------------------------------------------------------
# Lit les CSV écrits par run_fgsea_safe : <base_out>/<dataset>/<cn>/<cn>_<coll>_fgsea.csv
run_nes_concordance <- function(
    base_out,
    datasets    = c("fibro", "iPSC", "iNeurons"),
    contrasts   = c("control_vs_FA", "control_vs_FA_edited", "FA_edited_vs_FA"),
    collections = c("H", "C5_BP", "C5_CC", "C5_MF")
) {
  out_dir <- file.path(base_out, "CONCORDANCE")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  pairs <- utils::combn(datasets, 2, simplify = FALSE)
  rows  <- list()

  read_fg <- function(ds, cn, coll) {
    f <- file.path(base_out, ds, cn, paste0(cn, "_", coll, "_fgsea.csv"))
    if (!file.exists(f)) return(NULL)
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    d[, intersect(c("pathway", "NES", "padj"), colnames(d))]
  }

  for (cn in contrasts) for (coll in collections) for (pr in pairs) {
    a <- read_fg(pr[1], cn, coll); b <- read_fg(pr[2], cn, coll)
    if (is.null(a) || is.null(b)) next
    m <- merge(a, b, by = "pathway", suffixes = c("_a", "_b"))
    m <- m[is.finite(m$NES_a) & is.finite(m$NES_b), ]
    if (nrow(m) < 10) next
    ct  <- suppressWarnings(stats::cor.test(m$NES_a, m$NES_b, method = "spearman", exact = FALSE))
    sig <- m[(!is.na(m$padj_a) & m$padj_a < 0.05) | (!is.na(m$padj_b) & m$padj_b < 0.05), ]
    rows[[length(rows) + 1]] <- data.frame(
      contrast = cn, collection = coll, dataset_a = pr[1], dataset_b = pr[2],
      n_pathways = nrow(m), spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
      n_sig_either = nrow(sig),
      pct_same_sign_sig = if (nrow(sig) > 0)
        100 * mean(sign(sig$NES_a) == sign(sig$NES_b)) else NA_real_
    )
  }

  if (length(rows) > 0) {
    tab <- do.call(rbind, rows)
    utils::write.csv(tab, file.path(out_dir, "NES_concordance_summary.csv"), row.names = FALSE)
    message("  [concordance NES] ", nrow(tab), " comparaisons -> NES_concordance_summary.csv")
    invisible(tab)
  } else invisible(NULL)
}
