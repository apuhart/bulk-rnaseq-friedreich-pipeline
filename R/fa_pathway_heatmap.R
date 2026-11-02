# ============================================================
# fa_pathway_heatmap.R
#
# Heatmap du tableau produit par run_fa_pathway_sensitivity()
# (fa_pathway_sensitivity.R) : une ligne par gène du panel FA_PATHWAY_GENES
# (config.R), groupées par catégorie (facettes), une colonne par
# dataset x méthode. Couleur = LFC (control_vs_FA). Cadre noir + "*" =
# significatif après correction FDR restreinte au panel (colonne
# sig_panel du tableau long).
#
# Ne recalcule rien : prend directement le data.frame "long" déjà produit
# par run_fa_pathway_sensitivity() (en mémoire), pas besoin de relire le
# fichier Excel.
# ============================================================

plot_fa_pathway_heatmap <- function(
    df_long,
    out_dir,
    gene_list     = FA_PATHWAY_GENES,
    method_order  = c("DESeq2_baseline", "DESeq2_sex_adjustTRUE",
                       "DESeq2_filtre_souple", "limma_voom"),
    method_label  = c(DESeq2_baseline = "baseline", DESeq2_sex_adjustTRUE = "sex_adjust",
                       DESeq2_filtre_souple = "filtre souple", limma_voom = "limma-voom"),
    dataset_order = c("iPSC", "iNeurons"),
    file_name     = "heatmap_sensitivity.pdf"
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("  [INFO] ggplot2 indisponible : heatmap FA_PATHWAY_GENES ignoré.")
    return(invisible(NULL))
  }
  if (is.null(df_long) || nrow(df_long) == 0) {
    message("  [INFO] Heatmap FA_PATHWAY_GENES : aucune donnée à tracer.")
    return(invisible(NULL))
  }

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Table catégorie / symbole (mêmes catégories que config.R) ----
  gene_df <- do.call(rbind, lapply(names(gene_list), function(cat) {
    v <- gene_list[[cat]]
    data.frame(categorie = cat, symbole = unname(v), stringsAsFactors = FALSE)
  }))
  gene_df <- unique(gene_df)
  cat_order  <- names(gene_list)
  sym_to_cat <- setNames(gene_df$categorie, gene_df$symbole)

  # ---- Colonnes affichées (dataset x méthode, dans l'ordre voulu) ----
  df <- df_long[df_long$Methode %in% method_order & df_long$Dataset %in% dataset_order, ]
  if (nrow(df) == 0) {
    message("  [INFO] Heatmap FA_PATHWAY_GENES : aucune ligne après filtre dataset/méthode.")
    return(invisible(NULL))
  }

  df$col_label <- paste0(df$Dataset, " \u2014 ", method_label[df$Methode])
  col_levels <- unlist(lapply(dataset_order, function(d) paste0(d, " \u2014 ", method_label[method_order])))
  df$col_label <- factor(df$col_label, levels = col_levels)

  # ---- Ordre des gènes (catégorie puis alphabétique), pour l'axe Y ----
  gene_order <- gene_df$symbole[order(match(gene_df$categorie, cat_order), gene_df$symbole)]
  gene_order <- gene_order[gene_order %in% df$Symbole_HGNC]
  df$Symbole_HGNC <- factor(df$Symbole_HGNC, levels = rev(gene_order))
  df$categorie    <- factor(sym_to_cat[as.character(df$Symbole_HGNC)], levels = cat_order)

  # sig_panel peut contenir des NA (gène non testé par cette méthode) : on ne
  # marque comme "significatif" que les TRUE explicites.
  df_sig <- df[!is.na(df$sig_panel) & df$sig_panel == TRUE, ]

  vmax <- suppressWarnings(max(abs(df$LFC), na.rm = TRUE))
  if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  vmax <- max(vmax, 0.5)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_label, y = Symbole_HGNC, fill = LFC)) +
    ggplot2::geom_tile(color = "grey85") +
    ggplot2::geom_tile(
      data = df_sig, mapping = ggplot2::aes(x = col_label, y = Symbole_HGNC),
      fill = NA, color = "black", linewidth = 0.8, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = df_sig, mapping = ggplot2::aes(x = col_label, y = Symbole_HGNC),
      label = "*", color = "black", fontface = "bold", size = 5,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(categorie ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-vmax, vmax), na.value = "grey95",
      name = "log2FC\n(FA vs contr\u00f4le)"
    ) +
    ggplot2::labs(
      title = "FA_PATHWAY_GENES \u2014 LFC (control_vs_FA) par m\u00e9thode",
      subtitle = "Cadre noir + * = significatif (padj restreint au panel < seuil)",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y        = ggplot2::element_text(size = 7),
      panel.grid          = ggplot2::element_blank(),
      plot.title          = ggplot2::element_text(face = "bold", colour = "#1F4E5F"),
      plot.subtitle        = ggplot2::element_text(colour = "grey30"),
      strip.text.y.left     = ggplot2::element_text(angle = 0, size = 7, face = "bold"),
      strip.placement        = "outside",
      strip.background        = ggplot2::element_rect(fill = "#EAF2F0", color = NA),
      panel.spacing            = grid::unit(2, "pt"),
      legend.position            = "right"
    )

  n_genes <- length(unique(df$Symbole_HGNC))
  h <- max(7, n_genes * 0.16 + 1.5)

  out_file <- file.path(out_dir, file_name)
  grDevices::pdf(out_file, width = 9.5, height = h)
  print(p)
  grDevices::dev.off()

  message("  Heatmap FA_PATHWAY_GENES : \u00e9crit -> ", out_file)
  invisible(p)
}
