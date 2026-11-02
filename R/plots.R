# ============================================================
# R/functions/plots.R — All plotting helpers
# ============================================================

# ----------------------------
# Volcano plot
# ----------------------------

#' Standard volcano plot (all significant genes labelled)
plot_volcano <- function(
    res_df, dataset_name, contrast_name, direction_label,
    lfc_thr = 0.3, padj_thr = 0.05
) {
  df <- res_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(sig = padj < padj_thr & abs(log2FoldChange) > lfc_thr)
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    coord_cartesian(
      xlim = c(min(df$log2FoldChange) - 3, max(df$log2FoldChange) + 3)
    ) +
    geom_text_repel(
      data = subset(df, padj < 0.1 & abs(log2FoldChange) > lfc_thr),
      aes(label = gene), size = 5, max.overlaps = 50
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x     = "log2 Fold Change",
      y     = "-log10 adjusted p-value",
      color = "Significant",
      title = paste("Volcano -", dataset_name, contrast_name, "\n", direction_label)
    )
}


#' Focused volcano plot — labels only genes of interest
plot_volcano_focus <- function(
    res_df, dataset_name, contrast_name, direction_label,
    genes_of_interest,
    lfc_thr = 0.3, padj_thr = 0.05
) {
  df <- res_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      sig         = padj < padj_thr & abs(log2FoldChange) > lfc_thr,
      is_interest = gene %in% genes_of_interest
    )
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    geom_text_repel(
      data     = subset(df, is_interest),
      aes(label = gene),
      fontface = "bold", size = 4, max.overlaps = 50
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x     = "log2 Fold Change",
      y     = "-log10 adjusted p-value",
      color = "Significant",
      title = paste("Volcano (genes of interest) -", dataset_name, contrast_name, "\n", direction_label)
    )
}


# ----------------------------
# PCA plots
# ----------------------------

#' 2D PCA with condition colour + sex shape
plot_pca_2d <- function(vsd, dataset_name, contrast_name) {
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pca_data$sex <- colData(vsd)$sex[match(pca_data$name, rownames(colData(vsd)))]
  percentVar   <- round(100 * attr(pca_data, "percentVar"))
  
  ggplot(pca_data, aes(PC1, PC2)) +
    geom_point(aes(color = condition, shape = sex), size = 4) +
    geom_text_repel(aes(label = name), size = 4) +
    scale_color_manual(values = c(
      "control"   = "steelblue",
      "FA"        = "firebrick",
      "FA_edited" = "darkorange"
    )) +
    scale_shape_manual(values = c("M" = 17, "F" = 16)) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    ggtitle(paste("PCA -", dataset_name, contrast_name)) +
    theme_minimal(base_size = 14)
}


#' 3D interactive PCA (plotly) — returns a plotly object
plot_pca_3d <- function(vsd, dataset_name, contrast_name) {
  mat      <- t(assay(vsd))
  pca      <- prcomp(mat, scale. = TRUE)
  pca_data <- as.data.frame(pca$x[, 1:3])
  pca_data$condition <- colData(vsd)$condition
  pca_data$sex       <- colData(vsd)$sex
  pca_data$name      <- rownames(pca_data)
  pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  
  color_map <- c("control" = "steelblue", "FA" = "firebrick", "FA_edited" = "darkorange")
  shape_map <- c("M" = "triangle-up", "F" = "circle")
  
  p <- plot_ly(
    data         = pca_data,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color        = ~condition, colors  = color_map,
    symbol       = ~sex,       symbols = shape_map,
    text         = ~name,
    type         = "scatter3d",
    mode         = "markers+text",
    textposition = "top center",
    marker       = list(size = 6)
  ) %>% plotly::layout(
    scene = list(
      xaxis = list(title = paste0("PC1: ", pct[1], "%")),
      yaxis = list(title = paste0("PC2: ", pct[2], "%")),
      zaxis = list(title = paste0("PC3: ", pct[3], "%"))
    ),
    title = paste("3D PCA -", dataset_name, contrast_name)
  )
  p
}


# ----------------------------
# Heatmaps
# ----------------------------

#' Heatmap of top N DEGs (by padj)
plot_heatmap_top_degs <- function(
    vsd_mat, res_df, colData, sample_order,
    dataset_name, contrast_name,
    top_n = 100, padj_thr = 0.05
) {
  top_df <- res_df %>%
    dplyr::filter(padj < padj_thr) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = top_n)
  
  top_genes <- top_df$GENEID
  
  # Keep only genes present in vsd_mat
  keep <- top_genes %in% rownames(vsd_mat)
  top_df    <- top_df[keep, , drop = FALSE]
  top_genes <- top_genes[keep]
  
  if (length(top_genes) < 2) {
    message("Moins de 2 DEGs à afficher : heatmap ignorée pour ", contrast_name)
    return(invisible(NULL))
  }
  
  mat <- vsd_mat[top_genes, sample_order, drop = FALSE]
  
  # ENSG -> symbole (GENENAME) pour l'affichage des lignes ; repli sur l'ENSG
  # si aucun symbole n'est connu, et dédoublonnage si deux gènes partagent le
  # même symbole (nécessaire : pheatmap exige des noms de lignes uniques).
  symb <- top_df$GENENAME[match(top_genes, top_df$GENEID)]
  no_symb <- is.na(symb) | symb == ""
  symb[no_symb] <- top_genes[no_symb]
  rownames(mat) <- make.unique(symb)
  
  pheatmap(
    mat,
    scale          = "row",
    annotation_col = as.data.frame(colData[sample_order, ]),
    main           = paste("Top", top_n, "DEGs -", dataset_name, contrast_name),
    fontsize_row   = 16,
    fontsize_col   = 14,
    fontsize       = 16,
    cluster_cols   = FALSE
  )
}




plot_heatmap_goi <- function(
    vsd_mat, colData, sample_order,
    genes_of_interest,
    gene_map,
    dataset_name, contrast_name
) {
  
  # convert SYMBOL -> ENSG uniquement pour les entrées qui ne sont pas déjà
  # des ENSG (l'appelant peut fournir soit des symboles, soit des ENSG déjà
  # convertis - ex. goi_ensg dans mainPatchSave.R). Auparavant, la fonction
  # reconvertissait systématiquement dès qu'un seul gène manquait dans
  # vsd_mat, y compris quand genes_of_interest contenait déjà des ENSG :
  # gene_map$GENENAME %in% <des ENSG> ne matchait jamais rien, ce qui vidait
  # entièrement la liste et produisait un PDF de 0 page.
  is_ensg <- grepl("^ENSG", genes_of_interest)
  if (!all(is_ensg)) {
    converted <- gene_map$GENEID[match(genes_of_interest[!is_ensg], gene_map$GENENAME)]
    genes_of_interest <- c(genes_of_interest[is_ensg], converted)
    genes_of_interest <- genes_of_interest[!is.na(genes_of_interest)]
  }
  
  goi_present <- genes_of_interest[
    genes_of_interest %in% rownames(vsd_mat)
  ]
  
  if (length(goi_present) == 0) {
    message("No genes of interest found in vsd for ", contrast_name)
    return(invisible(NULL))
  }
  
  mat <- vsd_mat[goi_present, sample_order, drop = FALSE]
  
  # ENSG -> symbole (GENENAME) pour l'affichage des lignes ; repli sur l'ENSG
  # si aucun symbole n'est connu, et dédoublonnage si nécessaire.
  symb <- gene_map$GENENAME[match(goi_present, gene_map$GENEID)]
  no_symb <- is.na(symb) | symb == ""
  symb[no_symb] <- goi_present[no_symb]
  rownames(mat) <- make.unique(symb)
  
  pheatmap(
    mat,
    scale = "row",
    annotation_col = as.data.frame(colData[sample_order, ]),
    main = paste("Genes of Interest -", dataset_name, contrast_name),
    fontsize_row = 16,
    fontsize_col = 14,
    cluster_cols = FALSE,
    cluster_rows = FALSE
  )
}
#' 
#' #' Heatmap of genes of interest
#' plot_heatmap_goi <- function(
    #'     vsd_mat, colData, sample_order,
#'     genes_of_interest,
#'     dataset_name, contrast_name
#' ) {
#'   
#'   browser()
#'   goi_present <- genes_of_interest[genes_of_interest %in% rownames(vsd_mat)]
#' 
#'   if (length(goi_present) == 0) {
#'     message("No genes of interest found in vsd for ", contrast_name)
#'     return(invisible(NULL))
#'   }
#' 
#'   pheatmap(
#'     vsd_mat[goi_present, sample_order, drop = FALSE],
#'     scale          = "row",
#'     annotation_col = as.data.frame(colData[sample_order, ]),
#'     main           = paste("Genes of Interest -", dataset_name, contrast_name),
#'     fontsize_row   = 16,
#'     fontsize_col   = 14,
#'     fontsize       = 16,
#'     cluster_cols   = FALSE,
#'     cluster_rows   = FALSE
#'   )
#' }


# ----------------------------
# Scatter (cross-contrast)
# ----------------------------

#' Build and save scatter plot comparing two contrasts
#'
#' @param comp1,comp2   Names of contrasts in all_results
#' @param all_results   Named list of DESeq2 result data.frames
#' @param out_dir       Root output directory for this scatter comparison
#' @param dataset_name  Dataset label for plot titles
#' 
#' 
#' 
#' 
# run_scatter <- function(comp1, comp2, all_results, out_dir, dataset_name) {
# 
#   if (!all(c(comp1, comp2) %in% names(all_results))) return(invisible(NULL))
# 
#   # ---- directory layout ----
#   comp_out      <- file.path(out_dir, paste0(comp1, "_VS_", comp2))
#   full_dir      <- file.path(comp_out, "FULL")
#   sig_dir       <- file.path(comp_out, "SIGNIFICANT")
#   quad_full_dir <- file.path(full_dir,  "QUAD")
#   quad_sig_dir  <- file.path(sig_dir,   "QUAD")
#   for (d in c(quad_full_dir, quad_sig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
# 
#   diag_slope <- ifelse(comp2 == "FA_vs_FA_edited", -1, 1)
# 
#   fc1 <- paste0("log2FC_", comp1)
#   fc2 <- paste0("log2FC_", comp2)
#   p1  <- paste0("padj_",   comp1)
#   p2  <- paste0("padj_",   comp2)
# 
#   # ---- merge ----
#   merged <- dplyr::full_join(
#     all_results[[comp1]] %>%
#       dplyr::select(GENEID, gene_1 = gene, log2FoldChange, padj) %>%
#       dplyr::rename(!!fc1 := log2FoldChange, !!p1 := padj),
#     all_results[[comp2]] %>%
#       dplyr::select(GENEID, gene_2 = gene, log2FoldChange, padj) %>%
#       dplyr::rename(!!fc2 := log2FoldChange, !!p2 := padj),
#     by = "GENEID"
#   ) %>%
#     dplyr::mutate(
#       gene        = dplyr::coalesce(gene_1, gene_2),
#       gene        = ifelse(is.na(gene) | gene == "", GENEID, gene),
#       dist_diag   = abs(.data[[fc2]] - diag_slope * .data[[fc1]]) / sqrt(1 + diag_slope^2),
#       significant = (.data[[p1]] < 0.05 | .data[[p2]] < 0.05)
#     )
# 
#   write.csv(
#     merged,
#     file.path(comp_out, paste0(dataset_name, "_", comp1, "_VS_", comp2, "_scatter_table.csv")),
#     row.names = FALSE
#   )
# 
#   # ---- inner helpers ----
#   make_plot <- function(df, filename, save_dir, label_type = "none") {
#     p <- ggplot(df, aes(x = .data[[fc1]], y = .data[[fc2]])) +
#       geom_point(
#         aes(fill = dist_diag, color = significant),
#         shape = 21, size = 2.3, alpha = 0.85, stroke = 0.4
#       ) +
#       scale_fill_viridis_c(name = "Distance to diagonal") +
#       scale_color_manual(values = c("TRUE" = "red", "FALSE" = "transparent"), guide = "none") +
#       geom_abline(slope = diag_slope, intercept = 0, linetype = "dashed", color = "grey60") +
#       geom_hline(yintercept = 0, linetype = "dotted") +
#       geom_vline(xintercept = 0, linetype = "dotted") +
#       theme_minimal(base_size = 14) +
#       labs(
#         x     = paste("Log2FC", comp1),
#         y     = paste("Log2FC", comp2),
#         title = paste(dataset_name, "-", comp1, "vs", comp2)
#       )
# 
#     if (label_type != "none") {
#       label_col <- ifelse(label_type == "ENSG", "GENEID", "gene")
#       label_df  <- df %>%
#         dplyr::mutate(score = abs(.data[[fc1]]) + abs(.data[[fc2]])) %>%
#         dplyr::arrange(dplyr::desc(score)) %>%
#         dplyr::slice_head(n = 50)
#       p <- p + geom_text_repel(data = label_df, aes(label = .data[[label_col]]),
#                                size = 3, max.overlaps = 50)
#     }
#     ggplot2::ggsave(file.path(save_dir, filename), p, width = 10, height = 8, bg = "white")
#   }
# 
#   plot_set <- function(df, suffix, save_dir) {
#     base <- paste0(dataset_name, "_", comp1, "_VS_", comp2, "_", suffix)
#     make_plot(df, paste0(base, ".png"),       save_dir, "none")
#     make_plot(df, paste0(base, "_ENSG.png"),  save_dir, "ENSG")
#     make_plot(df, paste0(base, "_genes.png"), save_dir, "gene")
#   }
# 
#   get_quadrants <- function(df) {
#     if (diag_slope == 1) {
#       list(
#         hd = df %>% dplyr::filter(.data[[fc1]] > 0 & .data[[fc2]] > 0),
#         bg = df %>% dplyr::filter(.data[[fc1]] < 0 & .data[[fc2]] < 0)
#       )
#     } else {
#       list(
#         hd = df %>% dplyr::filter(.data[[fc1]] < 0 & .data[[fc2]] > 0),
#         bg = df %>% dplyr::filter(.data[[fc1]] > 0 & .data[[fc2]] < 0)
#       )
#     }
#   }
# 
#   # ---- run plots ----
#   plot_set(merged, "FULL", full_dir)
#   quads <- get_quadrants(merged)
#   plot_set(quads$hd, "quad_HD", quad_full_dir)
#   plot_set(quads$bg, "quad_BG", quad_full_dir)
# 
#   sig_df    <- merged %>% dplyr::filter(.data[[p1]] < 0.05 | .data[[p2]] < 0.05)
#   plot_set(sig_df, "SIGNIFICANT", sig_dir)
#   quads_sig <- get_quadrants(sig_df)
#   plot_set(quads_sig$hd, "quad_HD", quad_sig_dir)
#   plot_set(quads_sig$bg, "quad_BG", quad_sig_dir)
# 
#   invisible(merged)
# }


# ============================================================
# run_scatter() — scatter log2FC entre 2 contrastes, pour un dataset donné
# Version enrichie : courbe de tendance (lm), stats (r, pente, R², p-value),
# labels ggrepel, export xlsx des gènes utilisés.
# À placer dans R/cross_dataset_comparison.R
# ============================================================

run_scatter <- function(comp1, comp2, all_results, out_dir, dataset_name,
                        padj_thr = 0.05, lfc_thr = 1,
                        top_n_labels = 100, export_xlsx = TRUE) {
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  df1 <- as.data.frame(all_results[[comp1]])
  df2 <- as.data.frame(all_results[[comp2]])
  
  # Récupère l'ID du gène : colonne "gene"/"GENEID", ou rownames
  get_gene_col <- function(df) {
    if ("GENEID" %in% colnames(df)) return(df$GENEID)
    if ("gene" %in% colnames(df))   return(df$gene)
    rownames(df)
  }
  df1$GENEID <- get_gene_col(df1)
  df2$GENEID <- get_gene_col(df2)
  
  # Récupère le symbole si disponible, sinon retombe sur l'ID
  get_symbol_col <- function(df) {
    if ("GENENAME" %in% colnames(df)) return(df$GENENAME)
    if ("gene" %in% colnames(df))     return(df$gene)
    df$GENEID
  }
  df1$GENENAME <- get_symbol_col(df1)
  
  df_a <- df1 %>%
    dplyr::select(GENEID, GENENAME, log2FoldChange, padj) %>%
    dplyr::rename(lfc_a = log2FoldChange, padj_a = padj)
  
  df_b <- df2 %>%
    dplyr::select(GENEID, log2FoldChange, padj) %>%
    dplyr::rename(lfc_b = log2FoldChange, padj_b = padj)
  
  df_merge <- dplyr::inner_join(df_a, df_b, by = "GENEID") %>%
    dplyr::filter(!is.na(lfc_a), !is.na(lfc_b)) %>%
    dplyr::mutate(
      sig_a = !is.na(padj_a) & padj_a < padj_thr & abs(lfc_a) >= lfc_thr,
      sig_b = !is.na(padj_b) & padj_b < padj_thr & abs(lfc_b) >= lfc_thr,
      sig_status = dplyr::case_when(
        sig_a & sig_b  ~ "Both",
        sig_a & !sig_b ~ paste0("Only ", comp1),
        sig_b & !sig_a ~ paste0("Only ", comp2),
        TRUE           ~ "NS"
      )
    )
  
  status_levels <- c("NS", paste0("Only ", comp1), paste0("Only ", comp2), "Both")
  df_merge$sig_status <- factor(df_merge$sig_status, levels = status_levels)
  
  color_values <- setNames(
    c("grey70", "#2166ac", "#d6604d", "#7b2d8b"),
    status_levels
  )
  
  if (nrow(df_merge) < 10) {
    message("  [SKIP] ", dataset_name, " : ", comp1, " vs ", comp2,
            " -- trop peu de gènes communs (n=", nrow(df_merge), ")")
    return(invisible(NULL))
  }
  
  # ---- Corrélation ----
  r_val  <- cor(df_merge$lfc_a, df_merge$lfc_b, use = "complete.obs")
  r_test <- cor.test(df_merge$lfc_a, df_merge$lfc_b)
  
  # ---- Régression linéaire : lfc_b ~ lfc_a ----
  fit       <- lm(lfc_b ~ lfc_a, data = df_merge)
  slope     <- coef(fit)[["lfc_a"]]
  r_squared <- summary(fit)$r.squared
  slope_p   <- summary(fit)$coefficients["lfc_a", "Pr(>|t|)"]
  pct_diff  <- (slope - 1) * 100
  
  direction_txt <- if (slope > 1) {
    sprintf("%s montre un effet ~%.0f%% plus fort que %s", comp2, pct_diff, comp1)
  } else {
    sprintf("%s montre un effet ~%.0f%% plus faible que %s", comp2, abs(pct_diff), comp1)
  }
  
  message(sprintf(
    "  %s | %s vs %s : n=%d | r=%.3f (p=%.3g) | pente=%.3f (R2=%.3f, p=%.3g)",
    dataset_name, comp1, comp2, nrow(df_merge), r_val, r_test$p.value,
    slope, r_squared, slope_p
  ))
  
  # ---- Labels : gènes sig. dans les deux + top N les plus éloignés de la régression ----
  df_merge <- df_merge %>%
    dplyr::mutate(
      predicted_b = predict(fit, newdata = df_merge),
      abs_resid   = abs(lfc_b - predicted_b)
    )
  
  
  
  # 
  # 
  # label_both <- df_merge %>% dplyr::filter(sig_status == "Both")
  # label_only <- df_merge %>%
  #   dplyr::filter(sig_status != "Both", sig_status != "NS") %>%
  #   dplyr::arrange(dplyr::desc(abs_resid)) %>%
  #   dplyr::slice_head(n = top_n_labels)
  # 
  # label_df <- dplyr::bind_rows(label_both, label_only) %>%
  #   dplyr::distinct(GENEID, .keep_all = TRUE)
  # 
  
  
  # ---- Labels ----
  
  df_merge <- df_merge %>%
    dplyr::mutate(
      predicted_b = predict(fit, newdata = df_merge),
      abs_resid   = abs(lfc_b - predicted_b)
    )
  
  # Gènes significatifs dans les deux
  label_both <- df_merge %>%
    dplyr::filter(sig_status == "Both")
  
  # Top N gènes significatifs dans un seul contraste
  label_only <- df_merge %>%
    dplyr::filter(
      sig_status != "Both",
      sig_status != "NS"
    ) %>%
    dplyr::arrange(dplyr::desc(abs_resid)) %>%
    dplyr::slice_head(n = top_n_labels)
  
  # Labels finaux
  label_df <- dplyr::bind_rows(
    label_both,
    label_only
  ) %>%
    dplyr::distinct(GENEID, .keep_all = TRUE)
  
  
  
  
  # ---- Plot ----
  p <- ggplot(df_merge, aes(x = lfc_a, y = lfc_b, color = sig_status)) +
    geom_point(alpha = 0.5, size = 1.1) +
    geom_smooth(
      aes(x = lfc_a, y = lfc_b), inherit.aes = FALSE, data = df_merge,
      method = "lm", se = TRUE, color = "black", linewidth = 0.8
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey30") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggrepel::geom_label_repel(
      data = label_df, aes(label = GENENAME),
      size = 2.6, max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_color_manual(values = color_values, name = NULL) +
    annotate(
      "label",
      x = min(df_merge$lfc_a, na.rm = TRUE),
      y = max(df_merge$lfc_b, na.rm = TRUE),
      hjust = 0, vjust = 1,
      label = sprintf(
        "Pearson r = %.2f (p = %.2g)\nPente = %.2f (R2 = %.2f, p = %.2g)\nnGenes = %d\n%s",
        r_val, r_test$p.value, slope, r_squared, slope_p, nrow(df_merge), direction_txt
      ),
      size = 3, fill = "white", alpha = 0.85
    ) +
    theme_bw(base_size = 11) +
    labs(
      title    = paste0(dataset_name, " : ", comp1, " vs ", comp2),
      subtitle = "Ligne pointillée = y=x | Ligne pleine = régression observée",
      x = paste0("log2FC (", comp1, ")"),
      y = paste0("log2FC (", comp2, ")")
    ) +
    theme(legend.position = "top")
  
  out_file <- file.path(out_dir, paste0(dataset_name, "_", comp1, "_vs_", comp2, "_scatter.pdf"))
  pdf(out_file, width = 10, height = 9)
  print(p)
  dev.off()
  message("  Scatter sauvegardé : ", out_file)
  
  # ---- Export xlsx des gènes utilisés ----
  if (export_xlsx) {
    xlsx_file <- file.path(out_dir, paste0("genes_", dataset_name, "_", comp1, "_vs_", comp2, ".xlsx"))
    writexl::write_xlsx(df_merge, xlsx_file)
    message("  Table gènes exportée : ", xlsx_file)
  }
  
  return(invisible(df_merge))
}



run_quadrant_plots <- function(
    data_df,
    out_dir,
    dataset_name,
    top_n = 50
) {
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  # ============================================================
  # Vérification
  # ============================================================
  
  if (!is.data.frame(data_df)) {
    
    stop(
      "run_quadrant_plots(): data_df n'est pas un data.frame. ",
      "Classe reçue : ",
      paste(class(data_df), collapse = ", ")
    )
    
  }
  
  required_cols <- c(
    "GENEID",
    "GENENAME",
    "lfc_FA",
    "lfc_EDITED",
    "rescue_status"
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(data_df)
  )
  
  if (length(missing_cols) > 0) {
    
    stop(
      "Colonnes manquantes dans data_df : ",
      paste(missing_cols, collapse = ", ")
    )
    
  }
  
  # ============================================================
  # Définition des quadrants
  # ============================================================
  
  data_df <- data_df %>%
    dplyr::mutate(
      
      quadrant = dplyr::case_when(
        
        lfc_FA >= 0 &
          lfc_EDITED >= 0 ~
          "Q1_FApos_EDITEDpos",
        
        lfc_FA < 0 &
          lfc_EDITED >= 0 ~
          "Q2_FAneg_EDITEDpos",
        
        lfc_FA < 0 &
          lfc_EDITED < 0 ~
          "Q3_FAneg_EDITEDneg",
        
        lfc_FA >= 0 &
          lfc_EDITED < 0 ~
          "Q4_FApos_EDITEDneg",
        
        TRUE ~ NA_character_
      )
      
    )
  
  # ============================================================
  # Noms des quadrants
  # ============================================================
  
  quadrant_info <- c(
    
    Q1_FApos_EDITEDpos =
      "Q1 : FA+ / FA_edited+",
    
    Q2_FAneg_EDITEDpos =
      "Q2 : FA- / FA_edited+",
    
    Q3_FAneg_EDITEDneg =
      "Q3 : FA- / FA_edited-",
    
    Q4_FApos_EDITEDneg =
      "Q4 : FA+ / FA_edited-"
    
  )
  
  # ============================================================
  # Couleurs
  # ============================================================
  
  color_values <- c(
    
    "Not FA-DEG"     = "grey75",
    "Rescued"        = "#2166ac",
    "Partial rescue" = "#67a9cf",
    "Not rescued"    = "#d6604d",
    "Reversed"       = "#7b2d8b",
    "Other"          = "grey50"
    
  )
  
  # ============================================================
  # Un plot par quadrant
  # ============================================================
  
  for (q in names(quadrant_info)) {
    
    quadrant_df <- data_df %>%
      
      dplyr::filter(
        quadrant == q
      )
    
    if (nrow(quadrant_df) == 0) {
      
      message(
        "  [SKIP] ",
        q,
        " : aucun gène"
      )
      
      next
      
    }
    
    # ----------------------------------------------------------
    # Labels
    #
    # On garde uniquement les gènes FA-DEG
    # ----------------------------------------------------------
    
    label_df <- quadrant_df %>%
      
      dplyr::filter(
        rescue_status != "Not FA-DEG"
      ) %>%
      
      dplyr::arrange(
        dplyr::desc(abs(lfc_FA))
      ) %>%
      
      dplyr::slice_head(
        n = top_n
      )
    
    # ----------------------------------------------------------
    # Plot
    # ----------------------------------------------------------
    
    p <- ggplot(
      
      quadrant_df,
      
      aes(
        x = lfc_FA,
        y = lfc_EDITED,
        color = rescue_status
      )
      
    ) +
      
      geom_point(
        alpha = 0.55,
        size = 1.3
      ) +
      
      geom_label_repel(
        
        data = label_df,
        
        aes(
          label = GENENAME
        ),
        
        size = 2.5,
        
        max.overlaps = Inf,
        
        box.padding = 0.35,
        
        point.padding = 0.15,
        
        force = 2,
        
        min.segment.length = 0,
        
        show.legend = FALSE
        
      ) +
      
      geom_abline(
        
        slope = 1,
        
        intercept = 0,
        
        linetype = "dotted",
        
        color = "grey30"
        
      ) +
      
      geom_hline(
        
        yintercept = 0,
        
        linetype = "dashed",
        
        color = "grey60"
        
      ) +
      
      geom_vline(
        
        xintercept = 0,
        
        linetype = "dashed",
        
        color = "grey60"
        
      ) +
      
      scale_color_manual(
        
        values = color_values,
        
        name = NULL,
        
        drop = FALSE
        
      ) +
      
      theme_bw(
        base_size = 11
      ) +
      
      labs(
        
        title = paste0(
          dataset_name,
          " : ",
          quadrant_info[[q]]
        ),
        
        subtitle = paste0(
          "n = ",
          nrow(quadrant_df),
          " genes | ",
          min(top_n, nrow(label_df)),
          " labels"
        ),
        
        x = "log2FC (FA vs Control)",
        
        y = "log2FC (FA_edited vs Control)"
        
      ) +
      
      theme(
        legend.position = "top"
      )
    
    # ----------------------------------------------------------
    # Export
    # ----------------------------------------------------------
    
    out_file <- file.path(
      
      out_dir,
      
      paste0(
        dataset_name,
        "_",
        q,
        "_scatter.pdf"
      )
      
    )
    
    pdf(
      out_file,
      width = 10,
      height = 9
    )
    
    print(p)
    
    dev.off()
    
    message(
      "  Quadrant scatter sauvegardé : ",
      out_file,
      " | n = ",
      nrow(quadrant_df),
      " | labels = ",
      nrow(label_df)
    )
    
  }
  
  invisible(data_df)
  
}


# 
# run_rescue_scatter <- function(
    #     all_results,
#     out_dir,
#     dataset_name,
#     padj_thr = 0.05,
#     lfc_thr = 1,
#     top_n_labels = 30,
#     export_xlsx = TRUE
# ) {
#   
#   library(ggplot2)
#   library(ggrepel)
#   library(dplyr)
#   
#   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#   
#   # ============================================================
#   # Vérification des contrastes
#   # ============================================================
#   
#   required <- c("control_vs_FA", "control_vs_FA_edited")
#   
#   if (!all(required %in% names(all_results))) {
#     message(
#       "  [SKIP] ", dataset_name,
#       " : contrastes nécessaires absents."
#     )
#     return(invisible(NULL))
#   }
#   
#   df1 <- as.data.frame(all_results[["control_vs_FA"]])
#   df2 <- as.data.frame(all_results[["control_vs_FA_edited"]])
#   
#   # ============================================================
#   # Récupération gene ID / symbol
#   # ============================================================
#   
#   get_gene_col <- function(df) {
#     
#     if ("GENEID" %in% colnames(df)) return(df$GENEID)
#     if ("gene" %in% colnames(df))   return(df$gene)
#     
#     rownames(df)
#   }
#   
#   get_symbol_col <- function(df) {
#     
#     if ("GENENAME" %in% colnames(df)) return(df$GENENAME)
#     if ("gene" %in% colnames(df))     return(df$gene)
#     
#     rownames(df)
#   }
#   
#   df1$GENEID   <- get_gene_col(df1)
#   df1$GENENAME <- get_symbol_col(df1)
#   
#   df2$GENEID <- get_gene_col(df2)
#   
#   # ============================================================
#   # Merge
#   # ============================================================
#   
#   df_a <- df1 %>%
#     dplyr::select(
#       GENEID,
#       GENENAME,
#       log2FoldChange,
#       padj
#     ) %>%
#     dplyr::rename(
#       lfc_FA  = log2FoldChange,
#       padj_FA = padj
#     )
#   
#   df_b <- df2 %>%
#     dplyr::select(
#       GENEID,
#       log2FoldChange,
#       padj
#     ) %>%
#     dplyr::rename(
#       lfc_EDITED  = log2FoldChange,
#       padj_EDITED = padj
#     )
#   
#   df <- dplyr::inner_join(
#     df_a,
#     df_b,
#     by = "GENEID"
#   ) %>%
#     dplyr::filter(
#       !is.na(lfc_FA),
#       !is.na(lfc_EDITED)
#     )
#   
#   # ============================================================
#   # Significance
#   # ============================================================
#   
#   df <- df %>%
#     dplyr::mutate(
#       
#       sig_FA =
#         !is.na(padj_FA) &
#         padj_FA < padj_thr &
#         abs(lfc_FA) >= lfc_thr,
#       
#       sig_EDITED =
#         !is.na(padj_EDITED) &
#         padj_EDITED < padj_thr &
#         abs(lfc_EDITED) >= lfc_thr,
#       
#       # ========================================================
#       # Rescue score
#       #
#       # 0   = aucun rescue
#       # 1   = rescue complet
#       # >1  = dépassement / reversal
#       # ========================================================
#       
#       rescue_score =
#         1 - abs(lfc_EDITED) / pmax(abs(lfc_FA), 1e-10),
#       
#       # ========================================================
#       # Catégories
#       # ========================================================
#       
#       rescue_status = dplyr::case_when(
#         
#         # FA n'est pas réellement dérégulé
#         !sig_FA ~ "Not FA-DEG",
#         
#         # DEG dans FA, mais plus DEG après édition
#         sig_FA &
#           !sig_EDITED &
#           abs(lfc_EDITED) < abs(lfc_FA) ~ "Rescued",
#         
#         # Toujours DEG mais effet réduit
#         sig_FA &
#           sig_EDITED &
#           abs(lfc_EDITED) < abs(lfc_FA) ~ "Partial rescue",
#         
#         # Même direction mais effet similaire ou augmenté
#         sig_FA &
#           sig_EDITED &
#           sign(lfc_FA) == sign(lfc_EDITED) &
#           abs(lfc_EDITED) >= abs(lfc_FA) ~ "Not rescued",
#         
#         # Changement de direction
#         sig_FA &
#           sign(lfc_FA) != sign(lfc_EDITED) ~ "Reversed",
#         
#         TRUE ~ "Other"
#       )
#     )
#   
#   # ============================================================
#   # Facteur pour ordre légende
#   # ============================================================
#   
#   status_levels <- c(
#     "Not FA-DEG",
#     "Rescued",
#     "Partial rescue",
#     "Not rescued",
#     "Reversed",
#     "Other"
#   )
#   
#   df$rescue_status <- factor(
#     df$rescue_status,
#     levels = status_levels
#   )
#   
#   # ============================================================
#   # Couleurs
#   # ============================================================
#   
#   color_values <- c(
#     "Not FA-DEG"    = "grey75",
#     "Rescued"       = "#2166ac",
#     "Partial rescue"= "#67a9cf",
#     "Not rescued"   = "#d6604d",
#     "Reversed"      = "#7b2d8b",
#     "Other"         = "grey50"
#   )
#   
#   # ============================================================
#   # Corrélation
#   # ============================================================
#   
#   r_val <- cor(
#     df$lfc_FA,
#     df$lfc_EDITED,
#     use = "complete.obs"
#   )
#   
#   r_test <- cor.test(
#     df$lfc_FA,
#     df$lfc_EDITED
#   )
#   
#   # ============================================================
#   # Régression
#   # ============================================================
#   
#   fit <- lm(
#     lfc_EDITED ~ lfc_FA,
#     data = df
#   )
#   
#   slope <- coef(fit)[["lfc_FA"]]
#   
#   r_squared <- summary(fit)$r.squared
#   
#   slope_p <- summary(fit)$coefficients[
#     "lfc_FA",
#     "Pr(>|t|)"
#   ]
#   
#   # ============================================================
#   # Pourcentage de rescue
#   # ============================================================
#   
#   fa_deg <- df %>%
#     dplyr::filter(sig_FA)
#   
#   n_FA_DEG <- nrow(fa_deg)
#   
#   n_rescued <- sum(
#     fa_deg$rescue_status == "Rescued",
#     na.rm = TRUE
#   )
#   
#   n_partial <- sum(
#     fa_deg$rescue_status == "Partial rescue",
#     na.rm = TRUE
#   )
#   
#   n_not_rescued <- sum(
#     fa_deg$rescue_status == "Not rescued",
#     na.rm = TRUE
#   )
#   
#   n_reversed <- sum(
#     fa_deg$rescue_status == "Reversed",
#     na.rm = TRUE
#   )
#   
#   pct_rescued <- ifelse(
#     n_FA_DEG > 0,
#     100 * n_rescued / n_FA_DEG,
#     NA
#   )
#   
#   pct_rescue_or_partial <- ifelse(
#     n_FA_DEG > 0,
#     100 * (n_rescued + n_partial) / n_FA_DEG,
#     NA
#   )
#   
#   message(
#     sprintf(
#       "  %s | FA-DEGs=%d | Rescued=%d (%.1f%%) | Partial=%d | Not rescued=%d | Reversed=%d",
#       dataset_name,
#       n_FA_DEG,
#       n_rescued,
#       pct_rescued,
#       n_partial,
#       n_not_rescued,
#       n_reversed
#     )
#   )
#   # 
#   # # ============================================================
#   # # Labels
#   # # ============================================================
#   # 
#   # # Priorité :
#   # # 1. rescued
#   # # 2. reversed
#   # # 3. partial
#   # # 4. top effet FA
#   # 
#   # label_df <- df %>%
#   #   dplyr::filter(
#   #     rescue_status %in%
#   #       c(
#   #         "Rescued",
#   #         "Partial rescue",
#   #         "Reversed",
#   #         "Not rescued"
#   #       )
#   #   ) %>%
#   #   dplyr::mutate(
#   #     label_priority = dplyr::case_when(
#   #       rescue_status == "Rescued" ~ 1,
#   #       rescue_status == "Reversed" ~ 2,
#   #       rescue_status == "Partial rescue" ~ 3,
#   #       TRUE ~ 4
#   #     )
#   #   ) %>%
#   #   dplyr::arrange(
#   #     label_priority,
#   #     dplyr::desc(abs(lfc_FA))
#   #   ) %>%
#   #   dplyr::slice_head(
#   #     n = top_n_labels
#   #   )
#   # 
#   
#   
#   # ============================================================
#   # Labels
#   #
#   # 1. Top 30 de CHAQUE catégorie
#   # 2. Top 50 par QUADRANT
#   # 3. Suppression des doublons
#   # ============================================================
#   
#   N_CATEGORY <- 30
#   N_QUADRANT <- 50
#   
#   
#   # ------------------------------------------------------------
#   # 1. TOP 30 DE CHAQUE CATÉGORIE
#   # ------------------------------------------------------------
#   
#   label_category <- df %>%
#     dplyr::filter(
#       rescue_status %in% c(
#         "Rescued",
#         "Partial rescue",
#         "Not rescued",
#         "Reversed"
#       )
#     ) %>%
#     dplyr::group_by(rescue_status) %>%
#     dplyr::arrange(
#       dplyr::desc(abs(lfc_FA)),
#       .by_group = TRUE
#     ) %>%
#     dplyr::slice_head(n = N_CATEGORY) %>%
#     dplyr::ungroup()
#   
#   
#   # ------------------------------------------------------------
#   # 2. Définition des quadrants
#   # ------------------------------------------------------------
#   
#   df <- df %>%
#     dplyr::mutate(
#       
#       quadrant = dplyr::case_when(
#         
#         lfc_FA >= 0 & lfc_EDITED >= 0 ~
#           "FA + / Edited +",
#         
#         lfc_FA < 0 & lfc_EDITED < 0 ~
#           "FA - / Edited -",
#         
#         lfc_FA >= 0 & lfc_EDITED < 0 ~
#           "FA + / Edited -",
#         
#         lfc_FA < 0 & lfc_EDITED >= 0 ~
#           "FA - / Edited +",
#         
#         TRUE ~ "Other"
#       )
#     )
#   
#   
#   # ------------------------------------------------------------
#   # 3. TOP 50 PAR QUADRANT
#   # ------------------------------------------------------------
#   
#   label_quadrant <- df %>%
#     dplyr::filter(
#       rescue_status != "Not FA-DEG",
#       rescue_status != "Other",
#       quadrant != "Other"
#     ) %>%
#     dplyr::group_by(quadrant) %>%
#     dplyr::arrange(
#       dplyr::desc(abs(lfc_FA)),
#       .by_group = TRUE
#     ) %>%
#     dplyr::slice_head(n = N_QUADRANT) %>%
#     dplyr::ungroup()
#   
#   
#   # ------------------------------------------------------------
#   # 4. Combiner
#   # ------------------------------------------------------------
#   
#   label_df <- dplyr::bind_rows(
#     label_category,
#     label_quadrant
#   ) %>%
#     dplyr::distinct(
#       GENEID,
#       .keep_all = TRUE
#     )
#   
#   
#   # ------------------------------------------------------------
#   # 5. Vérification
#   # ------------------------------------------------------------
#   
#   message(
#     "Labels sélectionnés : ",
#     nrow(label_df)
#   )
#   
#   print(
#     table(
#       label_df$rescue_status,
#       label_df$quadrant
#     )
#   )
#   
#   
#   # ============================================================
#   # Plot
#   # ============================================================
#   
#   p <- ggplot(
#     df,
#     aes(
#       x = lfc_FA,
#       y = lfc_EDITED,
#       color = rescue_status
#     )
#   ) +
#     
#     geom_point(
#       alpha = 0.5,
#       size = 1.1
#     ) +
#     
#     # Régression observée
#     geom_smooth(
#       aes(
#         x = lfc_FA,
#         y = lfc_EDITED
#       ),
#       inherit.aes = FALSE,
#       data = df,
#       method = "lm",
#       se = TRUE,
#       color = "black",
#       linewidth = 0.8
#     ) +
#     
#     # y = x
#     geom_abline(
#       slope = 1,
#       intercept = 0,
#       linetype = "dotted",
#       color = "grey30"
#     ) +
#     
#     # Axes zéro
#     geom_hline(
#       yintercept = 0,
#       linetype = "dashed",
#       color = "grey60"
#     ) +
#     
#     geom_vline(
#       xintercept = 0,
#       linetype = "dashed",
#       color = "grey60"
#     ) +
#     
#     # Labels
#     ggrepel::geom_label_repel(
#       data = label_df,
#       aes(label = GENENAME),
#       size = 2.6,
#       max.overlaps = Inf,
#       show.legend = FALSE
#     ) +
#     
#     scale_color_manual(
#       values = color_values,
#       name = NULL,
#       drop = FALSE
#     ) +
#     
#     annotate(
#       "label",
#       x = min(df$lfc_FA, na.rm = TRUE),
#       y = max(df$lfc_EDITED, na.rm = TRUE),
#       hjust = 0,
#       vjust = 1,
#       label = sprintf(
#         paste0(
#           "Pearson r = %.2f (p = %.2g)\n",
#           "Pente = %.2f (R² = %.2f, p = %.2g)\n",
#           "FA-DEGs = %d\n",
#           "Rescued = %d (%.1f%%)\n",
#           "Partial = %d\n",
#           "Not rescued = %d\n",
#           "Reversed = %d"
#         ),
#         r_val,
#         r_test$p.value,
#         slope,
#         r_squared,
#         slope_p,
#         n_FA_DEG,
#         n_rescued,
#         pct_rescued,
#         n_partial,
#         n_not_rescued,
#         n_reversed
#       ),
#       size = 3,
#       fill = "white",
#       alpha = 0.85
#     ) +
#     
#     theme_bw(base_size = 11) +
#     
#     labs(
#       title = paste0(
#         dataset_name,
#         " : disease effect and rescue"
#       ),
#       subtitle =
#         "FA vs Control compared with FA_edited vs Control",
#       x = "log2FC (FA vs Control)",
#       y = "log2FC (FA_edited vs Control)"
#     ) +
#     
#     theme(
#       legend.position = "top"
#     )
#   
#   # ============================================================
#   # Export PDF
#   # ============================================================
#   
#   out_file <- file.path(
#     out_dir,
#     paste0(
#       dataset_name,
#       "_FA_vs_FAedited_rescue_scatter.pdf"
#     )
#   )
#   
#   pdf(
#     out_file,
#     width = 10,
#     height = 9
#   )
#   
#   print(p)
#   
#   dev.off()
#   
#   message(
#     "  Rescue scatter sauvegardé : ",
#     out_file
#   )
#   
#   # ============================================================
#   # Export Excel
#   # ============================================================
#   
#   if (export_xlsx) {
#     
#     xlsx_file <- file.path(
#       out_dir,
#       paste0(
#         "genes_",
#         dataset_name,
#         "_FA_vs_FAedited_rescue.xlsx"
#       )
#     )
#     
#     writexl::write_xlsx(
#       df,
#       xlsx_file
#     )
#     
#     message(
#       "  Table rescue sauvegardée : ",
#       xlsx_file
#     )
#   }
#   
#   
#   # ============================================================
#   # Quadrant plots
#   # ============================================================
#   
#   
#   run_quadrant_plots(
#     data_df = df,
#     out_dir = file.path(out_dir, "QUADRANTS"),
#     dataset_name = dataset_name,
#     top_n = 50
#   )
#   
#   
#   
#   return(
#     invisible(df)
#   )
# }
# 




run_rescue_scatter <- function(
    all_results,
    out_dir,
    dataset_name,
    padj_thr = 0.05,
    lfc_thr = 1,
    top_n_labels = 30,
    export_xlsx = TRUE,
    
    # ------------------------------------------------------------
    # GO / GSEA
    # ------------------------------------------------------------
    run_GO = TRUE,
    run_GSEA = TRUE,
    organism_db = "org.Hs.eg.db",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    top_n_GO = 20,
    top_n_GSEA = 20,     # ------------------------------------------------------------
    # Friedreich Ataxia focused analysis
    # ------------------------------------------------------------
    run_FA_focused = TRUE,
    FA_padj_cutoff = 0.05,
    FA_keywords = c(
      # Friedreich ataxia / frataxin
      "frataxin",
      "FXN",
      
      # mitochondria
      "mitochond",
      "mitochondrial",
      "mitochondrion",
      
      # respiration / OXPHOS
      "oxidative phosphorylation",
      "respiratory",
      "electron transport",
      "electron transfer",
      "respiratory chain",
      "ATP synthesis",
      "ATP metabolic",
      
      # TCA / energy metabolism
      "tricarboxylic acid",
      "TCA cycle",
      "citric acid cycle",
      "energy metabolism",
      "cellular respiration",
      
      # iron / Fe-S
      "iron",
      "iron-sulfur",
      "iron sulfur",
      "Fe-S",
      "heme",
      
      # oxidative stress / ROS
      "oxidative stress",
      "oxidation-reduction",
      "reactive oxygen",
      "ROS",
      "redox",
      "antioxidant",
      
      # mitochondrial quality control
      "mitophagy",
      "autophagy",
      "mitochondrial quality",
      "mitochondrial dynamics",
      "mitochondrial fission",
      "mitochondrial fusion",
      
      # protein homeostasis
      "proteostasis",
      "protein folding",
      "unfolded protein",
      "endoplasmic reticulum stress",
      
      # neuronal disease relevance
      "neuron",
      "neuronal",
      "neuro",
      "axon",
      "axonal",
      "nerve",
      "synaptic",
      "synapse"
    )
) {
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  # ============================================================
  # Packages GO / GSEA
  # ============================================================
  
  if (run_GO || run_GSEA) {
    
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      warning(
        "clusterProfiler n'est pas installé. ",
        "GO/GSEA désactivées."
      )
      run_GO <- FALSE
      run_GSEA <- FALSE
    }
    
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      warning(
        "AnnotationDbi n'est pas installé. ",
        "GO/GSEA désactivées."
      )
      run_GO <- FALSE
      run_GSEA <- FALSE
    }
  }
  
  
  # ============================================================
  # Vérification des contrastes
  # ============================================================
  
  required <- c(
    "control_vs_FA",
    "control_vs_FA_edited"
  )
  
  if (!all(required %in% names(all_results))) {
    
    message(
      "  [SKIP] ", dataset_name,
      " : contrastes nécessaires absents."
    )
    
    return(invisible(NULL))
  }
  
  
  df1 <- as.data.frame(
    all_results[["control_vs_FA"]]
  )
  
  df2 <- as.data.frame(
    all_results[["control_vs_FA_edited"]]
  )
  
  
  # ============================================================
  # Récupération gene ID / symbol
  # ============================================================
  
  get_gene_col <- function(df) {
    
    if ("GENEID" %in% colnames(df)) {
      return(df$GENEID)
    }
    
    if ("gene" %in% colnames(df)) {
      return(df$gene)
    }
    
    rownames(df)
  }
  
  
  get_symbol_col <- function(df) {
    
    if ("GENENAME" %in% colnames(df)) {
      return(df$GENENAME)
    }
    
    if ("gene" %in% colnames(df)) {
      return(df$gene)
    }
    
    rownames(df)
  }
  
  
  df1$GENEID <- get_gene_col(df1)
  df1$GENENAME <- get_symbol_col(df1)
  
  df2$GENEID <- get_gene_col(df2)
  
  
  # ============================================================
  # Merge
  # ============================================================
  
  df_a <- df1 %>%
    dplyr::select(
      GENEID,
      GENENAME,
      log2FoldChange,
      padj
    ) %>%
    dplyr::rename(
      lfc_FA = log2FoldChange,
      padj_FA = padj
    )
  
  
  df_b <- df2 %>%
    dplyr::select(
      GENEID,
      log2FoldChange,
      padj
    ) %>%
    dplyr::rename(
      lfc_EDITED = log2FoldChange,
      padj_EDITED = padj
    )
  
  
  df <- dplyr::inner_join(
    df_a,
    df_b,
    by = "GENEID"
  ) %>%
    dplyr::filter(
      !is.na(lfc_FA),
      !is.na(lfc_EDITED)
    )
  
  
  # ============================================================
  # Significance
  # ============================================================
  
  df <- df %>%
    dplyr::mutate(
      
      sig_FA =
        !is.na(padj_FA) &
        padj_FA < padj_thr &
        abs(lfc_FA) >= lfc_thr,
      
      sig_EDITED =
        !is.na(padj_EDITED) &
        padj_EDITED < padj_thr &
        abs(lfc_EDITED) >= lfc_thr,
      
      # ========================================================
      # Rescue score
      #
      # +1 = rescue complet théorique
      #  0  = effet inchangé
      # <0  = effet augmenté
      #
      # ATTENTION :
      # le score est basé sur la magnitude absolue.
      # Le reversal est traité séparément.
      # ========================================================
      
      rescue_score =
        1 -
        abs(lfc_EDITED) /
        pmax(abs(lfc_FA), 1e-10),
      
      
      # ========================================================
      # Delta de magnitude
      #
      # > 0 = diminution de l'effet
      # < 0 = augmentation de l'effet
      # ========================================================
      
      delta_abs_LFC =
        abs(lfc_FA) -
        abs(lfc_EDITED),
      
      
      # ========================================================
      # Catégories
      # ========================================================
      
      rescue_status = dplyr::case_when(
        
        !sig_FA ~
          "Not FA-DEG",
        
        sig_FA &
          !sig_EDITED &
          abs(lfc_EDITED) < abs(lfc_FA) ~
          "Rescued",
        
        sig_FA &
          sig_EDITED &
          abs(lfc_EDITED) < abs(lfc_FA) ~
          "Partial rescue",
        
        sig_FA &
          sig_EDITED &
          sign(lfc_FA) == sign(lfc_EDITED) &
          abs(lfc_EDITED) >= abs(lfc_FA) ~
          "Not rescued",
        
        sig_FA &
          sign(lfc_FA) != sign(lfc_EDITED) ~
          "Reversed",
        
        TRUE ~
          "Other"
      )
    )
  
  
  # ============================================================
  # Ordre des catégories
  # ============================================================
  
  status_levels <- c(
    "Not FA-DEG",
    "Rescued",
    "Partial rescue",
    "Not rescued",
    "Reversed",
    "Other"
  )
  
  
  df$rescue_status <- factor(
    df$rescue_status,
    levels = status_levels
  )
  
  
  # ============================================================
  # Couleurs
  # ============================================================
  
  color_values <- c(
    "Not FA-DEG" =
      "grey75",
    
    "Rescued" =
      "#2166ac",
    
    "Partial rescue" =
      "#67a9cf",
    
    "Not rescued" =
      "#d6604d",
    
    "Reversed" =
      "#7b2d8b",
    
    "Other" =
      "grey50"
  )
  
  
  # ============================================================
  # Corrélation
  # ============================================================
  
  r_val <- cor(
    df$lfc_FA,
    df$lfc_EDITED,
    use = "complete.obs"
  )
  
  
  r_test <- cor.test(
    df$lfc_FA,
    df$lfc_EDITED
  )
  
  
  # ============================================================
  # Régression
  # ============================================================
  
  fit <- lm(
    lfc_EDITED ~ lfc_FA,
    data = df
  )
  
  slope <- coef(fit)[["lfc_FA"]]
  
  r_squared <-
    summary(fit)$r.squared
  
  slope_p <-
    summary(fit)$coefficients[
      "lfc_FA",
      "Pr(>|t|)"
    ]
  
  
  # ============================================================
  # Statistiques rescue
  # ============================================================
  
  fa_deg <- df %>%
    dplyr::filter(sig_FA)
  
  
  n_FA_DEG <- nrow(fa_deg)
  
  
  n_rescued <- sum(
    fa_deg$rescue_status == "Rescued",
    na.rm = TRUE
  )
  
  
  n_partial <- sum(
    fa_deg$rescue_status == "Partial rescue",
    na.rm = TRUE
  )
  
  
  n_not_rescued <- sum(
    fa_deg$rescue_status == "Not rescued",
    na.rm = TRUE
  )
  
  
  n_reversed <- sum(
    fa_deg$rescue_status == "Reversed",
    na.rm = TRUE
  )
  
  
  pct_rescued <- ifelse(
    n_FA_DEG > 0,
    100 * n_rescued / n_FA_DEG,
    NA
  )
  
  
  pct_rescue_or_partial <- ifelse(
    n_FA_DEG > 0,
    100 *
      (n_rescued + n_partial) /
      n_FA_DEG,
    NA
  )
  
  
  message(
    sprintf(
      paste0(
        "  %s | FA-DEGs=%d | ",
        "Rescued=%d (%.1f%%) | ",
        "Partial=%d | ",
        "Not rescued=%d | ",
        "Reversed=%d"
      ),
      dataset_name,
      n_FA_DEG,
      n_rescued,
      pct_rescued,
      n_partial,
      n_not_rescued,
      n_reversed
    )
  )
  
  
  # ============================================================
  # QUADRANTS
  # ============================================================
  
  df <- df %>%
    dplyr::mutate(
      
      quadrant = dplyr::case_when(
        
        lfc_FA >= 0 &
          lfc_EDITED >= 0 ~
          "FA + / Edited +",
        
        lfc_FA < 0 &
          lfc_EDITED < 0 ~
          "FA - / Edited -",
        
        lfc_FA >= 0 &
          lfc_EDITED < 0 ~
          "FA + / Edited -",
        
        lfc_FA < 0 &
          lfc_EDITED >= 0 ~
          "FA - / Edited +",
        
        TRUE ~
          "Other"
      )
    )
  
  
  # ============================================================
  # LABELS
  # ============================================================
  
  # N_CATEGORY <- 30
  # N_QUADRANT <- 50
  N_CATEGORY <- 50
  N_QUADRANT <- 100
  
  label_category <- df %>%
    dplyr::filter(
      rescue_status %in%
        c(
          "Rescued",
          "Partial rescue",
          "Not rescued",
          "Reversed"
        )
    ) %>%
    dplyr::group_by(rescue_status) %>%
    dplyr::arrange(
      dplyr::desc(abs(lfc_FA)),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(
      n = N_CATEGORY
    ) %>%
    dplyr::ungroup()
  
  
  label_quadrant <- df %>%
    dplyr::filter(
      rescue_status != "Not FA-DEG",
      rescue_status != "Other",
      quadrant != "Other"
    ) %>%
    dplyr::group_by(quadrant) %>%
    dplyr::arrange(
      dplyr::desc(abs(lfc_FA)),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(
      n = N_QUADRANT
    ) %>%
    dplyr::ungroup()
  
  
  label_df <- dplyr::bind_rows(
    label_category,
    label_quadrant
  ) %>%
    dplyr::distinct(
      GENEID,
      .keep_all = TRUE
    )
  
  
  message(
    "Labels sélectionnés : ",
    nrow(label_df)
  )
  
  
  # ============================================================
  # SCATTER PRINCIPAL
  # ============================================================
  
  p <- ggplot(
    df,
    aes(
      x = lfc_FA,
      y = lfc_EDITED,
      color = rescue_status
    )
  ) +
    
    geom_point(
      alpha = 0.5,
      size = 1.1
    ) +
    
    geom_smooth(
      aes(
        x = lfc_FA,
        y = lfc_EDITED
      ),
      inherit.aes = FALSE,
      data = df,
      method = "lm",
      se = TRUE,
      color = "black",
      linewidth = 0.8
    ) +
    
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted",
      color = "grey30"
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    
    ggrepel::geom_label_repel(
      data = label_df,
      aes(
        label = GENENAME
      ),
      size = 2.6,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    scale_color_manual(
      values = color_values,
      name = NULL,
      drop = FALSE
    ) +
    
    annotate(
      "label",
      x = min(
        df$lfc_FA,
        na.rm = TRUE
      ),
      y = max(
        df$lfc_EDITED,
        na.rm = TRUE
      ),
      hjust = 0,
      vjust = 1,
      label = sprintf(
        paste0(
          "Pearson r = %.2f (p = %.2g)\n",
          "Slope = %.2f (R² = %.2f, p = %.2g)\n",
          "FA-DEGs = %d\n",
          "Rescued = %d (%.1f%%)\n",
          "Partial = %d\n",
          "Not rescued = %d\n",
          "Reversed = %d"
        ),
        r_val,
        r_test$p.value,
        slope,
        r_squared,
        slope_p,
        n_FA_DEG,
        n_rescued,
        pct_rescued,
        n_partial,
        n_not_rescued,
        n_reversed
      ),
      size = 3,
      fill = "white",
      alpha = 0.85
    ) +
    
    theme_bw(
      base_size = 11
    ) +
    
    labs(
      title = paste0(
        dataset_name,
        " : disease effect and rescue"
      ),
      subtitle =
        "FA vs Control compared with FA_edited vs Control",
      x = "log2FC (FA vs Control)",
      y = "log2FC (FA_edited vs Control)"
    ) +
    
    theme(
      legend.position = "top"
    )
  
  
  # ============================================================
  # Export scatter
  # ============================================================
  
  out_file <- file.path(
    out_dir,
    paste0(
      dataset_name,
      "_FA_vs_FAedited_rescue_scatter.pdf"
    )
  )
  
  
  pdf(
    out_file,
    width = 10,
    height = 9
  )
  
  print(p)
  
  dev.off()
  
  
  message(
    "  Rescue scatter sauvegardé : ",
    out_file
  )
  
  
  # ============================================================
  # ============================================================
  # GO ENRICHMENT
  # ============================================================
  # ============================================================
  
  GO_results <- list()
  
  
  if (run_GO) {
    
    message(
      "\n  ========================================\n",
      "  GO ENRICHMENT : ",
      dataset_name,
      "\n  ========================================"
    )
    
    
    # ----------------------------------------------------------
    # Charger la base d'annotation
    # ----------------------------------------------------------
    
    db_loaded <- tryCatch(
      
      get(
        organism_db,
        envir = asNamespace(
          "clusterProfiler"
        )
      ),
      
      error = function(e) {
        
        tryCatch(
          
          get(
            organism_db,
            envir = .GlobalEnv
          ),
          
          error = function(e2) NULL
        )
      }
    )
    
    
    # ----------------------------------------------------------
    # Méthode plus robuste :
    # charger le package explicitement
    # ----------------------------------------------------------
    
    if (is.null(db_loaded)) {
      
      tryCatch(
        
        suppressPackageStartupMessages(
          library(
            organism_db,
            character.only = TRUE
          )
        ),
        
        error = function(e) {
          
          message(
            "  [WARNING] Impossible de charger ",
            organism_db,
            ". GO désactivée."
          )
          
          run_GO <<- FALSE
        }
      )
    }
    
    
    if (run_GO) {
      
      # --------------------------------------------------------
      # Obtenir l'objet OrgDb
      # --------------------------------------------------------
      
      orgdb <- tryCatch(
        
        get(
          organism_db
        ),
        
        error = function(e) NULL
      )
      
      
      if (is.null(orgdb)) {
        
        message(
          "  [WARNING] Objet ",
          organism_db,
          " introuvable. GO désactivée."
        )
        
        run_GO <- FALSE
        
      }
    }
    
    
    if (run_GO) {
      
      # --------------------------------------------------------
      # Identifier le type de GENEID
      #
      # Si ce sont déjà des ENTREZ IDs :
      # conserver keyType = ENTREZID.
      #
      # Sinon on essaie ENSEMBL / SYMBOL.
      # --------------------------------------------------------
      
      id_candidates <- c(
        "ENSEMBL",
        "SYMBOL",
        "ENTREZID"
      )
      
      
      keytypes_available <-
        AnnotationDbi::keytypes(orgdb)
      
      
      key_type <- NULL
      
      
      for (candidate in id_candidates) {
        
        if (candidate %in% keytypes_available) {
          
          # Vérification grossière du format
          test_ids <- head(
            unique(
              as.character(
                df$GENEID
              )
            ),
            100
          )
          
          mapped_test <- tryCatch(
            
            AnnotationDbi::select(
              orgdb,
              keys = test_ids,
              columns = candidate,
              keytype = candidate
            ),
            
            error = function(e) NULL
          )
          
          if (
            !is.null(mapped_test) &&
            nrow(mapped_test) > 0 &&
            sum(
              !is.na(
                mapped_test[[candidate]]
              )
            ) > 0
          ) {
            
            key_type <- candidate
            
            break
          }
        }
      }
      
      
      if (is.null(key_type)) {
        
        message(
          "  [WARNING] Impossible d'identifier le type ",
          "de GENEID pour ",
          organism_db,
          "."
        )
        
      } else {
        
        message(
          "  Type d'identifiant utilisé pour GO : ",
          key_type
        )
        
        
        # ------------------------------------------------------
        # Tous les gènes testés = background
        # ------------------------------------------------------
        
        all_ids <- unique(
          as.character(
            df$GENEID
          )
        )
        
        
        # ------------------------------------------------------
        # Mapping vers ENTREZ
        # ------------------------------------------------------
        
        id_map <- tryCatch(
          
          AnnotationDbi::select(
            orgdb,
            keys = all_ids,
            columns = c(
              key_type,
              "ENTREZID"
            ),
            keytype = key_type
          ),
          
          error = function(e) NULL
        )
        
        
        if (
          is.null(id_map) ||
          !"ENTREZID" %in%
          colnames(id_map)
        ) {
          
          message(
            "  [WARNING] Mapping ENTREZID impossible."
          )
          
        } else {
          
          id_map <- id_map %>%
            dplyr::filter(
              !is.na(ENTREZID)
            ) %>%
            dplyr::distinct(
              .data[[key_type]],
              ENTREZID
            )
          
          
          # ----------------------------------------------------
          # Fonction GO
          # ----------------------------------------------------
          
          run_one_GO <- function(
    gene_ids,
    status_name
          ) {
            
            genes <- unique(
              as.character(
                gene_ids
              )
            )
            
            
            genes_entrez <-
              id_map %>%
              dplyr::filter(
                .data[[key_type]] %in% genes
              ) %>%
              dplyr::pull(
                ENTREZID
              ) %>%
              unique()
            
            
            universe_entrez <-
              id_map %>%
              dplyr::pull(
                ENTREZID
              ) %>%
              unique()
            
            
            message(
              "    ",
              status_name,
              " : ",
              length(genes),
              " gènes | ",
              length(genes_entrez),
              " ENTREZ IDs"
            )
            
            
            if (
              length(genes_entrez) < 5
            ) {
              
              message(
                "    [SKIP] Trop peu de gènes."
              )
              
              return(NULL)
            }
            
            
            ego <- tryCatch(
              
              clusterProfiler::enrichGO(
                
                gene = genes_entrez,
                
                universe = universe_entrez,
                
                OrgDb = orgdb,
                
                keyType = "ENTREZID",
                
                ont = ont,
                
                pAdjustMethod = "BH",
                
                pvalueCutoff =
                  pvalueCutoff,
                
                qvalueCutoff =
                  qvalueCutoff,
                
                minGSSize =
                  minGSSize,
                
                maxGSSize =
                  maxGSSize,
                
                readable = TRUE
              ),
              
              error = function(e) {
                
                message(
                  "    [ERROR] enrichGO : ",
                  e$message
                )
                
                NULL
              }
            )
            
            
            if (
              is.null(ego)
            ) {
              
              return(NULL)
            }
            
            
            ego_df <-
              as.data.frame(ego)
            
            fa_keywords <- paste(
              c(
                "mitochond",
                "oxidative",
                "respiratory",
                "electron transport",
                "iron",
                "iron ion",
                "oxidation-reduction",
                "reactive oxygen",
                "oxidative stress",
                "autophagy",
                "proteostasis",
                "protein folding",
                "neuron",
                "neuronal",
                "axon",
                "neuro",
                "energy metabolism",
                "metabolic process"
              ),
              collapse = "|"
            )
            
            ego_df <- ego_df %>%
              dplyr::mutate(
                FA_relevant =
                  grepl(
                    fa_keywords,
                    Description,
                    ignore.case = TRUE
                  )
              )
            
            if (
              nrow(ego_df) == 0
            ) {
              
              message(
                "    Aucun terme GO significatif."
              )
              
              return(NULL)
            }
            
            
            # --------------------------------------------------
            # Export PDF
            # --------------------------------------------------
            
            go_plot <- clusterProfiler::dotplot(
              ego,
              showCategory = top_n_GO
            ) +
              
              ggtitle(
                paste0(
                  dataset_name,
                  " - GO ",
                  ont,
                  " - ",
                  status_name
                )
              )
            
            
            go_file <- file.path(
              out_dir,
              paste0(
                dataset_name,
                "_GO_",
                gsub(
                  " ",
                  "_",
                  status_name
                ),
                ".pdf"
              )
            )
            
            
            pdf(
              go_file,
              width = 11,
              height = 8
            )
            
            print(go_plot)
            
            dev.off()
            
            
            # ----------------------------------------------------------
            # Plot GO : uniquement les termes FA-relevant
            # ----------------------------------------------------------
            
            ego_fa <- ego_df %>%
              dplyr::filter(
                FA_relevant
              )
            
            
            if (nrow(ego_fa) > 0) {
              
              ego_fa_plot <- ego_fa %>%
                dplyr::arrange(
                  p.adjust
                ) %>%
                dplyr::slice_head(
                  n = top_n_GO
                )
              
              # 
              # go_fa_plot <- ggplot(
              #   ego_fa_plot,
              #   aes(
              #     x = reorder(
              #       Description,
              #       -log10(p.adjust)
              #     ),
              #     y = -log10(p.adjust),
              #     size = Count
              #   )
              # ) +
              #   
              #   geom_point() +
              #   
              #   coord_flip() +
              #   
              #   theme_bw(
              #     base_size = 11
              #   ) +
              #   
              #   labs(
              #     title = paste0(
              #       dataset_name,
              #       " - GO FA-relevant - ",
              #       status_name
              #     ),
              #     x = NULL,
              #     y = "-log10 adjusted p-value",
              #     size = "Gene count"
              #   )
              # 
              
              
              # ----------------------------------------------------------
              # Préparation du GeneRatio
              # ----------------------------------------------------------
              
              ego_fa_plot <- ego_fa %>%
                dplyr::arrange(
                  p.adjust
                ) %>%
                dplyr::slice_head(
                  n = top_n_GO
                ) %>%
                dplyr::mutate(
                  GeneRatio_numeric =
                    sapply(
                      GeneRatio,
                      function(x) {
                        parts <- strsplit(
                          as.character(x),
                          "/"
                        )[[1]]
                        
                        as.numeric(parts[1]) /
                          as.numeric(parts[2])
                      }
                    )
                )
              
              
              # ----------------------------------------------------------
              # Plot GO FA-relevant
              # ----------------------------------------------------------
              
              go_fa_plot <- ggplot(
                ego_fa_plot,
                aes(
                  x = GeneRatio_numeric,
                  y = reorder(
                    Description,
                    GeneRatio_numeric
                  ),
                  size = Count,
                  color = p.adjust
                )
              ) +
                
                geom_point() +
                
                scale_color_gradient(
                  low = "red",
                  high = "blue",
                  name = "Adjusted\np-value"
                ) +
                
                theme_bw(
                  base_size = 11
                ) +
                
                labs(
                  title = paste0(
                    dataset_name,
                    " - GO FA-relevant - ",
                    status_name
                  ),
                  x = "Gene Ratio",
                  y = NULL,
                  size = "Gene Count"
                )
              
              
              
              go_fa_file <- file.path(
                out_dir,
                paste0(
                  dataset_name,
                  "_GO_FA_relevant_",
                  gsub(
                    " ",
                    "_",
                    status_name
                  ),
                  ".pdf"
                )
              )
              
              
              pdf(
                go_fa_file,
                width = 11,
                height = 8
              )
              
              print(
                go_fa_plot
              )
              
              dev.off()
              
              
              message(
                "    GO FA-relevant sauvegardé : ",
                go_fa_file
              )
              
            } else {
              
              message(
                "    Aucun terme GO FA-relevant."
              )
            }
            
            
            # --------------------------------------------------
            # Export CSV
            # --------------------------------------------------
            
            csv_file <- file.path(
              out_dir,
              paste0(
                dataset_name,
                "_GO_",
                gsub(
                  " ",
                  "_",
                  status_name
                ),
                ".csv"
              )
            )
            
            
            write.csv(
              ego_df,
              csv_file,
              row.names = FALSE
            )
            
            
            return(
              ego_df
            )
          }
          
          
          # ----------------------------------------------------
          # RESCUED
          # ----------------------------------------------------
          
          GO_results$Rescued <-
            run_one_GO(
              df %>%
                dplyr::filter(
                  rescue_status ==
                    "Rescued"
                ) %>%
                dplyr::pull(
                  GENEID
                ),
              "Rescued"
            )
          
          
          # ----------------------------------------------------
          # NOT RESCUED
          # ----------------------------------------------------
          
          GO_results$Not_rescued <-
            run_one_GO(
              df %>%
                dplyr::filter(
                  rescue_status ==
                    "Not rescued"
                ) %>%
                dplyr::pull(
                  GENEID
                ),
              "Not rescued"
            )
          
          
          # ----------------------------------------------------
          # PARTIAL RESCUE
          # ----------------------------------------------------
          
          GO_results$Partial_rescue <-
            run_one_GO(
              df %>%
                dplyr::filter(
                  rescue_status ==
                    "Partial rescue"
                ) %>%
                dplyr::pull(
                  GENEID
                ),
              "Partial rescue"
            )
          
          
          # ----------------------------------------------------
          # REVERSED
          # ----------------------------------------------------
          
          GO_results$Reversed <-
            run_one_GO(
              df %>%
                dplyr::filter(
                  rescue_status ==
                    "Reversed"
                ) %>%
                dplyr::pull(
                  GENEID
                ),
              "Reversed"
            )
        }
      }
    }
  }
  
  # 
  # # ============================================================
  # # ============================================================
  # # GSEA / FGSEA
  # # ============================================================
  # # ============================================================
  # 
  # GSEA_result <- NULL
  # 
  # 
  # if (run_GSEA) {
  #   
  #   message(
  #     "\n  ========================================\n",
  #     "  GSEA RESCUE : ",
  #     dataset_name,
  #     "\n  ========================================"
  #   )
  #   
  #   
  #   # ----------------------------------------------------------
  #   # IMPORTANT :
  #   #
  #   # On classe les gènes selon :
  #   #
  #   # delta_abs_LFC =
  #   # abs(LFC_FA) - abs(LFC_EDITED)
  #   #
  #   # > 0 = l'effet FA diminue après édition
  #   #       => rescue
  #   #
  #   # < 0 = l'effet FA est maintenu / augmente
  #   #       => not rescued
  #   #
  #   # Cela permet de ne pas dépendre uniquement
  #   # des seuils DEG.
  #   # ----------------------------------------------------------
  #     # gsea_df <- df %>%
  #   #   dplyr::filter(
  #   #     !is.na(delta_abs_LFC),
  #   #     !is.na(GENEID)
  #   #   ) %>%  
  #   gsea_df <- df %>%
  #     dplyr::filter(
  #       # sig_FA,
  #       !is.na(delta_abs_LFC),
  #       !is.na(GENEID)
  #     )%>%
  #     dplyr::group_by(
  #       GENEID
  #     ) %>%
  #     dplyr::slice_max(
  #       order_by = abs(delta_abs_LFC),
  #       n = 1,
  #       with_ties = FALSE
  #     ) %>%
  #     dplyr::ungroup()
  #   
  #   
  #   gene_rank <- gsea_df$delta_abs_LFC
  #   
  #   names(gene_rank) <-
  #     as.character(
  #       gsea_df$GENEID
  #     )
  #   
  #   
  #   gene_rank <- sort(
  #     gene_rank,
  #     decreasing = TRUE
  #   )
  #   
  #   
  #   # ----------------------------------------------------------
  #   # Retirer NA / doublons
  #   # ----------------------------------------------------------
  #   
  #   gene_rank <-
  #     gene_rank[
  #       !is.na(gene_rank)
  #     ]
  #   
  #   
  #   gene_rank <-
  #     gene_rank[
  #       !duplicated(
  #         names(gene_rank)
  #       )
  #     ]
  #   
  #   
  #   # ----------------------------------------------------------
  #   # Mapping vers ENTREZ
  #   # ----------------------------------------------------------
  #   # 
  #   # if (
  #   #   exists("id_map") &&
  #   #   !is.null(id_map)
  #   # ) {
  #   #   
  #   #   rank_map <- id_map %>%
  #   #     dplyr::filter(
  #   #       .data[[key_type]] %in%
  #   #         names(gene_rank)
  #   #     ) %>%
  #   #     dplyr::select(
  #   #       gene_id =
  #   #         .data[[key_type]],
  #   #       ENTREZID
  #   #     ) %>%
  #   #     dplyr::filter(
  #   #       !is.na(ENTREZID)
  #   #     ) %>%
  #   #     dplyr::distinct(
  #   #       gene_id,
  #   #       .keep_all = TRUE
  #   #     )
  #   #   
  #   #   
  #   #   rank_df <- data.frame(
  #   #     gene_id =
  #   #       names(gene_rank),
  #   #     score =
  #   #       as.numeric(gene_rank),
  #   #     stringsAsFactors = FALSE
  #   #   ) %>%
  #   #     dplyr::inner_join(
  #   #       rank_map,
  #   #       by = "gene_id"
  #   #     )
  #   #   
  #   #   
  #   #   # --------------------------------------------------------
  #   #   # Si plusieurs IDs correspondent au même ENTREZID :
  #   #   # garder le score absolu maximal
  #   #   # --------------------------------------------------------
  #   #   
  #   #   rank_df <- rank_df %>%
  #   #     dplyr::group_by(
  #   #       ENTREZID
  #   #     ) %>%
  #   #     dplyr::slice_max(
  #   #       order_by = abs(score),
  #   #       n = 1,
  #   #       with_ties = FALSE
  #   #     ) %>%
  #   #     dplyr::ungroup()
  #   #   
  #   #   
  #   #   geneList <- rank_df$score
  #   #   
  #   #   names(geneList) <-
  #   #     as.character(
  #   #       rank_df$ENTREZID
  #   #     )
  #   #   
  #   #   
  #   #   geneList <- sort(
  #   #     geneList,
  #   #     decreasing = TRUE
  #   #   )
  #   #   
  #   #   
  #   #   
  #   #   GSEA_rescue <- clusterProfiler::gseGO(
  #   #     geneList = geneList,
  #   #     OrgDb = orgdb,
  #   #     keyType = "ENTREZID",
  #   #     ont = ont,
  #   #     minGSSize = minGSSize,
  #   #     maxGSSize = maxGSSize,
  #   #     pvalueCutoff = pvalueCutoff,
  #   #     pAdjustMethod = "BH",
  #   #     verbose = FALSE
  #   #   )
  #   #   
  #   #   geneList_persist <- -geneList
  #   #   
  #   #   geneList_persist <- sort(
  #   #     geneList_persist,
  #   #     decreasing = TRUE
  #   #   )
  #   #   
  #   #   
  #   #   GSEA_persist <- clusterProfiler::gseGO(
  #   #     geneList = geneList_persist,
  #   #     OrgDb = orgdb,
  #   #     keyType = "ENTREZID",
  #   #     ont = ont,
  #   #     minGSSize = minGSSize,
  #   #     maxGSSize = maxGSSize,
  #   #     pvalueCutoff = pvalueCutoff,
  #   #     pAdjustMethod = "BH",
  #   #     verbose = FALSE
  #   #   )
  #   #   
  #   #   
  #   #   
  #   #   
  #   #   # --------------------------------------------------------
  #   #   # GSEA clusterProfiler
  #   #   #
  #   #   # Cela utilise une approche GSEA complète.
  #   #   # clusterProfiler utilise fgsea en interne
  #   #   # dans les versions modernes.
  #   #   # --------------------------------------------------------
  #   #   
  #   #   if (
  #   #     length(geneList) >= 10
  #   #   ) {
  #   #     
  #   #     GSEA_result <- tryCatch(
  #   #       
  #   #       clusterProfiler::gseGO(
  #   #         
  #   #         geneList =
  #   #           geneList,
  #   #         
  #   #         OrgDb =
  #   #           orgdb,
  #   #         
  #   #         keyType =
  #   #           "ENTREZID",
  #   #         
  #   #         ont =
  #   #           ont,
  #   #         
  #   #         minGSSize =
  #   #           minGSSize,
  #   #         
  #   #         maxGSSize =
  #   #           maxGSSize,
  #   #         
  #   #         pvalueCutoff =
  #   #           pvalueCutoff,
  #   #         
  #   #         pAdjustMethod =
  #   #           "BH",
  #   #         
  #   #         verbose =
  #   #           FALSE
  #   #       ),
  #   #       
  #   #       error = function(e) {
  #   #         
  #   #         message(
  #   #           "  [ERROR] GSEA : ",
  #   #           e$message
  #   #         )
  #   #         
  #   #         NULL
  #   #       }
  #   #     )
  #   #     
  #   #     
  #   #     if (
  #   #       !is.null(GSEA_result)
  #   #     ) {
  #   #       
  #   #       GSEA_df <-
  #   #         as.data.frame(
  #   #           GSEA_result
  #   #         )
  #   #       
  #   #       
  #   #       if (
  #   #         nrow(GSEA_df) > 0
  #   #       ) {
  #   #         
  #   #         # --------------------------------------------------
  #   #         # Export CSV
  #   #         # --------------------------------------------------
  #   #         
  #   #         gsea_csv <- file.path(
  #   #           out_dir,
  #   #           paste0(
  #   #             dataset_name,
  #   #             "_GSEA_rescue_GO.csv"
  #   #           )
  #   #         )
  #   #         
  #   #         
  #   #         write.csv(
  #   #           GSEA_df,
  #   #           gsea_csv,
  #   #           row.names = FALSE
  #   #         )
  #   #         
  #   #         
  #   #         # --------------------------------------------------
  #   #         # Dotplot
  #   #         # --------------------------------------------------
  #   #         
  #   #         gsea_plot <-
  #   #           clusterProfiler::dotplot(
  #   #             GSEA_result,
  #   #             showCategory =
  #   #               top_n_GSEA
  #   #           ) +
  #   #           
  #   #           ggtitle(
  #   #             paste0(
  #   #               dataset_name,
  #   #               " : GSEA rescue score"
  #   #             )
  #   #           )
  #   #         
  #   #         
  #   #         gsea_pdf <- file.path(
  #   #           out_dir,
  #   #           paste0(
  #   #             dataset_name,
  #   #             "_GSEA_rescue_GO.pdf"
  #   #           )
  #   #         )
  #   #         
  #   #         
  #   #         pdf(
  #   #           gsea_pdf,
  #   #           width = 11,
  #   #           height = 8
  #   #         )
  #   #         
  #   #         print(
  #   #           gsea_plot
  #   #         )
  #   #         
  #   #         dev.off()
  #   #         
  #   #         
  #   #         message(
  #   #           "  GSEA terminée : ",
  #   #           nrow(GSEA_df),
  #   #           " termes."
  #   #         )
  #   #         
  #   #       } else {
  #   #         
  #   #         message(
  #   #           "  Aucun terme GSEA significatif."
  #   #         )
  #   #       }
  #   #     }
  #   #     
  #   #   } else {
  #   #     
  #   #     message(
  #   #       "  [SKIP] Pas assez de gènes pour GSEA."
  #   #     )
  #   #   }
  #   # }
  #   # 
  #   # 
  #   # 
  #   # 
  #   # 
  #   # 
  #   # 
  #   
  #   
  #   # ----------------------------------------------------------
  #   # Mapping vers ENTREZ
  #   # ----------------------------------------------------------
  #   
  #   if (
  #     exists("id_map") &&
  #     !is.null(id_map) &&
  #     exists("key_type") &&
  #     !is.null(key_type)
  #   ) {
  #     
  #     rank_map <- id_map %>%
  #       dplyr::filter(
  #         .data[[key_type]] %in% names(gene_rank)
  #       ) %>%
  #       dplyr::select(
  #         gene_id = .data[[key_type]],
  #         ENTREZID
  #       ) %>%
  #       dplyr::filter(
  #         !is.na(ENTREZID)
  #       ) %>%
  #       dplyr::distinct(
  #         gene_id,
  #         ENTREZID
  #       )
  #     
  #     # --------------------------------------------------------
  #     # Construire le tableau de ranking
  #     # --------------------------------------------------------
  #     
  #     rank_df <- data.frame(
  #       gene_id = names(gene_rank),
  #       score = as.numeric(gene_rank),
  #       stringsAsFactors = FALSE
  #     ) %>%
  #       dplyr::inner_join(
  #         rank_map,
  #         by = "gene_id"
  #       )
  #     
  #     # --------------------------------------------------------
  #     # Un seul score par ENTREZID
  #     #
  #     # Si plusieurs IDs correspondent au même ENTREZID,
  #     # on garde celui ayant le score absolu maximal.
  #     # --------------------------------------------------------
  #     
  #     rank_df <- rank_df %>%
  #       dplyr::filter(
  #         !is.na(score),
  #         is.finite(score)
  #       ) %>%
  #       dplyr::group_by(
  #         ENTREZID
  #       ) %>%
  #       dplyr::slice_max(
  #         order_by = abs(score),
  #         n = 1,
  #         with_ties = FALSE
  #       ) %>%
  #       dplyr::ungroup()
  #     
  #     # --------------------------------------------------------
  #     # Vérifications avant GSEA
  #     # --------------------------------------------------------
  #     
  #     message(
  #       "  GSEA : ",
  #       nrow(rank_df),
  #       " gènes avec ENTREZID."
  #     )
  #     
  #     if (nrow(rank_df) < 10) {
  #       
  #       message(
  #         "  [SKIP] Pas assez de gènes mappés pour GSEA."
  #       )
  #       
  #     } else {
  #       
  #       geneList <- rank_df$score
  #       
  #       names(geneList) <-
  #         as.character(rank_df$ENTREZID)
  #       
  #       geneList <- geneList[
  #         !is.na(geneList) &
  #           is.finite(geneList)
  #       ]
  #       
  #       geneList <- geneList[
  #         !duplicated(names(geneList))
  #       ]
  #       
  #       geneList <- sort(
  #         geneList,
  #         decreasing = TRUE
  #       )
  #       
  #       # ------------------------------------------------------
  #       # Vérification supplémentaire
  #       # ------------------------------------------------------
  #       
  #       if (
  #         length(geneList) < 10 ||
  #         length(unique(names(geneList))) < 10
  #       ) {
  #         
  #         message(
  #           "  [SKIP] Ranking GSEA insuffisant."
  #         )
  #         
  #       } else {
  #         
  #         message(
  #           "  Lancement GSEA sur ",
  #           length(geneList),
  #           " gènes."
  #         )
  #         
  #         # ----------------------------------------------------
  #         # GSEA
  #         # ----------------------------------------------------
  #         
  #         GSEA_result <- tryCatch(
  #           
  #           {
  #             
  #             suppressWarnings(
  #               
  #               clusterProfiler::gseGO(
  #                 
  #                 geneList =
  #                   geneList,
  #                 
  #                 OrgDb =
  #                   orgdb,
  #                 
  #                 keyType =
  #                   "ENTREZID",
  #                 
  #                 ont =
  #                   ont,
  #                 
  #                 minGSSize =
  #                   minGSSize,
  #                 
  #                 maxGSSize =
  #                   maxGSSize,
  #                 
  #                 pvalueCutoff =
  #                   pvalueCutoff,
  #                 
  #                 pAdjustMethod =
  #                   "BH",
  #                 
  #                 verbose =
  #                   FALSE
  #               )
  #             )
  #           },
  #           
  #           error = function(e) {
  #             
  #             message(
  #               "  [ERROR] GSEA : ",
  #               conditionMessage(e)
  #             )
  #             
  #             NULL
  #           }
  #         )
  #         
  #         
  #         # ----------------------------------------------------
  #         # Résultats GSEA
  #         # ----------------------------------------------------
  #         
  #         if (!is.null(GSEA_result)) {
  #           
  #           GSEA_df <-
  #             as.data.frame(
  #               GSEA_result
  #             )
  #           
  #           if (nrow(GSEA_df) > 0) {
  #             
  #             # ------------------------------------------------
  #             # Export CSV
  #             # ------------------------------------------------
  #             
  #             gsea_csv <- file.path(
  #               out_dir,
  #               paste0(
  #                 dataset_name,
  #                 "_GSEA_rescue_GO.csv"
  #               )
  #             )
  #             
  #             write.csv(
  #               GSEA_df,
  #               gsea_csv,
  #               row.names = FALSE
  #             )
  #             
  #             
  #             # ------------------------------------------------
  #             # Dotplot
  #             # ------------------------------------------------
  #             
  #             gsea_plot <-
  #               clusterProfiler::dotplot(
  #                 GSEA_result,
  #                 showCategory = top_n_GSEA
  #               ) +
  #               ggtitle(
  #                 paste0(
  #                   dataset_name,
  #                   " : GSEA rescue score"
  #                 )
  #               )
  #             
  #             
  #             gsea_pdf <- file.path(
  #               out_dir,
  #               paste0(
  #                 dataset_name,
  #                 "_GSEA_rescue_GO.pdf"
  #               )
  #             )
  #             
  #             
  #             pdf(
  #               gsea_pdf,
  #               width = 11,
  #               height = 8
  #             )
  #             
  #             print(
  #               gsea_plot
  #             )
  #             
  #             dev.off()
  #             
  #             
  #             message(
  #               "  GSEA terminée : ",
  #               nrow(GSEA_df),
  #               " termes."
  #             )
  #             
  #           } else {
  #             
  #             message(
  #               "  Aucun terme GSEA significatif."
  #             )
  #           }
  #           
  #         } else {
  #           
  #           message(
  #             "  GSEA non disponible pour ",
  #             dataset_name,
  #             "."
  #           )
  #         }
  #       }
  #     }
  #     
  #   } else {
  #     
  #     message(
  #       "  [SKIP] Mapping ENTREZ non disponible pour GSEA."
  #     )
  #   }
  #   
  # }
  # 
  
  
  # ============================================================
  # ============================================================
  # GSEA MULTI-BDD
  # ============================================================
  # ============================================================
  
  GSEA_result <- list()
  
  if (run_GSEA) {
    
    message(
      "\n  ========================================\n",
      "  GSEA MULTI-PATHWAY : ", dataset_name,
      "\n  ========================================"
    )
    
    # ------------------------------------------------------------
    # Vérification des packages
    # ------------------------------------------------------------
    
    required_gsea <- c(
      "clusterProfiler",
      "AnnotationDbi"
    )
    
    missing_gsea <- required_gsea[
      !sapply(
        required_gsea,
        requireNamespace,
        quietly = TRUE
      )
    ]
    
    if (length(missing_gsea) > 0) {
      
      message(
        "  [SKIP] Packages manquants : ",
        paste(missing_gsea, collapse = ", ")
      )
      
    } else {
      
      # ==========================================================
      # 1. Charger OrgDb
      # ==========================================================
      
      orgdb <- NULL
      
      tryCatch({
        
        if (!requireNamespace(
          organism_db,
          quietly = TRUE
        )) {
          
          message(
            "  [WARNING] Package ",
            organism_db,
            " non installé."
          )
          
        } else {
          
          suppressPackageStartupMessages(
            library(
              organism_db,
              character.only = TRUE
            )
          )
          
          orgdb <- get(
            organism_db
          )
        }
        
      }, error = function(e) {
        
        message(
          "  [WARNING] Impossible de charger ",
          organism_db,
          " : ",
          e$message
        )
      })
      
      
      # ==========================================================
      # 2. Déterminer le type d'identifiant
      # ==========================================================
      
      key_type <- NULL
      
      if (!is.null(orgdb)) {
        
        available_keys <-
          AnnotationDbi::keytypes(
            orgdb
          )
        
        test_ids <- unique(
          as.character(
            df$GENEID
          )
        )
        
        test_ids <- test_ids[
          !is.na(test_ids) &
            test_ids != ""
        ]
        
        test_ids <- head(
          test_ids,
          500
        )
        
        for (candidate in c(
          "ENSEMBL",
          "SYMBOL",
          "ENTREZID"
        )) {
          
          if (
            candidate %in%
            available_keys
          ) {
            
            tmp <- tryCatch(
              
              AnnotationDbi::select(
                orgdb,
                keys = test_ids,
                columns = candidate,
                keytype = candidate
              ),
              
              error = function(e)
                NULL
            )
            
            if (
              !is.null(tmp) &&
              nrow(tmp) > 0 &&
              candidate %in%
              colnames(tmp)
            ) {
              
              n_valid <- sum(
                !is.na(
                  tmp[[candidate]]
                )
              )
              
              if (n_valid > 0) {
                
                key_type <- candidate
                
                break
              }
            }
          }
        }
      }
      
      
      # ==========================================================
      # 3. Mapping vers ENTREZ
      # ==========================================================
      
      id_map <- NULL
      
      if (
        !is.null(orgdb) &&
        !is.null(key_type)
      ) {
        
        message(
          "  Identifiant détecté : ",
          key_type
        )
        
        id_map <- tryCatch(
          
          AnnotationDbi::select(
            orgdb,
            keys =
              unique(
                as.character(
                  df$GENEID
                )
              ),
            columns = c(
              key_type,
              "ENTREZID"
            ),
            keytype = key_type
          ),
          
          error = function(e) {
            
            message(
              "  [WARNING] Mapping impossible : ",
              e$message
            )
            
            NULL
          }
        )
        
        
        if (!is.null(id_map)) {
          
          # ------------------------------------------------------
          # IMPORTANT :
          # on standardise les noms de colonnes ici
          # ------------------------------------------------------
          
          id_map <- id_map %>%
            dplyr::rename(
              gene_id = dplyr::all_of(
                key_type
              )
            ) %>%
            dplyr::filter(
              !is.na(gene_id),
              !is.na(ENTREZID),
              ENTREZID != ""
            ) %>%
            dplyr::distinct(
              gene_id,
              ENTREZID
            )
          
          message(
            "  Mapping : ",
            nrow(id_map),
            " correspondances gene → ENTREZID"
          )
        }
      }
      
      
      # ==========================================================
      # 4. Fonction pour construire un ranking
      # ==========================================================
      
      make_gene_rank <- function(
    score_column
      ) {
        
        tmp <- df %>%
          dplyr::select(
            GENEID,
            dplyr::all_of(
              score_column
            )
          ) %>%
          dplyr::filter(
            !is.na(GENEID),
            !is.na(
              .data[[score_column]]
            )
          ) %>%
          dplyr::rename(
            score = dplyr::all_of(
              score_column
            )
          ) %>%
          dplyr::group_by(
            GENEID
          ) %>%
          dplyr::slice_max(
            order_by = abs(score),
            n = 1,
            with_ties = FALSE
          ) %>%
          dplyr::ungroup()
        
        
        if (
          is.null(id_map) ||
          nrow(id_map) == 0
        ) {
          
          return(NULL)
        }
        
        
        # --------------------------------------------------------
        # Mapping
        # --------------------------------------------------------
        
        tmp <- tmp %>%
          dplyr::inner_join(
            id_map,
            by = c(
              "GENEID" =
                "gene_id"
            )
          )
        
        
        # --------------------------------------------------------
        # Un ENTREZID = une seule valeur
        # garder le score le plus extrême
        # --------------------------------------------------------
        
        tmp <- tmp %>%
          dplyr::group_by(
            ENTREZID
          ) %>%
          dplyr::slice_max(
            order_by = abs(score),
            n = 1,
            with_ties = FALSE
          ) %>%
          dplyr::ungroup()
        
        
        geneList <- tmp$score
        
        names(geneList) <-
          as.character(
            tmp$ENTREZID
          )
        
        
        geneList <- geneList[
          !is.na(geneList)
        ]
        
        
        geneList <- sort(
          geneList,
          decreasing = TRUE
        )
        
        
        return(
          geneList
        )
      }
      
      
      # ==========================================================
      # 5. Fonction GSEA générique
      # ==========================================================
      
      run_one_GSEA <- function(
    geneList,
    db_name,
    gsea_name,
    TERM2GENE = NULL,
    TERM2NAME = NULL
      ) {
        
        if (
          is.null(geneList) ||
          length(geneList) < 10
        ) {
          
          message(
            "    [SKIP] ",
            gsea_name,
            " : pas assez de gènes."
          )
          
          return(NULL)
        }
        
        
        message(
          "    ",
          db_name,
          " GSEA : ",
          gsea_name,
          " (",
          length(geneList),
          " genes)"
        )
        
        
        res <- tryCatch({
          
          if (
            db_name == "GO"
          ) {
            
            clusterProfiler::gseGO(
              
              geneList =
                geneList,
              
              OrgDb =
                orgdb,
              
              keyType =
                "ENTREZID",
              
              ont =
                ont,
              
              minGSSize =
                minGSSize,
              
              maxGSSize =
                maxGSSize,
              
              pvalueCutoff =
                1,
              
              pAdjustMethod =
                "BH",
              
              verbose =
                FALSE
            )
            
          } else {
            
            clusterProfiler::GSEA(
              
              geneList =
                geneList,
              
              TERM2GENE =
                TERM2GENE,
              
              TERM2NAME =
                TERM2NAME,
              
              minGSSize =
                minGSSize,
              
              maxGSSize =
                maxGSSize,
              
              pvalueCutoff =
                1,
              
              pAdjustMethod =
                "BH",
              
              verbose =
                FALSE
            )
          }
          
        }, error = function(e) {
          
          message(
            "    [ERROR] ",
            db_name,
            " / ",
            gsea_name,
            " : ",
            e$message
          )
          
          NULL
        })
        
        
        if (is.null(res)) {
          
          return(NULL)
        }
        
        
        res_df <- as.data.frame(
          res
        )
        
        
        if (
          nrow(res_df) == 0
        ) {
          
          message(
            "    Aucun pathway trouvé."
          )
          
          return(NULL)
        }
        
        
        # --------------------------------------------------------
        # Ajouter une colonne pour identifier l'analyse
        # --------------------------------------------------------
        
        res_df$GSEA_type <-
          gsea_name
        
        res_df$Database <-
          db_name
        
        
        # --------------------------------------------------------
        # Export CSV
        # --------------------------------------------------------
        
        safe_name <- gsub(
          "[^A-Za-z0-9_]+",
          "_",
          gsea_name
        )
        
        
        csv_file <- file.path(
          out_dir,
          paste0(
            dataset_name,
            "_GSEA_",
            safe_name,
            ".csv"
          )
        )
        
        
        write.csv(
          res_df,
          csv_file,
          row.names = FALSE
        )
        
        
        # --------------------------------------------------------
        # Garder les pathways significatifs
        # --------------------------------------------------------
        
        sig_df <- res_df %>%
          dplyr::filter(
            !is.na(p.adjust),
            p.adjust < 0.05
          )
        
        
        message(
          "    ",
          gsea_name,
          " : ",
          nrow(res_df),
          " pathways testés | ",
          nrow(sig_df),
          " significatifs"
        )
        
        
        # --------------------------------------------------------
        # Plot
        # --------------------------------------------------------
        
        if (
          nrow(sig_df) > 0
        ) {
          
          plot_res <- res
          
          plot_file <- file.path(
            out_dir,
            paste0(
              dataset_name,
              "_GSEA_",
              safe_name,
              ".pdf"
            )
          )
          
          
          pdf(
            plot_file,
            width = 11,
            height = 8
          )
          
          
          print(
            clusterProfiler::dotplot(
              plot_res,
              showCategory =
                min(
                  top_n_GSEA,
                  nrow(sig_df)
                )
            ) +
              ggtitle(
                paste0(
                  dataset_name,
                  " - ",
                  db_name,
                  " GSEA - ",
                  gsea_name
                )
              )
          )
          
          
          dev.off()
        }
        
        
        return(
          res
        )
      }
      
      
      # ==========================================================
      # 6. Construire les 3 rankings
      # ==========================================================
      
      message(
        "\n  Construction des rankings..."
      )
      
      
      # ----------------------------------------------------------
      # RESCUE
      #
      # positif = diminution de magnitude
      # ----------------------------------------------------------
      
      gene_rank_rescue <-
        make_gene_rank(
          "delta_abs_LFC"
        )
      
      
      # ----------------------------------------------------------
      # DIRECTIONAL
      #
      # FA - EDITED
      #
      # positif = déplacement vers 0 / direction rescue
      # ----------------------------------------------------------
      
      df$directional_score <-
        df$lfc_FA -
        df$lfc_EDITED
      
      
      gene_rank_directional <-
        make_gene_rank(
          "directional_score"
        )
      
      
      # ----------------------------------------------------------
      # FA REFERENCE
      #
      # ranking classique de la maladie
      # ----------------------------------------------------------
      
      gene_rank_FA <-
        make_gene_rank(
          "lfc_FA"
        )
      
      
      # ==========================================================
      # 7. GO
      # ==========================================================
      
      if (
        !is.null(orgdb) &&
        !is.null(gene_rank_rescue)
      ) {
        
        GSEA_result$GO_Rescue <-
          run_one_GSEA(
            gene_rank_rescue,
            "GO",
            "Rescue"
          )
        
        
        GSEA_result$GO_Directional <-
          run_one_GSEA(
            gene_rank_directional,
            "GO",
            "Directional"
          )
        
        
        GSEA_result$GO_FA_reference <-
          run_one_GSEA(
            gene_rank_FA,
            "GO",
            "FA reference"
          )
      }
      
      
      # ==========================================================
      # 8. KEGG
      # ==========================================================
      
      if (
        !is.null(gene_rank_rescue)
      ) {
        
        GSEA_result$KEGG_Rescue <-
          run_one_GSEA(
            gene_rank_rescue,
            "KEGG",
            "Rescue",
            TERM2GENE = NULL
          )
        
        # --------------------------------------------------------
        # KEGG nécessite une fonction dédiée
        # --------------------------------------------------------
        
        run_KEGG <- function(
    geneList,
    name
        ) {
          
          if (
            is.null(geneList) ||
            length(geneList) < 10
          ) {
            
            return(NULL)
          }
          
          
          message(
            "    KEGG GSEA : ",
            name
          )
          
          
          res <- tryCatch(
            
            clusterProfiler::gseKEGG(
              
              geneList =
                geneList,
              
              organism =
                "hsa",
              
              minGSSize =
                minGSSize,
              
              maxGSSize =
                maxGSSize,
              
              pvalueCutoff =
                1,
              
              pAdjustMethod =
                "BH",
              
              verbose =
                FALSE
            ),
            
            error = function(e) {
              
              message(
                "    [ERROR] KEGG : ",
                e$message
              )
              
              NULL
            }
          )
          
          
          if (
            is.null(res)
          ) {
            
            return(NULL)
          }
          
          
          res_df <-
            as.data.frame(res)
          
          
          res_df$GSEA_type <-
            name
          
          res_df$Database <-
            "KEGG"
          
          
          safe_name <- gsub(
            "[^A-Za-z0-9_]+",
            "_",
            name
          )
          
          
          write.csv(
            res_df,
            file.path(
              out_dir,
              paste0(
                dataset_name,
                "_GSEA_KEGG_",
                safe_name,
                ".csv"
              )
            ),
            row.names = FALSE
          )
          
          
          sig_df <- res_df %>%
            dplyr::filter(
              !is.na(p.adjust),
              p.adjust < 0.05
            )
          
          
          message(
            "    KEGG ",
            name,
            " : ",
            nrow(sig_df),
            " pathways significatifs"
          )
          
          
          if (
            nrow(sig_df) > 0
          ) {
            
            pdf(
              file.path(
                out_dir,
                paste0(
                  dataset_name,
                  "_GSEA_KEGG_",
                  safe_name,
                  ".pdf"
                )
              ),
              width = 11,
              height = 8
            )
            
            
            print(
              clusterProfiler::dotplot(
                res,
                showCategory =
                  min(
                    top_n_GSEA,
                    nrow(sig_df)
                  )
              ) +
                ggtitle(
                  paste0(
                    dataset_name,
                    " - KEGG GSEA - ",
                    name
                  )
                )
            )
            
            
            dev.off()
          }
          
          
          return(res)
        }
        
        
        GSEA_result$KEGG_Rescue <-
          run_KEGG(
            gene_rank_rescue,
            "Rescue"
          )
        
        
        GSEA_result$KEGG_Directional <-
          run_KEGG(
            gene_rank_directional,
            "Directional"
          )
        
        
        GSEA_result$KEGG_FA_reference <-
          run_KEGG(
            gene_rank_FA,
            "FA reference"
          )
      }
      
      
      # ==========================================================
      # 9. REACTOME
      # ==========================================================
      
      if (
        requireNamespace(
          "ReactomePA",
          quietly = TRUE
        )
      ) {
        
        run_Reactome <- function(
    geneList,
    name
        ) {
          
          if (
            is.null(geneList) ||
            length(geneList) < 10
          ) {
            
            return(NULL)
          }
          
          
          message(
            "    Reactome GSEA : ",
            name
          )
          
          
          res <- tryCatch(
            
            ReactomePA::gsePathway(
              
              geneList =
                geneList,
              
              organism =
                "human",
              
              minGSSize =
                minGSSize,
              
              maxGSSize =
                maxGSSize,
              
              pvalueCutoff =
                1,
              
              pAdjustMethod =
                "BH",
              
              verbose =
                FALSE
            ),
            
            error = function(e) {
              
              message(
                "    [ERROR] Reactome : ",
                e$message
              )
              
              NULL
            }
          )
          
          
          if (
            is.null(res)
          ) {
            
            return(NULL)
          }
          
          
          res_df <-
            as.data.frame(res)
          
          
          res_df$GSEA_type <-
            name
          
          res_df$Database <-
            "Reactome"
          
          
          safe_name <- gsub(
            "[^A-Za-z0-9_]+",
            "_",
            name
          )
          
          
          write.csv(
            res_df,
            file.path(
              out_dir,
              paste0(
                dataset_name,
                "_GSEA_Reactome_",
                safe_name,
                ".csv"
              )
            ),
            row.names = FALSE
          )
          
          
          sig_df <- res_df %>%
            dplyr::filter(
              !is.na(p.adjust),
              p.adjust < 0.05
            )
          
          
          message(
            "    Reactome ",
            name,
            " : ",
            nrow(sig_df),
            " pathways significatifs"
          )
          
          
          if (
            nrow(sig_df) > 0
          ) {
            
            pdf(
              file.path(
                out_dir,
                paste0(
                  dataset_name,
                  "_GSEA_Reactome_",
                  safe_name,
                  ".pdf"
                )
              ),
              width = 11,
              height = 8
            )
            
            
            print(
              clusterProfiler::dotplot(
                res,
                showCategory =
                  min(
                    top_n_GSEA,
                    nrow(sig_df)
                  )
              ) +
                ggtitle(
                  paste0(
                    dataset_name,
                    " - Reactome GSEA - ",
                    name
                  )
                )
            )
            
            
            dev.off()
          }
          
          
          return(res)
        }
        
        
        GSEA_result$Reactome_Rescue <-
          run_Reactome(
            gene_rank_rescue,
            "Rescue"
          )
        
        
        GSEA_result$Reactome_Directional <-
          run_Reactome(
            gene_rank_directional,
            "Directional"
          )
        
        
        GSEA_result$Reactome_FA_reference <-
          run_Reactome(
            gene_rank_FA,
            "FA reference"
          )
        
      } else {
        
        message(
          "  [INFO] ReactomePA non installé."
        )
      }
      
      
      # ==========================================================
      # 10. HALLMARK
      # ==========================================================
      
      if (
        requireNamespace(
          "msigdbr",
          quietly = TRUE
        )
      ) {
        
        message(
          "\n  Chargement MSigDB Hallmark..."
        )
        
        
        hallmark <- tryCatch(
          
          msigdbr::msigdbr(
            species = "Homo sapiens",
            category = "H"
          ),
          
          error = function(e) {
            
            message(
              "  [ERROR] msigdbr : ",
              e$message
            )
            
            NULL
          }
        )
        
        
        if (
          !is.null(hallmark)
        ) {
          
          hallmark_term_col <-
            if (
              "gs_name" %in%
              colnames(hallmark)
            ) {
              "gs_name"
            } else {
              "gs_name"
            }
          
          
          hallmark_gene_col <-
            if (
              "ncbi_gene" %in%
              colnames(hallmark)
            ) {
              "ncbi_gene"
            } else if (
              "entrez_gene" %in%
              colnames(hallmark)
            ) {
              "entrez_gene"
            } else {
              NULL
            }
          
          
          if (
            is.null(hallmark_gene_col)
          ) {
            
            message(
              "  [WARNING] Colonne ENTREZ absente de msigdbr."
            )
            
          } else {
            
            hallmark_t2g <-
              hallmark %>%
              dplyr::select(
                TERM =
                  dplyr::all_of(
                    hallmark_term_col
                  ),
                ENTREZID =
                  dplyr::all_of(
                    hallmark_gene_col
                  )
              ) %>%
              dplyr::filter(
                !is.na(ENTREZID)
              ) %>%
              dplyr::distinct()
            
            
            GSEA_result$Hallmark_Rescue <-
              run_one_GSEA(
                gene_rank_rescue,
                "Hallmark",
                "Rescue",
                TERM2GENE =
                  hallmark_t2g
              )
            
            
            GSEA_result$Hallmark_Directional <-
              run_one_GSEA(
                gene_rank_directional,
                "Hallmark",
                "Directional",
                TERM2GENE =
                  hallmark_t2g
              )
            
            
            GSEA_result$Hallmark_FA_reference <-
              run_one_GSEA(
                gene_rank_FA,
                "Hallmark",
                "FA reference",
                TERM2GENE =
                  hallmark_t2g
              )
          }
        }
        
      } else {
        
        message(
          "  [INFO] msigdbr non installé."
        )
      }
    }
  }
  
  # ============================================================
  # FRIEDREICH ATAXIA FOCUSED ANALYSIS
  # ============================================================
  
  FA_focused_results <- list()
  
  
  if (
    run_GSEA &&
    run_FA_focused &&
    length(GSEA_result) > 0
  ) {
    
    message(
      "\n  ========================================\n",
      "  FRIEDREICH ATAXIA FOCUSED PATHWAYS : ",
      dataset_name,
      "\n  ========================================"
    )
    
    
    # ----------------------------------------------------------
    # Construire regex robuste
    # ----------------------------------------------------------
    
    FA_regex <- paste(
      FA_keywords,
      collapse = "|"
    )
    
    
    # ----------------------------------------------------------
    # Parcourir tous les GSEA
    # ----------------------------------------------------------
    
    for (
      analysis_name in names(GSEA_result)
    ) {
      
      obj <- GSEA_result[[analysis_name]]
      
      
      if (is.null(obj)) {
        next
      }
      
      
      obj_df <- tryCatch(
        
        as.data.frame(obj),
        
        error = function(e) {
          
          message(
            "  [WARNING] Impossible de convertir ",
            analysis_name,
            " : ",
            e$message
          )
          
          NULL
        }
      )
      
      
      if (
        is.null(obj_df) ||
        nrow(obj_df) == 0
      ) {
        next
      }
      
      
      # --------------------------------------------------------
      # Vérifier la colonne Description
      # --------------------------------------------------------
      
      if (
        !"Description" %in%
        colnames(obj_df)
      ) {
        
        message(
          "  [WARNING] Pas de colonne Description pour ",
          analysis_name
        )
        
        next
      }
      
      
      # --------------------------------------------------------
      # Identifier les pathways FA-relevant
      # --------------------------------------------------------
      
      
      obj_df <- obj_df %>%
        dplyr::mutate(
          
          FA_relevant =
            grepl(
              FA_regex,
              Description,
              ignore.case = TRUE
            ),
          
          FA_keyword_match =
            ifelse(
              FA_relevant,
              regmatches(
                Description,
                regexpr(
                  FA_regex,
                  Description,
                  ignore.case = TRUE
                )
              ),
              NA_character_
            )
        )
      
      # --------------------------------------------------------
      # Pathways significatifs
      # --------------------------------------------------------
      
      sig_df <- obj_df %>%
        dplyr::filter(
          !is.na(p.adjust),
          p.adjust < FA_padj_cutoff
        )
      
      
      # --------------------------------------------------------
      # Pathways FA-relevant ET significatifs
      # --------------------------------------------------------
      
      fa_df <- sig_df %>%
        dplyr::filter(
          FA_relevant
        ) %>%
        dplyr::arrange(
          p.adjust
        )
      
      
      # --------------------------------------------------------
      # Export tous les pathways FA-relevant
      # --------------------------------------------------------
      
      FA_focused_results[[
        analysis_name
      ]] <- fa_df
      
      
      safe_name <- gsub(
        "[^A-Za-z0-9_]+",
        "_",
        analysis_name
      )
      
      
      if (
        nrow(fa_df) > 0
      ) {
        
        fa_csv <- file.path(
          out_dir,
          paste0(
            dataset_name,
            "_FA_focused_",
            safe_name,
            ".csv"
          )
        )
        
        
        write.csv(
          fa_df,
          fa_csv,
          row.names = FALSE
        )
        
        
        message(
          "  ",
          analysis_name,
          " : ",
          nrow(fa_df),
          " pathways FA-relevant parmi ",
          nrow(sig_df),
          " significatifs"
        )
        
        
        # ------------------------------------------------------
        # Plot FA-focused
        # ------------------------------------------------------
        
        plot_df <- fa_df %>%
          dplyr::slice_head(
            n = top_n_GSEA
          )
        
        
        if (
          nrow(plot_df) > 0
        ) {
          
          # Gestion robuste de NES
          if (
            "NES" %in% colnames(plot_df)
          ) {
            
            # ------------------------------------------------------------
            # Taille du pathway
            #
            # GSEA peut utiliser :
            # - Count
            # - setSize
            # - size
            #
            # On prend la première disponible.
            # ------------------------------------------------------------
            
            size_col <- NULL
            
            if ("Count" %in% colnames(plot_df)) {
              size_col <- "Count"
              
            } else if ("setSize" %in% colnames(plot_df)) {
              size_col <- "setSize"
              
            } else if ("size" %in% colnames(plot_df)) {
              size_col <- "size"
            }
            
            
            # ------------------------------------------------------------
            # Si aucune colonne de taille n'existe,
            # on utilise une taille constante.
            # ------------------------------------------------------------
            
            if (!is.null(size_col)) {
              
              plot_df$PathwaySize <-
                as.numeric(
                  plot_df[[size_col]]
                )
              
              size_label <- "Pathway size"
              
            } else {
              
              plot_df$PathwaySize <- 5
              
              size_label <- NULL
            }
            
            
            fa_plot <- ggplot(
              plot_df,
              aes(
                x = NES,
                y = reorder(
                  Description,
                  NES
                ),
                size = PathwaySize
              )
            ) +
              
              geom_point() +
              
              geom_vline(
                xintercept = 0,
                linetype = "dashed"
              ) +
              
              theme_bw(
                base_size = 11
              ) +
              
              labs(
                title = paste0(
                  dataset_name,
                  " - FA-focused pathways - ",
                  analysis_name
                ),
                subtitle =
                  "Significant pathways relevant to Friedreich ataxia biology",
                x = "Normalized Enrichment Score (NES)",
                y = NULL,
                size = size_label
              )
            
            
            fa_pdf <- file.path(
              out_dir,
              paste0(
                dataset_name,
                "_FA_focused_",
                safe_name,
                ".pdf"
              )
            )
            
            
            pdf(
              fa_pdf,
              width = 12,
              height = 8
            )
            
            print(
              fa_plot
            )
            
            dev.off()
          }
        }
        
      } else {
        
        message(
          "  ",
          analysis_name,
          " : aucun pathway FA-relevant significatif."
        )
      }
    }
    
    
    # ==========================================================
    # Tableau global FA-focused
    # ==========================================================
    
    non_empty <- FA_focused_results[
      sapply(
        FA_focused_results,
        function(x)
          !is.null(x) &&
          is.data.frame(x) &&
          nrow(x) > 0
      )
    ]
    
    
    if (
      length(non_empty) > 0
    ) {
      
      FA_focused_table <-
        dplyr::bind_rows(
          non_empty,
          .id = "GSEA_analysis"
        )
      
      
      # --------------------------------------------------------
      # Export global
      # --------------------------------------------------------
      
      FA_global_csv <- file.path(
        out_dir,
        paste0(
          dataset_name,
          "_FA_focused_ALL_GSEA.csv"
        )
      )
      
      
      write.csv(
        FA_focused_table,
        FA_global_csv,
        row.names = FALSE
      )
      
      
      message(
        "\n  FA-focused global : ",
        nrow(FA_focused_table),
        " pathways."
      )
      
      
      # --------------------------------------------------------
      # Excel
      # --------------------------------------------------------
      
      if (export_xlsx) {
        
        # Une feuille séparée pour éviter
        # les problèmes de tailles différentes
        
        FA_excel_sheets <- list(
          FA_focused_ALL =
            FA_focused_table
        )
        
        
        for (
          nm in names(non_empty)
        ) {
          
          sheet_name <- paste0(
            "FA_",
            nm
          )
          
          sheet_name <- gsub(
            "[^A-Za-z0-9_]",
            "_",
            sheet_name
          )
          
          sheet_name <- substr(
            sheet_name,
            1,
            31
          )
          
          FA_excel_sheets[[
            sheet_name
          ]] <-
            non_empty[[nm]]
        }
        
        
        FA_excel_file <- file.path(
          out_dir,
          paste0(
            dataset_name,
            "_FA_focused_GSEA.xlsx"
          )
        )
        
        
        writexl::write_xlsx(
          FA_excel_sheets,
          FA_excel_file
        )
        
        
        message(
          "  Excel FA-focused sauvegardé : ",
          FA_excel_file
        )
      }
      
    } else {
      
      message(
        "\n  Aucun pathway FA-relevant significatif ",
        "dans les GSEA."
      )
    }
  }
  # ============================================================
  # Export Excel
  # ============================================================
  
  if (export_xlsx) {
    
    xlsx_file <- file.path(
      out_dir,
      paste0(
        "genes_",
        dataset_name,
        "_FA_vs_FAedited_rescue.xlsx"
      )
    )
    
    
    # ----------------------------------------------------------
    # Sheet principal
    # ----------------------------------------------------------
    
    sheets <- list(
      rescue_genes = df
    )
    
    
    # ----------------------------------------------------------
    # GO
    # ----------------------------------------------------------
    
    if (
      length(GO_results) > 0
    ) {
      
      for (nm in names(GO_results)) {
        
        obj <- GO_results[[nm]]
        
        if (
          !is.null(obj) &&
          is.data.frame(obj) &&
          nrow(obj) > 0
        ) {
          
          sheet_name <- paste0(
            "GO_",
            nm
          )
          
          # Excel limite les noms de sheets à 31 caractères
          sheet_name <- substr(
            sheet_name,
            1,
            31
          )
          
          sheets[[sheet_name]] <- obj
        }
      }
    }
    
    
    # ----------------------------------------------------------
    # GSEA
    #
    # IMPORTANT :
    # chaque GSEA est exportée dans une feuille séparée
    # car les résultats ont des nombres de lignes différents.
    # ----------------------------------------------------------
    
    if (
      length(GSEA_result) > 0
    ) {
      
      for (nm in names(GSEA_result)) {
        
        obj <- GSEA_result[[nm]]
        
        if (
          !is.null(obj)
        ) {
          
          obj_df <- tryCatch(
            
            as.data.frame(obj),
            
            error = function(e) {
              
              message(
                "  [WARNING] Impossible d'exporter ",
                nm,
                " : ",
                e$message
              )
              
              NULL
            }
          )
          
          
          if (
            !is.null(obj_df) &&
            is.data.frame(obj_df) &&
            nrow(obj_df) > 0
          ) {
            
            sheet_name <- paste0(
              "GSEA_",
              nm
            )
            
            sheet_name <- substr(
              sheet_name,
              1,
              31
            )
            
            sheets[[sheet_name]] <-
              obj_df
          }
        }
      }
    }
    
    
    # ----------------------------------------------------------
    # Résumé GSEA
    # ----------------------------------------------------------
    
    gsea_summary <- data.frame()
    
    
    if (
      length(GSEA_result) > 0
    ) {
      
      summary_list <- list()
      
      
      for (nm in names(GSEA_result)) {
        
        obj <- GSEA_result[[nm]]
        
        
        if (
          !is.null(obj)
        ) {
          
          obj_df <- tryCatch(
            as.data.frame(obj),
            error = function(e) NULL
          )
          
          
          if (
            !is.null(obj_df) &&
            nrow(obj_df) > 0
          ) {
            
            n_sig <- sum(
              !is.na(obj_df$p.adjust) &
                obj_df$p.adjust < 0.05
            )
            
            
            summary_list[[nm]] <-
              data.frame(
                
                Analysis =
                  nm,
                
                N_pathways =
                  nrow(obj_df),
                
                N_significant =
                  n_sig,
                
                stringsAsFactors =
                  FALSE
              )
          }
        }
      }
      
      
      if (
        length(summary_list) > 0
      ) {
        
        gsea_summary <-
          dplyr::bind_rows(
            summary_list
          )
        
        
        sheets$GSEA_summary <-
          gsea_summary
      }
    }
    
    
    # ----------------------------------------------------------
    # Écriture Excel
    # ----------------------------------------------------------
    
    writexl::write_xlsx(
      sheets,
      xlsx_file
    )
    
    
    message(
      "  Excel sauvegardé : ",
      xlsx_file
    )
  }
  
  
  # ============================================================
  # Résumé des quadrants
  # ============================================================
  
  quadrant_summary <- df %>%
    dplyr::count(
      quadrant,
      rescue_status
    )
  
  
  write.csv(
    quadrant_summary,
    file.path(
      out_dir,
      paste0(
        dataset_name,
        "_quadrant_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ============================================================
  # Retour
  # ============================================================
  
  invisible(
    list(
      genes = df,
      GO = GO_results,
      GSEA = GSEA_result,
      quadrant_summary = quadrant_summary,
      scatter = p
    )
  )
}