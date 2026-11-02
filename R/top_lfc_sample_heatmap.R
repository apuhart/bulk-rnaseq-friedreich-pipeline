# ============================================================
# top_lfc_sample_heatmap.R
#
# Même design visuel que fa_pathway_sample_heatmap.R (ggplot2 : tuiles,
# z-score par gène, facettes par condition) mais SANS panel de gènes
# pré-défini : sélectionne les gènes SIGNIFICATIFS (padj < seuil) avec le
# |LFC| le plus fort pour CHAQUE comparaison, et écrit le résultat
# directement dans le dossier de cette comparaison
# (results_v4/<dataset>/<contraste>/), à côté des autres sorties du
# contraste (volcano, heatmap_topDEGs, etc.).
#
# Facette par direction (Up/Down) au lieu de pathway, puisqu'il n'y a pas
# de catégorie biologique ici -- juste "plus haut" ou "plus bas" dans la
# condition B par rapport à la condition A du nom du contraste
# (<condA>_vs_<condB>).
# ============================================================

# ------------------------------------------------------------
# Calcule le vst pour TOUT le dataset une seule fois (réutilisé pour tous
# les contrastes de ce dataset -- évite de refaire tximport/vst à chaque
# comparaison).
# ------------------------------------------------------------
.top_lfc_compute_dataset_vst <- function(
    dataset_name, tx2gene,
    min_counts  = get0("MIN_COUNTS",  ifnotfound = 10),
    min_samples = get0("MIN_SAMPLES", ifnotfound = 3)
) {
  
  dataset <- DATASETS[[dataset_name]]
  if (is.null(dataset)) {
    message("  [WARN] dataset '", dataset_name, "' introuvable dans DATASETS.")
    return(NULL)
  }
  
  samples <- dataset$samples
  files   <- file.path(dataset$path, samples, "quant.sf")
  names(files) <- samples
  if (!all(file.exists(files))) {
    message("  [WARN] quant.sf manquants pour ", dataset_name, " -> vst ignoré.")
    return(NULL)
  }
  
  txi <- tximport::tximport(
    files, type = "salmon",
    tx2gene = tx2gene[, c("TXNAME", "GENEID")],
    ignoreTxVersion = TRUE
  )
  
  colData <- data.frame(
    sample    = samples,
    condition = factor(dataset$conditions),
    row.names = samples
  )
  
  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = colData, design = ~ 1)
  dds <- filter_low_counts(dds, min_counts, min_samples)
  dds <- DESeq2::estimateSizeFactors(dds)
  vsd <- tryCatch(
    DESeq2::vst(dds, blind = TRUE),
    error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  )
  mat_full <- SummarizedExperiment::assay(vsd)
  rownames(mat_full) <- sub("\\..*", "", rownames(mat_full))
  
  list(mat = mat_full, colData = colData)
}

# ------------------------------------------------------------
# Une carte pour UN contraste (dataset + comparaison précise)
# ------------------------------------------------------------
plot_top_lfc_sample_heatmap <- function(
    dataset_name, contrast_name, res_df, out_dir,
    vst_data  = NULL,
    tx2gene   = NULL,
    top_n     = 30,
    padj_thr  = get0("PADJ_THRESHOLD", ifnotfound = 0.05),
    lfc_thr   = get0("LFC_THRESHOLD",  ifnotfound = 1),
    file_name = "heatmap_echantillons_topLFC.pdf"
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("  [INFO] ggplot2 indisponible : heatmap top-LFC ignoré.")
    return(invisible(NULL))
  }
  if (is.null(res_df) || nrow(res_df) == 0) {
    message("  [INFO] Heatmap top-LFC (", dataset_name, "/", contrast_name, ") : pas de résultats.")
    return(invisible(NULL))
  }
  
  if (is.null(vst_data)) {
    if (is.null(tx2gene)) {
      message("  [WARN] Heatmap top-LFC : ni vst_data ni tx2gene fournis -- ignoré.")
      return(invisible(NULL))
    }
    vst_data <- .top_lfc_compute_dataset_vst(dataset_name, tx2gene)
  }
  if (is.null(vst_data)) return(invisible(NULL))
  
  # ---- Sélection des gènes : significatifs, triés par |LFC| décroissant ----
  top_df <- res_df[!is.na(res_df$padj) & res_df$padj < padj_thr &
                     !is.na(res_df$log2FoldChange) & abs(res_df$log2FoldChange) >= lfc_thr, ]
  if (nrow(top_df) == 0) {
    message("  [INFO] Heatmap top-LFC (", dataset_name, "/", contrast_name,
            ") : aucun gène significatif -- ignoré.")
    return(invisible(NULL))
  }
  top_df <- top_df[order(-abs(top_df$log2FoldChange)), ]
  top_df <- top_df[seq_len(min(top_n, nrow(top_df))), ]
  top_df <- top_df[top_df$GENEID %in% rownames(vst_data$mat), ]
  if (nrow(top_df) < 2) {
    message("  [INFO] Heatmap top-LFC (", dataset_name, "/", contrast_name,
            ") : moins de 2 gènes exploitables -- ignoré.")
    return(invisible(NULL))
  }
  
  # ---- Échantillons de CE contraste uniquement (pas tout le dataset) ----
  cond_parts <- strsplit(contrast_name, "_vs_")[[1]]
  if (length(cond_parts) != 2) {
    message("  [WARN] Nom de contraste '", contrast_name, "' non reconnu -- ignoré.")
    return(invisible(NULL))
  }
  condA <- cond_parts[1]; condB <- cond_parts[2]
  colData <- vst_data$colData
  keep_samples <- colData$sample[colData$condition %in% c(condA, condB)]
  if (length(keep_samples) < 2) {
    message("  [WARN] Heatmap top-LFC (", dataset_name, "/", contrast_name,
            ") : moins de 2 échantillons trouvés pour ce contraste -- ignoré.")
    return(invisible(NULL))
  }
  
  mat <- vst_data$mat[top_df$GENEID, keep_samples, drop = FALSE]
  rownames(mat) <- top_df$GENENAME
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  
  direction <- setNames(ifelse(top_df$log2FoldChange > 0,
                               paste0("Plus haut en ", condB),
                               paste0("Plus bas en ", condB)),
                        top_df$GENENAME)
  lfc_val   <- setNames(round(top_df$log2FoldChange, 2), top_df$GENENAME)
  gene_order <- top_df$GENENAME[order(direction[top_df$GENENAME], -abs(top_df$log2FoldChange))]
  
  df_long <- as.data.frame(as.table(mat_z))
  colnames(df_long) <- c("gene", "sample", "zscore")
  df_long$gene   <- as.character(df_long$gene)
  df_long$sample <- as.character(df_long$sample)
  df_long$condition <- colData$condition[match(df_long$sample, colData$sample)]
  df_long$direction <- factor(direction[df_long$gene], levels = unique(direction[gene_order]))
  df_long$gene      <- factor(df_long$gene, levels = rev(gene_order))
  
  sample_order <- keep_samples[order(match(colData$condition[match(keep_samples, colData$sample)],
                                           c(condA, condB)), keep_samples)]
  df_long$sample <- factor(df_long$sample, levels = sample_order)
  
  labels_gene <- setNames(paste0(top_df$GENENAME, "  (", ifelse(lfc_val > 0, "+", ""), lfc_val, ")"),
                          top_df$GENENAME)
  
  vmax <- suppressWarnings(max(abs(df_long$zscore), na.rm = TRUE))
  if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = sample, y = gene, fill = zscore)) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::facet_grid(direction ~ condition, scales = "free", space = "free", switch = "y") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-vmax, vmax), na.value = "grey95",
      name = "z-score\n(expression)"
    ) +
    ggplot2::scale_y_discrete(labels = labels_gene) +
    ggplot2::labs(
      title = paste0(dataset_name, " / ", contrast_name, " \u2014 top ", nrow(top_df),
                     " DEGs significatifs (|LFC| le plus fort)"),
      subtitle = paste0("padj<", padj_thr, ", |LFC|\u2265", lfc_thr,
                        "\nLFC affich\u00e9 entre parenth\u00e8ses = log2FoldChange (", condB, " vs ", condA, ")"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y       = ggplot2::element_text(size = 7.5),
      panel.grid         = ggplot2::element_blank(),
      plot.title          = ggplot2::element_text(face = "bold", colour = "#1F4E5F",
                                                  margin = ggplot2::margin(b = 4)),
      plot.subtitle        = ggplot2::element_text(colour = "grey30", size = 9,
                                                   margin = ggplot2::margin(b = 10)),
      plot.margin            = ggplot2::margin(t = 12, r = 12, b = 8, l = 8),
      strip.text.y.left      = ggplot2::element_text(angle = 0, size = 7.5, face = "bold"),
      strip.text.x             = ggplot2::element_text(size = 8, face = "bold"),
      strip.placement            = "outside",
      strip.background            = ggplot2::element_rect(fill = "#EAF2F0", color = NA),
      panel.spacing                = grid::unit(2, "pt"),
      legend.position                = "right"
    )
  
  n_genes   <- length(unique(df_long$gene))
  n_samples <- length(unique(df_long$sample))
  h <- max(6, n_genes * 0.19 + 2.3)
  w <- max(9.5, n_samples * 0.45 + 3.5)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, file_name)
  
  grDevices::pdf(out_file, width = w, height = h)
  print(p)
  grDevices::dev.off()
  
  message("  Heatmap top-LFC (", dataset_name, "/", contrast_name, ") : \u00e9crit -> ", out_file)
  invisible(p)
}

# ------------------------------------------------------------
# Lance plot_top_lfc_sample_heatmap() pour TOUS les contrastes de TOUS les
# datasets disponibles dans all_dataset_results, en calculant le vst une
# seule fois par dataset (réutilisé pour ses différents contrastes).
# Écrit chaque heatmap dans results_v4/<dataset>/<contraste>/.
# ------------------------------------------------------------
run_top_lfc_sample_heatmaps <- function(
    all_dataset_results, tx2gene, base_out,
    top_n    = 30,
    padj_thr = get0("PADJ_THRESHOLD", ifnotfound = 0.05),
    lfc_thr  = get0("LFC_THRESHOLD",  ifnotfound = 1)
) {
  for (ds in names(all_dataset_results)) {
    
    vst_data <- tryCatch(
      .top_lfc_compute_dataset_vst(ds, tx2gene),
      error = function(e) { message("  [WARN] vst (", ds, ") \u00e9chou\u00e9 : ", conditionMessage(e)); NULL }
    )
    if (is.null(vst_data)) next
    
    for (cn in names(all_dataset_results[[ds]])) {
      tryCatch(
        plot_top_lfc_sample_heatmap(
          dataset_name  = ds, contrast_name = cn,
          res_df        = all_dataset_results[[ds]][[cn]],
          out_dir       = file.path(base_out, ds, cn),
          vst_data      = vst_data,
          top_n         = top_n, padj_thr = padj_thr, lfc_thr = lfc_thr
        ),
        error = function(e) message("  [WARN] heatmap top-LFC (", ds, "/", cn,
                                    ") \u00e9chou\u00e9e : ", conditionMessage(e))
      )
    }
  }
  invisible(NULL)
}