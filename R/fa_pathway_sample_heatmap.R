# ============================================================
# fa_pathway_sample_heatmap.R
#
# Même design visuel que fa_pathway_heatmap.R (ggplot2 : tuiles, facettes
# par catégorie/pathway, échelle de couleur divergente, astérisque pour la
# significativité) -- mais les colonnes sont ici les ÉCHANTILLONS individuels
# (facettés par condition) au lieu de dataset x méthode. Une carte par type
# cellulaire (fibro, iPSC, iNeurons), méthode standard (vst, design ~1,
# comme le QC global) -- pas de comparaison de méthodes ici.
# ============================================================

# ------------------------------------------------------------
# Une carte pour un seul dataset
# ------------------------------------------------------------
plot_fa_pathway_sample_heatmap <- function(
    dataset_name, tx2gene, gene_list, out_dir,
    all_dataset_results = NULL,
    contrast_sig  = "control_vs_FA",
    padj_thr      = get0("PADJ_THRESHOLD", ifnotfound = 0.05),
    lfc_thr       = get0("LFC_THRESHOLD",  ifnotfound = 1),
    file_name     = NULL,
    min_counts    = get0("MIN_COUNTS",  ifnotfound = 10),
    min_samples   = get0("MIN_SAMPLES", ifnotfound = 3)
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("  [INFO] ggplot2 indisponible : heatmap par échantillon ignoré.")
    return(invisible(NULL))
  }
  
  dataset <- DATASETS[[dataset_name]]
  if (is.null(dataset)) {
    message("  [WARN] dataset '", dataset_name, "' introuvable dans DATASETS.")
    return(invisible(NULL))
  }
  
  samples <- dataset$samples
  files   <- file.path(dataset$path, samples, "quant.sf")
  names(files) <- samples
  if (!all(file.exists(files))) {
    message("  [WARN] quant.sf manquants pour ", dataset_name,
            " -> heatmap par échantillon ignoré : ",
            paste(samples[!file.exists(files)], collapse = ", "))
    return(invisible(NULL))
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
  
  # ---- Gènes du panel présents dans ce dataset (+ catégorie/pathway) ----
  gm <- unique(tx2gene[, c("GENEID", "GENENAME")])
  gene_df <- do.call(rbind, lapply(names(gene_list), function(cat) {
    v <- gene_list[[cat]]
    data.frame(categorie = cat, symbole = unname(v), stringsAsFactors = FALSE)
  }))
  gene_df <- unique(gene_df)
  gene_df$ensg <- gm$GENEID[match(gene_df$symbole, gm$GENENAME)]
  gene_df <- gene_df[!is.na(gene_df$ensg) & gene_df$ensg %in% rownames(mat_full), ]
  
  if (nrow(gene_df) < 2) {
    message("  [WARN] Heatmap par échantillon (", dataset_name,
            ") : moins de 2 gènes du panel présents après filtrage -- ignoré.")
    return(invisible(NULL))
  }
  
  mat <- mat_full[gene_df$ensg, , drop = FALSE]
  rownames(mat) <- gene_df$symbole
  
  # z-score par gène (ligne) à travers les échantillons -- équivalent de
  # scale="row" en pheatmap, mais tracé en ggplot2 pour le même look que
  # fa_pathway_heatmap.R.
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  
  # ---- Significativité (contraste contrast_sig de CE dataset, si dispo) ----
  # Reprend all_dataset_results déjà calculé -- aucun recalcul DESeq2.
  sig_lookup <- NULL
  if (!is.null(all_dataset_results) &&
      !is.null(all_dataset_results[[dataset_name]]) &&
      !is.null(all_dataset_results[[dataset_name]][[contrast_sig]])) {
    res_df <- all_dataset_results[[dataset_name]][[contrast_sig]]
    sig_lookup <- setNames(
      !is.na(res_df$padj) & res_df$padj < padj_thr & abs(res_df$log2FoldChange) >= lfc_thr,
      res_df$GENENAME
    )
  }
  gene_df$significatif <- if (!is.null(sig_lookup)) {
    sig_lookup[gene_df$symbole] %in% TRUE
  } else {
    rep(FALSE, nrow(gene_df))
  }
  
  cat_order  <- names(gene_list)
  gene_order <- gene_df$symbole[order(match(gene_df$categorie, cat_order), gene_df$symbole)]
  
  # ---- Passage en format long pour ggplot2 ----
  df_long <- as.data.frame(as.table(mat_z))
  colnames(df_long) <- c("gene", "sample", "zscore")
  df_long$gene   <- as.character(df_long$gene)
  df_long$sample <- as.character(df_long$sample)
  df_long$condition    <- colData$condition[match(df_long$sample, colData$sample)]
  df_long$categorie    <- factor(gene_df$categorie[match(df_long$gene, gene_df$symbole)],
                                 levels = cat_order[cat_order %in% gene_df$categorie])
  df_long$significatif <- gene_df$significatif[match(df_long$gene, gene_df$symbole)]
  
  df_long$gene <- factor(df_long$gene, levels = rev(gene_order))
  
  sample_order <- samples[order(match(dataset$conditions, unique(dataset$conditions)), samples)]
  df_long$sample <- factor(df_long$sample, levels = sample_order)
  
  # Astérisque directement sur le label du gène (propriété de la ligne
  # entière, pas d'une cellule précise -- donc marqué sur l'axe, pas sur
  # chaque tuile).
  labels_gene <- setNames(
    ifelse(gene_df$significatif, paste0(gene_df$symbole, "  *"), gene_df$symbole),
    gene_df$symbole
  )
  
  vmax <- suppressWarnings(max(abs(df_long$zscore), na.rm = TRUE))
  if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = sample, y = gene, fill = zscore)) +
    ggplot2::geom_tile(color = "grey90") +
    ggplot2::facet_grid(categorie ~ condition, scales = "free", space = "free", switch = "y") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-vmax, vmax), na.value = "grey95",
      name = "z-score\n(expression)"
    ) +
    ggplot2::scale_y_discrete(labels = labels_gene) +
    ggplot2::labs(
      title = paste0(dataset_name, " \u2014 expression du panel par \u00e9chantillon"),
      subtitle = paste0("M\u00e9thode standard (vst)\n* = significatif dans ", contrast_sig,
                        " (padj<", padj_thr, ", |LFC|\u2265", lfc_thr, ")"),
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
      strip.text.y.left      = ggplot2::element_text(angle = 0, size = 7, face = "bold"),
      strip.text.x             = ggplot2::element_text(size = 8, face = "bold"),
      strip.placement            = "outside",
      strip.background            = ggplot2::element_rect(fill = "#EAF2F0", color = NA),
      panel.spacing                = grid::unit(2, "pt"),
      legend.position                = "right"
    )
  
  n_genes   <- length(unique(df_long$gene))
  n_samples <- length(unique(df_long$sample))
  h <- max(6, n_genes * 0.18 + 2.3)
  w <- max(9.5, n_samples * 0.45 + 3.5)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(file_name)) file_name <- paste0("heatmap_echantillons_", dataset_name, ".pdf")
  out_file <- file.path(out_dir, file_name)
  
  grDevices::pdf(out_file, width = w, height = h)
  print(p)
  grDevices::dev.off()
  
  message("  Heatmap par \u00e9chantillon (", dataset_name, ") : \u00e9crit -> ", out_file)
  invisible(p)
}

# ------------------------------------------------------------
# Lance plot_fa_pathway_sample_heatmap() pour plusieurs datasets d'un coup ;
# chaque dataset est protégé individuellement (l'échec de l'un n'empêche
# pas les autres de se générer).
# ------------------------------------------------------------
run_fa_pathway_sample_heatmaps <- function(
    tx2gene, gene_list, out_dir,
    all_dataset_results = NULL,
    datasets = c("fibro", "iPSC", "iNeurons")
) {
  for (ds in datasets) {
    tryCatch(
      plot_fa_pathway_sample_heatmap(
        dataset_name = ds, tx2gene = tx2gene, gene_list = gene_list, out_dir = out_dir,
        all_dataset_results = all_dataset_results
      ),
      error = function(e) message("  [WARN] heatmap par échantillon (", ds,
                                  ") échouée : ", conditionMessage(e))
    )
  }
  invisible(NULL)
}