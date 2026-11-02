# ============================================================
# cross_dataset_ipsc_vs_neurons.R
# Comparaison iPSC vs iNeurons à condition fixée
# (Control, FA, FA_edited)
# ============================================================

run_ipsc_vs_neurons <- function(all_vsd_results,
                                all_dataset_results,
                                base_out,
                                name_a          = "iPSC",
                                name_b          = "iNeurons",
                                conditions      = c("control", "FA", "FA_edited"),
                                padj_thr        = 0.05,
                                lfc_thr         = 1,
                                genes_highlight = NULL,
                                top_n_heatmap   = 50,
                                top_n_volcano   = 60) {
  
  out_dir <- file.path(base_out, "iPSC_vs_iNeurons")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_cond_results <- list()
  
  # ============================================================
  # Boucle par condition
  # ============================================================
  for (cond in conditions) {
    
    message("\n--- iPSC vs iNeurons : condition = ", cond, " ---")
    
    cond_dir <- file.path(out_dir, cond)
    dir.create(cond_dir, showWarnings = FALSE, recursive = TRUE)
    
    # ---- Samples ----
    samples_a <- get_samples_for_condition(DATASETS[[name_a]], cond)
    samples_b <- get_samples_for_condition(DATASETS[[name_b]], cond)
    
    if (length(samples_a) == 0) {
      message("  [SKIP] Pas de samples ", name_a, " pour condition ", cond)
      next
    }
    if (length(samples_b) == 0) {
      message("  [SKIP] Pas de samples ", name_b, " pour condition ", cond)
      next
    }
    
    # ---- Fichiers quant.sf ----
    files_a <- file.path(DATASETS[[name_a]]$path, samples_a, "quant.sf")
    files_b <- file.path(DATASETS[[name_b]]$path, samples_b, "quant.sf")
    names(files_a) <- paste0(name_a, "_", samples_a)
    names(files_b) <- paste0(name_b, "_", samples_b)
    
    files_all <- c(files_a, files_b)
    
    missing <- files_all[!file.exists(files_all)]
    if (length(missing) > 0) {
      message("  [WARN] Fichiers manquants : ",
              paste(basename(dirname(missing)), collapse = ", "))
      files_all <- files_all[file.exists(files_all)]
    }
    
    if (length(files_all) < 4) {
      message("  [SKIP] Pas assez de samples (n=", length(files_all), ")")
      next
    }
    
    # ---- tximport ----
    txi <- tximport::tximport(
      files_all,
      type            = "salmon",
      tx2gene         = tx2gene,
      ignoreTxVersion = TRUE
    )
    
    # ---- Métadonnées ----
    valid_a <- names(files_all)[names(files_all) %in% names(files_a)]
    valid_b <- names(files_all)[names(files_all) %in% names(files_b)]
    
    cell_type        <- c(rep(name_a, length(valid_a)),
                          rep(name_b, length(valid_b)))
    names(cell_type) <- c(valid_a, valid_b)
    cell_type        <- cell_type[names(files_all)]
    
    donor_a          <- get_donor_for_samples(
      DATASETS[[name_a]],
      sub(paste0("^", name_a, "_"), "", valid_a))
    donor_b          <- get_donor_for_samples(
      DATASETS[[name_b]],
      sub(paste0("^", name_b, "_"), "", valid_b))
    names(donor_a)   <- valid_a
    names(donor_b)   <- valid_b
    donor            <- c(donor_a, donor_b)[names(files_all)]
    
    col_data <- data.frame(
      cell_type = factor(cell_type, levels = c(name_a, name_b)),
      donor     = factor(donor),
      row.names = names(files_all)
    )
    
    # ---- Design : ignorer donor si c'est juste Male/Female ----
    use_donor <- length(unique(donor)) > 1 &&
      !all(donor %in% c("Male", "Female")) &&
      !any(is.na(donor)) &&
      length(unique(donor)) < length(donor)
    
    design_formula <- if (use_donor) ~ donor + cell_type else ~ cell_type
    message("  Design : ", deparse(design_formula))
    
    tryCatch({
      
      # ---- DESeq2 ----
      dds <- DESeq2::DESeqDataSetFromTximport(
        txi     = txi,
        colData = col_data,
        design  = design_formula
      )
      
      dds <- dds[rowSums(DESeq2::counts(dds) >= MIN_COUNTS) >= MIN_SAMPLES, ]
      dds <- DESeq2::DESeq(dds, quiet = FALSE)
      
      res <- DESeq2::results(
        dds,
        contrast             = c("cell_type", name_a, name_b),
        independentFiltering = TRUE
      )
      
      vsd <- DESeq2::vst(dds, blind = FALSE)
      
      res_df <- as.data.frame(res) %>%
        tibble::rownames_to_column("GENEID") %>%
        dplyr::left_join(gene_map, by = "GENEID") %>%
        dplyr::filter(!is.na(padj)) %>%
        dplyr::arrange(padj) %>%
        dplyr::mutate(
          sig       = ifelse(padj < padj_thr & abs(log2FoldChange) >= lfc_thr,
                             "Significant", "NS"),
          direction = dplyr::case_when(
            padj < padj_thr & log2FoldChange >  lfc_thr ~ paste0("Up in ", name_a),
            padj < padj_thr & log2FoldChange < -lfc_thr ~ paste0("Up in ", name_b),
            TRUE ~ "NS"
          ),
          condition = cond
        )
      
      all_cond_results[[cond]] <- res_df
      
      message("  Gènes significatifs : ", sum(res_df$sig == "Significant"))
      
      # ---- Export DEG table ----
      writexl::write_xlsx(
        res_df,
        file.path(cond_dir, paste0("DEG_iPSC_vs_iNeurons_", cond, ".xlsx"))
      )
      
      # ---- Volcano standard ----
      pdf(file.path(cond_dir,
                    paste0("volcano_iPSC_vs_iNeurons_", cond, ".pdf")),
          width = 12, height = 9)
      print(plot_volcano_ipsc_neurons(
        res_df,
        title           = paste0("iPSC vs iNeurons — ", cond),
        name_a          = name_a,
        name_b          = name_b,
        lfc_thr         = lfc_thr,
        padj_thr        = padj_thr,
        top_n           = top_n_volcano,
        genes_highlight = genes_highlight
      ))
      dev.off()
      
      # ---- Volcano GOI ----
      if (!is.null(genes_highlight) && length(genes_highlight) > 0) {
        pdf(file.path(cond_dir,
                      paste0("volcano_GOI_iPSC_vs_iNeurons_", cond, ".pdf")),
            width = 12, height = 9)
        print(plot_volcano_ipsc_neurons(
          res_df,
          title           = paste0("iPSC vs iNeurons — ", cond, " (GOI)"),
          name_a          = name_a,
          name_b          = name_b,
          lfc_thr         = lfc_thr,
          padj_thr        = padj_thr,
          top_n           = 0,
          genes_highlight = genes_highlight,
          force_label     = TRUE
        ))
        dev.off()
      }
      
      # ---- Heatmap top DEGs ----
      vsd_mat   <- SummarizedExperiment::assay(vsd)
      top_genes <- res_df %>%
        dplyr::filter(sig == "Significant") %>%
        dplyr::slice_min(padj, n = top_n_heatmap) %>%
        dplyr::pull(GENEID)
      
      if (length(top_genes) >= 2) {
        
        mat_heat <- vsd_mat[rownames(vsd_mat) %in% top_genes,
                            rownames(col_data), drop = FALSE]
        
        # Nommage propre : SYMBOL si dispo, sinon ENSG
        syms <- gene_map$GENENAME[match(rownames(mat_heat), gene_map$GENEID)]
        rownames(mat_heat) <- ifelse(is.na(syms) | syms == "",
                                     rownames(mat_heat), syms)
        
        # Dédupliquer : garder le LFC max par SYMBOL
        mat_heat <- dedup_by_max_lfc(mat_heat)
        
        mat_heat <- t(scale(t(mat_heat)))
        mat_heat[!is.finite(mat_heat)] <- 0   # lignes constantes
        
        ann_col <- data.frame(
          CellType = col_data[rownames(col_data), "cell_type"],
          Donor    = col_data[rownames(col_data), "donor"],
          row.names = rownames(col_data)
        )
        
        pdf(file.path(cond_dir,
                      paste0("heatmap_topDEGs_iPSC_vs_iNeurons_", cond, ".pdf")),
            width = 10, height = 14)
        pheatmap::pheatmap(
          mat_heat,
          annotation_col    = ann_col,
          annotation_colors = list(
            CellType = setNames(c("#4393c3", "#d6604d"), c(name_a, name_b))
          ),
          show_colnames     = FALSE,
          fontsize_row      = 7,
          cluster_rows      = nrow(mat_heat) >= 2,   # hclust exige >= 2 lignes
          clustering_method = "ward.D2",
          color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
          main  = paste0("Top ", top_n_heatmap,
                         " DEGs — iPSC vs iNeurons — ", cond)
        )
        dev.off()
      }
      
      # ---- PCA ----
      pdf(file.path(cond_dir,
                    paste0("PCA_iPSC_vs_iNeurons_", cond, ".pdf")))
      tryCatch({
        pca_data <- DESeq2::plotPCA(vsd, intgroup = "cell_type",
                                    returnData = TRUE)
        pct_var  <- round(100 * attr(pca_data, "percentVar"))
        print(
          ggplot2::ggplot(pca_data,
                          ggplot2::aes(PC1, PC2, color = cell_type,
                                       label = name)) +
            ggplot2::geom_point(size = 3) +
            ggrepel::geom_label_repel(size = 3, show.legend = FALSE) +
            ggplot2::scale_color_manual(
              values = setNames(c("#4393c3", "#d6604d"), c(name_a, name_b))
            ) +
            ggplot2::theme_bw(base_size = 11) +
            ggplot2::labs(
              title = paste0("PCA — iPSC vs iNeurons — ", cond),
              x     = paste0("PC1 (", pct_var[1], "%)"),
              y     = paste0("PC2 (", pct_var[2], "%)")
            )
        )
      }, error = function(e) message("  [WARN] PCA skipped : ", e$message))
      dev.off()
      
      # ---- fgsea ----
      ranks <- make_ranks(res_df %>% dplyr::rename(SYMBOL = GENENAME))
      run_all_fgsea(
        ranks             = ranks,
        out_dir           = cond_dir,
        label             = paste0("iPSC_vs_iNeurons_", cond),
        collections       = FGSEA_COLLECTIONS,
        neuronal_keywords = NEURONAL_KEYWORDS,
        direction_label   = paste0(name_b, " vs ", name_a),
        nperm             = FGSEA_NPERM,
        top_n             = FGSEA_TOP_N,
        padj_cut          = FGSEA_PADJ_CUT
      )
      
    }, error = function(e) {
      message("  [ERROR] ", cond, " : ", e$message)
    })
  }
  
  # ============================================================
  # Scatter cross-conditions
  # ============================================================
  if (length(all_cond_results) >= 2) {
    
    message("\n--- Scatter cross-conditions iPSC vs Neurons ---")
    
    scatter_dir <- file.path(out_dir, "SCATTER_cross_conditions")
    dir.create(scatter_dir, showWarnings = FALSE, recursive = TRUE)
    
    cond_names <- names(all_cond_results)
    
    for (i in seq_len(length(cond_names) - 1)) {
      for (j in seq(i + 1, length(cond_names))) {
        
        ca <- cond_names[i]
        cb <- cond_names[j]
        
        df_a <- all_cond_results[[ca]] %>%
          dplyr::select(GENEID, GENENAME, log2FoldChange, padj, sig) %>%
          dplyr::rename(lfc_a = log2FoldChange, padj_a = padj, sig_a = sig)
        
        df_b <- all_cond_results[[cb]] %>%
          dplyr::select(GENEID, log2FoldChange, padj, sig) %>%
          dplyr::rename(lfc_b = log2FoldChange, padj_b = padj, sig_b = sig)
        
        df_merge <- dplyr::inner_join(df_a, df_b, by = "GENEID") %>%
          dplyr::mutate(
            sig_both = dplyr::case_when(
              sig_a == "Significant" & sig_b == "Significant" ~ "Both",
              sig_a == "Significant" ~ paste0("Only ", ca),
              sig_b == "Significant" ~ paste0("Only ", cb),
              TRUE ~ "NS"
            )
          )
        
        r_val    <- round(cor(df_merge$lfc_a, df_merge$lfc_b,
                              use = "complete.obs"), 3)
        top_both <- df_merge %>%
          dplyr::filter(sig_both == "Both") %>%
          dplyr::mutate(dist = sqrt(lfc_a^2 + lfc_b^2)) %>%
          dplyr::slice_max(dist, n = 15)
        
        pdf(file.path(scatter_dir,
                      paste0("scatter_", ca, "_vs_", cb, ".pdf")),
            width = 10, height = 9)
        print(
          ggplot2::ggplot(df_merge,
                          ggplot2::aes(x = lfc_a, y = lfc_b,
                                       color = sig_both)) +
            ggplot2::geom_point(alpha = 0.4, size = 0.8) +
            ggrepel::geom_label_repel(
              data         = top_both,
              ggplot2::aes(label = GENENAME),
              size         = 2.5,
              max.overlaps = 15,
              show.legend  = FALSE
            ) +
            ggplot2::scale_color_manual(values = setNames(
              c("#7b2d8b", "#2166ac", "#d6604d", "grey70"),
              c("Both", paste0("Only ", ca), paste0("Only ", cb), "NS")
            )) +
            ggplot2::geom_hline(yintercept = 0,
                                linetype = "dashed", color = "grey50") +
            ggplot2::geom_vline(xintercept = 0,
                                linetype = "dashed", color = "grey50") +
            ggplot2::geom_abline(slope = 1, intercept = 0,
                                 linetype = "dotted", color = "grey30") +
            ggplot2::theme_bw(base_size = 11) +
            ggplot2::labs(
              title    = paste0("iPSC vs iNeurons : ", ca, " vs ", cb),
              subtitle = paste0("r = ", r_val,
                                " — n gènes = ", nrow(df_merge)),
              x        = paste0("log2FC iPSC/iNeurons — ", ca),
              y        = paste0("log2FC iPSC/iNeurons — ", cb),
              color    = NULL
            ) +
            ggplot2::theme(legend.position = "top")
        )
        dev.off()
      }
    }
  }
  
  # ============================================================
  # Heatmap récapitulative log2FC × conditions
  # ============================================================
  if (length(all_cond_results) >= 2) {
    
    message("\n--- Heatmap récapitulative log2FC ---")
    
    # Gènes significatifs dans au moins une condition
    sig_genes_all <- unique(unlist(lapply(all_cond_results, function(df) {
      df$GENEID[df$sig == "Significant"]
    })))
    
    if (length(sig_genes_all) > 1) {
      
      # Matrice LFC : gènes × conditions
      lfc_mat <- do.call(cbind, lapply(names(all_cond_results), function(cond) {
        df  <- all_cond_results[[cond]]
        lfc <- setNames(df$log2FoldChange, df$GENEID)
        lfc[sig_genes_all]
      }))
      colnames(lfc_mat) <- names(all_cond_results)
      
      # Diagnostic mapping
      symbols    <- gene_map$GENENAME[match(sig_genes_all, gene_map$GENEID)]
      n_total    <- length(sig_genes_all)
      n_no_match <- sum(is.na(symbols) | symbols == "")
      n_dup      <- sum(duplicated(symbols[!is.na(symbols) & symbols != ""]))
      message("Gènes significatifs total : ", n_total)
      message("Sans correspondance SYMBOL : ", n_no_match,
              " (", round(n_no_match / n_total * 100, 1), "%)")
      message("SYMBOL dupliqués : ", n_dup)
      
      # Nommage : SYMBOL si dispo, sinon ENSG
      rownames(lfc_mat) <- ifelse(is.na(symbols) | symbols == "",
                                  sig_genes_all, symbols)
      
      # Dédupliquer : pour chaque SYMBOL, garder la ligne
      # avec le LFC le plus extrême dans au moins une condition
      lfc_mat <- dedup_by_max_lfc(lfc_mat)
      
      # Sélection top 100 par variance inter-conditions
      if (nrow(lfc_mat) > 100) {
        row_var <- apply(lfc_mat, 1, var, na.rm = TRUE)
        lfc_mat <- lfc_mat[order(row_var, decreasing = TRUE)[1:100], ]
      }
      
      pdf(file.path(out_dir, "heatmap_summary_LFC_by_condition.pdf"),
          width = 8, height = 14)
      pheatmap::pheatmap(
        lfc_mat,
        color             = colorRampPalette(
          c("#2166ac", "white", "#d6604d"))(100),
        breaks            = seq(-4, 4, length.out = 101),
        cluster_rows      = nrow(lfc_mat) >= 2,   # hclust exige >= 2 lignes
        clustering_method = "ward.D2",
        fontsize_row      = 6,
        main              = paste0(
          "log2FC (iPSC / iNeurons) — top DEGs\n",
          "par condition (Control / FA / FA_edited)"
        ),
        na_col = "grey80"
      )
      dev.off()
    }
  }
  
  # ============================================================
  # Export global
  # ============================================================
  if (length(all_cond_results) > 0) {
    writexl::write_xlsx(
      all_cond_results,
      file.path(out_dir, "DEG_iPSC_vs_iNeurons_ALL_conditions.xlsx")
    )
    message("\n✔ Résultats exportés dans : ", out_dir)
  }
  
  return(invisible(all_cond_results))
}


# ============================================================
# Helpers
# ============================================================

# Dédupliquer une matrice dont les rownames sont des SYMBOL :
# pour chaque SYMBOL dupliqué, garder la ligne avec le max |LFC|
dedup_by_max_lfc <- function(mat) {
  # max absolu sur toutes les colonnes pour chaque ligne
  row_max_abs <- apply(mat, 1, function(x) max(abs(x), na.rm = TRUE))
  df <- data.frame(
    SYMBOL      = rownames(mat),
    row_max_abs = row_max_abs,
    row_idx     = seq_len(nrow(mat)),
    stringsAsFactors = FALSE
  )
  # Pour chaque SYMBOL, garder l'index avec le max
  keep_idx <- df %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(row_max_abs, n = 1, with_ties = FALSE) %>%
    dplyr::pull(row_idx)
  
  mat[keep_idx, , drop = FALSE]
}

get_samples_for_condition <- function(dataset, condition) {
  idx <- which(dataset$conditions == condition)
  dataset$samples[idx]
}

get_donor_for_samples <- function(dataset, sample_names) {
  idx    <- match(sample_names, dataset$samples)
  donors <- dataset$sex[idx]
  setNames(donors, sample_names)
}


# ============================================================
# Volcano plot
# ============================================================

plot_volcano_ipsc_neurons <- function(res_df,
                                      title           = "",
                                      name_a          = "iPSC",
                                      name_b          = "iNeurons",
                                      lfc_thr         = 1,
                                      padj_thr        = 0.05,
                                      top_n           = 20,
                                      genes_highlight = NULL,
                                      force_label     = FALSE) {
  
  res_df <- res_df %>%
    dplyr::mutate(log10p = -log10(pmax(padj, 1e-300)))
  
  top_auto <- if (top_n > 0) {
    res_df %>%
      dplyr::filter(sig == "Significant") %>%
      dplyr::slice_min(padj, n = top_n)
  } else {
    res_df[0, ]
  }
  
  goi_df <- if (!is.null(genes_highlight) && length(genes_highlight) > 0) {
    res_df %>% dplyr::filter(GENENAME %in% genes_highlight)
  } else {
    res_df[0, ]
  }
  
  label_df <- dplyr::bind_rows(top_auto, goi_df) %>%
    dplyr::distinct(GENEID, .keep_all = TRUE)
  
  ggplot2::ggplot(res_df,
                  ggplot2::aes(x = log2FoldChange, y = log10p,
                               color = direction)) +
    ggplot2::geom_point(alpha = 0.45, size = 1.0) +
    ggrepel::geom_label_repel(
      data         = label_df,
      ggplot2::aes(label = GENENAME),
      size         = 2.8,
      max.overlaps = 25,
      segment.size = 0.3,
      show.legend  = FALSE
    ) +
    ggplot2::scale_color_manual(values = setNames(
      c("grey70", "#4393c3", "#d6604d"),
      c("NS", paste0("Up in ", name_a), paste0("Up in ", name_b))
    )) +
    ggplot2::geom_hline(yintercept = -log10(padj_thr),
                        linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr),
                        linetype = "dashed", color = "grey40") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0("Bleu = enrichi dans ", name_a,
                        " | Rouge = enrichi dans ", name_b,
                        " | padj < ", padj_thr, ", |LFC| ≥ ", lfc_thr),
      x        = paste0("log2FC (", name_a, " / ", name_b, ")"),
      y        = "-log10(padj)",
      color    = NULL
    ) +
    ggplot2::theme(legend.position = "top")
}