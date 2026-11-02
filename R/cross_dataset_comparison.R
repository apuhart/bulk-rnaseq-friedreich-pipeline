# ============================================================
# cross_dataset_comparison.R
# Compare DEGs between two datasets (e.g. iPSC vs iNeurons)
# for matching contrasts: control, FA_edited, FA
#
# Source this file AFTER running all datasets, then call:
#   run_cross_dataset_comparison(all_dataset_results, BASE_OUT)
# ============================================================

# ----------------------------
# Dependencies (already loaded in main pipeline)
# Required additionally:
#   ggplot2, ggrepel, VennDiagram, grid, gridExtra
# ----------------------------

suppressPackageStartupMessages({
  for (pkg in c("ggplot2", "ggrepel", "VennDiagram", "grid", "gridExtra", "scales", "dplyr", "tidyr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = c(CRAN = "https://cloud.r-project.org"), quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
})


# ============================================================
# SECTION A — Collect results across all datasets
# ============================================================

#' Collect DEG result data frames from all datasets into a nested list
#'
#' Call run_dataset() with `return_results = TRUE` (see patched run_dataset below),
#' OR build this list manually:
#'
#'   all_dataset_results <- list(
#'     iPSC     = list(control_vs_FA = df1, control_vs_FA_edited = df2, FA_edited_vs_FA = df3),
#'     iNeurons = list(control_vs_FA = df4, ...)
#'   )
#'
#' Each df must contain columns: GENEID, GENENAME, log2FoldChange, padj, baseMean
collect_dataset_results <- function() {
  # This function is a placeholder documenting the expected structure.
  # See the patched run_dataset() at the bottom of this file.
  stop("Use the patched run_dataset() and collect results via `all_dataset_results`.")
}


# ============================================================
# SECTION B — Core comparison helpers
# ============================================================

#' Merge two DEG result data frames on gene ID
#'
#' @param df_a  DEG data.frame from dataset A (must have GENEID, log2FoldChange, padj)
#' @param df_b  DEG data.frame from dataset B
#' @param name_a  Label for dataset A (e.g. "iPSC")
#' @param name_b  Label for dataset B (e.g. "iNeurons")
#' @return Merged data.frame with suffixed columns
merge_deg_results <- function(df_a, df_b, name_a, name_b) {
  cols <- c("GENEID", "GENENAME", "log2FoldChange", "padj", "baseMean")
  
  # Keep only columns that exist
  cols_a <- intersect(cols, colnames(df_a))
  cols_b <- intersect(cols, colnames(df_b))
  
  merged <- dplyr::inner_join(
    df_a[, cols_a],
    df_b[, cols_b],
    by   = c("GENEID", "GENENAME"),
    suffix = c(paste0(".", name_a), paste0(".", name_b))
  )
  merged
}


#' Classify genes by significance in each dataset
#'
#' Returns a factor with levels:
#'   "Both", "Only_A", "Only_B", "Neither"
#'
classify_deg_overlap <- function(merged, name_a, name_b,
                                 padj_thr = 0.05, lfc_thr = 0) {
  
  sig_a <- !is.na(merged[[paste0("padj.", name_a)]]) &
    merged[[paste0("padj.", name_a)]] < padj_thr &
    abs(merged[[paste0("log2FoldChange.", name_a)]]) > lfc_thr
  
  sig_b <- !is.na(merged[[paste0("padj.", name_b)]]) &
    merged[[paste0("padj.", name_b)]] < padj_thr &
    abs(merged[[paste0("log2FoldChange.", name_b)]]) > lfc_thr
  
  dplyr::case_when(
    sig_a & sig_b  ~ "Both",
    sig_a & !sig_b ~ paste0("Only ", name_a),
    !sig_a & sig_b ~ paste0("Only ", name_b),
    TRUE           ~ "Neither"
  )
}


# ============================================================
# SECTION C — Plot: LFC scatter with overlap coloring
# ============================================================

#' Scatter plot of log2FoldChange: dataset A vs dataset B
#'
#' Points colored by overlap category; top genes labeled.
#'



plot_lfc_scatter <- function(merged, name_a, name_b,
                             contrast_name,
                             padj_thr  = 0.05,
                             lfc_thr   = 0,
                             top_label = 20,
                             genes_highlight = NULL) {
  
  lfc_col_a <- paste0("log2FoldChange.", name_a)
  lfc_col_b <- paste0("log2FoldChange.", name_b)
  
  merged$overlap <- classify_deg_overlap(merged, name_a, name_b, padj_thr, lfc_thr)
  
  # Pearson correlation on finite values
  finite_mask <- is.finite(merged[[lfc_col_a]]) & is.finite(merged[[lfc_col_b]])
  r_val <- cor(merged[[lfc_col_a]][finite_mask],
               merged[[lfc_col_b]][finite_mask],
               method = "pearson")
  
  # Color palette
  lvls <- c("Both",
            paste0("Only ", name_a),
            paste0("Only ", name_b),
            "Neither")

  pal <- c(
    "Both" = "#E41A1C",
    "Neither" = "grey70"
  )
  
  pal[paste0("Only ", name_a)] <- "#377EB8"
  pal[paste0("Only ", name_b)] <- "#FF7F00"
  
  
  # Pick top genes to label: prioritise "Both", then most extreme LFC in A
  label_df <- merged[merged$overlap != "Neither", ]
  label_df <- label_df[order(
    label_df$overlap == "Both",
    abs(label_df[[lfc_col_a]]),
    decreasing = TRUE
  ), ]
  label_df <- utils::head(label_df, top_label)
  
  # Add user-defined highlights
  if (!is.null(genes_highlight)) {
    extra <- merged[merged$GENENAME %in% genes_highlight & !merged$GENEID %in% label_df$GENEID, ]
    label_df <- rbind(label_df, extra)
  }
  
  
  
  
  
  # Define alpha levels
  alpha_vals <- c(
    "Both" = 0.8,
    "Neither" = 0.1
  )
  
  alpha_vals[paste0("Only ", name_a)] <- 0.8
  alpha_vals[paste0("Only ", name_b)] <- 0.8
  

  
  lfc_label_thr <- 5  # à ajuster
  
  extreme <- abs(merged[[lfc_col_a]]) > lfc_label_thr |
    abs(merged[[lfc_col_b]]) > lfc_label_thr
  
  label_df <- merged[
    (merged$overlap != "Neither" & extreme) |
      merged$overlap == "Both",
  ]
  
  
  ggplot2::ggplot(merged, ggplot2::aes(
    x     = .data[[lfc_col_a]],
    y     = .data[[lfc_col_b]],
    color = overlap,
    alpha = overlap   # <- map alpha
  )) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 1.2) +
    ggrepel::geom_text_repel(
      data          = label_df,
      ggplot2::aes(label = GENENAME),
      size          = 2.8,
      max.overlaps  = Inf,   # <- allow all labels
      segment.color = "grey40",
      show.legend   = FALSE
    ) +
    ggplot2::scale_color_manual(values = pal, name = "Significant in") +
    ggplot2::scale_alpha_manual(values = alpha_vals, guide = "none") 
    
    
    
  #   
  # ggplot2::ggplot(merged, ggplot2::aes(
  #   x     = .data[[lfc_col_a]],
  #   y     = .data[[lfc_col_b]],
  #   color = overlap
  # )) +
  #   ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  #   ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  #   ggplot2::geom_point(alpha = 0.5, size = 1.2) +
  #   ggrepel::geom_text_repel(
  #     data          = label_df,
  #     ggplot2::aes(label = GENENAME),
  #     size          = 2.8,
  #     max.overlaps  = 30,
  #     segment.color = "grey40",
  #     show.legend   = FALSE
  #   ) +
  #   ggplot2::scale_color_manual(values = pal, name = "Significant in") +
  #   ggplot2::labs(
  #     title    = paste0("LFC Scatter: ", name_a, " vs ", name_b),
  #     subtitle = paste0("Contrast: ", contrast_name,
  #                       "  |  Pearson r = ", round(r_val, 3),
  #                       "  |  padj < ", padj_thr),
  #     x        = paste0("log2FC (", name_a, ")"),
  #     y        = paste0("log2FC (", name_b, ")")
  #   ) +
  #   ggplot2::theme_bw(base_size = 13) +
  #   ggplot2::theme(legend.position = "bottom")
}


# ============================================================
# SECTION D — Plot: Venn diagram of significant DEGs
# ============================================================

#' Venn diagram of significant DEG gene sets
#'
plot_deg_venn <- function(merged, name_a, name_b,
                          contrast_name,
                          padj_thr = 0.05,
                          lfc_thr  = 0) {
  
  sig_a <- merged$GENEID[
    !is.na(merged[[paste0("padj.", name_a)]]) &
      merged[[paste0("padj.", name_a)]] < padj_thr &
      abs(merged[[paste0("log2FoldChange.", name_a)]]) > lfc_thr
  ]
  sig_b <- merged$GENEID[
    !is.na(merged[[paste0("padj.", name_b)]]) &
      merged[[paste0("padj.", name_b)]] < padj_thr &
      abs(merged[[paste0("log2FoldChange.", name_b)]]) > lfc_thr
  ]
  
  n_only_a  <- length(setdiff(sig_a, sig_b))
  n_only_b  <- length(setdiff(sig_b, sig_a))
  n_both    <- length(intersect(sig_a, sig_b))
  
  # Use base graphics via VennDiagram package
  vp <- VennDiagram::draw.pairwise.venn(
    area1         = length(sig_a),
    area2         = length(sig_b),
    cross.area    = n_both,
    category      = c(name_a, name_b),
    fill          = c("#377EB8", "#FF7F00"),
    alpha         = 0.4,
    cex           = 1.4,
    cat.cex       = 1.4,
    cat.fontface  = "bold",
    print.mode    = "raw",
    sigdigs       = 3,
    ind           = FALSE     # return grob, don't draw yet
  )
  
  grid::grid.newpage()
  grid::grid.draw(vp)
  grid::grid.text(
    paste0("Contrast: ", contrast_name, "  |  padj < ", padj_thr),
    x    = 0.5, y = 0.05,
    gp   = grid::gpar(fontsize = 11, col = "grey30")
  )
  
  invisible(list(only_a = n_only_a, only_b = n_only_b, both = n_both))
}


# ============================================================
# SECTION E — Plot: Heatmap of shared significant DEGs
# ============================================================

#' Heatmap of genes significant in BOTH datasets
#'
#' Requires the normalised (vsd) matrices from each dataset.
#' If vsd matrices are not available, falls back to a dot-plot of LFC values.
#'
plot_shared_deg_heatmap <- function(merged,
                                    name_a, name_b,
                                    contrast_name,
                                    padj_thr  = 0.05,
                                    lfc_thr   = 0,
                                    top_n     = 50,
                                    vsd_mat_a = NULL,
                                    vsd_mat_b = NULL,
                                    coldata_a = NULL,
                                    coldata_b = NULL) {
  
  
  # ---- Shared significant genes ----
  lfc_a   <- paste0("log2FoldChange.", name_a)
  lfc_b   <- paste0("log2FoldChange.", name_b)
  padj_a  <- paste0("padj.", name_a)
  padj_b  <- paste0("padj.", name_b)
  
  both <- merged[
    !is.na(merged[[padj_a]]) & merged[[padj_a]] < padj_thr &
      !is.na(merged[[padj_b]]) & merged[[padj_b]] < padj_thr &
      abs(merged[[lfc_a]]) > lfc_thr &
      abs(merged[[lfc_b]]) > lfc_thr,
  ]
  
  if (nrow(both) == 0) {
    message("  [Heatmap] No shared DEGs found for contrast: ", contrast_name)
    return(invisible(NULL))
  }
  
  # Rank by mean |LFC| across both datasets
  both$mean_abs_lfc <- (abs(both[[lfc_a]]) + abs(both[[lfc_b]])) / 2
  both <- both[order(both$mean_abs_lfc, decreasing = TRUE), ]
  both <- utils::head(both, top_n)
  
  message("  [Heatmap] Shared DEGs: ", nrow(both), " (showing top ", min(nrow(both), top_n), ")")
  
  # ---- If vsd matrices provided: expression heatmap ----
  if (!is.null(vsd_mat_a) && !is.null(vsd_mat_b)) {
    
    # Subset to shared genes available in both matrices
    genes_in_a <- intersect(both$GENEID, rownames(vsd_mat_a))
    genes_in_b <- intersect(both$GENEID, rownames(vsd_mat_b))
    genes_use  <- intersect(genes_in_a, genes_in_b)
    
    if (length(genes_use) == 0) {
      message("  [Heatmap] No shared genes found in vsd matrices; falling back to LFC dot-plot.")
      vsd_mat_a <- NULL  # trigger fallback
    } else {
      mat_a  <- vsd_mat_a[genes_use, , drop = FALSE]
      mat_b  <- vsd_mat_b[genes_use, , drop = FALSE]
      
      # z-score per row within each dataset
      # (lignes constantes -> NaN : remplacées par 0 pour ne pas casser le clustering)
      zscore <- function(m) { z <- t(scale(t(m))); z[!is.finite(z)] <- 0; z }
      mat_a  <- zscore(mat_a)
      mat_b  <- zscore(mat_b)
      
      colnames(mat_a) <- paste0(name_a, "_", colnames(mat_a))
      colnames(mat_b) <- paste0(name_b, "_", colnames(mat_b))
      
      # Combine matrices side-by-side
      combined <- cbind(mat_a, mat_b)
      
      # Row labels: gene symbol
      sym_map      <- setNames(both$GENENAME, both$GENEID)
      row_lab      <- unname(sym_map[rownames(combined)])
      row_lab[is.na(row_lab) | row_lab == ""] <- rownames(combined)[is.na(row_lab) | row_lab == ""]
      rownames(combined) <- make.unique(row_lab)   # pheatmap refuse les noms de lignes dupliqués
      
      # Column annotation
      ann_df <- data.frame(
        Dataset   = c(rep(name_a, ncol(mat_a)), rep(name_b, ncol(mat_b))),
        row.names = colnames(combined)
      )
      ann_colors <- list(Dataset = setNames(c("#377EB8", "#FF7F00"), c(name_a, name_b)))
      
      # Draw heatmap
      pheatmap::pheatmap(
        combined,
        annotation_col   = ann_df,
        annotation_colors = ann_colors,
        cluster_cols     = FALSE,
        cluster_rows     = nrow(combined) >= 2,   # hclust exige >= 2 gènes (1 seul DEG partagé -> pas de clustering)
        show_colnames    = TRUE,
        fontsize_row     = max(5, 10 - nrow(combined) %/% 10),
        fontsize_col     = 8,
        scale            = "none",
        main             = paste0("Shared DEGs (z-scored VST)\n",
                                  contrast_name, "  |  ",
                                  name_a, " || ", name_b),
        color            = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100)
      )
      
      
      return(invisible(genes_use))
    }
  }
  
  # ---- Fallback: LFC dot-plot ----
  # Disambiguate duplicate gene symbols (same symbol, different ENSG IDs)
  sym_count  <- table(both$GENENAME)
  both$label <- ifelse(
    sym_count[both$GENENAME] > 1,
    paste0(both$GENENAME, " (", both$GENEID, ")"),
    both$GENENAME
  )
  
  plot_df <- data.frame(
    Gene  = both$label,
    LFC_A = both[[lfc_a]],
    LFC_B = both[[lfc_b]]
  )
  # Reshape to long
  plot_long <- tidyr::pivot_longer(
    plot_df,
    cols      = c(LFC_A, LFC_B),
    names_to  = "Dataset",
    values_to = "LFC"
  )
  plot_long$Dataset <- ifelse(plot_long$Dataset == "LFC_A", name_a, name_b)
  # unique() guards against any remaining duplicates after disambiguation
  plot_long$Gene    <- factor(plot_long$Gene,
                              levels = rev(unique(plot_df$Gene)))
  
  ggplot2::ggplot(plot_long, ggplot2::aes(x = LFC, y = Gene, color = Dataset)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c("#377EB8", "#FF7F00")) +
    ggplot2::labs(
      title    = paste0("Shared DEGs — LFC: ", name_a, " vs ", name_b),
      subtitle = paste0("Contrast: ", contrast_name, "  |  padj < ", padj_thr),
      x        = "log2FoldChange",
      y        = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}


# ============================================================
# SECTION F — Summary table of shared DEGs
# ============================================================

#' Write a CSV of shared significant DEGs with LFC and padj from both datasets
#'
# save_shared_deg_table <- function(merged, name_a, name_b,
#                                   contrast_name,
#                                   out_dir,
#                                   padj_thr = 0.05,
#                                   lfc_thr  = 0) {
#   
#   padj_a <- paste0("padj.", name_a)
#   padj_b <- paste0("padj.", name_b)
#   lfc_a  <- paste0("log2FoldChange.", name_a)
#   lfc_b  <- paste0("log2FoldChange.", name_b)
#   
#   both <- merged[
#     !is.na(merged[[padj_a]]) & merged[[padj_a]] < padj_thr &
#       !is.na(merged[[padj_b]]) & merged[[padj_b]] < padj_thr &
#       abs(merged[[lfc_a]]) > lfc_thr &
#       abs(merged[[lfc_b]]) > lfc_thr,
#   ]
#   
#   both$mean_abs_lfc <- (abs(both[[lfc_a]]) + abs(both[[lfc_b]])) / 2
#   both$direction_concordant <- sign(both[[lfc_a]]) == sign(both[[lfc_b]])
#   both <- both[order(both$mean_abs_lfc, decreasing = TRUE), ]
#   
#   out_file <- file.path(out_dir, paste0(contrast_name, "_shared_DEGs_",
#                                         name_a, "_vs_", name_b, ".csv"))
#   utils::write.csv(both, out_file, row.names = FALSE)
#   message("  Saved shared DEG table → ", out_file)
#   invisible(both)
# }
# 



save_shared_deg_table <- function(merged, name_a, name_b,
                                  contrast_name,
                                  out_dir,
                                  padj_thr = 0.05,
                                  lfc_thr  = 0) {
  
  padj_a <- paste0("padj.", name_a)
  padj_b <- paste0("padj.", name_b)
  lfc_a  <- paste0("log2FoldChange.", name_a)
  lfc_b  <- paste0("log2FoldChange.", name_b)
  
  sig_a <- !is.na(merged[[padj_a]]) & merged[[padj_a]] < padj_thr &
    abs(merged[[lfc_a]]) > lfc_thr
  sig_b <- !is.na(merged[[padj_b]]) & merged[[padj_b]] < padj_thr &
    abs(merged[[lfc_b]]) > lfc_thr
  
  # Catégorie de chaque gène
  merged$overlap_category <- dplyr::case_when(
    sig_a & sig_b                         ~ "Both",
    sig_a & !sig_b                        ~ paste0("Only_", name_a),
    !sig_a & sig_b                        ~ paste0("Only_", name_b),
    TRUE                                  ~ "Neither"
  )
  
  # Garder uniquement les gènes sig dans au moins un dataset
  out <- merged[merged$overlap_category != "Neither", ]
  
  # mean_abs_lfc et concordance (NA pour les "Only")
  out$mean_abs_lfc <- (abs(out[[lfc_a]]) + abs(out[[lfc_b]])) / 2
  out$direction_concordant <- ifelse(
    out$overlap_category == "Both",
    sign(out[[lfc_a]]) == sign(out[[lfc_b]]),
    NA
  )
  
  out <- out[order(out$overlap_category, -out$mean_abs_lfc), ]
  
  out_file <- file.path(out_dir, paste0(contrast_name, "_shared_DEGs_",
                                        name_a, "_vs_", name_b, ".csv"))
  utils::write.csv(out, out_file, row.names = FALSE)
  message("  Saved DEG table (Both + Only) → ", out_file)
  
  # Retourner uniquement les "Both" pour le calcul du % concordant dans le summary
  invisible(out[out$overlap_category == "Both", ])
}



# ============================================================
# SECTION G — Master comparison runner
# ============================================================

#' Run all cross-dataset comparisons
#'
#' @param all_dataset_results Named list of lists:
#'   list(
#'     iPSC     = list(control_vs_FA = df, control_vs_FA_edited = df, FA_edited_vs_FA = df),
#'     iNeurons = list(...)
#'   )
#'   Each df needs: GENEID, GENENAME, log2FoldChange, padj
#'
#' @param base_out  Root output directory (BASE_OUT in your config)
#' @param name_a    Name of the first dataset  (default: first element)
#' @param name_b    Name of the second dataset (default: second element)
#' @param padj_thr  Adjusted p-value threshold
#' @param lfc_thr   Absolute LFC threshold
#' @param top_label Number of genes to label on scatter plots
#' @param vsd_list  Optional: named list of vsd matrices (same names as dataset results)
#' @param coldata_list Optional: named list of colData data.frames
#'
run_cross_dataset_comparison <- function(all_dataset_results,
                                         base_out,
                                         name_a       = names(all_dataset_results)[1],
                                         name_b       = names(all_dataset_results)[2],
                                         padj_thr     = 0.05,
                                         lfc_thr      = 0,
                                         top_label    = 20,
                                         genes_highlight = NULL,   # e.g. GENES_OF_INTEREST
                                         vsd_list     = NULL,
                                         coldata_list = NULL) {
  
  message("\n==========================================")
  message("CROSS-DATASET COMPARISON: ", name_a, " vs ", name_b)
  message("==========================================\n")
  
  results_a <- all_dataset_results[[name_a]]
  results_b <- all_dataset_results[[name_b]]
  
  if (is.null(results_a)) stop("Dataset '", name_a, "' not found in all_dataset_results")
  if (is.null(results_b)) stop("Dataset '", name_b, "' not found in all_dataset_results")
  
  # Contrasts present in both datasets
  shared_contrasts <- intersect(names(results_a), names(results_b))
  
  if (length(shared_contrasts) == 0) {
    warning("No shared contrasts between ", name_a, " and ", name_b)
    return(invisible(NULL))
  }
  
  # Output directory
  cross_dir <- file.path(base_out, paste0("CROSS_", name_a, "_vs_", name_b))
  dir.create(cross_dir, showWarnings = FALSE, recursive = TRUE)
  
  summary_rows <- list()
  
  for (cn in shared_contrasts) {
    
    message("\n--- Cross-contrast: ", cn, " ---")
    
    cn_dir <- file.path(cross_dir, cn)
    dir.create(cn_dir, showWarnings = FALSE, recursive = TRUE)
    
    df_a <- results_a[[cn]]
    df_b <- results_b[[cn]]
    
    # Sanity check
    if (is.null(df_a) || nrow(df_a) == 0) { message("  Skipping: empty results for ", name_a); next }
    if (is.null(df_b) || nrow(df_b) == 0) { message("  Skipping: empty results for ", name_b); next }
    
    # ── Merge ──
    merged <- merge_deg_results(df_a, df_b, name_a, name_b)
    
    if (nrow(merged) == 0) {
      message("  No overlapping genes between datasets for contrast: ", cn)
      next
    }
    
    message("  Genes in common: ", nrow(merged))
    
    # ── LFC Scatter ──
    p_scatter <- plot_lfc_scatter(
      merged, name_a, name_b, cn,
      padj_thr        = padj_thr,
      lfc_thr         = lfc_thr,
      top_label       = top_label,
      genes_highlight = genes_highlight
    )
    ggplot2::ggsave(
      file.path(cn_dir, paste0(cn, "_LFC_scatter_", name_a, "_vs_", name_b, ".pdf")),
      plot   = p_scatter,
      width  = 10, height = 9
    )
    
    # ── Venn ──
    pdf(file.path(cn_dir, paste0(cn, "_Venn_", name_a, "_vs_", name_b, ".pdf")),
        width = 7, height = 6)
    venn_counts <- plot_deg_venn(merged, name_a, name_b, cn, padj_thr, lfc_thr)
    dev.off()
    
    # ── Heatmap / LFC dot-plot ──
    vsd_a      <- if (!is.null(vsd_list)) vsd_list[[name_a]][[cn]] else NULL
    vsd_b      <- if (!is.null(vsd_list)) vsd_list[[name_b]][[cn]] else NULL
    coldata_a  <- if (!is.null(coldata_list)) coldata_list[[name_a]][[cn]] else NULL
    coldata_b  <- if (!is.null(coldata_list)) coldata_list[[name_b]][[cn]] else NULL
    
    pdf(file.path(cn_dir, paste0(cn, "_shared_heatmap_", name_a, "_vs_", name_b, ".pdf")),
        width = 12, height = 14)
    # tryCatch : une heatmap impossible ne doit plus interrompre toute la comparaison,
    # et dev.off() est toujours exécuté
    tryCatch({
      hm <- plot_shared_deg_heatmap(
        merged, name_a, name_b, cn,
        padj_thr  = padj_thr,
        lfc_thr   = lfc_thr,
        top_n     = 50,
        vsd_mat_a = vsd_a,
        vsd_mat_b = vsd_b,
        coldata_a = coldata_a,
        coldata_b = coldata_b
      )
      if (inherits(hm, "ggplot")) print(hm)     # repli "dot-plot" : renvoie un ggplot à afficher
    }, error = function(e) message("  [WARN] heatmap partagée ignorée (", cn, ") : ", conditionMessage(e)))
    dev.off()
    
    # ── CSV table ──
    shared_tbl <- save_shared_deg_table(
      merged, name_a, name_b, cn, cn_dir, padj_thr, lfc_thr
    )
    
    # ── Collect summary ──
    summary_rows[[cn]] <- data.frame(
      contrast     = cn,
      n_genes_common = nrow(merged),
      n_sig_in_A   = venn_counts$only_a + venn_counts$both,
      n_sig_in_B   = venn_counts$only_b + venn_counts$both,
      n_shared_sig = venn_counts$both,
      pct_concordant = if (nrow(shared_tbl) > 0)
        round(100 * mean(shared_tbl$direction_concordant), 1) else NA
    )
  }
  
  # ── Summary table ──
  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)
    out_summary <- file.path(cross_dir, paste0("SUMMARY_", name_a, "_vs_", name_b, ".csv"))
    utils::write.csv(summary_df, out_summary, row.names = FALSE)
    message("\n  Summary table → ", out_summary)
    print(summary_df)
  }
  
  message("\n✔ Cross-dataset comparison done: ", name_a, " vs ", name_b)
  invisible(summary_rows)
}


# ============================================================
# SECTION H — Patched run_dataset() that captures results
# ============================================================
# Replace the run_dataset() in main.R with this version,
# OR add the two lines marked [ADD] to your existing run_dataset().
#
# This version returns all_results invisibly so you can
# accumulate them in all_dataset_results.
#
# Usage in main.R:
#
#   all_dataset_results <- list()
#
#   for (ds_name in names(DATASETS)) {
#     all_dataset_results[[ds_name]] <- run_dataset(ds_name, DATASETS[[ds_name]])
#   }
#
#   # Then, after the loop:
#   run_cross_dataset_comparison(
#     all_dataset_results  = all_dataset_results,
#     base_out             = BASE_OUT,
#     name_a               = "iPSC",       # first dataset name
#     name_b               = "iNeurons",   # second dataset name
#     padj_thr             = PADJ_THRESHOLD,
#     lfc_thr              = LFC_THRESHOLD,
#     genes_highlight      = GENES_OF_INTEREST
#   )
#
# CHANGES vs original run_dataset():
#   - added `return(invisible(all_results))` at the end   [ADD]
# ============================================================

# run_dataset <- function(dataset_name, dataset) {
#   ...
#   # (everything unchanged until the very end)
#   ...
#   message("\nFinished dataset: ", dataset_name)
#   return(invisible(all_results))   # [ADD THIS LINE]
# }