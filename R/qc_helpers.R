# ============================================================
# R/qc_helpers.R — Contrôles qualité
# ============================================================
#  - run_dataset_qc()          : PCA globale (tous les échantillons), corrélations,
#                                 contrôle du sexe, marqueurs de différenciation
#  - plot_deseq_diagnostics()  : dispersions, MA-plot, histogramme des p-values
# ============================================================

# Petit utilitaire : ouvre un PDF, évalue expr, ferme toujours le device
with_pdf <- function(file, expr, width = 7, height = 6) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  tryCatch(expr, error = function(e)
    message("  [WARN] ", basename(file), " : ", conditionMessage(e)))
  invisible(NULL)
}


# ------------------------------------------------------------
# QC global d'un dataset (tous les échantillons ensemble)
# ------------------------------------------------------------
#' @param dataset_name  nom du dataset
#' @param dataset       élément de DATASETS (samples, conditions, sex, path, [line])
#' @param tx2gene       data.frame TXNAME, GENEID, GENENAME
#' @param out_dir       dossier du dataset ; les sorties vont dans <out_dir>/QC
#' @param markers       liste nommée de marqueurs (DIFF_MARKERS)
run_dataset_qc <- function(
    dataset_name, dataset, tx2gene, out_dir,
    min_counts  = get0("MIN_COUNTS",  ifnotfound = 10),
    min_samples = get0("MIN_SAMPLES", ifnotfound = 3),
    markers     = get0("DIFF_MARKERS", ifnotfound = list())
) {

  message("\n--- QC global : ", dataset_name, " ---")
  qc_dir <- file.path(out_dir, "QC")
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

  samples <- dataset$samples
  files   <- file.path(dataset$path, samples, "quant.sf")
  names(files) <- samples

  if (!all(file.exists(files))) {
    message("  [WARN] quant.sf manquants -> QC ignoré : ",
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
    sex       = factor(dataset$sex),
    line      = factor(if (!is.null(dataset$line)) dataset$line else samples),
    row.names = samples
  )

  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = colData, design = ~ 1)
  dds <- filter_low_counts(dds, min_counts, min_samples)
  dds <- DESeq2::estimateSizeFactors(dds)

  vsd <- tryCatch(
    DESeq2::vst(dds, blind = TRUE),
    error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  )
  mat <- SummarizedExperiment::assay(vsd)

  # ---------------- PCA globale ----------------
  rv  <- apply(mat, 1, stats::var)
  sel <- order(rv, decreasing = TRUE)[seq_len(min(500, nrow(mat)))]
  pca <- stats::prcomp(t(mat[sel, ]))
  pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  n_pc <- ncol(pca$x)

  scores <- data.frame(pca$x[, seq_len(min(4, n_pc)), drop = FALSE])
  scores$sample    <- rownames(scores)
  scores$condition <- colData[scores$sample, "condition"]
  scores$sex       <- colData[scores$sample, "sex"]
  scores$line      <- colData[scores$sample, "line"]

  # lignes pointillées entre échantillons de la même lignée (paires isogéniques)
  multi_line <- names(which(table(scores$line) > 1))
  scores_pair <- scores[scores$line %in% multi_line, ]

  plot_pc <- function(px, py) {
    xs <- paste0("PC", px); ys <- paste0("PC", py)
    g <- ggplot2::ggplot(scores, ggplot2::aes(x = .data[[xs]], y = .data[[ys]]))
    if (nrow(scores_pair) > 0) {
      g <- g + ggplot2::geom_path(
        data = scores_pair, ggplot2::aes(group = line),
        colour = "grey60", linetype = "dashed"
      )
    }
    g +
      ggplot2::geom_point(ggplot2::aes(colour = condition, shape = sex), size = 4) +
      ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 3, max.overlaps = Inf) +
      ggplot2::labs(
        title    = paste0("PCA globale — ", dataset_name),
        subtitle = "top 500 genes variables (VST) ; pointilles = meme lignee",
        x = paste0(xs, " (", pct[px], "%)"),
        y = paste0(ys, " (", pct[py], "%)")
      ) +
      ggplot2::theme_bw(base_size = 12)
  }

  with_pdf(file.path(qc_dir, paste0(dataset_name, "_PCA_global.pdf")), {
    print(plot_pc(1, 2))
    if (n_pc >= 3) print(plot_pc(2, 3))
    if (n_pc >= 4) print(plot_pc(3, 4))
  })
  utils::write.csv(scores, file.path(qc_dir, paste0(dataset_name, "_PCA_global_scores.csv")),
                   row.names = FALSE)

  # ---------------- Corrélation entre échantillons ----------------
  ann <- data.frame(condition = colData$condition, line = colData$line,
                    sex = colData$sex, row.names = rownames(colData))
  with_pdf(file.path(qc_dir, paste0(dataset_name, "_sample_correlation.pdf")), {
    pheatmap::pheatmap(
      stats::cor(mat), annotation_col = ann, annotation_row = ann,
      main = paste0("Correlation entre echantillons (VST) — ", dataset_name)
    )
  }, width = 9, height = 8)

  # ---------------- Contrôle du sexe ----------------
  sx <- check_sex_expression(dds, tx2gene)
  if (nrow(sx) > 0) {
    sx_out <- data.frame(gene = rownames(sx), round(sx, 2), check.names = FALSE)
    utils::write.csv(sx_out, file.path(qc_dir, paste0(dataset_name, "_sex_check.csv")),
                     row.names = FALSE)
    if (nrow(sx) >= 2) {
      with_pdf(file.path(qc_dir, paste0(dataset_name, "_sex_check.pdf")), {
        ann_s <- data.frame(sex_declare = colData$sex, row.names = rownames(colData))
        pheatmap::pheatmap(
          log2(sx + 1), annotation_col = ann_s,
          cluster_rows = FALSE, cluster_cols = FALSE,
          main = paste0("Genes lies au sexe (log2 comptages normalises + 1) — ", dataset_name)
        )
      }, width = 8, height = 5)
    }
  } else {
    message("  [INFO] aucun gène du sexe retrouvé (filtre ou annotation)")
  }

  # ---------------- Marqueurs de différenciation ----------------
  if (length(markers) > 0) {
    gm <- unique(tx2gene[, c("GENEID", "GENENAME")])
    nc <- DESeq2::counts(dds, normalized = TRUE)
    rownames(nc) <- sub("\\..*", "", rownames(nc))

    mk <- data.frame(
      group = rep(names(markers), lengths(markers)),
      gene  = unlist(markers, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    mk <- mk[!duplicated(mk$gene), ]
    mk$GENEID <- gm$GENEID[match(mk$gene, gm$GENENAME)]
    absent <- mk$gene[is.na(mk$GENEID) | !(mk$GENEID %in% rownames(nc))]
    mk <- mk[!is.na(mk$GENEID) & mk$GENEID %in% rownames(nc), ]

    if (length(absent) > 0)
      message("  Marqueurs non exprimés/filtrés : ", paste(absent, collapse = ", "))

    if (nrow(mk) >= 2) {
      m_log <- log2(nc[mk$GENEID, , drop = FALSE] + 1)
      rownames(m_log) <- mk$gene

      utils::write.csv(
        data.frame(group = mk$group, gene = mk$gene, round(m_log, 3), check.names = FALSE),
        file.path(qc_dir, paste0(dataset_name, "_differentiation_markers.csv")),
        row.names = FALSE
      )

      ann_row <- data.frame(group = mk$group, row.names = mk$gene)
      gaps <- cumsum(rle(mk$group)$lengths)
      gaps <- gaps[-length(gaps)]

      with_pdf(file.path(qc_dir, paste0(dataset_name, "_differentiation_markers.pdf")), {
        pheatmap::pheatmap(
          m_log, annotation_col = ann[, c("condition", "line"), drop = FALSE],
          annotation_row = ann_row,
          cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = gaps,
          main = paste0("Marqueurs (log2 comptages norm. + 1) — ", dataset_name)
        )
        pheatmap::pheatmap(
          m_log - rowMeans(m_log), annotation_col = ann[, c("condition", "line"), drop = FALSE],
          annotation_row = ann_row,
          cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = gaps,
          main = paste0("Marqueurs (centres par gene) — ", dataset_name)
        )
      }, width = 9, height = 10)

      # ---------------- Scores de différenciation par échantillon ----------------
      # Score d'un groupe = moyenne des log2(comptages normalisés + 1) de ses marqueurs.
      # Diagnostic (pas une correction) : repère les échantillons atypiques et permet de
      # vérifier si un signal (ex. EMT) suit la composition cellulaire.
      grp <- unique(mk$group)
      sc  <- vapply(grp, function(g) colMeans(m_log[mk$group == g, , drop = FALSE]),
                    numeric(ncol(m_log)))
      if (is.null(dim(sc))) sc <- matrix(sc, nrow = ncol(m_log), dimnames = list(colnames(m_log), grp))
      rownames(sc) <- colnames(m_log)

      # Échantillon "atypique" pour un groupe si z robuste (médiane / MAD) > 2,5 en valeur absolue
      # ET écart à la médiane >= 0,5 unité log2 (évite de signaler du bruit quand la MAD est minuscule)
      rz <- apply(sc, 2, function(x) {
        m <- stats::mad(x)
        if (is.na(m) || m == 0) rep(0, length(x)) else (x - stats::median(x)) / m
      })
      dev <- apply(sc, 2, function(x) x - stats::median(x))
      rownames(rz) <- rownames(sc); rownames(dev) <- rownames(sc)
      is_out <- abs(rz) > 2.5 & abs(dev) >= 0.5
      flags <- apply(is_out, 1, function(z) paste(colnames(is_out)[z], collapse = ";"))

      score_df <- data.frame(
        sample    = rownames(sc),
        condition = colData[rownames(sc), "condition"],
        line      = colData[rownames(sc), "line"],
        round(sc, 2),
        atypique_pour = unname(flags[rownames(sc)]),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      utils::write.csv(score_df, file.path(qc_dir, paste0(dataset_name, "_differentiation_scores.csv")),
                       row.names = FALSE)

      # corrélation (Spearman) de chaque score avec les 3 premières PC de la PCA globale
      npc <- min(3, n_pc)
      pc_cor <- t(vapply(grp, function(g)
        vapply(seq_len(npc), function(k)
          suppressWarnings(stats::cor(sc[, g], pca$x[rownames(sc), k], method = "spearman")),
          numeric(1)), numeric(npc)))
      if (npc == 1) pc_cor <- matrix(pc_cor, ncol = 1)
      colnames(pc_cor) <- paste0("PC", seq_len(npc)); rownames(pc_cor) <- grp
      utils::write.csv(data.frame(groupe = grp, round(pc_cor, 2)),
                       file.path(qc_dir, paste0(dataset_name, "_score_vs_PC.csv")), row.names = FALSE)

      long <- data.frame(
        sample = rep(rownames(sc), times = ncol(sc)),
        groupe = rep(colnames(sc), each = nrow(sc)),
        score  = as.vector(sc),
        condition = rep(score_df$condition, times = ncol(sc)),
        stringsAsFactors = FALSE
      )
      g <- ggplot2::ggplot(long, ggplot2::aes(condition, score, colour = condition)) +
        ggplot2::geom_point(size = 3) +
        ggrepel::geom_text_repel(ggplot2::aes(label = sample), size = 2.5, max.overlaps = Inf,
                                 show.legend = FALSE) +
        ggplot2::facet_wrap(~ groupe, scales = "free_y") +
        ggplot2::labs(title = paste0("Scores de differenciation (moyenne log2 comptages norm. + 1) — ", dataset_name),
                      x = NULL, y = "score") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
      with_pdf(file.path(qc_dir, paste0(dataset_name, "_differentiation_scores.pdf")),
               print(g), width = 11, height = 8)
    }
  }

  message("  QC écrit dans ", qc_dir)
  invisible(list(dds = dds, vsd = vsd, pca = pca))
}


# ------------------------------------------------------------
# Diagnostics DESeq2 par contraste
# ------------------------------------------------------------
#' @param dds          DESeqDataSet (après DESeq())
#' @param res          objet results() NON shrunk (res$res dans run_deseq2_contrast)
#' @param contrast_dir dossier de sortie du contraste
#' @param cn           nom du contraste
plot_deseq_diagnostics <- function(dds, res, contrast_dir, cn, alpha = 0.05) {

  with_pdf(file.path(contrast_dir, paste0(cn, "_diag_dispersion.pdf")), {
    DESeq2::plotDispEsts(dds, main = paste0("Estimation des dispersions — ", cn))
  })

  with_pdf(file.path(contrast_dir, paste0(cn, "_diag_MA.pdf")), {
    DESeq2::plotMA(res, alpha = alpha, ylim = c(-5, 5),
                   main = paste0("MA-plot (LFC non shrunk) — ", cn))
  })

  with_pdf(file.path(contrast_dir, paste0(cn, "_diag_pvalue_hist.pdf")), {
    p <- res$pvalue[!is.na(res$pvalue) & res$baseMean > 1]
    graphics::hist(p, breaks = 50, col = "grey70", border = "white",
                   xlab = "p-value", main = paste0("Distribution des p-values — ", cn))
  })

  invisible(NULL)
}
