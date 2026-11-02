# ============================================================
# R/functions/deseq2_helpers.R — DESeq2 wrappers
# ============================================================

#' Run DESeq2 for a pairwise contrast
#'
#' @param files_pair   Named character vector of quant.sf paths
#' @param samples_pair Character vector of sample names
#' @param conds_pair   Character vector of conditions (same order as samples)
#' @param donor_pair   Character vector of sex covariates
#' @param line_pair    Character vector of patient-line IDs (same order as samples).
#'                     If every line is present in BOTH conditions (isogenic pairs,
#'                     e.g. FA vs FA_edited), design becomes ~ line + condition.
#'                     Otherwise: ~ sex + condition (sex_adjust = TRUE) or ~ condition (FALSE).
#' @param tx2gene      tx2gene data.frame
#' @param condA        Reference condition (denominator)
#' @param condB        Numerator condition
#' @param min_counts   Minimum raw count per sample for a gene to be "expressed"
#' @param min_samples  Number of samples that must reach min_counts (gene kept if
#'                     >= min_samples samples have >= min_counts)
#' @param alpha        Significance level passed to results() (independent filtering)
#' @param sex_adjust   TRUE : ~ sex + condition. FALSE : ~ condition et exclusion de exclude_genes
#'                     (contrastes non appariés uniquement ; sans effet si le design est apparié)
#' @param exclude_genes ENSG (sans version) à retirer avant le test quand sex_adjust = FALSE
#'                     (gènes chrX / chrY / XIST, cf. get_sex_chr_genes())
#' @return list: dds, res, resLFC, res_df (colonnes historiques, pour les figures),
#'         res_df_full (res_df + lfc_mle = LFC NON shrunk), vsd, colData, design, paired,
#'         n_excluded (nombre de gènes exclus)
#' 
#' 
#' 



run_deseq2_contrast <- function(
    files_pair, samples_pair, conds_pair, donor_pair,
    tx2gene, condA, condB,
    min_counts  = 10,
    min_samples = 3,
    line_pair   = NULL,
    alpha       = 0.05,
    sex_adjust  = TRUE,
    exclude_genes = NULL
) {
  
  txi <- tximport(
    files_pair,
    type       = "salmon",
    tx2gene    = tx2gene[, c("TXNAME", "GENEID")],
    ignoreTxVersion = TRUE
  )
  
  colData <- data.frame(
    sample    = samples_pair,
    sex       = factor(donor_pair),
    condition = factor(conds_pair, levels = c(condA, condB)),
    row.names = samples_pair
  )
  if (!is.null(line_pair)) colData$line <- factor(line_pair)
  
  # --- Design apparié ? ---
  # Apparié si chaque lignée est présente dans les DEUX conditions du contraste
  # (ex. F4259 et F4259ed pour FA_edited_vs_FA). Pour control_vs_FA les lignées
  # sont toutes différentes -> pas d'appariement possible.
  paired <- FALSE
  if (!is.null(line_pair)) {
    tab    <- table(colData$line, colData$condition)
    paired <- all(rowSums(tab > 0) == 2)
  }
  
  if (paired) {
    # line contient déjà sex (1 sexe par lignée) -> on ne met pas sex en plus
    mm_test <- model.matrix(~ line + condition, data = colData)
    if (qr(mm_test)$rank < ncol(mm_test))
      stop("Design apparié non de rang plein pour ", condA, " vs ", condB)
    design_formula <- ~ line + condition
    message(condA, " vs ", condB, " -> design: ", deparse(design_formula),
            "  (apparié par lignée)")
  } else if (sex_adjust) {
    # --- Sex confondu avec condition pour CE contraste ? ---
    mm_test   <- model.matrix(~ sex + condition, data = colData)
    full_rank <- qr(mm_test)$rank == ncol(mm_test)
    design_formula <- if (full_rank) ~ sex + condition else ~ condition
    message(
      condA, " vs ", condB, " -> design: ", deparse(design_formula),
      if (!full_rank) "  (sex retiré : confondu avec condition dans ce sous-groupe)" else ""
    )
  } else {
    design_formula <- ~ condition
    message(condA, " vs ", condB, " -> design: ", deparse(design_formula),
            "  (sex non ajusté ; gènes chrX/chrY/XIST exclus du test)")
  }
  
  dds <- DESeqDataSetFromTximport(txi, colData = colData, design = design_formula)
  
  # Exclusion des gènes liés au sexe (seulement contrastes non appariés non ajustés sur le sexe)
  n_excluded <- 0L
  if (!paired && !sex_adjust && length(exclude_genes) > 0) {
    drop <- sub("\\..*", "", rownames(dds)) %in% exclude_genes
    n_excluded <- sum(drop)
    message("  Gènes chrX/chrY/XIST exclus : ", n_excluded)
    dds <- dds[!drop, ]
  }
  
  dds <- filter_low_counts(dds, min_counts, min_samples)
  dds <- DESeq(dds)
  
  
  
  
  res <- results(dds, contrast = c("condition", condB, condA), alpha = alpha)
  
  # --- lfcShrink corrigé : apeglm veut "coef", pas "contrast" ---
  coef_name <- paste0("condition_", condB, "_vs_", condA)
  resLFC <- tryCatch(
    if (coef_name %in% resultsNames(dds)) {
      lfcShrink(dds, coef = coef_name, type = "apeglm")
    } else {
      stop("coef '", coef_name, "' absent de resultsNames(dds): ",
           paste(resultsNames(dds), collapse = ", "))
    },
    error = function(e) {
      message("apeglm lfcShrink a échoué (", e$message, ") — fallback 'normal'")
      lfcShrink(dds, contrast = c("condition", condB, condA), type = "normal")
    }
  )
  
  # --- Fusion : LFC/lfcSE shrunk (resLFC) + pvalue/padj/stat du test (res) ---
  # Les deux objets partagent les mêmes rownames dans le même ordre (même dds),
  # mais on aligne quand même explicitement par rowname par sécurité.
  stopifnot(all(rownames(res) == rownames(resLFC)))
  
  res_shrunk_df <- as.data.frame(resLFC)[, c("log2FoldChange", "lfcSE")]
  res_test_df   <- as.data.frame(res)[, c("baseMean", "stat", "pvalue", "padj")]
  lfc_mle_vec   <- res$log2FoldChange   # LFC NON shrunk (MLE)
  
  res_combined <- cbind(res_test_df, res_shrunk_df)
  res_combined <- res_combined[, c("baseMean", "log2FoldChange", "lfcSE",
                                   "stat", "pvalue", "padj")]
  
  # --- Annotate ---
  res_df <- res_combined %>%
    as.data.frame() %>%
    tibble::rownames_to_column("GENEID") %>%
    dplyr::mutate(GENEID = sub("\\..*", "", GENEID)) %>%
    annotate_results(tx2gene) %>%
    dplyr::arrange(padj) %>%
    dplyr::distinct(GENEID, .keep_all = TRUE) %>%
    dplyr::select(GENEID, GENENAME, gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
  
  # Version étendue : + LFC non shrunk (apeglm écrase certains effets réels vers 0 ;
  # à utiliser pour les seuils de LFC et le ratio de rescue). res_df reste inchangé
  # pour ne pas perturber les figures existantes.
  mle_map <- stats::setNames(lfc_mle_vec, sub("\\..*", "", rownames(res)))
  res_df_full <- res_df
  res_df_full$lfc_mle <- unname(mle_map[res_df_full$GENEID])
  res_df_full <- res_df_full[, c("GENEID", "GENENAME", "gene", "baseMean", "log2FoldChange",
                                 "lfc_mle", "lfcSE", "stat", "pvalue", "padj")]
  
  # --- VST ---
  vsd <- tryCatch(
    vst(dds),
    error = function(e) { message("vst failed: ", e$message); NULL }
  )
  
  list(dds = dds, res = res, resLFC = resLFC, res_df = res_df, res_df_full = res_df_full,
       vsd = vsd, colData = colData, design = design_formula, paired = paired,
       n_excluded = n_excluded)
}


#' Filter low-count genes
#'
#' Keeps genes with >= min_counts raw counts in >= min_samples samples.
#'
#' @param dds          DESeqDataSet
#' @param min_counts   minimum raw count per sample
#' @param min_samples  minimum number of samples reaching min_counts
#' @return filtered DESeqDataSet
filter_low_counts <- function(dds, min_counts = 10, min_samples = 3) {
  keep <- rowSums(counts(dds) >= min_counts) >= min_samples
  message("  Filtre comptages : ", sum(keep), " / ", length(keep),
          " gènes conservés (>= ", min_counts, " comptages dans >= ",
          min_samples, " échantillons)")
  dds[keep, ]
}


#' Build ranked gene list for fgsea from a DESeq2 results data.frame
#'
#' Classement sur la statistique de Wald (`stat`). Les symboles dupliqués
#' (plusieurs ENSG pour un même nom) sont dédoublonnés en gardant le |stat| maximal,
#' sinon fgsea reçoit des noms non uniques.
#'
#' @param res_df  data.frame with columns gene and stat
#' @return Named numeric vector sorted decreasingly
make_ranks <- function(res_df) {
  df <- data.frame(gene = res_df$gene, stat = res_df$stat, stringsAsFactors = FALSE)
  df <- df[!is.na(df$gene) & !is.na(df$stat), ]
  df <- df[order(-abs(df$stat)), ]
  n_dup <- sum(duplicated(df$gene))
  if (n_dup > 0) message("  make_ranks : ", n_dup, " symboles dupliqués écartés (|stat| max conservé)")
  df <- df[!duplicated(df$gene), ]
  ranks <- stats::setNames(df$stat, df$gene)
  sort(ranks, decreasing = TRUE)
}


#' Check which genes of interest are present/absent in DESeq2 results
#'
#' @param res_df            data.frame with gene column
#' @param genes_of_interest character vector of gene symbols
#' @return list: present, missing
check_genes_of_interest <- function(res_df, genes_of_interest) {
  present <- intersect(genes_of_interest, res_df$gene)
  missing <- setdiff(genes_of_interest, res_df$gene)
  cat("Genes of interest — total:", length(genes_of_interest),
      "| present:", length(present),
      "| missing:", length(missing), "\n")
  list(present = present, missing = missing)
}


#' Extract sex-chromosome marker expression for QC
#'
#' Les gènes sont retrouvés par SYMBOLE via tx2gene (l'ancienne version utilisait
#' des ENSG codés en dur, dont plusieurs mal attribués).
#'
#' @param dds       DESeqDataSet after estimateSizeFactors
#' @param tx2gene   data.frame TXNAME, GENEID, GENENAME
#' @param sex_genes symboles à extraire
#' @return matrix of normalized counts (rows = gene symbols)
check_sex_expression <- function(
    dds, tx2gene,
    sex_genes = c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")
) {
  gm <- unique(tx2gene[, c("GENEID", "GENENAME")])
  ids <- gm$GENEID[gm$GENENAME %in% sex_genes]
  norm_counts <- counts(dds, normalized = TRUE)
  rownames(norm_counts) <- sub("\\..*", "", rownames(norm_counts))
  out <- norm_counts[rownames(norm_counts) %in% ids, , drop = FALSE]
  rownames(out) <- gm$GENENAME[match(rownames(out), gm$GENEID)]
  out
}
