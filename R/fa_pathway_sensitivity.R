# ============================================================
# fa_pathway_sensitivity.R
#
# Analyse de sensibilité : recalcule control_vs_FA (iPSC, iNeurons) avec 4
# méthodes différentes et compare, pour chaque gène de FA_PATHWAY_GENES
# (config.R), s'il ressort significatif ou non selon la méthode :
#
#   1. DESeq2_baseline       -> réglages actuels du pipeline (référence)
#   2. DESeq2_sex_adjustTRUE -> sexe ajouté au design (~ sex + condition)
#   3. DESeq2_filtre_souple  -> filtre de comptage assoupli (moins de gènes
#                                exclus avant le test)
#   4. limma_voom            -> méthode alternative à DESeq2 (mêmes comptages
#                                bruts et même design que la baseline)
#
# Pour chaque méthode, deux corrections de multi-tests sont rapportées :
#   - padj_genome : correction BH standard sur TOUT le génome testé (ce que
#     produit nativement DESeq2/limma)
#   - padj_panel  : correction BH recalculée UNIQUEMENT sur les gènes du
#     panel testés par cette méthode (pratique standard et défendable pour
#     un petit panel de gènes pré-spécifiés par hypothèse biologique).
#
# Sortie : results_v4/FA_PATHWAY_GENES/genes_candidats_FA_pathway_sensitivity.xlsx
# ============================================================

# ------------------------------------------------------------
# Ré-exécute control_vs_FA pour un dataset avec des réglages donnés.
# Réutilise run_deseq2_contrast() (R/deseq2_helpers.R), donc identique en
# tout point à ce que fait mainPatchSave.R, seuls sex_adjust/min_counts/
# min_samples changent.
# ------------------------------------------------------------
.fa_sens_run_deseq2 <- function(dataset_name, sex_adjust, min_counts, min_samples,
                                tx2gene, padj_thr) {
  
  dataset    <- DATASETS[[dataset_name]]
  conditions <- dataset$conditions
  cond_levels <- unique(conditions)
  if (!all(c("control", "FA") %in% cond_levels)) return(NULL)
  
  idx <- which(conditions %in% c("control", "FA"))
  samples   <- dataset$samples
  sex       <- dataset$sex
  line      <- dataset$line
  base_path <- dataset$path
  
  files <- file.path(base_path, samples[idx], "quant.sf")
  names(files) <- samples[idx]
  
  excl <- if (!sex_adjust) sexchr_ids else NULL
  
  run_deseq2_contrast(
    files_pair   = files,
    samples_pair = samples[idx],
    conds_pair   = conditions[idx],
    donor_pair   = sex[idx],
    line_pair    = line[idx],
    tx2gene      = tx2gene,
    condA        = "control",
    condB        = "FA",
    min_counts   = min_counts,
    min_samples  = min_samples,
    alpha        = padj_thr,
    sex_adjust   = sex_adjust,
    exclude_genes = excl
  )
}

# ------------------------------------------------------------
# Construit une table standardisée (GENEID, GENENAME, LFC, pvalue) en
# limma-voom, à partir du dds "baseline" (mêmes comptages bruts, même
# design) pour une comparaison apples-to-apples avec DESeq2_baseline.
# Retourne aussi l'objet voom/design/coef (réutilisés par les tests de
# groupe roast/camera, pour ne pas refaire le calcul).
# ------------------------------------------------------------
.fa_sens_run_limma <- function(dds_baseline, tx2gene) {
  
  if (is.null(dds_baseline)) return(NULL)
  if (!requireNamespace("limma", quietly = TRUE) || !requireNamespace("edgeR", quietly = TRUE)) {
    message("  [INFO] limma/edgeR non installés : méthode limma_voom ignorée.")
    return(NULL)
  }
  
  counts_mat <- DESeq2::counts(dds_baseline, normalized = FALSE)
  cd         <- as.data.frame(SummarizedExperiment::colData(dds_baseline))
  design_f   <- DESeq2::design(dds_baseline)
  
  dge <- edgeR::DGEList(counts = counts_mat)
  dge <- edgeR::calcNormFactors(dge)
  
  design_mat <- stats::model.matrix(design_f, data = cd)
  coef_name  <- grep("^condition", colnames(design_mat), value = TRUE)
  if (length(coef_name) == 0) {
    message("  [WARN] limma_voom : coefficient 'condition' introuvable dans le design.")
    return(NULL)
  }
  coef_name <- coef_name[1]
  
  v    <- limma::voom(dge, design_mat)
  fit  <- limma::lmFit(v, design_mat)
  fit  <- limma::eBayes(fit)
  tt   <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  
  gm <- unique(tx2gene[, c("GENEID", "GENENAME")])
  ensg <- sub("\\..*", "", rownames(tt))
  
  table_std <- data.frame(
    GENEID   = ensg,
    GENENAME = gm$GENENAME[match(ensg, gm$GENEID)],
    LFC      = tt$logFC,
    pvalue   = tt$P.Value,
    stringsAsFactors = FALSE
  )
  
  list(table = table_std, voom = v, design = design_mat, coef = coef_name)
}

# ------------------------------------------------------------
# Test de groupe sur le panel FA_PATHWAY_GENES entier (pas gène par gène) :
#   - roast (self-contained) : "ce groupe de gènes bouge-t-il collectivement,
#     indépendamment de ce qui se passe hors du groupe ?" -- c'est le test
#     qui peut détecter un signal quand chaque gène pris seul est trop
#     faible pour être significatif individuellement.
#   - camera (compétitif) : "ce groupe de gènes bouge-t-il PLUS que le reste
#     du génome ?", en tenant compte de la corrélation inter-gènes -- plus
#     conservateur, complémentaire de roast.
# S'appuie sur le même objet voom/design que limma_voom (aucun recalcul).
# ------------------------------------------------------------
.fa_sens_run_geneset_tests <- function(limma_obj, tx2gene, panel_symbols, set_name = "FA_pathway_panel") {
  
  if (is.null(limma_obj)) return(NULL)
  
  v          <- limma_obj$voom
  design_mat <- limma_obj$design
  coef_name  <- limma_obj$coef
  
  gm         <- unique(tx2gene[, c("GENEID", "GENENAME")])
  panel_ensg <- unique(gm$GENEID[gm$GENENAME %in% panel_symbols])
  ensg_voom  <- sub("\\..*", "", rownames(v))
  idx        <- which(ensg_voom %in% panel_ensg)
  
  if (length(idx) < 3) {
    message("  [WARN] Test de groupe : trop peu de gènes du panel présents (", length(idx), ") -- ignoré.")
    return(NULL)
  }
  
  gene_index <- setNames(list(idx), set_name)
  
  roast_res <- tryCatch(
    limma::mroast(v, index = gene_index, design = design_mat, contrast = coef_name, nrot = 9999),
    error = function(e) { message("  [WARN] roast échoué : ", conditionMessage(e)); NULL }
  )
  camera_res <- tryCatch(
    limma::camera(v, index = gene_index, design = design_mat, contrast = coef_name),
    error = function(e) { message("  [WARN] camera échoué : ", conditionMessage(e)); NULL }
  )
  
  list(n_genes = length(idx), roast = roast_res, camera = camera_res)
}

# ------------------------------------------------------------
# Accesseur défensif : renvoie toujours une valeur de longueur 1 (NA si la
# colonne n'existe pas, si le data.frame est NULL/vide, ou si l'accès
# échoue pour une raison quelconque). Évite le classique "arguments imply
# differing number of rows" quand une colonne attendue (ex. PValue.Mixed)
# n'existe pas sous ce nom exact selon la version de limma installée.
# ------------------------------------------------------------
.fa_sens_safe_val <- function(df, col) {
  out <- tryCatch({
    if (is.null(df) || !(col %in% colnames(df)) || nrow(df) == 0) NA
    else {
      v <- df[[col]][1]
      if (length(v) == 0) NA else v
    }
  }, error = function(e) NA)
  out
}

# ------------------------------------------------------------
# Standardise une sortie DESeq2 (res_df) vers (GENEID, GENENAME, LFC, pvalue)
# ------------------------------------------------------------
.fa_sens_std_deseq2 <- function(res_df) {
  if (is.null(res_df)) return(NULL)
  data.frame(
    GENEID   = res_df$GENEID,
    GENENAME = res_df$GENENAME,
    LFC      = res_df$log2FoldChange,
    pvalue   = res_df$pvalue,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Fonction principale
# ------------------------------------------------------------
run_fa_pathway_sensitivity <- function(
    base_out,
    tx2gene,
    gene_list           = FA_PATHWAY_GENES,
    datasets_to_test    = c("iPSC", "iNeurons"),
    padj_thr            = PADJ_THRESHOLD,
    lfc_thr             = LFC_THRESHOLD,
    min_counts_relaxed  = 5,
    min_samples_relaxed = 2,
    out_subdir          = "FA_PATHWAY_GENES",
    file_name            = "genes_candidats_FA_pathway_sensitivity.xlsx"
) {
  
  out_dir <- file.path(base_out, out_subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  gene_df <- do.call(rbind, lapply(names(gene_list), function(cat) {
    v <- gene_list[[cat]]
    data.frame(categorie = cat, libelle = names(v), symbole = unname(v),
               stringsAsFactors = FALSE)
  }))
  panel_symbols <- unique(gene_df$symbole)
  
  all_rows    <- list()
  summary_rows <- list()
  geneset_rows <- list()
  
  for (ds in datasets_to_test) {
    
    message("\n--- Sensibilité FA_PATHWAY_GENES : ", ds, " (control_vs_FA) ---")
    
    dataset_sex_adjust <- isTRUE(DATASETS[[ds]]$sex_adjust)
    
    bas <- tryCatch(
      .fa_sens_run_deseq2(ds, sex_adjust = dataset_sex_adjust,
                          min_counts = MIN_COUNTS, min_samples = MIN_SAMPLES,
                          tx2gene = tx2gene, padj_thr = padj_thr),
      error = function(e) { message("  [WARN] baseline échouée : ", conditionMessage(e)); NULL }
    )
    sx <- tryCatch(
      .fa_sens_run_deseq2(ds, sex_adjust = TRUE,
                          min_counts = MIN_COUNTS, min_samples = MIN_SAMPLES,
                          tx2gene = tx2gene, padj_thr = padj_thr),
      error = function(e) { message("  [WARN] sex_adjust=TRUE échouée : ", conditionMessage(e)); NULL }
    )
    rf <- tryCatch(
      .fa_sens_run_deseq2(ds, sex_adjust = dataset_sex_adjust,
                          min_counts = min_counts_relaxed, min_samples = min_samples_relaxed,
                          tx2gene = tx2gene, padj_thr = padj_thr),
      error = function(e) { message("  [WARN] filtre assoupli échoué : ", conditionMessage(e)); NULL }
    )
    lv <- tryCatch(
      .fa_sens_run_limma(if (!is.null(bas)) bas$dds else NULL, tx2gene),
      error = function(e) { message("  [WARN] limma-voom échoué : ", conditionMessage(e)); NULL }
    )
    
    # ---- Test de groupe (roast self-contained + camera compétitif) ----
    # Réutilise le même objet voom/design que limma_voom -- pas de recalcul.
    # Entièrement protégé par tryCatch : une erreur ici (ex. nom de colonne
    # différent selon la version de limma) ne doit jamais faire perdre les
    # résultats gène-par-gène déjà calculés pour les autres méthodes.
    tryCatch({
      gs <- .fa_sens_run_geneset_tests(lv, tx2gene, panel_symbols)
      if (!is.null(gs)) {
        roast_row  <- if (!is.null(gs$roast))  as.data.frame(gs$roast)  else NULL
        camera_row <- if (!is.null(gs$camera)) as.data.frame(gs$camera) else NULL
        geneset_rows[[length(geneset_rows) + 1]] <- data.frame(
          Dataset             = ds,
          N_genes_testes      = gs$n_genes,
          roast_Direction     = .fa_sens_safe_val(roast_row, "Direction"),
          roast_PValue        = .fa_sens_safe_val(roast_row, "PValue"),
          roast_FDR           = .fa_sens_safe_val(roast_row, "FDR"),
          roast_PValue_Mixed  = .fa_sens_safe_val(roast_row, "PValue.Mixed"),
          camera_PValue       = .fa_sens_safe_val(camera_row, "PValue"),
          camera_FDR          = .fa_sens_safe_val(camera_row, "FDR"),
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      message("  [WARN] test de groupe échoué (ignoré, le reste continue) : ", conditionMessage(e))
    })
    
    methods <- list(
      DESeq2_baseline       = .fa_sens_std_deseq2(if (!is.null(bas)) bas$res_df_full else NULL),
      DESeq2_sex_adjustTRUE = .fa_sens_std_deseq2(if (!is.null(sx))  sx$res_df_full  else NULL),
      DESeq2_filtre_souple  = .fa_sens_std_deseq2(if (!is.null(rf))  rf$res_df_full  else NULL),
      limma_voom            = if (!is.null(lv)) lv$table else NULL
    )
    
    for (meth_name in names(methods)) {
      
      tbl <- methods[[meth_name]]
      if (is.null(tbl)) next
      
      # Correction BH sur TOUT le génome testé par cette méthode (déjà standard
      # pour DESeq2 -> padj ; recalculée ici de façon identique pour limma afin
      # que les deux soient comparables).
      tbl$padj_genome <- stats::p.adjust(tbl$pvalue, method = "BH")
      
      # Correction BH restreinte au panel : seulement les gènes du panel
      # effectivement testés (pvalue non NA) par cette méthode.
      panel_rows <- tbl[!is.na(tbl$GENENAME) & tbl$GENENAME %in% panel_symbols & !is.na(tbl$pvalue), ]
      panel_rows$padj_panel <- stats::p.adjust(panel_rows$pvalue, method = "BH")
      
      n_sig_genome <- sum(panel_rows$padj_genome < padj_thr & abs(panel_rows$LFC) >= lfc_thr, na.rm = TRUE)
      n_sig_panel  <- sum(panel_rows$padj_panel  < padj_thr & abs(panel_rows$LFC) >= lfc_thr, na.rm = TRUE)
      
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Dataset = ds, Methode = meth_name,
        Genes_panel_testes = nrow(panel_rows),
        Sig_padj_genome    = n_sig_genome,
        Sig_padj_panel     = n_sig_panel,
        stringsAsFactors = FALSE
      )
      
      for (i in seq_len(nrow(gene_df))) {
        sym <- gene_df$symbole[i]
        hit <- panel_rows[panel_rows$GENENAME == sym, ]
        if (nrow(hit) == 0) next
        hit <- hit[1, ]
        all_rows[[length(all_rows) + 1]] <- data.frame(
          Categorie    = gene_df$categorie[i],
          Nom_fourni   = gene_df$libelle[i],
          Symbole_HGNC = sym,
          Dataset      = ds,
          Methode      = meth_name,
          LFC          = round(hit$LFC, 3),
          pvalue       = hit$pvalue,
          padj_genome  = hit$padj_genome,
          padj_panel   = hit$padj_panel,
          sig_genome   = !is.na(hit$padj_genome) & hit$padj_genome < padj_thr & abs(hit$LFC) >= lfc_thr,
          sig_panel    = !is.na(hit$padj_panel)  & hit$padj_panel  < padj_thr & abs(hit$LFC) >= lfc_thr,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(all_rows) == 0 && length(geneset_rows) == 0) {
    message("  [WARN] FA_PATHWAY_GENES sensitivity : aucun résultat produit.")
    return(invisible(NULL))
  }
  if (length(all_rows) == 0) {
    message("  [WARN] Aucun résultat gène-par-gène (seul le test de groupe a produit des résultats).")
  }
  
  df_long <- do.call(rbind, all_rows)
  df_summary <- do.call(rbind, summary_rows)
  df_geneset <- if (length(geneset_rows) > 0) do.call(rbind, geneset_rows) else NULL
  
  legende <- data.frame(
    Cle = c("Objectif", "DESeq2_baseline", "DESeq2_sex_adjustTRUE", "DESeq2_filtre_souple",
            "limma_voom", "padj_genome", "padj_panel",
            "Test de groupe - roast", "Test de groupe - camera"),
    Valeur = c(
      "Compare 4 methodes statistiques sur control_vs_FA (iPSC, iNeurons) pour le panel FA_PATHWAY_GENES.",
      paste0("Reglages actuels du pipeline (min_counts=", MIN_COUNTS, ", min_samples=", MIN_SAMPLES, ", sex_adjust = valeur declaree du dataset dans config.R)."),
      "Meme filtre que la baseline, mais design force a ~ sex + condition (sexe ajoute comme covariable).",
      paste0("Meme design que la baseline, mais filtre de comptage assoupli (min_counts=", min_counts_relaxed, ", min_samples=", min_samples_relaxed, ")."),
      "Meme comptages bruts et meme design que DESeq2_baseline, mais analyse par limma-voom (edgeR::calcNormFactors + limma::voom/lmFit/eBayes) au lieu de DESeq2.",
      "Correction Benjamini-Hochberg standard sur TOUT le genome teste par cette methode (equivalent au padj habituel).",
      "Correction Benjamini-Hochberg recalculee UNIQUEMENT sur les genes du panel (75 genes) effectivement testes par cette methode -- plus adaptee a un panel pre-specifie.",
      "Test AUTO-CONTENU (limma::mroast, 9999 rotations) : teste si le panel ENTIER (les genes du panel entre eux uniquement) devie collectivement de zero, meme si aucun gene individuel n'est significatif seul. roast_PValue teste la direction dominante (roast_Direction : Up/Down/Mixed) ; roast_PValue_Mixed teste un effet de groupe non-directionnel.",
      "Test COMPETITIF (limma::camera) : teste si le panel bouge PLUS que le reste du genome, en tenant compte de la correlation inter-genes -- plus conservateur que roast, complementaire."
    ),
    stringsAsFactors = FALSE
  )
  
  out_file <- file.path(out_dir, file_name)
  sheets <- list("Comparaison methodes" = df_long, "Resume" = df_summary)
  if (!is.null(df_geneset)) sheets[["Test de groupe"]] <- df_geneset
  sheets[["Legende"]] <- legende
  
  writexl::write_xlsx(sheets, out_file)
  
  message("\n  FA_PATHWAY_GENES sensitivity : classeur ecrit -> ", out_file)
  print(df_summary)
  if (!is.null(df_geneset)) {
    message("\n  Test de groupe (panel entier) :")
    print(df_geneset)
  }
  
  invisible(list(long = df_long, summary = df_summary, geneset = df_geneset))
}
