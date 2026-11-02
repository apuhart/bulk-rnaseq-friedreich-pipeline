options(repos = c(CRAN = "https://cloud.r-project.org"))

# ============================================================
# main.R — RNA-seq pipeline orchestrator
# ============================================================

# ----------------------------
# 0. Sources
# ----------------------------
# Racine du projet = dossier de ce script (fonctionne avec Rscript, source() et RStudio).
# Si besoin, forcer : PROJECT_ROOT <- "/chemin/vers/Bulk-RNAseq"
get_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) == 1 && nzchar(f)) return(dirname(normalizePath(f)))          # Rscript
  for (i in seq_len(sys.nframe())) {                                          # source()
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p)))                          # RStudio
  }
  getwd()
}
PROJECT_ROOT <- get_project_root()
message("PROJECT_ROOT : ", PROJECT_ROOT)

source(file.path(PROJECT_ROOT, "config/packages.R"))
check_packages()   # stoppe avec un message clair si un paquet requis manque
source(file.path(PROJECT_ROOT, "config/config.R"))
source(file.path(PROJECT_ROOT, "config/biomart.R"))
source(file.path(PROJECT_ROOT, "R/deseq2_helpers.R"))
source(file.path(PROJECT_ROOT, "R/fgsea_helpers.R"))
source(file.path(PROJECT_ROOT, "R/plots.R"))
source(file.path(PROJECT_ROOT, "R/cross_dataset_comparison.R"))
source(file.path(PROJECT_ROOT, "R/cross_dataset_ipsc_vs_neurons.R"))
source(file.path(PROJECT_ROOT, "R/cross_dataset_disease_effect_scatter.R"))
source(file.path(PROJECT_ROOT, "R/qc_helpers.R"))
source(file.path(PROJECT_ROOT, "R/extra_analyses.R"))
source(file.path(PROJECT_ROOT, "R/fa_pathway_genes.R"))
source(file.path(PROJECT_ROOT, "R/fa_pathway_sensitivity.R"))
source(file.path(PROJECT_ROOT, "R/fa_pathway_heatmap.R"))
source(file.path(PROJECT_ROOT, "R/fa_pathway_sample_heatmap.R"))
source(file.path(PROJECT_ROOT, "R/top_lfc_sample_heatmap.R"))

# ----------------------------
# 1. Global setup
# ----------------------------
library(biomaRt)
library(tximport)
library(DESeq2)
library(plotly)
library(pheatmap)
library(msigdbr)
library(fgsea)



tx2gene <- get_tx2gene()


# ou "uswest", "asia"
# gene mapping (ENSG <-> SYMBOL)
gene_map <- tx2gene %>%
  dplyr::select(GENEID, GENENAME) %>%
  dplyr::distinct()

# convert GOI symbols -> ENSG
goi_ensg <- gene_map$GENEID[
  gene_map$GENENAME %in% GENES_OF_INTEREST
]

# Gènes chrX / chrY / XIST : exclus des contrastes non appariés sans ajustement sur le sexe
sexchr_ids <- if (any(vapply(DATASETS, function(d) identical(d$sex_adjust, FALSE), logical(1)))) {
  get_sex_chr_genes(tx2gene)
} else character(0)

dir.create(BASE_OUT, showWarnings = FALSE, recursive = TRUE)

# Traçabilité : copie de la config + versions des paquets utilisées pour ce run
file.copy(file.path(PROJECT_ROOT, "config/config.R"),
          file.path(BASE_OUT, "config_snapshot.R"), overwrite = TRUE)
writeLines(capture.output(sessionInfo()), file.path(BASE_OUT, "sessionInfo.txt"))

# ----------------------------
# 2. Helper: run one dataset
# ----------------------------
run_dataset <- function(dataset_name, dataset) {
  
  message("\n==============================")
  message("DATASET: ", dataset_name)
  message("==============================\n")
  
  out_dir <- file.path(BASE_OUT, dataset_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  samples    <- dataset$samples
  conditions <- dataset$conditions
  sex        <- dataset$sex
  line       <- dataset$line   # NULL si non défini -> pas d'appariement
  sex_adjust <- if (is.null(dataset$sex_adjust)) TRUE else isTRUE(dataset$sex_adjust)
  base_path  <- dataset$path
  
  # ----------------------------
  # 3. Define contrasts
  # ----------------------------
  cond_levels <- unique(conditions)
  contrasts   <- list()
  
  if (all(c("control", "FA") %in% cond_levels))
    contrasts[["control_vs_FA"]] <- c("control", "FA")
  
  if (all(c("control", "FA_edited") %in% cond_levels))
    contrasts[["control_vs_FA_edited"]] <- c("control", "FA_edited")
  
  if (all(c("FA_edited", "FA") %in% cond_levels))
    contrasts[["FA_edited_vs_FA"]] <- c("FA_edited", "FA")
  
  all_results      <- list()
  all_results_full <- list()   # idem + lfc_mle (LFC non shrunk) : analyses complémentaires
  vsd_store        <- list()   # <-- collects vsd per contrast
  
  # ----------------------------
  # 3b. QC global du dataset (PCA tous échantillons, sexe, marqueurs de différenciation)
  # ----------------------------
  tryCatch(
    run_dataset_qc(dataset_name, dataset, tx2gene, out_dir),
    error = function(e) message("  [WARN] QC global échoué : ", conditionMessage(e))
  )
  
  # ----------------------------
  # 4. Run DESeq2 per contrast
  # ----------------------------
  for (cn in names(contrasts)) {
    
    contrast_dir <- file.path(out_dir, cn)
    dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)
    
    message("\n--- Contrast: ", cn, " ---")
    
    condA <- contrasts[[cn]][1]
    condB <- contrasts[[cn]][2]
    idx   <- which(conditions %in% c(condA, condB))
    
    files <- file.path(base_path, samples[idx], "quant.sf")
    names(files) <- samples[idx]
    
    res <- run_deseq2_contrast(
      files_pair  = files,
      samples_pair = samples[idx],
      conds_pair   = conditions[idx],
      donor_pair   = sex[idx],
      line_pair    = line[idx],
      tx2gene      = tx2gene,
      condA        = condA,
      condB        = condB,
      min_counts   = MIN_COUNTS,
      min_samples  = MIN_SAMPLES,
      alpha        = PADJ_THRESHOLD,
      sex_adjust   = sex_adjust,
      exclude_genes = sexchr_ids
    )
    
    # ----- Sauvegardes : tables de résultats, design utilisé, dds -----
    utils::write.csv(res$res_df_full, file.path(contrast_dir, paste0(cn, "_DEG_results.csv")),
                     row.names = FALSE)
    writeLines(
      paste0("design: ", paste(deparse(res$design), collapse = ""),
             " | apparié: ", res$paired,
             " | filtre: >= ", MIN_COUNTS, " comptages dans >= ", MIN_SAMPLES, " échantillons",
             " | alpha: ", PADJ_THRESHOLD,
             " | sex_adjust: ", sex_adjust,
             " | gènes chrX/chrY/XIST exclus: ", res$n_excluded,
             " | LFC: log2FoldChange = shrunk (apeglm), lfc_mle = non shrunk"),
      file.path(contrast_dir, paste0(cn, "_design.txt"))
    )
    if (isTRUE(SAVE_DDS)) saveRDS(res$dds, file.path(contrast_dir, paste0(cn, "_dds.rds")))
    
    # ----- Diagnostics DESeq2 (dispersions, MA-plot, histogramme des p-values) -----
    plot_deseq_diagnostics(res$dds, res$res, contrast_dir, cn, alpha = PADJ_THRESHOLD)
    
    all_results[[cn]]      <- res$res_df
    all_results_full[[cn]] <- res$res_df_full
    vsd_store[[cn]]   <- res$vsd   # save vsd object for cross-dataset heatmap
    
    # ----------------------------
    # 5. PCA
    # ----------------------------
    if (!is.null(res$vsd)) {
      pdf(file.path(contrast_dir, paste0(cn, "_PCA_2D.pdf")))
      print(plot_pca_2d(res$vsd, dataset_name, cn))
      dev.off()
      
      p3d <- plot_pca_3d(res$vsd, dataset_name, cn)
      htmlwidgets::saveWidget(
        p3d,
        file.path(contrast_dir, paste0(cn, "_PCA_3D.html")),
        selfcontained = TRUE
      )
    }
    
    # ----------------------------
    # 6. Volcano plots
    # ----------------------------
    pdf(file.path(contrast_dir, paste0(cn, "_volcano.pdf")), width = 15, height = 12)
    print(plot_volcano(
      res$res_df, dataset_name, cn,
      direction_label = paste(condB, "vs", condA),
      lfc_thr  = LFC_THRESHOLD,
      padj_thr = PADJ_THRESHOLD
    ))
    dev.off()
    
    pdf(file.path(contrast_dir, paste0(cn, "_volcano_GOI.pdf")), width = 15, height = 12)
    print(plot_volcano_focus(
      res$res_df, dataset_name, cn,
      direction_label   = paste(condB, "vs", condA),
      genes_of_interest = GENES_OF_INTEREST,
      lfc_thr  = LFC_THRESHOLD,
      padj_thr = PADJ_THRESHOLD
    ))
    dev.off()
    
    # ----------------------------
    # 7. Heatmaps
    # ----------------------------
    sample_order <- rownames(res$colData)
    
    pdf(file.path(contrast_dir, paste0(cn, "_heatmap_topDEGs.pdf")), width = 10, height = 20)
    print(plot_heatmap_top_degs(
      vsd_mat      = assay(res$vsd),
      res_df       = res$res_df,
      colData      = res$colData,
      sample_order = sample_order,
      dataset_name = dataset_name,
      contrast_name = cn,
      top_n        = HEATMAP_TOP_N,
      padj_thr     = PADJ_THRESHOLD
    ))
    dev.off()
    
    pdf(file.path(contrast_dir, paste0(cn, "_heatmap_GOI.pdf")), width = 10, height = 20)
    print(plot_heatmap_goi(
      vsd_mat           = assay(res$vsd),
      colData           = res$colData,
      sample_order      = sample_order,
      gene_map          = gene_map,
      genes_of_interest = goi_ensg,
      dataset_name      = dataset_name,
      contrast_name     = cn
    ))
    dev.off()
    
    # ----------------------------
    # 8. fgsea
    # ----------------------------
    ranks <- make_ranks(res$res_df)
    
    run_all_fgsea(
      ranks             = ranks,
      out_dir           = contrast_dir,
      label             = cn,
      collections       = FGSEA_COLLECTIONS,
      neuronal_keywords = NEURONAL_KEYWORDS,
      direction_label   = paste(condB, "vs", condA),
      top_n             = FGSEA_TOP_N,
      padj_cut          = FGSEA_PADJ_CUT
    )
  }
  
  
  
  
  
  
  
  
  
  
  
  # ----------------------------
  
  # 9. FXN expression / DE results
  
  # ----------------------------
  
  cnames <- names(all_results)
  
  for (cn in cnames) {
    
    res_df <- all_results[[cn]]
    
    # Recherche de FXN par symbole
    
    fxn <- res_df[
      res_df$GENENAME == "FXN" |
        rownames(res_df) == "FXN",
      ,
      drop = FALSE
    ]
    
    cat("\n====================================\n")
    cat("FXN —", dataset_name, "—", cn, "\n")
    cat("====================================\n")
    
    if (nrow(fxn) == 0) {
      cat("FXN not found in this result table.\n")
    } else {
      print(fxn)
    }
  }
  
  # ----------------------------
  
  # Cross-contrast scatter
  
  # ----------------------------
  
  if (length(cnames) >= 2) {
    for (i in 1:(length(cnames) - 1)) {
      for (j in (i + 1):length(cnames)) {
        
        run_scatter(
          comp1        = cnames[i],
          comp2        = cnames[j],
          all_results  = all_results,
          out_dir      = file.path(out_dir, "SCATTER"),
          dataset_name = dataset_name
        )
        
      }
      
    }
  }
  
  
  
  # ============================================================
  # Rescue scatter
  # ============================================================
  
  run_rescue_scatter(
    all_results = all_results,
    out_dir = file.path(out_dir, "RESCUE_SCATTER"),
    dataset_name = dataset_name,
    padj_thr = PADJ_THRESHOLD,
    lfc_thr = LFC_THRESHOLD,
    top_n_labels = 100,
    export_xlsx = TRUE
  )
  
  # Nombre de DEGs par contraste (LFC shrunk vs non shrunk)
  tryCatch(
    summarise_deg_counts(all_results_full, out_dir, dataset_name,
                         padj_thr = PADJ_THRESHOLD, lfc_thr = LFC_THRESHOLD),
    error = function(e) message("  [WARN] résumé DEGs échoué : ", conditionMessage(e))
  )
  
  # Ratio de rescue (LFC non shrunk) : LFC(FA vs edited) / LFC(FA vs control) pour les FA-DEGs
  tryCatch(
    run_rescue_ratio(
      all_results  = all_results_full,
      out_dir      = file.path(out_dir, "RESCUE_RATIO"),
      dataset_name = dataset_name,
      padj_thr     = PADJ_THRESHOLD,
      lfc_thr      = LFC_THRESHOLD
    ),
    error = function(e) message("  [WARN] rescue ratio échoué : ", conditionMessage(e))
  )
  
  # DEGs appariés : effet FA-spécifique ou effet de clone (éditée) ?
  tryCatch(
    classify_paired_degs(
      all_results  = all_results_full,
      out_dir      = file.path(out_dir, "PAIRED_DEG_CLASSIFICATION"),
      dataset_name = dataset_name,
      padj_thr     = PADJ_THRESHOLD
    ),
    error = function(e) message("  [WARN] classification des DEGs appariés échouée : ", conditionMessage(e))
  )
  
  
  #   
  #   # ----------------------------
  #   # 9. Cross-contrast scatter
  #   # ----------------------------
  #   cnames <- names(all_results)
  #   
  #   if (length(cnames) >= 2) {
  #     for (i in 1:(length(cnames) - 1)) {
  #       for (j in (i + 1):length(cnames)) {
  #         run_scatter(
  #           comp1       = cnames[i],
  #           comp2       = cnames[j],
  #           all_results = all_results,
  #           out_dir     = file.path(out_dir, "SCATTER"),
  #           dataset_name = dataset_name
  #         )
  #       }
  #     }
  #   }
  #   
  
  
  message("\nFinished dataset: ", dataset_name)
  
  # Return both DEG results and vsd matrices
  return(invisible(list(deg = all_results, vsd = vsd_store)))
}





# ----------------------------
# 10. Run all datasets
# ----------------------------
all_dataset_results <- list()
all_vsd_results     <- list()

for (ds_name in names(DATASETS)) {
  out <- run_dataset(ds_name, DATASETS[[ds_name]])
  all_dataset_results[[ds_name]] <- out$deg
  all_vsd_results[[ds_name]]     <- lapply(out$vsd, function(v) assay(v))
  
}

# Étape protégée : une erreur dans une comparaison inter-datasets ne doit pas stopper le reste
# (elle ferme aussi les PDF laissés ouverts).
safe_step <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    message("\n  [ERREUR] ", label, " : ", conditionMessage(e))
    while (!is.null(grDevices::dev.list())) grDevices::dev.off()
    invisible(NULL)
  })
}

# ----------------------------
# 11. Cross-dataset comparison
# ----------------------------
DATASET_A <- "iPSC"       # must match keys in DATASETS (config.R)
DATASET_B <- "iNeurons"   # must match keys in DATASETS (config.R)

safe_step("comparaison iPSC vs iNeurons", run_cross_dataset_comparison(
  all_dataset_results = all_dataset_results,
  base_out            = BASE_OUT,
  name_a              = DATASET_A,
  name_b              = DATASET_B,
  padj_thr            = PADJ_THRESHOLD,
  lfc_thr             = LFC_THRESHOLD,
  genes_highlight     = GENES_OF_INTEREST,
  vsd_list            = all_vsd_results    # enables z-scored expression heatmaps
))


# ----------------------------
# 12. iPSC vs iNeurons par condition
# ----------------------------

safe_step("iPSC vs iNeurons par condition", run_ipsc_vs_neurons(
  all_vsd_results     = all_vsd_results,
  all_dataset_results = all_dataset_results,
  base_out            = BASE_OUT,
  name_a              = "iPSC",       # doit matcher les clés de DATASETS
  name_b              = "iNeurons",
  conditions          = c("control", "FA", "FA_edited"),
  padj_thr            = PADJ_THRESHOLD,
  lfc_thr             = LFC_THRESHOLD,
  genes_highlight     = GENES_OF_INTEREST
))




safe_step("scatter effet maladie", run_disease_effect_scatter(
  all_dataset_results = all_dataset_results,
  base_out            = BASE_OUT,
  name_a              = "iPSC",
  name_b              = "iNeurons",
  contrasts           = c("control_vs_FA", "control_vs_FA_edited", "FA_edited_vs_FA"),
  padj_thr            = PADJ_THRESHOLD, lfc_thr         = LFC_THRESHOLD
))


# ----------------------------
# 13. Concordance sans seuil entre datasets (stat de Wald, NES fgsea)
# ----------------------------
tryCatch(
  run_stat_concordance(all_dataset_results, BASE_OUT,
                       datasets = names(all_dataset_results)),
  error = function(e) message("  [WARN] concordance (stat) échouée : ", conditionMessage(e))
)
tryCatch(
  run_nes_concordance(BASE_OUT, datasets = names(all_dataset_results)),
  error = function(e) message("  [WARN] concordance (NES) échouée : ", conditionMessage(e))
)

# ----------------------------
# 14. Gènes candidats de la voie Fe-S / frataxine (liste curée, config.R)
# ----------------------------
safe_step("genes candidats FA pathway", run_fa_pathway_gene_summary(
  all_dataset_results = all_dataset_results,
  base_out            = BASE_OUT
))

# ----------------------------
# 15. Sensibilité méthodologique sur le panel FA_PATHWAY_GENES
#     (DESeq2 baseline / sex_adjust=TRUE / filtre assoupli / limma-voom)
# ----------------------------
sens_result <- safe_step("sensibilite FA pathway (methodes)", run_fa_pathway_sensitivity(
  base_out = BASE_OUT,
  tx2gene  = tx2gene
))

# ----------------------------
# 16. Heatmap LFC + significativité du panel FA_PATHWAY_GENES
#     (à partir du tableau produit à l'étape 15, sans recalcul)
# ----------------------------
safe_step("heatmap FA pathway (sensibilite)", plot_fa_pathway_heatmap(
  df_long = if (!is.null(sens_result)) sens_result$long else NULL,
  out_dir = file.path(BASE_OUT, "FA_PATHWAY_GENES")
))

# ----------------------------
# 17-19. Même trio (tableau / sensibilité / heatmap), mais pour
#        CANDIDATE_GENES_PATHWAYS (config.R) : les gènes identifiés comme
#        intéressants par l'analyse rescue/classification appariée, pas la
#        liste Fe-S fournie par l'utilisateur. Écrit dans un sous-dossier
#        séparé (CANDIDATE_GENES/) pour ne rien écraser de FA_PATHWAY_GENES/.
# ----------------------------
safe_step("genes candidats (rescue) - tableau", run_fa_pathway_gene_summary(
  all_dataset_results = all_dataset_results,
  base_out            = BASE_OUT,
  gene_list           = CANDIDATE_GENES_PATHWAYS,
  out_subdir          = "CANDIDATE_GENES",
  file_name           = "genes_candidats_rescue.xlsx"
))

candidate_sens_result <- safe_step("genes candidats (rescue) - sensibilite", run_fa_pathway_sensitivity(
  base_out   = BASE_OUT,
  tx2gene    = tx2gene,
  gene_list  = CANDIDATE_GENES_PATHWAYS,
  out_subdir = "CANDIDATE_GENES",
  file_name  = "genes_candidats_rescue_sensitivity.xlsx"
))

safe_step("genes candidats (rescue) - heatmap", plot_fa_pathway_heatmap(
  df_long   = if (!is.null(candidate_sens_result)) candidate_sens_result$long else NULL,
  out_dir   = file.path(BASE_OUT, "CANDIDATE_GENES"),
  gene_list = CANDIDATE_GENES_PATHWAYS,
  file_name = "heatmap_candidats_rescue.pdf"
))

message("\n✔ ALL ANALYSES COMPLETED (résultats dans ", BASE_OUT, ")")

# ----------------------------
# 20-21. Heatmaps d'expression PAR ÉCHANTILLON (méthode standard, pas de
#        comparaison de méthodes), une carte par type cellulaire
#        (fibro, iPSC, iNeurons), pour les deux listes de gènes.
# ----------------------------
safe_step("heatmaps par echantillon - FA_PATHWAY_GENES", run_fa_pathway_sample_heatmaps(
  tx2gene             = tx2gene,
  gene_list           = FA_PATHWAY_GENES,
  out_dir             = file.path(BASE_OUT, "FA_PATHWAY_GENES"),
  all_dataset_results = all_dataset_results
))

safe_step("heatmaps par echantillon - CANDIDATE_GENES", run_fa_pathway_sample_heatmaps(
  tx2gene             = tx2gene,
  gene_list           = CANDIDATE_GENES_PATHWAYS,
  out_dir             = file.path(BASE_OUT, "CANDIDATE_GENES"),
  all_dataset_results = all_dataset_results
))

message("\n✔ TOUT EST TERMINÉ (y compris heatmaps par échantillon)")

# ----------------------------
# 22. Heatmaps des DEGs significatifs tries par |LFC| le plus fort, un par
#     comparaison, deposes dans le dossier de cette comparaison
#     (results_v4/<dataset>/<contraste>/) -- pas un panel pre-defini.
# ----------------------------
safe_step("heatmaps top-LFC par comparaison", run_top_lfc_sample_heatmaps(
  all_dataset_results = all_dataset_results,
  tx2gene             = tx2gene,
  base_out            = BASE_OUT,
  top_n               = 50
))

message("\n✔ VRAIMENT TOUT EST TERMIN\u00c9")