# ============================================================
# fa_pathway_genes.R
#
# Croise la liste curée FA_PATHWAY_GENES (config.R : machinerie Fe-S
# mitochondriale et cytosolique, protéines Fe-S nucléaires/ADN, OXPHOS/TCA,
# ISR, métabolisme sérine/1-carbone, AMPK, ferroptose, NRF2, fer) avec les
# résultats DESeq2 de TOUS les datasets et contrastes disponibles, et exporte
# un classeur Excel dans results_v4/FA_PATHWAY_GENES/.
#
# Contexte : demandé à la main le 2026-09-22 à partir d'une liste de gènes
# fournie par l'utilisateur pour l'ataxie de Friedreich ; intégré ici au
# pipeline pour se régénérer automatiquement à chaque run.
# ============================================================

run_fa_pathway_gene_summary <- function(
    all_dataset_results,
    base_out,
    gene_list = FA_PATHWAY_GENES,
    contrasts = c("control_vs_FA", "control_vs_FA_edited", "FA_edited_vs_FA"),
    padj_thr  = PADJ_THRESHOLD,
    lfc_thr   = LFC_THRESHOLD,
    out_subdir = "FA_PATHWAY_GENES",
    file_name  = "genes_candidats_FA_pathway.xlsx"
) {
  
  out_dir <- file.path(base_out, out_subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------------------------------------------------
  # Aplatit la liste de gènes en data.frame (categorie, libelle, symbole)
  # ------------------------------------------------------------
  gene_df <- do.call(rbind, lapply(names(gene_list), function(cat) {
    v <- gene_list[[cat]]
    data.frame(
      categorie = cat,
      libelle   = names(v),
      symbole   = unname(v),
      stringsAsFactors = FALSE
    )
  }))
  
  # ------------------------------------------------------------
  # Repère un gène (par symbole HGNC) dans une table de résultats DESeq2
  # ------------------------------------------------------------
  lookup_gene <- function(res_df, symbol) {
    empty <- list(ensg = NA_character_, lfc = NA_real_, padj = NA_real_, sig = NA)
    if (is.null(res_df) || !("GENENAME" %in% colnames(res_df))) return(empty)
    hit <- res_df[!is.na(res_df$GENENAME) & res_df$GENENAME == symbol, , drop = FALSE]
    if (nrow(hit) == 0) return(empty)
    hit <- hit[1, ]
    sig <- !is.na(hit$padj) & hit$padj < padj_thr & abs(hit$log2FoldChange) >= lfc_thr
    list(ensg = hit$GENEID, lfc = hit$log2FoldChange, padj = hit$padj, sig = isTRUE(sig))
  }
  
  # Statut rescue simplifié (même logique que run_rescue_scatter() dans plots.R),
  # calculé uniquement quand control_vs_FA ET control_vs_FA_edited existent tous
  # les deux pour un dataset.
  rescue_label <- function(fa_info, ed_info) {
    if (!isTRUE(fa_info$sig)) return(NULL)
    if (is.na(ed_info$lfc)) return("Not testable (FA_edited)")
    if (!isTRUE(ed_info$sig) && abs(ed_info$lfc) < abs(fa_info$lfc)) return("Rescued")
    if (isTRUE(ed_info$sig)  && abs(ed_info$lfc) < abs(fa_info$lfc)) return("Partial rescue")
    if (sign(fa_info$lfc) != sign(ed_info$lfc)) return("Reversed")
    "Not rescued"
  }
  
  datasets <- names(all_dataset_results)
  rows <- vector("list", nrow(gene_df))
  missing_genes <- character(0)
  
  for (i in seq_len(nrow(gene_df))) {
    
    sym <- gene_df$symbole[i]
    row <- list(
      Categorie    = gene_df$categorie[i],
      Nom_fourni   = gene_df$libelle[i],
      Symbole_HGNC = sym
    )
    
    ensg_found    <- NA_character_
    rescue_status <- character(0)
    
    for (ds in datasets) {
      
      ds_results <- all_dataset_results[[ds]]
      
      for (cn in contrasts) {
        info <- lookup_gene(ds_results[[cn]], sym)
        if (!is.na(info$ensg)) ensg_found <- info$ensg
        
        row[[paste0(ds, "_", cn, "_LFC")]]  <- if (is.na(info$lfc)) NA_real_ else round(info$lfc, 3)
        row[[paste0(ds, "_", cn, "_padj")]] <- info$padj
        row[[paste0(ds, "_", cn, "_sig")]]  <- info$sig
      }
      
      fa_info <- lookup_gene(ds_results[["control_vs_FA"]], sym)
      ed_info <- lookup_gene(ds_results[["control_vs_FA_edited"]], sym)
      st <- rescue_label(fa_info, ed_info)
      if (!is.null(st)) rescue_status <- c(rescue_status, paste0(ds, ": ", st))
    }
    
    row[["ENSG"]]          <- ensg_found
    row[["Rescue_status"]] <- if (length(rescue_status) == 0) "-" else paste(rescue_status, collapse = " | ")
    row[["Note"]]          <- if (is.na(ensg_found)) "Gene non retrouve (verifier symbole HGNC / filtre DESeq2)" else ""
    
    if (is.na(ensg_found)) missing_genes <- c(missing_genes, sym)
    
    rows[[i]] <- row
  }
  
  # ------------------------------------------------------------
  # Assemble le tableau final (une ligne par gène, colonnes ordonnées)
  # ------------------------------------------------------------
  # Toutes les lignes doivent avoir les mêmes noms de colonnes (mêmes datasets
  # x contrastes) : on les uniformise avant le rbind.
  all_names <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) {
    missing <- setdiff(all_names, names(r))
    r[missing] <- NA
    r[all_names]
  })
  df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  
  front_cols <- c("Categorie", "Nom_fourni", "Symbole_HGNC", "ENSG")
  back_cols  <- c("Rescue_status", "Note")
  mid_cols   <- setdiff(names(df), c(front_cols, back_cols))
  df <- df[, c(front_cols, mid_cols, back_cols)]
  
  # ------------------------------------------------------------
  # Feuille légende
  # ------------------------------------------------------------
  sig_cols <- grep("_sig$", names(df), value = TRUE)
  any_sig  <- if (length(sig_cols) > 0) apply(df[, sig_cols, drop = FALSE], 1, function(x) any(x %in% TRUE)) else rep(FALSE, nrow(df))
  n_sig    <- sum(any_sig)
  
  legende <- data.frame(
    Cle = c("Objectif", "Colonnes _LFC / _padj / _sig", "Colonne Rescue_status",
            "Colonne Note", "Constat", "Genes non retrouves"),
    Valeur = c(
      "Croise FA_PATHWAY_GENES (config.R) avec les resultats DESeq2 (control_vs_FA, control_vs_FA_edited, FA_edited_vs_FA) de tous les datasets disponibles dans ce run.",
      "Une triplette par dataset x contraste (ex. iPSC_control_vs_FA_LFC / _padj / _sig). NA si ce contraste n'existe pas pour ce dataset (ex. fibro n'a pas de lignee FA_edited).",
      "Calcule uniquement quand control_vs_FA ET control_vs_FA_edited existent pour un dataset donne (Rescued / Partial rescue / Not rescued / Reversed) -- meme logique que RESCUE_SCATTER/run_rescue_scatter().",
      "Signale les genes non retrouves sous ce symbole HGNC dans les tables DESeq2 (a verifier : alias different, ou filtre pour faible expression).",
      paste0(n_sig, "/", nrow(df), " genes de la liste sont un DEG significatif (padj < ", padj_thr, ", |LFC| >= ", lfc_thr, ") dans au moins un dataset/contraste de ce run."),
      if (length(missing_genes) == 0) "Aucun" else paste(unique(missing_genes), collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
  
  out_file <- file.path(out_dir, file_name)
  writexl::write_xlsx(
    list("Genes candidats FA" = df, "Legende" = legende),
    out_file
  )
  
  if (length(missing_genes) > 0) {
    message("  [INFO] FA_PATHWAY_GENES : ", length(unique(missing_genes)),
            " gene(s) non retrouve(s) : ", paste(unique(missing_genes), collapse = ", "))
  }
  message("  FA_PATHWAY_GENES : ", n_sig, "/", nrow(df),
          " genes significatifs -- classeur ecrit -> ", out_file)
  
  invisible(df)
}