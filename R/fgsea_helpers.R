# ============================================================
# R/functions/fgsea_helpers.R — fgsea wrappers
# ============================================================

#' Fetch MSigDB pathways, optionally filtered by neuronal/biological keywords
#'
#' @param category        MSigDB category (e.g. "H", "C5")
#' @param subcategory     MSigDB subcategory (e.g. "BP", NULL)
#' @param neuronal_only   If TRUE, filter pathway names by keyword list
#' @param keywords        Character vector of keywords for filtering
#' @return Named list of gene vectors
get_pathways <- function(
    category,
    subcategory    = NULL,
    neuronal_only  = FALSE,
    keywords       = NULL
) {
  # msigdbr >= 10 : arguments `collection` / `subcollection` (GO:BP, GO:CC, GO:MF, ...)
  # msigdbr <  10 : arguments `category`   / `subcategory`   (BP, CC, MF, ...)
  new_api <- utils::packageVersion("msigdbr") >= "10.0.0"

  msig <- if (new_api) {
    subcol <- subcategory
    if (!is.null(subcol) && subcol %in% c("BP", "CC", "MF")) subcol <- paste0("GO:", subcol)
    if (is.null(subcol)) {
      msigdbr::msigdbr(species = "Homo sapiens", collection = category)
    } else {
      msigdbr::msigdbr(species = "Homo sapiens", collection = category, subcollection = subcol)
    }
  } else {
    if (is.null(subcategory)) {
      msigdbr::msigdbr(species = "Homo sapiens", category = category)
    } else {
      msigdbr::msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory)
    }
  }

  if (neuronal_only && !is.null(keywords)) {
    # Mots-clés courts (<= 4 lettres) : cherchés en DÉBUT de mot ("gaba" -> GABAERGIC, "atp" -> ATPASE) ;
    # "ros", "isc", "upr", "chop" : mot ENTIER seulement (sinon "ros" matche PROSTATE, "isc" matche DISCRETE).
    whole <- c("ros", "isc", "upr", "chop")            # ambigus -> mot entier uniquement
    pats  <- ifelse(tolower(keywords) %in% whole,
                    paste0("\\b", keywords, "\\b"),
                    ifelse(nchar(keywords) <= 4, paste0("\\b", keywords), keywords))
    name_spaced <- gsub("_", " ", msig$gs_name)
    hit <- grepl(paste(pats, collapse = "|"), name_spaced, ignore.case = TRUE, perl = TRUE)
    msig_filt <- msig[hit, ]

    if (length(unique(msig_filt$gs_name)) < 10) {
      message("Too few filtered pathways — falling back to full set.")
      return(split(msig$gene_symbol, msig$gs_name))
    }
    message("Neuronal/keyword filter: ", length(unique(msig_filt$gs_name)), " pathways retained.")
    return(split(msig_filt$gene_symbol, msig_filt$gs_name))
  }

  split(msig$gene_symbol, msig$gs_name)
}


#' Run fgsea and save results + dotplot
#'
#' @param ranks          Named numeric vector (gene -> log2FC), sorted decreasingly
#' @param out_dir        Output directory
#' @param label          Label prefix for output files
#' @param category       MSigDB category
#' @param subcategory    MSigDB subcategory (or NULL)
#' @param top_n          Number of top pathways to show on dotplot
#' @param padj_cut       Adjusted p-value cut-off for dotplot display
#' @param nperm          DEPRECATED, ignored (fgseaMultilevel n'utilise plus de permutations)
#' @param min_size       Taille min des gene sets (défaut : FGSEA_MIN_SIZE ou 15)
#' @param max_size       Taille max des gene sets (défaut : FGSEA_MAX_SIZE ou 500)
#' @param seed           Graine aléatoire (défaut : FGSEA_SEED ou 42)
#' @param neuronal_only  Restrict to neuronal/biological keywords
#' @param keywords       Keywords for neuronal_only filter
#' @param direction_label  String describing log2FC direction
#' @return fgsea result data.frame (with leadingEdge as string), or NULL on error
run_fgsea_safe <- function(
    ranks,
    out_dir,
    label,
    category        = "H",
    subcategory     = NULL,
    top_n           = 30,
    padj_cut        = 0.3,
    nperm           = NULL,
    min_size        = get0("FGSEA_MIN_SIZE", ifnotfound = 15),
    max_size        = get0("FGSEA_MAX_SIZE", ifnotfound = 500),
    seed            = get0("FGSEA_SEED",     ifnotfound = 42),
    neuronal_only   = FALSE,
    keywords        = NULL,
    direction_label = ""
) {
  tryCatch({

    pathways <- get_pathways(
      category      = category,
      subcategory   = subcategory,
      neuronal_only = neuronal_only,
      keywords      = keywords
    )

    set.seed(seed)
    fg    <- fgsea::fgseaMultilevel(
      pathways = pathways,
      stats    = ranks,
      minSize  = min_size,
      maxSize  = max_size
    )
    fg    <- fg %>% dplyr::arrange(padj)
    fg_out <- fg %>%
      dplyr::mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

    # ---- file suffix ----
    coll_tag <- paste0(category, ifelse(is.null(subcategory), "", paste0("_", subcategory)))
    suffix   <- if (neuronal_only) "_GOI" else ""
    base     <- file.path(out_dir, paste0(label, "_", coll_tag, "_fgsea", suffix))

    write.csv(fg_out, paste0(base, ".csv"), row.names = FALSE)

    # ---- dotplot ----
    top_fg <- fg %>%
      dplyr::filter(!is.na(padj), padj < padj_cut) %>%
      dplyr::arrange(padj) %>%
      head(top_n) %>%
      dplyr::mutate(
        genes_in_pathway = sapply(leadingEdge, length),
        padj_cat = cut(
          padj,
          breaks = c(-Inf, 0.05, 0.15, 0.30, Inf),
          labels = c("< 0.05", "0.05-0.15", "0.15-0.30", "> 0.30")
        )
      )

    if (nrow(top_fg) == 0) {
      message("No pathways with padj < ", padj_cut, " for ", label, " (", coll_tag, ")")
      return(fg_out)
    }

    gsea_plot <- ggplot(
      top_fg,
      aes(x = reorder(pathway, NES), y = NES, size = genes_in_pathway, color = padj_cat)
    ) +
      geom_point(alpha = 0.9) +
      scale_color_manual(
        values = c(
          "< 0.05"    = "firebrick",
          "0.05-0.15" = "orange",
          "0.15-0.30" = "gold",
          "> 0.30"    = "grey70"
        ),
        name = "Adjusted p-value"
      ) +
      scale_size_continuous(range = c(3, 8)) +
      coord_flip() +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 12)
      ) +
      labs(
        x     = "Pathway",
        y     = "NES",
        size  = "Number of genes",
        title = paste0("GSEA (", coll_tag, ") - ", label,
                       ifelse(neuronal_only, " [keyword filter]", ""),
                       "\n", direction_label)
      )

    png(paste0(base, "_dotplot.png"), width = 2000, height = 1000)
    print(gsea_plot)
    dev.off()

    message("fgsea done: ", base)
    return(fg_out)

  }, error = function(e) {
    message("fgsea failed for ", label, " (", category,
            ifelse(is.null(subcategory), "", paste0("_", subcategory)), "): ", e$message)
    return(NULL)
  })
}


#' Run all configured fgsea collections for a contrast
#'
#' @param ranks             Named ranked vector
#' @param out_dir           Output directory
#' @param label             Label prefix
#' @param collections       List of lists with $category and $subcategory
#' @param neuronal_keywords Keywords for neuronal filter
#' @param direction_label   Direction label string
#' @param ...               Extra args passed to run_fgsea_safe
run_all_fgsea <- function(
    ranks, out_dir, label,
    collections, neuronal_keywords,
    direction_label = "", ...
) {
  for (col in collections) {
    # full gene set
    run_fgsea_safe(
      ranks           = ranks,
      out_dir         = out_dir,
      label           = label,
      category        = col$category,
      subcategory     = col$subcategory,
      neuronal_only   = FALSE,
      direction_label = direction_label,
      ...
    )
    # keyword-filtered
    run_fgsea_safe(
      ranks           = ranks,
      out_dir         = out_dir,
      label           = label,
      category        = col$category,
      subcategory     = col$subcategory,
      neuronal_only   = TRUE,
      keywords        = neuronal_keywords,
      direction_label = direction_label,
      ...
    )
  }
}
