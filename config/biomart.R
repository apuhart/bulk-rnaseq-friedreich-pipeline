# ============================================================
# config/biomart.R — Transcript-to-gene mapping
# ============================================================
# Deux sources possibles (dans cet ordre) :
#   1. TX2GENE_FILE : fichier nf-core `salmon_tx2gene.tsv`
#      -> même annotation que celle utilisée par Salmon (recommandé)
#   2. GTF Ensembl local (GTF_PATH), téléchargé si absent.
#      Le résultat est mis en cache (.rds) pour éviter de relire le GTF à chaque run.
#
# Retourne un data.frame : TXNAME, GENEID, GENENAME (sans numéros de version)

get_tx2gene <- function(
    tx2gene_file = get0("TX2GENE_FILE", ifnotfound = NULL),
    gtf_path     = get0("GTF_PATH",
                        ifnotfound = file.path(getwd(), "data", "Homo_sapiens.GRCh38.112.gtf.gz")),
    gtf_url      = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
) {

  # ---------- 1. fichier nf-core ----------
  if (!is.null(tx2gene_file) && nzchar(tx2gene_file) && file.exists(tx2gene_file)) {

    message("tx2gene : lecture de ", tx2gene_file)
    tx <- utils::read.delim(tx2gene_file, header = FALSE,
                            stringsAsFactors = FALSE, quote = "")
    if (ncol(tx) == 2) tx$V3 <- tx$V2          # pas de nom de gène -> ENSG
    tx <- tx[, 1:3]
    colnames(tx) <- c("TXNAME", "GENEID", "GENENAME")
    # en-tête éventuel
    if (tx$TXNAME[1] %in% c("transcript_id", "tx_id", "TXNAME", "tx_name")) tx <- tx[-1, ]

    tx$TXNAME <- sub("\\..*", "", tx$TXNAME)
    tx$GENEID <- sub("\\..*", "", tx$GENEID)
    tx <- dplyr::distinct(tx)

    message("tx2gene : ", nrow(tx), " transcrits, ",
            length(unique(tx$GENEID)), " gènes.")
    return(tx)
  }

  if (!is.null(tx2gene_file) && nzchar(tx2gene_file))
    warning("TX2GENE_FILE introuvable (", tx2gene_file, ") -> repli sur le GTF.")

  # ---------- 2. GTF local + cache ----------
  cache_file <- sub("\\.gtf(\\.gz)?$", "_tx2gene.rds", gtf_path)
  if (file.exists(cache_file)) {
    message("tx2gene : lecture du cache ", cache_file)
    return(readRDS(cache_file))
  }

  if (!file.exists(gtf_path)) {
    message("GTF introuvable localement, téléchargement depuis : ", gtf_url)
    dir.create(dirname(gtf_path), showWarnings = FALSE, recursive = TRUE)
    utils::download.file(gtf_url, gtf_path, mode = "wb")
  }

  message("Lecture du GTF : ", gtf_path, " (peut prendre quelques minutes) ...")

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    BiocManager::install("rtracklayer", update = FALSE, ask = FALSE)
  }

  gtf    <- rtracklayer::import(gtf_path)
  gtf_df <- as.data.frame(gtf)

  tx_df <- gtf_df[gtf_df$type == "transcript",
                  c("transcript_id", "gene_id", "gene_name")]
  colnames(tx_df) <- c("TXNAME", "GENEID", "GENENAME")

  tx_df$TXNAME <- sub("\\..*", "", tx_df$TXNAME)
  tx_df$GENEID <- sub("\\..*", "", tx_df$GENEID)
  tx_df <- dplyr::distinct(tx_df)

  message("tx2gene : ", nrow(tx_df), " transcrits, ",
          length(unique(tx_df$GENEID)), " gènes.")

  saveRDS(tx_df, cache_file)
  message("tx2gene : cache écrit dans ", cache_file)

  return(tx_df)
}

annotate_results <- function(res_df, tx2gene) {
  gene_map <- tx2gene %>%
    dplyr::select(GENEID, GENENAME) %>%
    dplyr::distinct()
  res_df %>%
    dplyr::left_join(gene_map, by = "GENEID") %>%
    dplyr::mutate(
      gene = dplyr::coalesce(GENENAME, GENEID),
      gene = ifelse(is.na(gene) | gene == "", GENEID, gene)
    )
}


# ============================================================
# Chromosome de chaque gène (pour exclure chrX / chrY / XIST des contrastes non ajustés sur le sexe)
# ============================================================
# Lit les lignes "gene" du GTF (rapide via zcat|grep sous Linux/macOS, sinon readLines) et met en cache
# un data.frame (GENEID, chr) dans data/<gtf>_gene_chr.rds.
get_gene_chr <- function(
    gtf_path = get0("GTF_PATH", ifnotfound = file.path(getwd(), "data", "Homo_sapiens.GRCh38.112.gtf.gz"))
) {
  cache_file <- sub("\\.gtf(\\.gz)?$", "_gene_chr.rds", gtf_path)
  if (file.exists(cache_file)) return(readRDS(cache_file))

  if (!file.exists(gtf_path)) {
    warning("GTF introuvable (", gtf_path, ") : chromosomes indisponibles.")
    return(NULL)
  }

  message("Lecture des chromosomes depuis le GTF (une seule fois, résultat mis en cache) ...")
  reader <- if (grepl("\\.gz$", gtf_path)) "zcat" else "cat"
  lines <- tryCatch(
    suppressWarnings(system(
      sprintf("%s '%s' | grep -P '\\tgene\\t' | cut -f1,9", reader, gtf_path), intern = TRUE
    )),
    error = function(e) character(0)
  )
  if (length(lines) == 0) {                       # repli portable (plus lent, plus de mémoire)
    con <- if (grepl("\\.gz$", gtf_path)) gzfile(gtf_path, "rt") else file(gtf_path, "rt")
    all_lines <- readLines(con); close(con)
    all_lines <- all_lines[!startsWith(all_lines, "#")]
    all_lines <- all_lines[grepl("\tgene\t", all_lines, fixed = TRUE)]
    parts <- strsplit(all_lines, "\t", fixed = TRUE)
    lines <- vapply(parts, function(p) paste(p[1], p[9], sep = "\t"), character(1))
  }

  parts <- strsplit(lines, "\t", fixed = TRUE)
  out <- data.frame(
    GENEID = vapply(parts, function(p) sub('.*gene_id "([^"]+)".*', "\\1", p[2]), character(1)),
    chr    = vapply(parts, `[`, character(1), 1),
    stringsAsFactors = FALSE
  )
  out$GENEID <- sub("\\..*", "", out$GENEID)
  out <- out[!duplicated(out$GENEID), ]

  saveRDS(out, cache_file)
  message("Chromosomes : ", nrow(out), " gènes (cache : ", cache_file, ")")
  out
}


#' Identifiants (ENSG sans version) des gènes à exclure quand on ne s'ajuste pas sur le sexe :
#' tous les gènes chrX et chrY + XIST (et TSIX). Sans GTF, seuls XIST/TSIX sont exclus (avertissement).
get_sex_chr_genes <- function(tx2gene, gtf_path = get0("GTF_PATH", ifnotfound = NULL)) {
  gm <- unique(tx2gene[, c("GENEID", "GENENAME")])
  ids_xist <- gm$GENEID[gm$GENENAME %in% c("XIST", "TSIX")]

  gc <- if (is.null(gtf_path)) get_gene_chr() else get_gene_chr(gtf_path)
  if (is.null(gc)) {
    warning("Chromosomes indisponibles : seuls XIST/TSIX sont exclus (les gènes chrX/chrY restent dans le test).")
    return(unique(ids_xist))
  }
  ids <- unique(c(gc$GENEID[gc$chr %in% c("X", "Y")], ids_xist))
  message("Gènes exclus des contrastes non ajustés sur le sexe (chrX, chrY, XIST) : ", length(ids))
  ids
}
