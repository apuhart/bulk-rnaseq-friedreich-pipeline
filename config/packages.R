# ============================================================
# R/functions/packages.R — Package installation & loading
# ============================================================

load_packages <- function() {

  bioc_pkgs <- c(
    "tximport", "DESeq2", "biomaRt", "apeglm",
    "fgsea", "limma", "edgeR"
  )

  cran_pkgs <- c(
    "dplyr", "ggplot2", "ggrepel", "pheatmap",
    "matrixStats", "msigdbr", "WGCNA",
    "plotly", "htmlwidgets", "tibble"
  )

  # Install BiocManager if needed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # Install missing Bioconductor packages
  missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE)
  }

  # Install missing CRAN packages
  missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_cran) > 0) {
    message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran)
  }

  # Load all
  all_pkgs <- c(bioc_pkgs, cran_pkgs)
  invisible(lapply(all_pkgs, library, character.only = TRUE))
  message("All packages loaded.")
}


# ============================================================
# Vérification des paquets nécessaires (n'attache RIEN : évite les masquages
# de fonctions dplyr::select, etc.)
#   check_packages()            : liste ce qui manque (stop si un paquet requis manque)
#   install_missing_packages()  : installe les paquets manquants (à lancer une fois)
# ============================================================
# Requis : le pipeline s'arrête s'ils manquent
REQUIRED_CRAN <- c("dplyr", "tibble", "ggplot2", "ggrepel", "pheatmap", "plotly",
                   "htmlwidgets", "msigdbr", "jsonlite", "VennDiagram")
REQUIRED_BIOC <- c("DESeq2", "tximport", "fgsea")
# Optionnels : simple avertissement (fonctionnalité dégradée ou ignorée si absents)
#  apeglm -> sinon lfcShrink retombe sur le type "normal" (moins précis)
#  rtracklayer -> seulement si tx2gene est lu depuis un GTF sans cache
OPTIONAL_PKGS <- c("apeglm", "rtracklayer", "writexl", "tidyr", "matrixStats",
                   "clusterProfiler", "AnnotationDbi", "ReactomePA", "readr", "rstudioapi",
                   "limma", "edgeR")

.missing_pkgs <- function(pkgs) pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

check_packages <- function() {
  miss_req <- .missing_pkgs(c(REQUIRED_CRAN, REQUIRED_BIOC))
  miss_opt <- .missing_pkgs(OPTIONAL_PKGS)
  if (length(miss_opt) > 0)
    message("Paquets optionnels absents (fonctionnalités dégradées ou ignorées) : ",
            paste(miss_opt, collapse = ", "))
  if (length(miss_req) > 0)
    stop("Paquets requis manquants : ", paste(miss_req, collapse = ", "),
         "\n  -> lance install_missing_packages() puis relance le pipeline.")
  invisible(TRUE)
}

install_missing_packages <- function(include_optional = TRUE) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  bioc <- c(REQUIRED_BIOC, if (include_optional) c("apeglm", "rtracklayer", "clusterProfiler",
                                                   "AnnotationDbi", "ReactomePA", "limma", "edgeR"))
  cran <- c(REQUIRED_CRAN, if (include_optional) c("writexl", "tidyr", "matrixStats", "readr",
                                                   "rstudioapi"))
  mb <- .missing_pkgs(bioc); mc <- .missing_pkgs(cran)
  if (length(mb) > 0) BiocManager::install(mb, ask = FALSE, update = FALSE)
  if (length(mc) > 0) install.packages(mc)
  invisible(TRUE)
}
