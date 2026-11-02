# =====================================================================
# 00_config.R : shared configuration for the whole pathway pipeline.
# EDIT ONLY THIS FILE (paths, cell types, thresholds, themes, reference pathways).
# Every numbered script sources it. Run order: 01 -> 02/03/04/05 ; 06 and 07 need DESeq2.
# =====================================================================

## ---- paths ----------------------------------------------------------
RES_DIR <- Sys.getenv("FA_RES_DIR", file.path(getwd(), "results_v4"))  # <- EDIT
OUT_DIR <- file.path(RES_DIR, "pathway_ABC")          # everything is written here

## ---- data layout ----------------------------------------------------
CELLS    <- c("iPSC", "iNeurons")                      # cell types with the three contrasts A, B, C
FIBRO    <- "fibro"                                    # only contrast A (set NULL to skip)
CELL_DDS <- "iNeurons"                                 # cell type used by 06_module_scores.R and 07_gene_heatmap.R
CON_A <- "control_vs_FA"; CON_B <- "control_vs_FA_edited"; CON_C <- "FA_edited_vs_FA"   # A = FA vs ctrl ; B = edited vs ctrl ; C = FA vs edited
COLLS    <- c("C5_BP", "C5_CC", "C5_MF", "H")          # fgsea collections

## ---- thresholds -----------------------------------------------------
PADJ     <- 0.05
SIZE_MIN <- 15; SIZE_MAX <- 300                        # pathway size window for candidates
EMT_NAME <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
DROP_EMT <- TRUE                                       # drop EMT/matrix-driven pathways from theme plots/heatmap (FALSE = hollow points)

## ---- themes (first matching keyword rule wins -> order matters) -------
THEMES <- list(
  "Ribosome / translation"        = "ribosom|translation|rrna|trna|ribonucleoprotein|amino acid activation",
  "Mitochondria / OXPHOS"         = "mitochond|respiratory|oxidative phosph|atp synth|electron transport|nadh|proton motive|iron sulfur|iron-sulfur|2 iron",
  "Synapse / neuronal"            = "synap|neuron|axon|dendrit|vesicle|exocyt|glutamate|gaba|nervous|neurotransmit|projection",
  "Muscle / cell fusion"          = "muscle|myo|syncytium|sarcomere|striated|contraction|fusion",
  "Matrix / adhesion"             = "collagen|matrix|adhesion|junction|integrin|basement|cell cell|glycosaminoglycan|heparin",
  "Splicing / transcription"      = "splic|snrnp|transcription|mrna|rna process|polymerase|chromatin|histone",
  "Lysosome / lipid / metabolism" = "lysosom|vacuol|lipid|cholesterol|steroid|sterol|isoprenoid|terpenoid|metabol|biosynth|glycosyl|carbohydrate",
  "Development / morphogenesis"   = "development|morphogenesis|differentiation|kidney|renal|bone|gland|mesenchym|ossification|mineraliz|glomerul",
  "Immune / stress signalling"    = "interferon|inflamm|immune|cytokine|stress|p53|kras|tnfa|nf-?kb|uv response|signaling",
  "Ion transport / receptors"     = "channel|chloride|ion transport|receptor activity|gated|depolariz|potassium|calcium")
THEME_COL <- c("Ribosome / translation" = "#1F77B4", "Mitochondria / OXPHOS" = "#D62728", "Synapse / neuronal" = "#2CA02C",
               "Muscle / cell fusion" = "#FF7F0E", "Matrix / adhesion" = "#8C564B", "Splicing / transcription" = "#9467BD",
               "Lysosome / lipid / metabolism" = "#BCBD22", "Development / morphogenesis" = "#17BECF",
               "Immune / stress signalling" = "#E377C2", "Ion transport / receptors" = "#7F7F7F", "Other" = "#C8C8C8")

## ---- reference pathways (those you highlight) used by 05_novel_pathways.R ----
REFS <- list(
  "FXN axis + RNA/translation" = c("GOBP_METALLO_SULFUR_CLUSTER_ASSEMBLY", "GOCC_IRON_SULFUR_CLUSTER_ASSEMBLY_COMPLEX",
      "GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY", "GOBP_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
      "GOBP_MITOCHONDRIAL_TRANSLATION", "GOBP_TRICARBOXYLIC_ACID_CYCLE", "GOBP_MRNA_PROCESSING", "GOBP_RNA_SPLICING", "GOBP_RIBOSOME_BIOGENESIS"),
  "Inflammation" = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE", "GOBP_CHEMOKINE_PRODUCTION"),
  "Neuronal" = c("GOBP_AXONAL_TRANSPORT", "GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS"),
  "Lipids" = c("GOBP_CHOLESTEROL_BIOSYNTHETIC_PROCESS_VIA_DESMOSTEROL"),
  "Lysosome / oxidative stress" = c("GOCC_LYSOSOMAL_LUMEN", "GOBP_REGULATION_OF_LYSOSOMAL_LUMEN_PH", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS"))

## ---- helpers shared by all scripts -----------------------------------
clean <- function(x) tolower(gsub("_", " ", sub("^HALLMARK_", "H: ", sub("^(GOBP|GOCC|GOMF)_", "", x))))
theme_of <- function(p) vapply(clean(p), function(s) { for (th in names(THEMES)) if (grepl(THEMES[[th]], s)) return(th); "Other" }, "")

# Expression (VST) + sample table of one cell type: control + FA (from CON_A dds) and FA + edited (from CON_C dds).
# Gene identifiers are converted to symbols with the DEG table (GENEID -> GENENAME).
# FAKE=1 (environment variable) only builds simulated data to test the layout (needs fake_genes).
load_expr_meta <- function(cell, fake_genes = NULL) {
  rf <- read.csv(file.path(RES_DIR, cell, CON_A, paste0(CON_A, "_DEG_results.csv")), stringsAsFactors = FALSE)
  if (Sys.getenv("FAKE") == "1") {
    set.seed(1); smp <- c("C6718", "C6719", "C7522", "F4259", "F4676", "F86_R", "F4259ed", "F4676ed", "F86ed")
    ids <- rf$GENEID[rf$GENENAME %in% fake_genes]
    expr <- matrix(rnorm(length(ids) * 9, 8, 1), length(ids), 9, dimnames = list(ids, smp))
    meta <- data.frame(sample = smp, condition = c(rep("control", 3), rep("FA", 3), rep("FA_edited", 3)))
  } else {
    suppressPackageStartupMessages(library(DESeq2))
    get <- function(con) { d <- readRDS(file.path(RES_DIR, cell, con, paste0(con, "_dds.rds"))); list(cnt = counts(d), meta = as.data.frame(colData(d))) }
    parts <- lapply(c(CON_A, CON_C), get)                 # A: control + FA ; C: FA + edited
    allg <- Reduce(intersect, lapply(parts, function(p) rownames(p$cnt)))
    cnt <- do.call(cbind, lapply(parts, function(p) p$cnt[allg, ]))
    meta <- do.call(rbind, lapply(parts, function(p) data.frame(sample = rownames(p$meta), condition = as.character(p$meta$condition))))
    keep <- !duplicated(meta$sample); cnt <- cnt[, keep]; meta <- meta[keep, ]
    dds <- DESeqDataSetFromMatrix(round(cnt), data.frame(row.names = meta$sample, condition = meta$condition), ~1)
    expr <- assay(vst(dds, blind = TRUE))
  }
  if (mean(rownames(expr) %in% rf$GENEID) > 0.5) {
    sym <- rf$GENENAME[match(rownames(expr), rf$GENEID)]; ok <- !is.na(sym) & sym != "" & !duplicated(sym)
    expr <- expr[ok, ]; rownames(expr) <- sym[ok] }
  meta$group <- ifelse(grepl("edit", meta$condition, ignore.case = TRUE), "FA_edited", ifelse(grepl("FA", meta$condition), "FA", "control"))
  meta$line  <- sub("(_R|ed)$", "", meta$sample)           # F86_R / F86ed -> F86
  meta$group <- factor(meta$group, levels = c("control", "FA", "FA_edited"))
  list(expr = expr, meta = meta)
}
