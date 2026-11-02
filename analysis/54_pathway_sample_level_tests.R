#!/usr/bin/env Rscript
# Replicate-aware pathway tests (complements fgsea; same sample sheet as module_scores.R)
#  1) module score per sample (mean z-scored VST over ALL genes of the set) + t-test on the 9 samples
#       A  FA vs control      Welch t-test (3 vs 3)
#       B  edited vs control  Welch t-test (3 vs 3)
#       C  FA vs edited       paired t-test, pairs = same line (3 pairs)
#  2) rotation test limma::fry on voom(counts) (accounts for gene-gene correlation and the small n)
#  3) lines_concordant = number of lines (out of 3) whose effect has the sign of the mean effect
# Output: <OUT_DIR>/pathway_sample_tests_<CELL>.csv   (effect = FA minus other group, like the heatmaps)
suppressPackageStartupMessages({library(DESeq2); library(limma); library(edgeR); library(msigdbr)})

## ----------------------------- CONFIG --------------------------------
RES_DIR <- Sys.getenv("RNA_V4_ROOT", file.path(getwd(), "results_v4"))   # <- EDIT
CELL    <- "iNeurons"           # "iNeurons" or "iPSC"
OUT_DIR <- file.path(RES_DIR, "pathway_ABC")
CON_A <- "control_vs_FA"; CON_C <- "FA_edited_vs_FA"
MIN_SIZE <- 15; MAX_SIZE <- 300
# pathways to test: all of Pathways_ABC_<CELL>_all.csv (FXN-related and others); set to a character vector to restrict
ONLY <- NULL
## ---------------------------------------------------------------------

# ---- counts + sample sheet (same logic as module_scores.R) ----
get <- function(con) { d <- readRDS(file.path(RES_DIR, CELL, con, paste0(con, "_dds.rds")))
  list(cnt = counts(d), meta = as.data.frame(colData(d))) }
parts <- lapply(c(CON_A, CON_C), get)
allg  <- Reduce(intersect, lapply(parts, function(p) rownames(p$cnt)))
cnt   <- do.call(cbind, lapply(parts, function(p) p$cnt[allg, ]))
meta  <- do.call(rbind, lapply(parts, function(p) { m <- p$meta; data.frame(sample = rownames(m), condition = as.character(m$condition)) }))
keep  <- !duplicated(meta$sample); cnt <- cnt[, keep]; meta <- meta[keep, ]
meta$group <- ifelse(grepl("edit", meta$condition, ignore.case = TRUE), "edited", ifelse(grepl("FA", meta$condition), "FA", "control"))
meta$line  <- sub("(_R|ed)$", "", meta$sample)
stopifnot(all(table(meta$group) == 3))

# ---- gene symbols ----
rf <- read.csv(file.path(RES_DIR, CELL, CON_A, paste0(CON_A, "_DEG_results.csv")), stringsAsFactors = FALSE)
toSym <- function(x) { s <- rf$GENENAME[match(x, rf$GENEID)]; s }
vst_m <- assay(vst(DESeqDataSetFromMatrix(round(cnt), data.frame(row.names = meta$sample, g = meta$group), ~1), blind = TRUE))
if (mean(rownames(vst_m) %in% rf$GENEID) > 0.5) rownames(cnt) <- rownames(vst_m) <- toSym(rownames(vst_m))
ok <- !is.na(rownames(vst_m)) & rownames(vst_m) != "" & !duplicated(rownames(vst_m)); vst_m <- vst_m[ok, ]; cnt <- cnt[ok, ]
z <- t(scale(t(vst_m[apply(vst_m, 1, sd) > 0, ])))

# ---- gene sets (full sets, not leading edges: avoids circularity) ----
mdb <- function(cat, sub = NULL) { if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = cat) else msigdbr(species = "Homo sapiens", collection = cat, subcollection = sub) }
m <- rbind(mdb("H"), mdb("C5", "GO:BP"), mdb("C5", "GO:CC"), mdb("C5", "GO:MF"))
sets <- split(m$gene_symbol, m$gs_name)
abc  <- read.csv(file.path(OUT_DIR, sprintf("Pathways_ABC_%s_all.csv", CELL)), stringsAsFactors = FALSE)
paths <- if (is.null(ONLY)) unique(abc$pathway) else ONLY
sets <- lapply(sets[intersect(paths, names(sets))], function(g) intersect(unique(g), rownames(z)))
sets <- sets[lengths(sets) >= MIN_SIZE & lengths(sets) <= MAX_SIZE]
cat(length(sets), "pathways tested\n")

# ---- 1) module scores + t-tests ----
S <- t(vapply(sets, function(g) colMeans(z[g, , drop = FALSE]), numeric(ncol(z)))); colnames(S) <- colnames(z)
idx <- function(g) meta$sample[meta$group == g]
FA <- idx("FA"); CT <- idx("control"); ED <- idx("edited")
# pair FA and edited by line
FA <- FA[order(meta$line[match(FA, meta$sample)])]; ED <- ED[order(meta$line[match(ED, meta$sample)])]
stopifnot(identical(meta$line[match(FA, meta$sample)], meta$line[match(ED, meta$sample)]))
tA <- {a <- FA; b <- CT;  t(apply(S, 1, function(s) { r <- t.test(s[a], s[b]); c(effect = mean(s[a]) - mean(s[b]), d = (mean(s[a]) - mean(s[b])) / sd(c(s[a] - mean(s[a]), s[b] - mean(s[b]))), p = r$p.value, conc = sum(sign(s[a] - mean(s[b])) == sign(mean(s[a]) - mean(s[b])))) }))}
tB <- {a <- ED; b <- CT;  t(apply(S, 1, function(s) { r <- t.test(s[a], s[b]); c(effect = mean(s[a]) - mean(s[b]), d = (mean(s[a]) - mean(s[b])) / sd(c(s[a] - mean(s[a]), s[b] - mean(s[b]))), p = r$p.value, conc = sum(sign(s[a] - mean(s[b])) == sign(mean(s[a]) - mean(s[b])))) }))}
tC <- t(apply(S, 1, function(s) { dlt <- s[FA] - s[ED]; r <- t.test(s[FA], s[ED], paired = TRUE)
  c(effect = mean(dlt), d = mean(dlt) / max(sd(dlt), 1e-6), p = r$p.value, conc = sum(sign(dlt) == sign(mean(dlt)))) }))

# ---- 2) fry on voom ----
dge <- calcNormFactors(DGEList(counts = round(cnt)))
fryrun <- function(samples, design, coef) {
  v <- voom(dge[, samples], design, plot = FALSE)
  fry(v, index = sets, design = design, contrast = coef) }
sA <- c(FA, CT); gA <- factor(meta$group[match(sA, meta$sample)], levels = c("control", "FA"));    dA <- model.matrix(~gA)
sB <- c(ED, CT); gB <- factor(meta$group[match(sB, meta$sample)], levels = c("control", "edited")); dB <- model.matrix(~gB)
sC <- c(FA, ED); gC <- factor(meta$group[match(sC, meta$sample)], levels = c("edited", "FA")); lC <- factor(meta$line[match(sC, meta$sample)]); dC <- model.matrix(~lC + gC)
fA <- fryrun(sA, dA, 2); fB <- fryrun(sB, dB, 2); fC <- fryrun(sC, dC, ncol(dC))      # FA minus edited

mk <- function(tab, f, con) data.frame(pathway = rownames(tab), contrast = con, effect = tab[, "effect"], d = tab[, "d"],
  p_score = tab[, "p"], lines_concordant = tab[, "conc"], fry_direction = f[rownames(tab), "Direction"], fry_p = f[rownames(tab), "PValue"],
  fry_fdr = f[rownames(tab), "FDR"], fry_mixed_fdr = f[rownames(tab), "FDR.Mixed"], stringsAsFactors = FALSE)
out <- rbind(mk(tA, fA, "A"), mk(tB, fB, "B"), mk(tC, fC, "C"))
out$fdr_score <- ave(out$p_score, out$contrast, FUN = function(p) p.adjust(p, "BH"))
out$size <- lengths(sets)[out$pathway]
write.csv(out, file.path(OUT_DIR, sprintf("pathway_sample_tests_%s.csv", CELL)), row.names = FALSE)
cat("Done ->", file.path(OUT_DIR, sprintf("pathway_sample_tests_%s.csv", CELL)), "\n")
