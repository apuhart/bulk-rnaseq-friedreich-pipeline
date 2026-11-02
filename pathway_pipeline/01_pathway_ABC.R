#!/usr/bin/env Rscript
# Genome-wide pathway classification across contrasts A / B / C
#   A = control_vs_FA        (FA vs ctrl)
#   B = control_vs_FA_edited (edited vs ctrl)
#   C = FA_edited_vs_FA      (FA vs edited)
# Input : <RES_DIR>/<cell>/<contrast>/<contrast>_<coll>_fgsea.csv  (cols: pathway, NES, padj, size, leadingEdge)
# Output: <OUT_DIR>/Pathways_ABC_<cell>_all.csv, _candidates.csv, _clusters.csv, shared_iPSC_iNeurons.csv
#         + Pathways_ABC.xlsx if openxlsx is installed.

# ---- shared configuration (00_config.R lives next to this script) ----
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()      # in RStudio: setwd() to the script folder first
source(file.path(HERE, "00_config.R"))

## ---- parameters specific to this script ------------------------------
JACCARD   <- 0.4                        # leading-edge overlap to merge pathways into one cluster
EMT_FRAC  <- 0.25                       # flag pathway if >= this fraction of its leading edge is in the Hallmark EMT leading edge
## ---------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

read_con <- function(cell, con) {
  do.call(rbind, lapply(COLLS, function(cl) {
    f <- file.path(RES_DIR, cell, con, sprintf("%s_%s_fgsea.csv", con, cl))
    if (!file.exists(f)) stop("Missing: ", f)
    d <- read.csv(f, stringsAsFactors = FALSE)
    data.frame(collection = cl, pathway = d$pathway, size = d$size, NES = d$NES,
               padj = d$padj, leadingEdge = d$leadingEdge, stringsAsFactors = FALSE)
  }))
}

classify <- function(m) {
  A <- m$padj_A < PADJ; B <- m$padj_B < PADJ; C <- m$padj_C < PADJ
  A[is.na(A)] <- FALSE; B[is.na(B)] <- FALSE; C[is.na(C)] <- FALSE
  sA <- sign(m$NES_A); sB <- sign(m$NES_B); sC <- sign(m$NES_C)
  k <- rep("n.s.", nrow(m))
  k[A | B | C] <- "Other"
  k[!A & B & C] <- "Edited-vs-ctrl shift"
  k[A & !B & !C & sC == sA] <- "Intermediate"
  k[A & C & sC != sA] <- "Discordant (A vs C)"
  k[A & B & !C & sB == sA] <- "Persistent"
  k[A & C & B & sB != sA & sC == sA] <- "Over-corrected"
  k[A & C & B & sB == sA & sC == sA] <- "Partly rescued"
  k[A & C & !B & sC == sA] <- "Rescued"
  k
}

jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
split_le <- function(x) strsplit(x, ";", fixed = TRUE)

cluster_greedy <- function(d) {            # d sorted by priority; returns cluster rep for each row
  les <- split_le(d$leadingEdge_A); reps <- integer(0); rep_of <- character(nrow(d))
  for (i in seq_len(nrow(d))) {
    hit <- NA
    for (r in reps) if (jac(les[[i]], les[[r]]) > JACCARD) { hit <- r; break }
    if (is.na(hit)) { reps <- c(reps, i); rep_of[i] <- d$pathway[i] } else rep_of[i] <- d$pathway[hit]
  }
  rep_of
}

res <- list()
for (cell in CELLS) {
  A <- read_con(cell, CON_A); B <- read_con(cell, CON_B); C <- read_con(cell, CON_C)
  names(A)[4:6] <- c("NES_A", "padj_A", "leadingEdge_A")
  names(B)[4:6] <- c("NES_B", "padj_B", "leadingEdge_B")
  names(C)[4:6] <- c("NES_C", "padj_C", "leadingEdge_C")
  m <- Reduce(function(x, y) merge(x, y, by = c("collection", "pathway"), all = FALSE),
              list(A, B[, -3], C[, -3]))           # drop duplicate size column
  m$class <- classify(m)
  m$FXN_dependent <- m$class %in% c("Rescued", "Partly rescued", "Over-corrected")

  # EMT / mesenchymal overlap (leading edge of A vs Hallmark EMT leading edge of A)
  emt <- split_le(m$leadingEdge_A[m$pathway == EMT_NAME])
  emt <- if (length(emt)) emt[[1]] else character(0)
  m$EMT_overlap <- vapply(split_le(m$leadingEdge_A), function(l) length(intersect(l, emt)) / max(1, length(l)), 0)
  m$EMT_flag <- m$EMT_overlap >= EMT_FRAC | m$pathway == EMT_NAME

  # candidates = FXN-dependent, size window, not EMT-driven
  cand <- m[m$FXN_dependent & m$size >= SIZE_MIN & m$size <= SIZE_MAX & !m$EMT_flag, ]
  cand <- cand[order(-abs(cand$NES_C)), ]
  cand$cluster <- if (nrow(cand)) cluster_greedy(cand) else character(0)
  m$cluster <- cand$cluster[match(paste(m$collection, m$pathway), paste(cand$collection, cand$pathway))]

  # one row per cluster (representative = largest |NES_C|), with members listed
  clus <- do.call(rbind, lapply(split(cand, cand$cluster), function(g) {
    g <- g[order(-abs(g$NES_C)), ]
    data.frame(cluster_rep = g$pathway[1], class = g$class[1], NES_A = g$NES_A[1], NES_B = g$NES_B[1],
               NES_C = g$NES_C[1], padj_C = g$padj_C[1], n_pathways = nrow(g),
               members = paste(head(g$pathway, 8), collapse = " | "), stringsAsFactors = FALSE)
  }))
  clus <- clus[order(-abs(clus$NES_C)), ]
  m$leadingEdge_B <- NULL            # leadingEdge_C is kept (used by novel_pathways.R)
  res[[cell]] <- list(all = m, cand = cand, clus = clus)
  cat(sprintf("\n== %s: %d pathways | %d FXN-dependent | %d after size/EMT filter | %d clusters\n",
              cell, nrow(m), sum(m$FXN_dependent), nrow(cand), nrow(clus)))
  print(table(m$class))
  write.csv(m[m$class != "n.s.", ], file.path(OUT_DIR, sprintf("Pathways_ABC_%s_all.csv", cell)), row.names = FALSE)
  write.csv(cand, file.path(OUT_DIR, sprintf("Pathways_ABC_%s_candidates.csv", cell)), row.names = FALSE)
  write.csv(clus, file.path(OUT_DIR, sprintf("Pathways_ABC_%s_clusters.csv", cell)), row.names = FALSE)
}

# shared across cell types: A significant with the same sign (iPSC & iNeurons), fibro A added as a column
a <- res[[CELLS[1]]]$all; b <- res[[CELLS[2]]]$all
s <- merge(a[, c("collection", "pathway", "size", "NES_A", "padj_A", "class", "EMT_flag")],
           b[, c("collection", "pathway", "NES_A", "padj_A", "class", "EMT_flag")],
           by = c("collection", "pathway"), suffixes = paste0("_", CELLS))
if (!is.null(FIBRO)) {
  fb <- read_con(FIBRO, CON_A)[, c("collection", "pathway", "NES", "padj")]
  names(fb)[3:4] <- c("NES_fibro", "padj_fibro")
  s <- merge(s, fb, by = c("collection", "pathway"), all.x = TRUE)
}
na <- paste0("NES_A_", CELLS); pa <- paste0("padj_A_", CELLS)
s <- s[s[[pa[1]]] < PADJ & s[[pa[2]]] < PADJ & sign(s[[na[1]]]) == sign(s[[na[2]]]), ]
s <- s[order(-(abs(s[[na[1]]]) + abs(s[[na[2]]]))), ]
write.csv(s, file.path(OUT_DIR, "Pathways_ABC_shared_iPSC_iNeurons.csv"), row.names = FALSE)
cat(sprintf("\nShared (A sig, same sign, %s & %s): %d\n", CELLS[1], CELLS[2], nrow(s)))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()
  for (cell in CELLS) {
    openxlsx::addWorksheet(wb, paste0(cell, "_clusters")); openxlsx::writeData(wb, paste0(cell, "_clusters"), res[[cell]]$clus)
    openxlsx::addWorksheet(wb, paste0(cell, "_candidates")); openxlsx::writeData(wb, paste0(cell, "_candidates"), res[[cell]]$cand[, setdiff(names(res[[cell]]$cand), c("leadingEdge_A", "leadingEdge_C"))])
    openxlsx::addWorksheet(wb, paste0(cell, "_all")); openxlsx::writeData(wb, paste0(cell, "_all"), res[[cell]]$all[res[[cell]]$all$class != "n.s.", setdiff(names(res[[cell]]$all), c("leadingEdge_A", "leadingEdge_C"))])
  }
  openxlsx::addWorksheet(wb, "shared"); openxlsx::writeData(wb, "shared", s)
  for (sh in names(wb)) openxlsx::freezePane(wb, sh, firstRow = TRUE)
  openxlsx::saveWorkbook(wb, file.path(OUT_DIR, "Pathways_ABC.xlsx"), overwrite = TRUE)
} else message("openxlsx not installed -> CSV only (install.packages('openxlsx') for the Excel file)")
cat("\nDone ->", OUT_DIR, "\n")
