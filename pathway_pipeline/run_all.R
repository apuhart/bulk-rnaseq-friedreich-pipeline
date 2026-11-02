#!/usr/bin/env Rscript
# Runs the whole pipeline in order. Edit 00_config.R first.   Usage:  Rscript run_all.R [--with-dds]
#   01 pathway classification (A/B/C)     -> tables + Pathways_ABC.xlsx
#   02 rescue map  /  03 theme plots  /  04 theme heatmap
#   05 new pathways following the reference trends -> New_pathways_following_reference_trends.xlsx + plots
#   06 module scores / 07 gene heatmap (need DESeq2 and the dds.rds files)  -> only with --with-dds
.a <- commandArgs(FALSE); .f <- sub("^--file=", "", .a[grepl("^--file=", .a)])
HERE <- if (length(.f)) dirname(normalizePath(.f)) else getwd()
steps <- c("01_pathway_ABC.R", "02_rescue_map.R", "03_theme_plots.R", "04_theme_heatmap.R", "05_novel_pathways.R")
if ("--with-dds" %in% commandArgs(TRUE)) steps <- c(steps, "06_module_scores.R", "07_gene_heatmap.R")
for (s in steps) {
  cat("\n>>>>>", s, "\n")
  status <- system2("Rscript", shQuote(file.path(HERE, s)))
  if (status != 0) stop("Step failed: ", s)
}
cat("\nAll done.\n")
