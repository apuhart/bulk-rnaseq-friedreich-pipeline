PATHWAY PIPELINE (FA multi-omic RNA-seq)
1. Edit 00_config.R (RES_DIR, cell types, thresholds, themes, reference pathways). Nothing else needs editing.
2. Rscript run_all.R            -> scripts 01-05 (need only the fgsea CSVs)
   Rscript run_all.R --with-dds -> + 06 module scores and 07 gene heatmap (need DESeq2 and the *_dds.rds files)
Scripts can also be run one by one (Rscript 05_novel_pathways.R), in numeric order the first time.
Outputs: <RES_DIR>/pathway_ABC/ (tables, figures; novel/ sub-folder for script 05).
Packages: ggplot2, ggrepel, patchwork, igraph, openxlsx (optional, Excel), DESeq2 (06-07 only).
