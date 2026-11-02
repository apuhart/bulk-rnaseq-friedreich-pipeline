#!/usr/bin/env bash
# Lance les analyses APRES le pipeline principal (run_bulk_pipeline.R), depuis la racine du dossier.
# Usage :  ./run_downstream.sh  <étape>      étapes : pathway | analysis | gsea_loo | figures | all
#   RES_DIR=/chemin/vers/results_v4 ./run_downstream.sh pathway     (par défaut : ./results_v4)
set -euo pipefail
cd "$(dirname "$0")"
export RES_DIR="${RES_DIR:-$PWD/results_v4}"
export RNA_V4_ROOT="$RES_DIR" FA_RES_DIR="$RES_DIR" RNA_V4_BASE="$RES_DIR/iNeurons" PROJECT_ROOT="$PWD"
[ -d "$RES_DIR" ] || { echo "results_v4 introuvable : $RES_DIR (lance d'abord Rscript run_bulk_pipeline.R, ou donne RES_DIR=...)"; exit 1; }
step=${1:-help}
case "$step" in
  pathway)   Rscript pathway_pipeline/run_all.R ;;                       # 01-05 ; ajouter --with-dds pour 06-07
  analysis)  for s in 40_gsea_rescue_classification 42_rescue_composition_split 51_rescue_tables_ipsc_ineurons 52_paired_rescue_allgenes 54_pathway_sample_level_tests; do
               echo ">>>>> $s"; Rscript "analysis/$s.R"; done ;;
  gsea_loo)  OUT_DIR="${OUT_DIR:-gsea_loo/out}" DATASETS="${DATASETS:-fibro,iPSC,iNeurons}" Rscript gsea_loo/01_gsea_loo.R ;;
  figures)   python3 gsea_loo/02_figures_and_table.py "${GSEA_DIR:-gsea_loo/out}" "${FIG_DIR:-plot_enrichment}" ;;
  all)       "$0" pathway; "$0" analysis; "$0" gsea_loo; "$0" figures ;;
  *)         sed -n '2,5p' "$0" ;;
esac
