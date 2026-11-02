#!/usr/bin/env bash
#
# pdf_to_png.sh — convertit récursivement tous les PDF d'un dossier en PNG.
# Utile pour importer des figures R (ggplot/pheatmap) dans Google Slides /
# PowerPoint, qui n'acceptent pas directement le PDF comme image.
#
# Chaque PDF d'une seule page devient "<nom>.png" (à côté du PDF d'origine).
# Un PDF multi-pages devient "<nom>-1.png", "<nom>-2.png", ...
#
# Usage :
#   ./pdf_to_png.sh [dossier_racine] [dpi]
#
# Exemples :
#   ./pdf_to_png.sh                      # dossier courant, 150 dpi
#   ./pdf_to_png.sh results_v4           # un dossier précis
#   ./pdf_to_png.sh results_v4 300       # haute résolution (impression)
#
# Prérequis : poppler-utils (pdftoppm, pdfinfo)
#   Ubuntu/Debian : sudo apt install poppler-utils
#   macOS         : brew install poppler

set -euo pipefail

ROOT="${1:-.}"
DPI="${2:-150}"

if ! command -v pdftoppm &> /dev/null || ! command -v pdfinfo &> /dev/null; then
  echo "Erreur : pdftoppm/pdfinfo introuvables." >&2
  echo "Installe poppler-utils (Ubuntu: sudo apt install poppler-utils | macOS: brew install poppler)." >&2
  exit 1
fi

if [ ! -d "$ROOT" ]; then
  echo "Erreur : le dossier '$ROOT' n'existe pas." >&2
  exit 1
fi

total=0
ok=0

while IFS= read -r -d '' pdf; do
  total=$((total + 1))
  dir=$(dirname "$pdf")
  base=$(basename "$pdf" .pdf)
  prefix="$dir/$base"

  npages=$(pdfinfo "$pdf" 2>/dev/null | awk -F': *' '/^Pages:/{print $2}') || true
  npages=${npages:-1}

  if [ -z "$npages" ] || [ "$npages" -le 0 ] 2>/dev/null; then
    echo "SAUTÉ (PDF vide/cassé, 0 page) : $pdf"
    continue
  fi

  if ! timeout 30 pdftoppm -png -r "$DPI" "$pdf" "$prefix" 2>/dev/null; then
    echo "ÉCHEC : $pdf"
    continue
  fi

  # Pour un PDF d'une seule page, pdftoppm produit "<prefix>-1.png" :
  # on le renomme en "<prefix>.png" pour un nom plus simple à réutiliser.
  if [ "$npages" -le 1 ] && [ -f "${prefix}-1.png" ]; then
    mv "${prefix}-1.png" "${prefix}.png"
  fi

  ok=$((ok + 1))
  echo "OK   : $pdf"
done < <(find "$ROOT" -type f -iname "*.pdf" -print0)

echo ""
echo "Terminé : $ok/$total PDF convertis en PNG (${DPI} dpi) dans '$ROOT'."
