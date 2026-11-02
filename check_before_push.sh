#!/usr/bin/env bash
# À lancer AVANT "git push" : liste ce qui sera publié et cherche chemins personnels / identifiants de lignées.
cd "$(dirname "$0")"
git add -A -n >/dev/null 2>&1 || { echo "Pas de dépôt git ici (fais d'abord : git init)"; exit 1; }
echo "== Fichiers qui seront versionnés :"; git add -A; git ls-files; echo
echo "== Fichiers lourds (> 1 Mo) :"; git ls-files -z | xargs -0 du -k 2>/dev/null | awk '$1>1024'; echo "(vide = OK)"
echo "== Chemins personnels :"; git grep -nI "apuhart\|/media/" -- . ':!check_before_push.sh' | cut -c1-140; echo "(vide = OK)"
echo "== Identifiants de lignées (C6718 C6719 C7522 F4259 F4676 F86) :"; git grep -nIE "C6718|C6719|C7522|F4259|F4676|F86\b" -- . ':!check_before_push.sh' | cut -c1-140; echo "(vide = OK)"
