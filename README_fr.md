# Pipeline bulk RNA-seq — notes en français

Ordre d'exécution (depuis la racine de ce dossier ; copie d'abord `config/config.example.R` en `config/config.R` si tu pars d'un clone git) :

| # | Commande | Rôle | Durée |
|---|---|---|---|
| 0 | `Rscript -e 'source("config/packages.R"); install_missing_packages()'` | installe les paquets R | une fois |
| 1 | `Rscript run_bulk_pipeline.R` | Salmon → tximport → DESeq2 → fgsea → QC, cross-datasets → `results_v4/` | long |
| 2 | `./run_downstream.sh pathway` | classification A/B/C, rescue map, thèmes, nouveaux pathways (`pathway_pipeline/`) | court |
| 3 | `./run_downstream.sh analysis` | scripts 40, 42, 51, 52, 54 (rescue, composition, tables, tests au niveau échantillon) | moyen |
| 4 | `./run_downstream.sh gsea_loo` | GSEA leave-one-line-out, 3 modèles, tous contrastes → `gsea_loo/out/` | long |
| 5 | `./run_downstream.sh figures` | figures + `Pathways_leave_one_line_out.xlsx` → `plot_enrichment/` | 1 min |

Si tes résultats existent déjà ailleurs : `RES_DIR=/chemin/results_v4 ./run_downstream.sh pathway`.
Les chemins Salmon de `config/config.R` se changent par variables d'environnement (`FA_SALMON_FIBRO`, `FA_SALMON_iPSC`, `FA_SALMON_iNeurons`) ou directement dans le fichier.
Le GTF (`data/Homo_sapiens.GRCh38.112.gtf.gz`, 54 Mo) n'est pas inclus : les caches `.rds` de `data/` suffisent tant qu'ils existent ; sinon copie le GTF dans `data/`.

## Contenu
- `run_bulk_pipeline.R` (ex `mainPatchSave.R`), `config/`, `R/` : pipeline principal.
- `pathway_pipeline/` : analyses de pathways configurées par `00_config.R`.
- `analysis/` : scripts 40, 42, 51, 52, 54 (issus du dossier scripts ChIP/intégration, mais purement bulk).
- `gsea_loo/` : GSEA leave-one-line-out + toutes les figures et le tableau xlsx.

## Convention de signe
A = FA vs control ; B = edited vs control ; C = FA vs edited (NES/log2FC > 0 = plus haut dans le premier groupe).
Attention : le dossier de contraste `FA_edited_vs_FA` porte un nom inversé par rapport au signe (log2FC = FA − édité).
