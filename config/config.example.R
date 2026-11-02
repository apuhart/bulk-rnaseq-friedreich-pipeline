# config.example.R — MODÈLE de configuration (identifiants anonymisés, chemins factices).
# Copie ce fichier en config/config.R puis remplace : chemins Salmon, noms d'échantillons, lignées, sexe.
# config/config.R est ignoré par git (.gitignore) : il contient tes vrais identifiants.

# ============================================================
# config.R — Central configuration for the RNA-seq pipeline
# ============================================================
# Edit this file before running any script.

# ----------------------------
# Racine du projet (définie par mainPatchSave.R ; repli = dossier courant)
# ----------------------------
if (!exists("PROJECT_ROOT")) PROJECT_ROOT <- getwd()

# ----------------------------
# Output directory
# results_v4 -> n'écrase pas les anciens résultats (results/, results_v2/) : comparaison possible
# ----------------------------
BASE_OUT <- file.path(PROJECT_ROOT, "results_v4")

# Sauvegarder l'objet DESeqDataSet de chaque contraste (RDS, quelques dizaines de Mo)
SAVE_DDS <- TRUE

# ----------------------------
# Ensembl mirror
# Options: "useast", "uswest", "asia", "www"
# ----------------------------
ENSEMBL_MIRROR <- "https://useast.ensembl.org"   # (n'est plus utilisé : tx2gene via fichier/GTF)

# ----------------------------
# tx2gene
# Priorité 1 : fichier nf-core (garantit la cohérence avec l'annotation utilisée par Salmon)
#   ex. ".../results_fibro/star_salmon/salmon_tx2gene.tsv"  (3 colonnes, sans en-tête)
# Priorité 2 : GTF Ensembl local (téléchargé si absent). Un cache .rds est créé à côté.
# ----------------------------
TX2GENE_FILE <- NULL
GTF_PATH     <- file.path(PROJECT_ROOT, "data", "Homo_sapiens.GRCh38.112.gtf.gz")
# ----------------------------
# DESeq2 / filtering
# ----------------------------
# Filtre : un gène est conservé s'il a >= MIN_COUNTS comptages bruts dans >= MIN_SAMPLES échantillons
# (avant : rowSums(counts) > 3, beaucoup trop permissif -> ~25 000 gènes, LFC aberrants)
MIN_COUNTS      <- 10
MIN_SAMPLES     <- 3      # = taille du plus petit groupe (3 réplicats)
CLASSIFY_RATIO  <- 0.5    # classification FA-spécifique / clone-édité : un effet est "proche du contrôle" s'il vaut < 50 % de la différence FA vs édité
LFC_THRESHOLD   <- 1    # |log2FC| threshold for volcano labelling
PADJ_THRESHOLD  <- 0.05   # adjusted p-value threshold (aussi utilisé comme alpha de results())

# ----------------------------
# fgsea
# ----------------------------
FGSEA_NPERM     <- 10000  # DÉPRÉCIÉ : ignoré (fgseaMultilevel n'utilise plus nperm) ; conservé pour compatibilité
FGSEA_MIN_SIZE  <- 15     # taille min d'un gene set (après intersection avec les gènes classés)
FGSEA_MAX_SIZE  <- 500    # taille max
FGSEA_SEED      <- 42     # reproductibilité
FGSEA_PADJ_CUT  <- 0.3    # pre-filter before dotplot
FGSEA_TOP_N     <- 30     # top N pathways to plot

# MSigDB collections to run  (list of lists: category + optional subcategory)
FGSEA_COLLECTIONS <- list(
  list(category = "H",  subcategory = NULL),
  list(category = "C5", subcategory = "BP"),
  list(category = "C5", subcategory = "CC"),
  list(category = "C5", subcategory = "MF")
)

# ----------------------------
# Heatmap
# ----------------------------
HEATMAP_TOP_N   <- 100    # top N DEGs (by padj) for heatmap

# ----------------------------
# Ajustement sur le sexe (contrastes NON appariés : control_vs_FA, control_vs_FA_edited)
#   sex_adjust = TRUE  : design ~ sex + condition
#   sex_adjust = FALSE : design ~ condition, et les gènes chrX / chrY / XIST sont EXCLUS
#                        du test (effets de sexe et d'inactivation de X, pas de FA)
# Le contraste apparié FA_edited_vs_FA (~ line + condition) n'est pas concerné.
# Pourquoi FALSE pour iPSC et iNeurons : avec n = 3 vs 3, la covariable sex coûte un degré de
# liberté sans corriger le décalage global des gènes X (testé sur ces données).
# Les gènes exclus sont listés dans <contraste>_design.txt.
# ----------------------------
# Datasets
# ----------------------------
DATASETS <- list(
  
  fibro = list(
    samples    = c("CTRL1", "CTRL2", "CTRL3", "FA1", "FA2", "FA3"),
    conditions = c("control", "control", "control", "FA", "FA", "FA"),
    sex        = c("M", "F", "F", "M", "M", "F"),
    # line = identifiant de la lignée/patient d'origine. Les lignées éditées
    # partagent la "line" de leur lignée FA parentale -> permet le design apparié.
    line       = c("CTRL1", "CTRL2", "CTRL3", "FA1", "FA2", "FA3"),
    sex_adjust = TRUE,    # 1M/2F contre 2M/1F : sexe équilibré entre groupes
    path       = Sys.getenv("FA_SALMON_FIBRO", "/path/to/FIBRO/star_salmon/")
  ),
  
  iNeurons = list(
    samples    = c("CTRL1", "CTRL2", "CTRL3", "FA1", "FA1ed", "FA2", "FA2ed", "FA3ed", "FA3_R"),
    conditions = c("control", "control", "control", "FA", "FA_edited", "FA", "FA_edited", "FA_edited", "FA"),
    # CTRL1 corrigé M -> F le 2026-09-22 : XIST tres eleve (3847) et 0 gene Y
    # dans ce dataset (voir iNeurons_sex_check.csv), incompatible avec "M".
    # Note : le meme echantillon "CTRL1" est XY dans le dataset fibro (voir
    # fibro_sex_check.csv) -> possible probleme d'identite d'echantillon entre
    # les fibroblastes et la lignee iPSC/iNeurons portant ce nom, a documenter.
    sex        = c("F", "F", "F", "M", "M", "M", "M", "F", "F"),
    # CTRL1 CTRL2 CTRL3 FA1 FA1ed FA2 FA2ed FA3ed FA3_R
    line       = c("CTRL1", "CTRL2", "CTRL3", "FA1", "FA1", "FA2", "FA2", "FA3", "FA3"),
    sex_adjust = FALSE,   # cf. note ci-dessus. QC : XIST/gènes Y de CTRL1 sont incohérents avec "M"
    path       = Sys.getenv("FA_SALMON_iNeurons", "/path/to/iNeurons/star_salmon/")
  ),
  
  iPSC = list(
    samples    = c("CTRL1", "CTRL2", "CTRL3", "ED86", "ED4259", "ED4676", "FA1", "FA2", "FA3"),
    conditions = c("control", "control", "control", "FA_edited", "FA_edited", "FA_edited", "FA", "FA", "FA"),
    # CTRL1 corrigé M -> F le 2026-09-22 : XIST tres eleve (3715) et 0 gene Y
    # dans ce dataset (voir iPSC_sex_check.csv), incompatible avec "M".
    # Note : le meme echantillon "CTRL1" est XY dans le dataset fibro (voir
    # fibro_sex_check.csv) -> possible probleme d'identite d'echantillon entre
    # les fibroblastes et la lignee iPSC/iNeurons portant ce nom, a documenter.
    sex        = c("F", "F", "F", "F", "M", "M", "M", "M", "F"),
    # CTRL1 CTRL2 CTRL3 ED86 ED4259 ED4676 FA1 FA2 FA3
    line       = c("CTRL1", "CTRL2", "CTRL3", "FA3", "FA1", "FA2", "FA1", "FA2", "FA3"),
    sex_adjust = FALSE,   # cf. note ci-dessus. QC : XIST/gènes Y de CTRL1 sont incohérents avec "M"
    path       = Sys.getenv("FA_SALMON_iPSC", "/path/to/iPSC/star_salmon/")
  )
)

# ----------------------------
# Marqueurs de différenciation (QC : pureté / maturité neuronale, identité cellulaire)
# ----------------------------
DIFF_MARKERS <- list(
  pluripotence = c("POU5F1", "NANOG", "SOX2", "LIN28A", "DPPA4"),
  progeniteurs = c("NES", "PAX6", "SOX1", "VIM", "HES1", "DCX"),
  neurones     = c("MAP2", "TUBB3", "RBFOX3", "SYP", "SNAP25", "STMN2", "NEFL", "GAP43"),
  glutamate    = c("SLC17A7", "SLC17A6", "GRIN1", "GRIA2"),
  gaba         = c("GAD1", "GAD2", "SLC32A1"),
  forebrain    = c("FOXG1", "EMX1", "SATB2", "TBR1"),
  gliaux       = c("GFAP", "S100B", "AQP4", "OLIG2", "SOX10"),
  mesenchyme   = c("COL1A1", "COL1A2", "FN1", "THY1", "ACTA2", "TAGLN")
)

# ----------------------------
# Genes of interest
# ----------------------------
GENES_OF_INTEREST <- unique(c(
  # Fe-S biogenesis — mitochondrial ISC
  "ISCU", "NFS1", "LYRM4", "NDUFAB1", "FXN", "FDX2", "FDXR", "HSPA9", "HSCB", "GLRX5",
  "ABCB7", "ISCA1", "ISCA2", "IBA57", "NFU1",
  # Fe-S biogenesis — CIA pathway
  "CFDP1", "NBP35", "CIAPIN1", "NDOR1", "NARFL", "CIAO1", "CIAO2B", "CIAO2A", "MMS19",
  # Fe-S client proteins — DNA metabolism
  "PRIM2", "DDX11", "DNA2", "BRIP1", "RTEL1", "ERCC2", "POLA1", "POLD1", "POLE",
  "MUTYH", "NTHL1", "KIF4A",
  # Fe-S client proteins — OXPHOS / metabolism
  "NDUFS1", "NDUFS7", "NDUFS8", "NDUFV1", "SDHB", "UQCRFS1", "ACO2", "LIAS", "FECH", "CISD1",
  "ABCE1", "ACO1", "GPAM", "DPYD",
  # ISR / ATF4 axis
  "ATF4", "ATF3", "ASNS", "DDIT3", "GDF15", "TRIB3",
  # One-carbon / serine metabolism
  "PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2", "MTHFD1L", "MTHFD2", "GSS",
  "ALDH1L1", "ALDH1L2",
  # AMPK / mTOR / ferroptosis
  "PRKAA1", "STK11", "MTOR", "GPX4", "SLC25A28", "IREB2", "SLC11A2",
  # NF-kB
  "RELA", "RELB", "NFKB1", "NFKB2", "REL",
  # cGAS-STING
  "TMEM173", "CGAS", "ZBP1", "RIPK1", "RIPK3", "STING1",
  # Cytokines / ECM
  "IL17A", "MMP1", "MMP2", "MMP9",
  # Interferon
  "IRF3", "IRF7", "IRF9", "STAT1", "STAT2", "STAT3",
  # Stress / transcription
  "TP53", "FOXO3", "HIF1A",
  # AP-1
  "JUN", "FOS", "ATF2",
  # Adaptive immunity
  "RORC", "BATF"
))

# ----------------------------
# Gènes de la voie Fe-S / stress lié à la frataxine (liste curée, fournie par
# l'utilisateur le 2026-09-22). Chaque vecteur de catégorie est nommé
# "libellé donné" = "symbole HGNC actuel utilisé dans tx2gene" : les alias
# historiques (HSC20, IOP1, CHLR1, XPD, CFD1, NBP35, POLA, POLE1, MitoNEET,
# IRP1, GPAT, AMPK, LKB1, NRF2, IRP2, DMT1, DDIT3/CHOP...) sont convertis vers
# le symbole HGNC actuel, mais le libellé d'origine reste affiché dans les
# exports (voir R/fa_pathway_genes.R -> run_fa_pathway_gene_summary()).
# ----------------------------
FA_PATHWAY_GENES <- list(
  "Fe-S biogenesis - machinerie mitochondriale" = c(
    "ISCU" = "ISCU", "NFS1" = "NFS1", "ISD11" = "LYRM4", "ACPM" = "NDUFAB1",
    "Frataxin (FXN)" = "FXN", "Ferredoxin (FDX2)" = "FDX2", "FDXR" = "FDXR",
    "HSPA9" = "HSPA9", "HSC20" = "HSCB", "GLRX5" = "GLRX5", "ABCB7" = "ABCB7",
    "ISCA1" = "ISCA1", "ISCA2" = "ISCA2", "IBA57" = "IBA57", "NFU1" = "NFU1"
  ),
  "Fe-S biogenesis - machinerie cytosolique (CIA)" = c(
    "CFD1" = "NUBP2", "NBP35" = "NUBP1", "CIAPIN1" = "CIAPIN1", "NDOR1" = "NDOR1",
    "IOP1" = "NARFL", "CIAO1" = "CIAO1", "CIAO2B" = "CIAO2B", "MMS19" = "MMS19",
    "CIAO2A" = "CIAO2A"
  ),
  "Proteine Fe-S nucleaire - metabolisme ADN" = c(
    "PRIM2" = "PRIM2", "CHLR1" = "DDX11", "DNA2" = "DNA2", "FANCJ" = "BRIP1",
    "RTEL1" = "RTEL1", "XPD" = "ERCC2", "POLA" = "POLA1", "POLD1" = "POLD1",
    "POLE1" = "POLE", "MUTYH" = "MUTYH", "NTHL1" = "NTHL1", "KIF4A" = "KIF4A"
  ),
  "Proteine Fe-S mitochondriale (OXPHOS / TCA / heme)" = c(
    "NDUFS8" = "NDUFS8", "NDUFS7" = "NDUFS7", "NDUFS1" = "NDUFS1", "NDUFV2" = "NDUFV2",
    "NDUV1" = "NDUFV1", "SDHB" = "SDHB", "UQCRFS1" = "UQCRFS1", "ACO2" = "ACO2",
    "LIAS" = "LIAS", "FECH" = "FECH", "MitoNEET" = "CISD1"
  ),
  "Proteine Fe-S cytosolique" = c(
    "ABCE1" = "ABCE1", "IRP1 (ACO1)" = "ACO1", "GPAT (PPAT)" = "PPAT", "DPYD" = "DPYD"
  ),
  "Integrated Stress Response (ISR)" = c(
    "ATF4" = "ATF4", "ATF3" = "ATF3", "ASNS" = "ASNS", "DDIT3 (CHOP)" = "DDIT3",
    "GDF15" = "GDF15", "TRIB3" = "TRIB3"
  ),
  "Metabolisme serine / 1-carbone" = c(
    "PHGDH" = "PHGDH", "PSAT1" = "PSAT1", "PSPH" = "PSPH", "SHMT1" = "SHMT1",
    "SHMT2" = "SHMT2", "MTHFD1L" = "MTHFD1L", "MTHFD2" = "MTHFD2", "GSS" = "GSS",
    "ALDH1L1" = "ALDH1L1", "ALDH1L2" = "ALDH1L2"
  ),
  "AMPK signaling" = c(
    "AMPK (PRKAA1)" = "PRKAA1", "LKB1 (STK11)" = "STK11", "MTOR" = "MTOR"
  ),
  "Ferroptose" = c("GPX4" = "GPX4"),
  "Voie NRF2" = c("NFE2L2 (NRF2)" = "NFE2L2"),
  "Metabolisme du fer" = c(
    "SLC25A28" = "SLC25A28", "IREB2 (IRP2)" = "IREB2", "SLC11A2 (DMT1)" = "SLC11A2"
  )
)

# ----------------------------
# Gènes candidats identifiés par l'analyse rescue / classification appariée
# (pas la liste Fe-S fournie par l'utilisateur -- ici des gènes trouvés
# empiriquement comme "rescued"/"FA_specifique" dans les données, puis
# regroupés par fonction biologique probable). Les catégories sont
# volontairement différentes selon le type cellulaire d'origine : les gènes
# trouvés en iPSC pointent vers le stress oxydant / signalisation calcique /
# jonction neuromusculaire, ceux trouvés en iNeurons vers l'identité gliale /
# la protection protéique / la synaptogenèse -- cf. discussion du 2026-09-22.
# ----------------------------
CANDIDATE_GENES_PATHWAYS <- list(
  "Controle positif" = c(
    "Frataxin (FXN)" = "FXN"
  ),
  "Stress oxydant / detoxification (evidence forte)" = c(
    "GSTO2 (iPSC)" = "GSTO2", "SOD3 (iNeurons)" = "SOD3"
  ),
  "Signalisation calcique mitochondriale" = c(
    "NFATC4 (iPSC)" = "NFATC4"
  ),
  "Chaperonne / protection mitochondriale" = c(
    "CRYAB (iNeurons)" = "CRYAB"
  ),
  "Jonction neuromusculaire" = c(
    "ALG14 (iPSC)" = "ALG14", "COL13A1 (iPSC)" = "COL13A1"
  ),
  "Identite gliale (cellules de Schwann)" = c(
    "SOX10 (iNeurons)" = "SOX10"
  ),
  "Synaptogenese astrocytaire" = c(
    "SPARCL1 (iNeurons)" = "SPARCL1"
  ),
  "Guidage axonal / neurodeveloppement" = c(
    "NTN1 (iPSC)" = "NTN1"
  ),
  "Facteurs de transcription developpementaux (evidence moderee)" = c(
    "TFAP2C (iPSC)" = "TFAP2C", "SOX21 (iPSC)" = "SOX21", "KMT2C (iPSC)" = "KMT2C"
  ),
  "Reponse au stress / IGF (evidence moderee)" = c(
    "IGFBP3 (iPSC)" = "IGFBP3"
  ),
  "Stress oxydant (evidence faible, generique)" = c(
    "H19 (iNeurons)" = "H19"
  ),
  "Signalisation vasculaire - PRUDENCE, possible biais EMT/composition cellulaire" = c(
    "EDNRA (iNeurons)" = "EDNRA", "EPAS1 (iNeurons)" = "EPAS1",
    "DDR2 (iNeurons)" = "DDR2", "ANXA1 (iNeurons)" = "ANXA1"
  ),
  "Autres genes rescued (iPSC) - fonction FA non etablie" = c(
    "FLT4 (iPSC)" = "FLT4", "PKP3 (iPSC)" = "PKP3", "IRX1 (iPSC)" = "IRX1",
    "RPSA2 (iPSC)" = "RPSA2", "PROCA1 (iPSC)" = "PROCA1", "TRHDE (iPSC)" = "TRHDE",
    "ZC3H6 (iPSC)" = "ZC3H6", "PAPLN (iPSC)" = "PAPLN", "ID4 (iPSC)" = "ID4",
    "AHNAK (iPSC)" = "AHNAK", "BCHE (iPSC)" = "BCHE", "QNG1 (iPSC)" = "QNG1"
  )
)

# ----------------------------
# Neuronal keyword filter (used when neuronal_only = TRUE in fgsea)
# ----------------------------
NEURONAL_KEYWORDS <- c(
  "neuro", "synap", "axon", "dend",
  "glutamate", "gaba", "calcium",
  "neuron", "neuronal", "transmission", "plasticity",
  "mitochond", "oxidative", "respiratory", "electron transport",
  "atp", "energy", "tca", "krebs",
  "iron", "heme", "sulfur", "fe-s", "isc", "cluster",
  "oxidative stress", "ros", "reactive oxygen",
  "glutathione", "redox", "detox",
  "metabolic", "lipid", "fatty acid", "amino acid",
  "one carbon", "folate", "serine",
  "stress response", "unfolded", "upr", "er stress",
  "atf4", "chop", "ddit3", "integrated stress",
  "inflamm", "nfkb", "cytokine", "interferon", "immune",
  "apoptosis", "autophagy", "ferroptosis"
)