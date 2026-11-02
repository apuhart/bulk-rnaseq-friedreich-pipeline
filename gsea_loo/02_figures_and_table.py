#!/usr/bin/env python3
"""02_figures_and_table.py — TOUTES les figures + le tableau xlsx, à partir des seuls CSV de 01_gsea_loo.R

Usage : python3 02_figures_and_table.py  GSEA_DIR  OUT_DIR
  GSEA_DIR : dossier contenant <fibro|iPSC|iNeurons>_<A|B|C>.csv (sortie de 01_gsea_loo.R)
  OUT_DIR  : dossier de sortie (plot_enrichment/ranking, robust, corrected, persistent, partage_entre_modeles.png,
             Pathways_leave_one_line_out.xlsx)
Tout vient de la même source : mêmes NES / padj / leading edge dans les figures et dans le tableau.
Filtres des figures : padj < 0.05, 15 <= size <= 300, pathways redondants fusionnés (Jaccard des leading edges > 0.6),
sauf ceux dont le leading edge contient FXN (toujours gardés)."""
import sys, os, re, textwrap
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils.dataframe import dataframe_to_rows

GSEA, OUT = sys.argv[1], sys.argv[2]
MODELS = {'Fibroblasts': ('fibro', 'A'), 'iNeurons': ('iNeurons', 'ABC'), 'iPSC': ('iPSC', 'ABC')}
CL = {'A': 'FA vs control', 'B': 'Edited vs control', 'C': 'FA vs edited'}
SIZE_MIN, SIZE_MAX, PADJ, JACC = 15, 300, 0.05, 0.6

# ---------------------------------------------------------------- thèmes / noms
FXN = re.compile(r"iron.sulfur|iron sulfur|metallo.sulfur|sulfur cluster|\bfe.s\b|respiratory chain|electron transport|oxidative phosphorylation|nadh dehydrogenase|ubiquinone|mitochondrial (translation|gene expression|ribosome|respiratory|electron|large ribosomal|small ribosomal|protein.containing)|\bheme\b|iron ion|cellular iron|aerobic respiration|atp synthesis|proton.transporting|tricarboxylic|citric acid|ferroptosis|respirasome|oxidoreductase complex")
CATS = [('Mitochondria / respiration', r"mitochond|respirat|oxidative phosph|electron transport|atp synth|oxidoreductase complex", '#D55E00'),
 ('Gene expression (RNA, ribosome)', r"\brna\b|splic|ribosom|translation|transcription|mrna|trna|rrna|nucle(ar|olus)|polymerase", '#0072B2'),
 ('Cell cycle / DNA repair', r"cell cycle|mitotic|chromosom|dna repair|dna damage|checkpoint|dna replication|e2f|g2m|spindle|kinetochore|mitosis|dna metabolic|double.strand", '#009E73'),
 ('Matrix / development', r"matrix|collagen|extracellular|adhesion|development|morphogenesis|differentiation|ossification|cartilage|skeletal|junction|epithelial|angiogenesis|blood vessel|encapsulating|connective|mesenchym|muscle|tissue|cell migration|motility|wound", '#CC79A7'),
 ('Immune / stress', r"immun|inflamm|cytokine|interferon|tnf|leukocyte|chemokine|response to|defense|apoptot|hypoxia", '#E69F00'),
 ('Neuronal / synapse', r"synap|neuro|axon|channel|dendri|gaba|glutamat|behavior|cognition", '#56B4E9'),
 ('Lipid / lysosome / metabolism', r"lipid|sterol|cholesterol|lysosom|metabolic|biosynth|vacuol|catabolic|transport|golgi|endoplasmic|vesicle", '#999999')]
CC = {k: c for k, _, c in CATS}; CC['Other'] = '#BBBBBB'
def plain(p): return re.sub(r'^(GOBP|GOCC|GOMF)_', '', p).replace('_', ' ').lower()
def nice(p):
    if p.startswith('HALLMARK_'): return 'Hallmark: ' + p[9:].replace('_', ' ').lower()
    q = plain(p); return q[0].upper() + q[1:]
def cat(p):
    n = plain(p)
    if re.search(r"\bt cell|lymphocyte|\bb cell|helper|leukocyte|interferon|cytokine", n): return 'Immune / stress'
    for k, rx, _ in CATS:
        if re.search(rx, n): return k
    return 'Other'

# ---------------------------------------------------------------- chargement + scores LOO
def load(ds, c):
    d = pd.read_csv(os.path.join(GSEA, f'{ds}_{c}.csv'))
    runs = [x.split('|')[1] for x in d.columns if x.startswith('NES|')]
    loo = [r for r in runs if r != 'all']
    s = np.sign(d['NES|all'])
    for r in loo: d['contrib|' + r.replace('without_', '')] = s * (d['NES|all'] - d['NES|' + r])      # perte de |NES| quand on retire la lignée
    nes = d[['NES|' + r for r in runs]].values
    d['worst_case_NES_signed'] = s * np.maximum((s.values[:, None] * nes).min(axis=1), 0)       # 0 si le signe s'inverse
    d['n_LOO_padj<0.05'] = (d[['padj|' + r for r in loo]] < PADJ).sum(axis=1)
    d['sig_all'] = d['padj|all'] < PADJ
    d['same_sign_all_LOO'] = np.all([np.sign(d['NES|' + r]) == s for r in loo], axis=0)
    d['nominal_all_LOO'] = np.all([d['pval|' + r] < 0.05 for r in loo], axis=0)
    d['robust'] = d.sig_all & d.same_sign_all_LOO & d.nominal_all_LOO
    cb = d[[x for x in d.columns if x.startswith('contrib|')]].clip(lower=0)
    tot = cb.sum(axis=1)
    d['driver_line'] = ['F' + x.split('|')[1] for x in cb.idxmax(axis=1)]
    d['driver_share'] = cb.max(axis=1) / tot.where(tot > 0)
    lost = np.any([(np.sign(d['NES|' + r]) != s) | (d['pval|' + r] >= 0.05) for r in loo], axis=0)
    d['single_line_driven'] = d.sig_all & lost            # définition de 53_gsea_loo_v3.R : le signal disparaît quand on retire UNE lignée
    d['FXN_in_leading_edge'] = d.leadingEdge.map(lambda v: ('FXN' in set(v.split(';'))) if isinstance(v, str) else np.nan)
    d['FXN_related_name'] = d.pathway.map(lambda p: bool(FXN.search(plain(p))))
    d['lset'] = d.leadingEdge.map(lambda v: set(v.split(';')) if isinstance(v, str) else set())
    d['fxg'] = d['FXN_in_leading_edge'].fillna(False).astype(bool); d['fxn'] = d['FXN_related_name']
    d['NES'] = d['NES|all']; d['padj'] = d['padj|all']
    return d

def figset(d):
    d = d[(d['size'] >= SIZE_MIN) & (d['size'] <= SIZE_MAX)].copy()
    d['proc'] = d.pathway.map(cat); d.loc[d.fxn, 'proc'] = 'Mitochondria / respiration'
    return d
def dedup(d):
    keep = []
    for i, r in d.iterrows():
        if (not r.fxg) and any(len(r.lset & k) / max(1, len(r.lset | k)) > JACC for _, k in keep if r.lset and k): continue
        keep.append((i, r.lset))
    return d.loc[[i for i, _ in keep]].reset_index(drop=True)

# ---------------------------------------------------------------- une barre-figure pour tout
def bars(d, score, sign, title, xlabel, fn, N=30, nfx=7, nle=12, best_line=False, rank_style=False):
    d = d.sort_values(score, ascending=False).reset_index(drop=True); n_all = len(d); d['rank'] = np.arange(1, n_all + 1)
    if n_all == 0: print('  (vide)', fn); return 0
    top = d.head(N); fx0 = d[(d.fxn | d.fxg) & ~d.index.isin(top.index)]
    fx = pd.concat([fx0[fx0.fxg], fx0[~fx0.fxg].head(nfx)]).sort_values(score, ascending=False).head(nle + nfx)
    rows = pd.concat([top, fx]); n1 = len(top); ypos = list(range(n1)) + [n1 + 0.8 + i for i in range(len(fx))]
    nleg = len(set(rows['proc'])) + 2; fh = 0.27 * (len(rows) + 2) + 1.2 + 0.2 * nleg; fig, ax = plt.subplots(figsize=(11, fh)); labs = []; xmax = rows[score].max()
    for y, (_, r) in zip(ypos, rows.iterrows()):
        ax.barh(y, r[score], color=CC[r.proc], edgecolor='k' if r.fxg else 'none', linewidth=2.2 if r.fxg else 0, height=.82, hatch='///' if r.fxg else None)
        ax.text(r[score] + 0.01 * xmax, y, (' ▲' if r[sign] > 0 else ' ▼') + (f"  rank {int(r['rank'])}/{n_all}" if y >= n1 else ''), va='center', fontsize=8, color='#444')
        labs.append(('★ ' if r.fxn else '') + nice(r.pathway)[:62] + ('  [FXN]' if r.fxg else ''))
    ax.set_yticks(ypos); ax.set_yticklabels(labs, fontsize=9)
    for t, (_, r) in zip(ax.get_yticklabels(), rows.iterrows()):
        if r.fxg: t.set_fontweight('bold'); t.set_color('#B00020')
        elif r.fxn: t.set_fontstyle('italic')
    ax.invert_yaxis()
    if len(fx): ax.axhline(n1 - 0.1 + 0.4, color='#888', ls=':', lw=1)
    if best_line:
        b = d[d.fxn | d.fxg].head(1)
        if len(b):
            v = b[score].iloc[0]; ax.axvline(v, color='k', ls='--', lw=1)
            ax.text(v, -0.75, ' best FXN pathway', fontsize=8, va='bottom', ha='left' if v < ax.get_xlim()[1] * 0.8 else 'right')
    ax.spines[['top', 'right']].set_visible(False); ax.set_xlabel(xlabel, fontsize=12)
    fig.suptitle(textwrap.fill(title, 90), fontsize=12, fontweight='bold', x=0.01, ha='left', y=0.995)
    leg = [Patch(color=v, label=k) for k, v in CC.items() if k in set(rows['proc'])] + [
        Patch(facecolor='white', edgecolor='none', label='★ = FXN-related term by name (Fe-S, iron, respiratory chain, mito translation)'),
        Patch(facecolor='#ddd', edgecolor='k', linewidth=2, hatch='///', label='hatched + [FXN] (red) = FXN gene is in the pathway leading edge')]
    fig.legend(handles=leg, loc='lower left', bbox_to_anchor=(0.25, 0.004), ncol=1, fontsize=8, frameon=False)
    os.makedirs(os.path.dirname(fn), exist_ok=True); nl = len(textwrap.fill(title, 90).split('\n')); plt.tight_layout(rect=[0, (0.2 * nleg + 0.1) / fh, 1, 1 - (0.02 + 0.2 * nl) / fh]); plt.savefig(fn, dpi=170); plt.close(); return n_all

# ---------------------------------------------------------------- chargement de tout
DATA = {m: {c: load(ds, c) for c in cs} for m, (ds, cs) in MODELS.items()}
for m in DATA:
    for c, d in DATA[m].items(): print(f'{m} {c}: {len(d)} pathways, {int(d.sig_all.sum())} significatives, {int(d.robust.sum())} robustes')

# ---------------------------------------------------------------- 1) ranking FDR / NES  + 2) robust
for m in DATA:
    for c, d in DATA[m].items():
        f = figset(d); f = f[f.padj < PADJ].copy(); f['fdr'] = -np.log10(f.padj.clip(lower=1e-300)); f['abs'] = f.NES.abs(); f['sgn'] = np.sign(f.NES)
        for met, lab in (('FDR', '−log10(FDR)'), ('NES', '|NES|')):
            sc = 'fdr' if met == 'FDR' else 'abs'; g = dedup(f.sort_values(sc, ascending=False))
            bars(g, sc, 'sgn', f'{m}: {CL[c]}, ranked by {lab}', lab, f'{OUT}/ranking/ranking_{met}_{m}_{c}.png', N=28, best_line=True)
        r = figset(d); r = r[(r.padj < PADJ) & (r.worst_case_NES_signed != 0)].copy(); r['w'] = r.worst_case_NES_signed.abs(); r['sgn'] = np.sign(r.NES)
        r = dedup(r.sort_values('w', ascending=False))
        bars(r, 'w', 'sgn', f'{m}: {CL[c]}, ranked by worst-case |NES| (leave-one-line-out)', 'worst-case |NES| over leave-one-line-out', f'{OUT}/robust/robust_{m}_{c}.png')

# ---------------------------------------------------------------- 3) corrected / persistent (iNeurons, iPSC)
for m in ('iPSC', 'iNeurons'):
    A, B, C = (figset(DATA[m][k]).set_index('pathway') for k in 'ABC')
    j = A[['NES', 'padj', 'size', 'lset', 'fxg', 'fxn', 'proc']].join(B[['NES', 'padj']], rsuffix='_B').join(C[['NES', 'padj', 'lset', 'fxg']], rsuffix='_C').dropna(subset=['NES', 'padj', 'NES_B', 'padj_B', 'NES_C', 'padj_C'])
    j['lset'] = j.apply(lambda r: r.lset | r.lset_C, axis=1); j['fxg'] = j.fxg | j.fxg_C
    j = j.reset_index()
    cor = j[(j.padj < PADJ) & (j.padj_C < PADJ) & (j.padj_B >= PADJ) & (np.sign(j.NES) == np.sign(j.NES_C))].copy(); cor['s'] = cor.NES_C.abs(); cor['sg'] = np.sign(cor.NES)
    bars(dedup(cor.sort_values('s', ascending=False)), 's', 'sg', f'{m}: pathways corrected by editing (A and C significant, B not), ranked by |NES| FA vs edited', '|NES| of FA vs edited (C)', f'{OUT}/corrected/corrected_{m}.png')
    per = j[(j.padj < PADJ) & (j.padj_B < PADJ) & (np.sign(j.NES) == np.sign(j.NES_B))].copy(); per['s'] = per.NES_B.abs(); per['sg'] = np.sign(per.NES)
    bars(dedup(per.sort_values('s', ascending=False)), 's', 'sg', f'{m}: pathways that persist after editing (A and B significant, same sign), ranked by |NES| edited vs control', '|NES| of edited vs control (B)', f'{OUT}/persistent/persistent_{m}.png')

# ---------------------------------------------------------------- 4) partage entre modèles (contraste A)
COL = {'Fibroblasts': '#B07A3A', 'iPSC': '#6B4FA3', 'iNeurons': '#2E8089'}; ms = list(COL); P = None
for m in ms:
    t = figset(DATA[m]['A'])[['pathway', 'NES', 'padj', 'lset', 'fxg']].set_index('pathway'); t.columns = [f'{c}|{m}' for c in t.columns]
    P = t if P is None else P.join(t, how='outer')
rows = []
for p, r in P.iterrows():
    sig = [m for m in ms if pd.notna(r[f'padj|{m}']) and r[f'padj|{m}'] < PADJ]
    for sg in (1, -1):
        s = [m for m in sig if np.sign(r[f'NES|{m}']) == sg]
        if len(s) >= 2:
            rows.append((p, sg, len(s), min(abs(r[f'NES|{m}']) for m in s), any(bool(r[f'fxg|{m}']) for m in ms if pd.notna(r[f'fxg|{m}'])),
                         set().union(*[r[f'lset|{m}'] for m in s if isinstance(r[f'lset|{m}'], set)])))
d = pd.DataFrame(rows, columns=['pathway', 'sgn', 'n', 'score', 'fxg', 'lset']).sort_values(['n', 'score'], ascending=False).reset_index(drop=True)
d = dedup(d); print('partage: n==3', int((d.n == 3).sum()), 'sur', len(d)); d = d.head(30)
fig, ax = plt.subplots(figsize=(11, 0.42 * len(d) + 2.2)); h = 0.26
for i, r in d.iterrows():
    for k, m in enumerate(ms):
        v = P.loc[r.pathway, f'NES|{m}']; pdj = P.loc[r.pathway, f'padj|{m}']
        if pd.isna(v): ax.text(0.02, i + (k - 1) * h, '–', va='center', fontsize=7, color='#999'); continue
        ax.barh(i + (k - 1) * h, v, height=h * 0.92, color=COL[m], alpha=1 if pdj < PADJ else 0.28)
ax.axvline(0, color='k', lw=.8); ax.set_yticks(range(len(d)))
ax.set_yticklabels([('★ ' if FXN.search(plain(p)) else '') + nice(p)[:60] + ('  [FXN]' if f else '') for p, f in zip(d.pathway, d.fxg)], fontsize=9)
for t, f in zip(ax.get_yticklabels(), d.fxg):
    if f: t.set_fontweight('bold'); t.set_color('#B00020')
ax.invert_yaxis(); ax.spines[['top', 'right']].set_visible(False)
ax.set_xlabel('NES, FA vs control (A)   ← lower in FA | higher in FA →', fontsize=11)
fig.suptitle('Pathways significant in the same direction in at least 2 of 3 models\nranked by number of models, then by the weakest |NES| among them', fontsize=11, fontweight='bold', x=0.01, ha='left', y=0.995)
leg = [Patch(color=COL[m], label=m) for m in ms] + [Patch(color='#999', alpha=.28, label='faded = not significant (padj ≥ 0.05) in that model'),
      Patch(facecolor='white', edgecolor='none', label='★ = FXN-related term by name (Fe-S, iron, respiratory chain, mito translation)'),
      Patch(facecolor='#ddd', edgecolor='none', label='[FXN] (red) = FXN gene in the leading edge of the pathway, in at least one model')]
ax.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.4, -0.07), ncol=1, fontsize=8, frameon=False)
plt.tight_layout(rect=[0, 0, 1, 0.96]); plt.savefig(f'{OUT}/partage_entre_modeles.png', dpi=170); plt.close()

# ---------------------------------------------------------------- 5) tableau xlsx
LEG = [('pathway', 'Nom du gène set (GO BP/CC/MF ou Hallmark)'), ('collection', "Collection d'origine : C5_BP, C5_CC, C5_MF (GO) ou H (Hallmark)"),
 ('size', 'Nombre de gènes du gène set présents dans le classement'), ('FXN_in_leading_edge', 'TRUE si FXN fait partie du leading edge (analyse complète)'),
 ('FXN_related_name', "TRUE si le nom évoque l'axe FXN (Fe-S, fer, chaîne respiratoire, traduction mitochondriale) ; = ★ des figures"),
 ('NES_all / padj_all', 'NES et padj avec toutes les lignées'), ('NES_without_X / padj_without_X', 'Idem après retrait de la lignée X (une colonne par lignée)'),
 ('worst_case_NES_signed', "Plus petit |NES| parmi l'analyse complète et tous les retraits, avec le signe du NES complet ; 0 si le signe s'inverse. Score des figures « robust »"),
 ('n_LOO_padj<0.05', 'Nombre de retraits où padj < 0,05'), ('sig_all', 'padj < 0,05 avec toutes les lignées'),
 ('same_sign_all_LOO', 'Même signe de NES dans tous les retraits'), ('robust', 'sig_all, même signe dans tous les retraits, et p nominale < 0,05 dans chaque retrait'),
 ('contrib_X', 'Perte de |NES| quand on retire la lignée X (NES complet − NES sans X, orienté dans le sens du NES complet)'),
 ('driver_line', 'Lignée dont le retrait fait perdre le plus de |NES|'), ('driver_share', 'Part de cette lignée dans la somme des pertes positives (1 = elle porte tout le signal)'),
 ('single_line_driven', 'TRUE si sig_all et si, en retirant au moins UNE lignée, le NES change de signe ou la p nominale passe à ≥ 0,05')]
wb = Workbook(); ws = wb.active; ws.title = 'LISEZ-MOI'
for t in ['Lisez-moi',
 'Une feuille par modèle et par comparaison (iNeurons / iPSC ; A, B, C). Fibroblastes : feuille Fibroblasts_A (FA vs control uniquement).',
 'A = FA vs control ; B = edited vs control ; C = FA vs edited. NES > 0 = plus haut dans le premier groupe cité.',
 'Toutes les feuilles et toutes les figures viennent des mêmes fichiers (01_gsea_loo.R) : mêmes NES, padj, tailles et leading edges partout.',
 'Les feuilles contiennent tous les pathways testés (GO BP, CC, MF + Hallmark ; 15 à 500 gènes), sans fusion des redondants, triés par padj_all. Les figures ajoutent : size 15-300, padj < 0,05 et fusion des redondants.',
 "Fibroblastes : contraste A uniquement (pas de lignées éditées), donc pas de figures corrected / persistent / partage pour ce modèle sauf la figure de partage (contraste A)."]: ws.append([t])
ws['A1'].font = Font(bold=True, size=14)
wl = wb.create_sheet('Légendes'); wl.append(['colonne', 'signification'])
for x in LEG: wl.append(list(x))
wsum = wb.create_sheet('Résumé'); wsum.append(['feuille', 'modèle', 'comparaison', 'n pathways', 'n sig_all', 'n robust', 'n single_line_driven'])
for m in ('iNeurons', 'iPSC', 'Fibroblasts'):
    for c, d in DATA[m].items():
        out = pd.DataFrame({'pathway': d.pathway, 'collection': d.collection, 'size': d['size'], 'FXN_in_leading_edge': d.FXN_in_leading_edge, 'FXN_related_name': d.FXN_related_name})
        for x in [c_ for c_ in d.columns if c_.startswith('NES|')]: out['NES_' + x.split('|')[1]] = d[x]
        for x in [c_ for c_ in d.columns if c_.startswith('padj|')]: out['padj_' + x.split('|')[1]] = d[x]
        for k in ('worst_case_NES_signed', 'n_LOO_padj<0.05', 'sig_all', 'same_sign_all_LOO', 'robust'): out[k] = d[k]
        for x in [c_ for c_ in d.columns if c_.startswith('contrib|')]: out['contrib_' + x.split('|')[1]] = d[x]
        for k in ('driver_line', 'driver_share', 'single_line_driven'): out[k] = d[k]
        out = out.sort_values('padj_all').reset_index(drop=True); name = f'{m}_{c}'
        w = wb.create_sheet(name)
        for r in dataframe_to_rows(out, index=False, header=True): w.append(r)
        for cell in w[1]: cell.font = Font(bold=True, color='FFFFFF'); cell.fill = PatternFill('solid', fgColor='2B4558')
        w.freeze_panes = 'B2'; w.auto_filter.ref = w.dimensions; w.column_dimensions['A'].width = 70
        wsum.append([name, m, CL[c], len(out), int(out.sig_all.sum()), int(out.robust.sum()), int(out.single_line_driven.sum())])
for w in (wl,):
    for cell in w[1]: cell.font = Font(bold=True, color='FFFFFF'); cell.fill = PatternFill('solid', fgColor='2B4558')
    w.column_dimensions['A'].width = 32; w.column_dimensions['B'].width = 110
    for row in w.iter_rows(min_row=2):
        row[0].font = Font(bold=True)
        for cell in row: cell.alignment = Alignment(wrap_text=True, vertical='top')
ws.column_dimensions['A'].width = 160
for row in ws.iter_rows(min_row=2):
    for cell in row: cell.alignment = Alignment(wrap_text=True, vertical='top')
# ordre des feuilles : LISEZ-MOI, Légendes, Résumé, iNeurons A-C, iPSC A-C, Fibroblasts_A
order = ['LISEZ-MOI', 'Légendes', 'Résumé'] + [f'{m}_{c}' for m in ('iNeurons', 'iPSC', 'Fibroblasts') for c in DATA[m]]
wb._sheets = [wb[n] for n in order]; os.makedirs(OUT, exist_ok=True); wb.save(f'{OUT}/Pathways_leave_one_line_out.xlsx'); print('OK ->', OUT)
