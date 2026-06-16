"""Build region (subsequence) tables for embedding.

A region is a slice of a parent protein with explicit (start, end) coordinates,
1-indexed and inclusive (UniProt convention). The output is a DataFrame ready
to feed into `embed_regions()` — one row per region, with the sliced sequence,
a stable item_id, and provenance columns.

Provenance: each function attaches a dict to `df.attrs['provenance']` describing
what produced the regions (universe, set, operation, args). Downstream caches
inherit and extend this — see `esm_embed.embed_regions(provenance=...)`.
"""

import re

import pandas as pd

from .definitions import resolve_dynamic, NAMED_SETS
from .interpro import load_interpro, db_prefixes


# UniProt feature columns: keyword + colloquial short name + boundary policy.
# Boundary policy: how to fill in unknown coordinates ('?' or missing in the
# UniProt cell). `fill_unknown_start` / `fill_unknown_end`: if not None, use
# that integer; if None, drop the entry.
UNIPROT_FEATURE_KEYWORDS = {
    'Signal peptide':       {'keyword': 'SIGNAL',   'short': ['sp', 'signal'],
                             'fill_unknown_start': 1, 'fill_unknown_end': None},
    'Transmembrane':        {'keyword': 'TRANSMEM', 'short': ['tm', 'transmembrane']},
    'Intramembrane':        {'keyword': 'INTRAMEM'},
    'Topological domain':   {'keyword': 'TOPO_DOM'},
    'Chain':                {'keyword': 'CHAIN'},
    'Propeptide':           {'keyword': 'PROPEP', 'short': ['pp', 'propeptide']},
    'Peptide':              {'keyword': 'PEPTIDE'},
    'Lipidation':           {'keyword': 'LIPID',
                             'short': ['lipid', 'lipidation', 'lipidated']},
    'Disulfide bond':       {'keyword': 'DISULFID', 'short': ['ssbond', 'disulfide']},
    'Glycosylation':        {'keyword': 'CARBOHYD',
                             'short': ['glyc', 'glycan', 'glycosylation', 'glycosylated']},
    'Modified residue':     {'keyword': 'MOD_RES',
                             'short': ['modres', 'modified', 'ptm']},
    'Initiator methionine': {'keyword': 'INIT_MET'},
    'Transit peptide':      {'keyword': 'TRANSIT', 'short': ['transit', 'tp'],
                             'fill_unknown_start': 1, 'fill_unknown_end': None},
    'Motif':                {'keyword': 'MOTIF',   'short': ['motif']},
}


def _shorts(cfg):
    s = cfg.get('short')
    if s is None: return ()
    return (s,) if isinstance(s, str) else tuple(s)


def _resolve_feature(name):
    """Accept full column name, KEYWORD, or any short form. Returns column name."""
    for col, cfg in UNIPROT_FEATURE_KEYWORDS.items():
        if name == col or name == cfg['keyword'] or name in _shorts(cfg):
            return col
    all_shorts = [s for cfg in UNIPROT_FEATURE_KEYWORDS.values() for s in _shorts(cfg)]
    raise KeyError(
        f'unknown feature {name!r}; have '
        f'columns={list(UNIPROT_FEATURE_KEYWORDS)}, '
        f'keywords={[c["keyword"] for c in UNIPROT_FEATURE_KEYWORDS.values()]}, '
        f'shorts={all_shorts}'
    )


def interpro_domain_regions(query, db=None, entry_type=None, dedup=True):
    """Return InterPro domain/site regions for every protein in `query`.

    Parameters
    ----------
    query : str | list | DataFrame
        Anything `resolve()` handles: a registered named-set name, a single
        UniProt ID, a comma-separated string, a list of IDs, a file path, or
        a DataFrame with an `Entry` column.
    db : str or None
        Filter to a member db: 'pfam', 'smart', 'cdd', 'gene3d', 'superfamily', etc.
    entry_type : str or None
        Filter to entry_type: 'Domain', 'Family', 'Repeat', 'Active_site', etc.
    dedup : bool

    Returns a DataFrame; provenance is attached as `df.attrs['provenance']`.
    """
    from .lookup import resolve_with_seqs
    if hasattr(query, 'columns') and 'Entry' in query.columns:
        query = query['Entry'].tolist()
    ids, seqs = resolve_with_seqs(query)
    seq_lookup = seqs
    id_set = set(ids)

    ipr = load_interpro(with_entry_types=(entry_type is not None))
    ipr = ipr[ipr.uniprot.isin(id_set)]

    if db is not None:
        prefixes = db_prefixes(db)
        ipr = ipr[ipr.member_db_acc.str.startswith(tuple(prefixes), na=False)]
    if entry_type is not None:
        ipr = ipr[ipr.entry_type == entry_type]
    if dedup:
        ipr = ipr.drop_duplicates(['uniprot', 'member_db_acc', 'start', 'end'])

    rows = []
    skipped = 0
    for _, r in ipr.iterrows():
        uid = r['uniprot']
        seq = seq_lookup.get(uid)
        if not seq:
            skipped += 1
            continue
        start, end = int(r['start']), int(r['end'])
        if start < 1 or end > len(seq) or end < start:
            skipped += 1
            continue
        rows.append({
            'item_id':     f'{uid}_{r["member_db_acc"]}_{start}_{end}',
            'uniprot':     uid,
            'source':      r.get('member_db', r['member_db_acc']),
            'source_id':   r['member_db_acc'],
            'interpro_id': r['interpro_id'],
            'entry_type':  r.get('entry_type', ''),
            'start':       start,
            'end':         end,
            'length':      end - start + 1,
            'sequence':    seq[start - 1:end],
        })
    out = pd.DataFrame(rows)
    if skipped:
        print(f'  domain_regions: dropped {skipped} regions (missing seq or out-of-bounds)')

    from .layout import attach_provenance, SUBGROUP_FROM_INTERPRO_DB
    subgroup = (SUBGROUP_FROM_INTERPRO_DB.get(db) if db else 'all_domains')
    attach_provenance(
        out,
        set=query if isinstance(query, str) else 'custom',
        operation='interpro_domain_regions',
        args={'db': db, 'entry_type': entry_type},
        source_category='interpro',
        subgroup=subgroup,
        kind=None,
    )
    return out


# Backcompat alias: old name kept so existing scripts keep working.
domain_regions = interpro_domain_regions


def pfam_regions(named_set, **kw):
    """Convenience: interpro_domain_regions filtered to Pfam."""
    return interpro_domain_regions(named_set, db='pfam', **kw)


def _parse_uniprot_feature_cell(cell, keyword):
    """Yield (start, end, qualifiers) for each `keyword` occurrence in cell."""
    if not isinstance(cell, str) or not cell:
        return
    pat = re.compile(
        rf'\b{keyword}\s+(\??\d+|\?)(?:\.\.(\??\d+|\?))?(.*?)(?=;\s*{keyword}\s|$)',
        re.DOTALL,
    )

    def _coord(token):
        if token is None:
            return None
        token = token.lstrip('?')
        return int(token) if token else None

    for m in pat.finditer(cell):
        start = _coord(m.group(1))
        end_tok = m.group(2)
        end = _coord(end_tok) if end_tok is not None else start
        yield start, end, m.group(3).strip()


def uniprot_regions(query, feature='Signal peptide', classifier=None, kind=None):
    """Return UniProt feature regions for every protein in `query`.

    `query` accepts anything `resolve()` handles (named-set name, single ID,
    list, comma-separated string, file path, or DataFrame with `Entry`).

    `feature` accepts the column name ('Signal peptide'), the keyword
    ('SIGNAL'), or the short form ('sp').

    Boundary policy from `UNIPROT_FEATURE_KEYWORDS[feature]`:
      `fill_unknown_start` / `fill_unknown_end`: if not None, fill in;
      if None, drop entries with that bound unknown.

    `classifier` (optional): callable note_str -> class_str. When provided, the
    output gets `class` and `note` columns. Used for features with multiple
    biologically-distinct subtypes (transit peptide kinds, motif kinds, etc.).
    `kind` (optional): if set, keep only rows whose classifier output equals it.
    """
    from .lookup import resolve as _resolve_basic
    from .uniprot import load_proteome
    feature = _resolve_feature(feature)
    cfg = UNIPROT_FEATURE_KEYWORDS[feature]
    keyword = cfg['keyword']
    fill_start = cfg.get('fill_unknown_start')
    fill_end   = cfg.get('fill_unknown_end')

    if hasattr(query, 'columns') and 'Entry' in query.columns:
        query = query['Entry'].tolist()
    if isinstance(query, str) and query in NAMED_SETS:
        _ids, _seqs, df = resolve_dynamic(query)
    else:
        ids = _resolve_basic(query)
        df = load_proteome()
        df = df[df['Entry'].isin(set(ids))]
    if feature not in df.columns:
        raise KeyError(f'feature column {feature!r} not in proteome DataFrame')

    rows = []
    skipped_unknown = 0
    skipped_bounds  = 0
    for _, r in df.iterrows():
        uid = r['Entry']
        seq = r['Sequence']
        cell = r[feature]
        for start, end, quals in _parse_uniprot_feature_cell(cell, keyword):
            if start is None:
                if fill_start is None:
                    skipped_unknown += 1
                    continue
                start = fill_start
            if end is None:
                if fill_end is None:
                    skipped_unknown += 1
                    continue
                end = fill_end
            if start < 1 or end > len(seq) or end < start:
                skipped_bounds += 1
                continue
            row = {
                'item_id':     f'{uid}_{keyword}_{start}_{end}',
                'uniprot':     uid,
                'source':      'uniprot',
                'source_id':   keyword,
                'interpro_id': '',
                'entry_type':  '',
                'start':       start,
                'end':         end,
                'length':      end - start + 1,
                'sequence':    seq[start - 1:end],
            }
            if classifier is not None:
                note_match = re.search(r'/note="([^"]+)"', quals)
                note = note_match.group(1) if note_match else ''
                # classifiers may optionally accept the sequence as a 2nd arg
                # (used e.g. by motif classifier to split KDEL-canonical vs
                # dibasic/KKxx based on the literal motif residues, since
                # UniProt's /note= text is identical for both).
                try:
                    cls = classifier(note, row['sequence'])
                except TypeError:
                    cls = classifier(note)
                if kind is not None and cls != kind:
                    continue
                row['class'] = cls
                row['note']  = note
            rows.append(row)
    out = pd.DataFrame(rows)
    if skipped_unknown:
        print(f'  uniprot_feature_regions: dropped {skipped_unknown} (unknown bound)')
    if skipped_bounds:
        print(f'  uniprot_feature_regions: dropped {skipped_bounds} (out-of-bounds)')

    from .layout import attach_provenance, SUBGROUP_FROM_UNIPROT_FEATURE
    attach_provenance(
        out,
        set=query if isinstance(query, str) else 'custom',
        operation='uniprot_regions',
        args={'feature': feature, 'kind': kind},
        source_category='uniprot_features',
        subgroup=SUBGROUP_FROM_UNIPROT_FEATURE.get(feature, feature.lower().replace(' ', '_')),
        kind=kind,
    )
    return out


# Backcompat alias.
uniprot_feature_regions = uniprot_regions


def uniprot_regions_remove(named_set, feature='Signal peptide'):
    """Return the COMPLEMENT of a UniProt feature, one region per parent.

    For each protein, the residues NOT in any annotated occurrence of the
    feature MUST form a single contiguous span. For N-terminal cleavable
    features (SP, PROPEP) this is `(end+1, length)` — unambiguous. For
    interior features (TM, etc.) the complement splits into multiple spans
    and the operation is ambiguous: in that case this function raises with
    examples, and you must define an explicit per-feature `remove_policy`
    in `UNIPROT_FEATURE_KEYWORDS` (or use a different operation).

    item_id is set to the parent UniProt accession (one row per parent).
    Use this for populating `human_<short>_removed/` universes.
    """
    feature = _resolve_feature(feature)
    if named_set not in NAMED_SETS:
        raise KeyError(f'unknown named set {named_set!r}; have {sorted(NAMED_SETS)}')

    cfg = UNIPROT_FEATURE_KEYWORDS[feature]
    keyword = cfg['keyword']
    fill_start = cfg.get('fill_unknown_start')
    fill_end   = cfg.get('fill_unknown_end')
    short      = (_shorts(cfg) or (keyword.lower(),))[0]

    _ids, _seqs, df = resolve_dynamic(named_set)
    if feature not in df.columns:
        raise KeyError(f'feature column {feature!r} not in proteome DataFrame')

    rows = []
    ambiguous = []   # proteins where complement isn't a single contiguous span
    for _, r in df.iterrows():
        uid = r['Entry']
        seq = r['Sequence']
        L = len(seq)
        feature_spans = []
        for start, end, _ in _parse_uniprot_feature_cell(r[feature], keyword):
            if start is None:
                start = fill_start
            if end is None:
                end = fill_end
            if start is None or end is None or start < 1 or end > L or end < start:
                continue
            feature_spans.append((start, end))

        if not feature_spans:
            continue

        # find all contiguous uncovered spans
        covered = [False] * (L + 1)
        for s, e in feature_spans:
            for i in range(s, e + 1):
                covered[i] = True
        spans = []
        i = 1
        while i <= L:
            if not covered[i]:
                j = i
                while j <= L and not covered[j]:
                    j += 1
                spans.append((i, j - 1))
                i = j
            else:
                i += 1

        if len(spans) == 0:
            continue
        if len(spans) > 1:
            ambiguous.append((uid, feature_spans, spans))
            continue

        s, e = spans[0]
        rows.append({
            'item_id':     uid,
            'uniprot':     uid,
            'source':      'uniprot',
            'source_id':   f'{short}_removed',
            'interpro_id': '',
            'entry_type':  '',
            'start':       s,
            'end':         e,
            'length':      e - s + 1,
            'sequence':    seq[s - 1:e],
        })

    if ambiguous:
        examples = '\n'.join(
            f'  {uid}: feature spans={fs}, complement spans={cs}'
            for uid, fs, cs in ambiguous[:5]
        )
        more = f'\n  ... and {len(ambiguous)-5} more' if len(ambiguous) > 5 else ''
        raise ValueError(
            f"uniprot_feature_regions_remove: feature {feature!r} produced "
            f"non-contiguous complement for {len(ambiguous)} proteins (operation is ambiguous "
            f"for interior features). Examples:\n{examples}{more}\n"
            f"Either define an explicit `remove_policy` for {feature!r} in "
            f"UNIPROT_FEATURE_KEYWORDS, or use a different operation."
        )

    out = pd.DataFrame(rows)

    from .layout import attach_provenance, SUBGROUP_FROM_UNIPROT_FEATURE
    sub = SUBGROUP_FROM_UNIPROT_FEATURE.get(feature, feature.lower().replace(' ', '_'))
    attach_provenance(
        out,
        set=named_set,
        operation='uniprot_regions_remove',
        args={'feature': feature},
        source_category='uniprot_features',
        subgroup=f'{sub}_removed',
        kind=None,
    )
    return out


# Backcompat alias.
uniprot_feature_regions_remove = uniprot_regions_remove


def signal_peptide_regions(named_set):
    """Convenience: uniprot_regions for the Signal peptide column."""
    return uniprot_regions(named_set, feature='Signal peptide')


def propeptide_regions(named_set):
    """Convenience: uniprot_regions for the Propeptide column."""
    return uniprot_regions(named_set, feature='Propeptide')


def _classify_transit(note):
    n = (note or '').lower()
    if 'mitochondrion' in n: return 'mts'
    if 'peroxisome' in n:    return 'pts2'
    if 'chloroplast' in n:   return 'chloroplast'
    return 'other'


def transit_peptide_regions(query, kind=None):
    """N-terminal transit peptide regions from the Transit peptide column.

    kind: None (all) | 'mts' | 'pts2' | 'chloroplast' | 'other'

    Human proteome currently carries 563 MTS (mitochondrial) + 3 PTS2 (peroxisomal)
    entries (no chloroplast). N-terminal cleavable signals analogous to SP.
    Note: MTS and PTS2 target different organelles via different machinery —
    don't combine their embeddings into a single mean. Use per-class subsets.
    """
    return uniprot_regions(query, feature='Transit peptide',
                           classifier=_classify_transit, kind=kind)


def mts_regions(query):  return transit_peptide_regions(query, kind='mts')
def pts2_regions(query): return transit_peptide_regions(query, kind='pts2')
mito_transit_regions = mts_regions  # backward-compat alias


def _classify_motif(note, sequence=None):
    """Coarse class for UniProt /note= text on MOTIF entries.

    Most classes are decided from /note= alone. ER-retention is special:
    UniProt's note text ('Prevents secretion from ER', 'ER retention motif',
    etc.) is identical for two mechanistically-distinct classes, so we split
    by the actual residues:

      - er_retention_kdel: short C-terminal soluble retention motifs (KDEL family,
        recognized by KDEL receptor). Pattern: 3-6 residues ending in L or F,
        no consecutive basic residues. E.g. KDEL, HDEL, REEL, RDEF.
      - er_retention_dibasic: KKxx / dibasic membrane-anchored retention motifs
        (recognized by COPI for retrograde retrieval). Pattern: contains
        consecutive K/R. E.g. RR, RKR, KKCS, KRKKK.

    Everything else with the ER-retention note falls into 'er_retention_other'.
    """
    n = (note or '').lower()
    if 'nuclear localization signal' in n or 'nuclear localization site' in n:
        return 'nls'
    if 'nuclear export signal' in n or 'nuclear export sequence' in n:
        return 'nes'
    if ('microbody targeting signal' in n
        or ('peroxisom' in n and ('targeting signal' in n
                                  or 'localization signal' in n))
        or 'pts1' in n or 'pts-1' in n):
        return 'pts1'
    if 'cell attachment' in n:
        seq = (sequence or '').upper()
        if seq == 'RGD' and 'atypical' not in n and 'transferrin' not in n:
            return 'cell_attachment_rgd'
        return 'cell_attachment_atypical'
    if 'pdz-binding' in n or 'pdz binding' in n:
        return 'pdz_binding'
    if 'sh3-binding' in n or 'sh3 binding' in n:
        return 'sh3_binding'
    if ('prevents secretion from er' in n
        or 'er retention' in n
        or 'endoplasmic reticulum reten' in n
        or 'retention in the endoplasmic reticulum' in n
        or 'kdel' in n):
        seq = (sequence or '').upper()
        if (3 <= len(seq) <= 6
            and seq.endswith(('L', 'F'))
            and not re.search(r'[KR]{2}', seq)):
            return 'er_retention_kdel'
        if re.search(r'[KR]{2}', seq):
            return 'er_retention_dibasic'
        return 'er_retention_other'
    return 'other'


def motif_regions(query, kind=None):
    """Short functional motifs from the UniProt Motif column.

    kind: None (all) | 'nls' | 'nes' | 'er_retention_kdel' | 'er_retention_dibasic'
          | 'er_retention_other' | 'other'

    Definitions live in `_classify_motif`. NLS / NES / ER-retention are
    surfaced as named classes today; PDZ-binding, SH3-binding, PTS1,
    cell-attachment, selectivity-filter, etc. are ready to be promoted out of
    'other' when needed.
    """
    return uniprot_regions(query, feature='Motif',
                           classifier=_classify_motif, kind=kind)


def nls_regions(query):                  return motif_regions(query, kind='nls')
def nes_regions(query):                  return motif_regions(query, kind='nes')
def pts1_regions(query):                 return motif_regions(query, kind='pts1')
def er_retention_kdel_regions(query):    return motif_regions(query, kind='er_retention_kdel')
def er_retention_dibasic_regions(query): return motif_regions(query, kind='er_retention_dibasic')

def cell_attachment_rgd_regions(query):      return motif_regions(query, kind='cell_attachment_rgd')
def cell_attachment_atypical_regions(query): return motif_regions(query, kind='cell_attachment_atypical')
rgd_regions = cell_attachment_rgd_regions   # household-name alias

def pdz_binding_regions(query): return motif_regions(query, kind='pdz_binding')
def sh3_binding_regions(query): return motif_regions(query, kind='sh3_binding')


def cell_attachment_regions(query, canonical=True):
    """Cell-attachment motifs (integrin-binding RGD + variants).

    canonical=True  (default): literal RGD sequence with standard "Cell
                                attachment site" note. ~128 integrin-binding
                                ECM motifs.
    canonical=False:           includes atypical (RGE / DGR / RGS) and the
                                TFRC-transferrin case.
    """
    df = motif_regions(query)
    if canonical:
        mask = df['class'] == 'cell_attachment_rgd'
    else:
        mask = df['class'].str.startswith('cell_attachment')
    return df[mask].reset_index(drop=True)


# ---- canonicality wrappers (compositions across columns/classes) ----

def er_retention_regions(query, canonical=True):
    """ER-retention motifs across both mechanisms.

    canonical=True  (default): KDEL family only — soluble, recognized by KDEL
                                receptor. The classical ER retention motif.
    canonical=False:           also includes COPI-mediated dibasic/KKxx motifs
                                and atypical residual entries (mechanistically
                                distinct from KDEL).
    """
    df = motif_regions(query)
    if canonical:
        mask = df['class'] == 'er_retention_kdel'
    else:
        mask = df['class'].str.startswith('er_retention')
    return df[mask].reset_index(drop=True)


def mitochondrial_targeting_regions(query):
    """Mitochondrial targeting signal (MTS). Alias for mts_regions.
    No canonicality flag — MTS is the only mechanism for matrix import
    (TIM/TOM machinery). Outer-membrane β-barrels and IMS proteins use
    different routes but are not annotated as Transit peptide.
    """
    return mts_regions(query)


def nuclear_localization_regions(query):
    """Nuclear localization signal (NLS). Alias for nls_regions."""
    return nls_regions(query)


# British-spelling alias
nuclear_localisation_regions = nuclear_localization_regions


def nuclear_export_regions(query):
    """Nuclear export signal (NES). Alias for nes_regions."""
    return nes_regions(query)


def peroxisome_targeting_regions(query, canonical=True):
    """Peroxisomal targeting signals across both mechanisms.

    canonical=True  (default): PTS1 only — C-terminal tripeptide, not cleaved,
                                recognized by PEX5. Dominant mechanism.
    canonical=False:           also includes PTS2 (N-terminal cleavable
                                nonapeptide, recognized by PEX7). PTS2 lives in
                                the Transit peptide column, distinct from PTS1
                                in Motif.
    """
    pts1 = pts1_regions(query)
    if canonical:
        return pts1
    pts2 = pts2_regions(query)
    return pd.concat([pts1, pts2], ignore_index=True)


# ---- site-level features (1-residue annotations) ----

def _classify_lipid(note):
    n = (note or '').lower()
    if 'gpi-anchor' in n: return 'gpi'
    if 'palmitoyl' in n: return 'palmitoyl'
    if 'myristoyl' in n: return 'myristoyl'
    if 'geranylgeranyl' in n or 'farnesyl' in n: return 'prenyl'
    return 'other'


def _classify_glycan(note):
    n = (note or '').lower()
    if n.startswith('n-linked'): return 'n_linked'
    if n.startswith('o-linked'): return 'o_linked'
    if n.startswith('c-linked'): return 'c_linked'
    if n.startswith('s-linked'): return 's_linked'
    return 'other'


def _classify_modres(note):
    """Classify a Modified residue /note into a coarse PTM class.

    Phospho-: Phosphoserine/threonine/tyrosine, etc.
    Acetyl: N-acetyl- (N-terminal) and N6-acetyllysine.
    Methyl: N-methyl-, N6-methyl-, asymmetric/symmetric dimethylarginine, etc.
    """
    n = (note or '').lower()
    if n.startswith('phospho'): return 'phospho'
    if 'acetyl' in n and 'hydroxyisobutyryl' not in n: return 'acetyl'
    if 'methyl' in n: return 'methyl'
    return 'other'


def _site_table(query, feature, classifier, kind=None):
    """Generic 1-residue site extractor with /note classification.

    Reads UniProt feature cells (semicolon-separated keyword entries), keeps
    only those with a known position, classifies each occurrence by /note,
    and optionally filters to a single class.
    """
    from .lookup import resolve as _resolve_basic
    from .uniprot import load_proteome
    feature = _resolve_feature(feature)
    cfg = UNIPROT_FEATURE_KEYWORDS[feature]
    keyword = cfg['keyword']

    if hasattr(query, 'columns') and 'Entry' in query.columns:
        query = query['Entry'].tolist()
    if isinstance(query, str) and query in NAMED_SETS:
        _ids, _seqs, df = resolve_dynamic(query)
    else:
        ids = _resolve_basic(query)
        df = load_proteome()
        df = df[df['Entry'].isin(set(ids))]

    rows = []
    skipped_unknown = skipped_isoform = 0
    pat_iso = re.compile(rf'\b{keyword}\s+[A-Z0-9]+-\d+:')
    for _, r in df.iterrows():
        cell = r[feature]
        if not isinstance(cell, str) or not cell:
            continue
        # Drop isoform-specific entries (e.g. "LIPID P14138-2:741") — these
        # reference non-canonical isoform sequences not in load_proteome.
        if pat_iso.search(cell):
            skipped_isoform += sum(1 for _ in pat_iso.finditer(cell))
        uid = r['Entry']
        seq = r['Sequence']
        for start, end, quals in _parse_uniprot_feature_cell(cell, keyword):
            if start is None or end is None:
                skipped_unknown += 1
                continue
            if start < 1 or end > len(seq) or end < start:
                continue
            note_match = re.search(r'/note="([^"]+)"', quals)
            note = note_match.group(1) if note_match else ''
            cls = classifier(note)
            if kind is not None and cls != kind:
                continue
            rows.append({
                'item_id':     f'{uid}_{keyword}_{start}_{end}',
                'uniprot':     uid,
                'source':      'uniprot',
                'source_id':   keyword,
                'interpro_id': '',
                'entry_type':  '',
                'start':       start,
                'end':         end,
                'length':      end - start + 1,
                'sequence':    seq[start - 1:end],
                'class':       cls,
                'note':        note,
            })
    out = pd.DataFrame(rows)
    if skipped_isoform:
        print(f'  {feature.lower()}_sites: skipped {skipped_isoform} isoform-only entries')
    if skipped_unknown:
        print(f'  {feature.lower()}_sites: dropped {skipped_unknown} (unknown bound)')
    from .layout import attach_provenance, SUBGROUP_FROM_UNIPROT_FEATURE
    attach_provenance(
        out,
        set=query if isinstance(query, str) else 'custom',
        operation=f'{feature.lower().replace(" ", "_")}_sites',
        args={'kind': kind},
        source_category='uniprot_features',
        subgroup=SUBGROUP_FROM_UNIPROT_FEATURE.get(feature, feature.lower().replace(' ', '_')),
        kind=kind,
    )
    return out


def lipidation_sites(query, kind=None):
    """1-residue lipid attachment sites from the Lipidation column.

    kind: None (all) | 'gpi' | 'palmitoyl' | 'myristoyl' | 'prenyl' | 'other'
    """
    return _site_table(query, 'Lipidation', _classify_lipid, kind=kind)


def gpi_sites(query):       return lipidation_sites(query, kind='gpi')
def palmitoyl_sites(query): return lipidation_sites(query, kind='palmitoyl')
def myristoyl_sites(query): return lipidation_sites(query, kind='myristoyl')
def prenyl_sites(query):    return lipidation_sites(query, kind='prenyl')


def glycosylation_sites(query, kind=None):
    """1-residue glycosylation sites from the Glycosylation column.

    kind: None (all) | 'n_linked' | 'o_linked' | 'c_linked' | 's_linked' | 'other'
    """
    return _site_table(query, 'Glycosylation', _classify_glycan, kind=kind)


def n_glycan_sites(query): return glycosylation_sites(query, kind='n_linked')
def o_glycan_sites(query): return glycosylation_sites(query, kind='o_linked')
def c_glycan_sites(query): return glycosylation_sites(query, kind='c_linked')


def modified_residue_sites(query, kind=None):
    """1-residue PTM sites from the Modified residue column.

    Three coarse classes are surfaced as `kind`:
      'phospho'  — Phosphoserine/threonine/tyrosine                 (~41k sites)
      'acetyl'   — N-terminal acetylation + N6-acetyllysine          (~6k)
      'methyl'   — N-/N6-/dimethyl-Lys/Arg                            (~2k)

    Everything else (succinyl-, lactoyl-, crotonyl-, glutaryl-, hydroxyproline,
    sulfo-, ADP-ribosyl-, etc.) currently falls into 'other'. Break out as
    additional classes if/when needed — see _classify_modres in regions.py.
    """
    return _site_table(query, 'Modified residue', _classify_modres, kind=kind)


def phospho_sites(query): return modified_residue_sites(query, kind='phospho')
def acetyl_sites(query):  return modified_residue_sites(query, kind='acetyl')
def methyl_sites(query):  return modified_residue_sites(query, kind='methyl')


def disulfide_sites(query):
    """1-residue Cys sites that participate in a disulfide bond.

    UniProt annotates each bond as `DISULFID a..b` (the two Cys partners). This
    function emits ONE ROW PER PARTICIPATING Cys (so each bond contributes two
    rows). Pair structure is biology, not preserved at the row level — recover
    it from UniProt's annotation if needed.
    """
    from .lookup import resolve as _resolve_basic
    from .uniprot import load_proteome

    if hasattr(query, 'columns') and 'Entry' in query.columns:
        query = query['Entry'].tolist()
    if isinstance(query, str) and query in NAMED_SETS:
        _ids, _seqs, df = resolve_dynamic(query)
    else:
        ids = _resolve_basic(query)
        df = load_proteome()
        df = df[df['Entry'].isin(set(ids))]

    rows = []
    skipped_isoform = skipped_unknown = 0
    pat_iso = re.compile(r'\bDISULFID\s+[A-Z0-9]+-\d+:')
    for _, r in df.iterrows():
        cell = r['Disulfide bond']
        if not isinstance(cell, str) or not cell:
            continue
        if pat_iso.search(cell):
            skipped_isoform += sum(1 for _ in pat_iso.finditer(cell))
        uid = r['Entry']
        seq = r['Sequence']
        for start, end, _q in _parse_uniprot_feature_cell(cell, 'DISULFID'):
            if start is None or end is None:
                skipped_unknown += 1
                continue
            for pos in (start, end):
                if pos < 1 or pos > len(seq):
                    continue
                rows.append({
                    'item_id':     f'{uid}_DISULFID_{pos}_{pos}',
                    'uniprot':     uid,
                    'source':      'uniprot',
                    'source_id':   'DISULFID',
                    'interpro_id': '',
                    'entry_type':  '',
                    'start':       pos,
                    'end':         pos,
                    'length':      1,
                    'sequence':    seq[pos - 1],
                })
    out = pd.DataFrame(rows).drop_duplicates(subset=['item_id']).reset_index(drop=True)
    if skipped_isoform:
        print(f'  disulfide_sites: skipped {skipped_isoform} isoform-only entries')
    if skipped_unknown:
        print(f'  disulfide_sites: dropped {skipped_unknown} (unknown bound)')
    from .layout import attach_provenance
    attach_provenance(
        out,
        set=query if isinstance(query, str) else 'custom',
        operation='disulfide_sites',
        source_category='uniprot_features',
        subgroup='disulfide',
        kind=None,
    )
    return out


def disulfide_bonds(query):
    """One row per disulfide bond, with both Cys positions and loop length.

    Bond-level companion to `disulfide_sites` (which is site-level for embedding).
    Use this for counts, loop-length distributions, or any question of the form
    'which Cys is paired with which'.
    """
    from .lookup import resolve as _resolve_basic
    from .uniprot import load_proteome

    if hasattr(query, 'columns') and 'Entry' in query.columns:
        query = query['Entry'].tolist()
    if isinstance(query, str) and query in NAMED_SETS:
        _ids, _seqs, df = resolve_dynamic(query)
    else:
        ids = _resolve_basic(query)
        df = load_proteome()
        df = df[df['Entry'].isin(set(ids))]

    rows = []
    skipped_isoform = skipped_unknown = 0
    pat_iso = re.compile(r'\bDISULFID\s+[A-Z0-9]+-\d+:')
    for _, r in df.iterrows():
        cell = r['Disulfide bond']
        if not isinstance(cell, str) or not cell:
            continue
        if pat_iso.search(cell):
            skipped_isoform += sum(1 for _ in pat_iso.finditer(cell))
        uid = r['Entry']
        L = len(r['Sequence'])
        for start, end, _q in _parse_uniprot_feature_cell(cell, 'DISULFID'):
            if start is None or end is None:
                skipped_unknown += 1
                continue
            cys1, cys2 = sorted((start, end))
            if cys1 < 1 or cys2 > L or cys1 == cys2:
                continue
            rows.append({
                'bond_id':     f'{uid}_{cys1}_{cys2}',
                'uniprot':     uid,
                'cys1':        cys1,
                'cys2':        cys2,
                'loop_length': cys2 - cys1 + 1,
                'kind':        'intra',
            })
    out = pd.DataFrame(rows).drop_duplicates(subset=['bond_id']).reset_index(drop=True)
    if skipped_isoform:
        print(f'  disulfide_bonds: skipped {skipped_isoform} isoform-only entries')
    if skipped_unknown:
        print(f'  disulfide_bonds: dropped {skipped_unknown} (unknown bound)')
    # Note: disulfide_bonds is analytical (one row per bond, two positions
    # per row) — not a feature for embedding. We still attach provenance for
    # consistency, but this universe path is not normally materialized.
    from .layout import attach_provenance
    attach_provenance(
        out,
        set=query if isinstance(query, str) else 'custom',
        operation='disulfide_bonds',
        source_category='uniprot_features',
        subgroup='disulfide_bonds',
        kind=None,
    )
    return out


def from_external_table(positions, *, study, window, anchor='center',
                        meta_cols=None, sequence_source='reviewed_human',
                        kind=None, drop_oob=True):
    """Build region rows from an external (uniprot, position, **meta) table.

    For each row in `positions`, windows the parent sequence around `position`
    and emits a standard region row. Designed for site-style external studies
    (cleavage sites, phospho sites, etc.) — see `protein_sets.studies`.

    Parameters
    ----------
    positions : DataFrame
        Must have columns 'uniprot' and 'position' (1-indexed, parent coords).
        Any extra columns listed in `meta_cols` are carried through onto the
        output and included in the item_id for uniqueness.
    study : str
        Study name. item_id is namespaced as f"{study}__...".
    window : int
        Window length (in residues). For ``anchor='center'``, must be odd.
    anchor : {'center', 'n_term', 'c_term'}
        How `position` maps to the window:
          'center' → start = position - (window-1)//2, end = start + window-1
          'n_term' → start = position, end = position + window-1
          'c_term' → end = position, start = position - window + 1
    meta_cols : list[str] or None
        Columns from `positions` to keep on the output (and use in item_id).
        Defaults to all non-(uniprot,position) columns.
    sequence_source : str
        Resolver query for the parent sequences. Default 'reviewed_human'.
    kind : str or None
        Provenance ``kind`` field; if None, inferred as f"{anchor}_w{window}".
    drop_oob : bool
        Drop rows whose window falls outside the parent sequence. If False,
        rows are kept but the sequence is truncated to the in-bounds slice.

    Returns
    -------
    DataFrame with columns:
        item_id, uniprot, sequence, start, end, length, study, position,
        + meta_cols.
    """
    from .lookup import resolve_with_seqs

    if anchor == 'center' and window % 2 == 0:
        raise ValueError(f'center anchor requires odd window, got {window}')
    if anchor not in ('center', 'n_term', 'c_term'):
        raise ValueError(f"anchor must be 'center', 'n_term', or 'c_term'; got {anchor!r}")

    if 'uniprot' not in positions.columns or 'position' not in positions.columns:
        raise ValueError("positions must have 'uniprot' and 'position' columns")

    if meta_cols is None:
        meta_cols = [c for c in positions.columns if c not in ('uniprot', 'position')]

    _, seqs = resolve_with_seqs(positions['uniprot'].drop_duplicates().tolist())

    rows = []
    skipped = 0
    for _, r in positions.iterrows():
        uid = r['uniprot']
        seq = seqs.get(uid)
        if not seq:
            skipped += 1
            continue
        pos = int(r['position'])
        if anchor == 'center':
            half = (window - 1) // 2
            start, end = pos - half, pos + half
        elif anchor == 'n_term':
            start, end = pos, pos + window - 1
        else:  # c_term
            start, end = pos - window + 1, pos

        if drop_oob and (start < 1 or end > len(seq)):
            skipped += 1
            continue
        s_clip, e_clip = max(start, 1), min(end, len(seq))
        if e_clip < s_clip:
            skipped += 1
            continue

        meta = {c: r[c] for c in meta_cols}
        meta_slug = '__'.join(str(meta[c]).replace('|', '-').replace(' ', '_')
                              for c in meta_cols) if meta_cols else ''
        item_id = f'{study}__{uid}__{pos}'
        if meta_slug:
            item_id = f'{item_id}__{meta_slug}'

        rows.append({
            'item_id':  item_id,
            'uniprot':  uid,
            'sequence': seq[s_clip - 1:e_clip],
            'start':    start,
            'end':      end,
            'length':   e_clip - s_clip + 1,
            'study':    study,
            'position': pos,
            **meta,
        })

    out = pd.DataFrame(rows)
    if skipped:
        print(f'  from_external_table: dropped {skipped} rows (missing seq or out-of-bounds)')

    from .layout import attach_provenance
    attach_provenance(
        out,
        set=study,
        operation='from_external_table',
        args={'window': window, 'anchor': anchor, 'sequence_source': sequence_source},
        source_category='studies',
        subgroup=study,
        kind=kind or f'{anchor}_w{window}',
    )
    return out
