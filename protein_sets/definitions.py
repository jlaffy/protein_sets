"""Named filter functions over a cached UniProt proteome.

A definition takes the proteome DataFrame and returns a filtered DataFrame.
The strict variant additionally requires `ECO:0000269` (experimental existence)
in the relevant feature column.
"""

from .uniprot import load_proteome, load_alt_isoforms
from .hpa import hpa_secretome_uniprots
from .cspa import cspa_uniprots
from .humantfs import humantfs_uniprots
from . import oma as _oma
from ._config import predictions_path

STRICT_CODE = 'ECO:0000269'


def _has_value(df, col):
    return df[col].notna() & (df[col] != '')


def reviewed_human(refresh=False):
    """All reviewed human proteins."""
    return load_proteome('UP000005640', reviewed=True, refresh=refresh)


def reviewed_human_alt_isoforms(refresh=False):
    """All reviewed human *alternative* (non-canonical) isoforms.

    One row per non-canonical isoform with full sequence; `Entry` carries the
    `-N` suffix. The canonical isoform of each parent lives in `reviewed_human`
    under the bare accession. Use `all_isoforms_of(accession)` to combine them.
    """
    return load_alt_isoforms('UP000005640', refresh=refresh)


def all_isoforms_of(accession, refresh=False):
    """Return canonical row + alt-isoform rows for a parent accession.

    Pulls the canonical (bare-accession) row from `reviewed_human` and the
    `-N` rows from `reviewed_human_alt_isoforms`. Returns one DataFrame.
    """
    import pandas as pd
    canon = reviewed_human(refresh=refresh)
    alts = reviewed_human_alt_isoforms(refresh=refresh)
    canon_rows = canon[canon['Entry'] == accession]
    alt_rows = alts[alts['Entry'].str.startswith(accession + '-')]
    return pd.concat([canon_rows, alt_rows], ignore_index=True, sort=False).fillna('')


# --- Alt-isoform "no TM" subset ---
#
# `alt_isoforms_no_tm` is defined by per-isoform tmbed predictions: alts whose
# tmbed-predicted topology lacks any TM segment (alpha or beta). Reads
# `predictions/tmbed/reviewed_human_alt_isoforms/protein_level.tsv`.
#
# Includes splicing-introduced soluble forms of canonical-TM parents
# ("lost TM") — i.e. the biologically interesting class. The protein_level.tsv
# carries `parent_canonical_has_tm` so callers can filter to lost_tm
# (parent had TM but alt no longer does) without needing a separate named set.

TMBED_ALT_ISOFORM_TSV = predictions_path(
    'tmbed', 'reviewed_human_alt_isoforms', 'protein_level.tsv'
)


def _load_tmbed_alt_predictions():
    """Load per-protein tmbed predictions for alt isoforms.

    Returns DataFrame with Entry, topology_class, has_tm, n_tm, has_signal_peptide,
    parent, parent_canonical_has_tm.
    Raises FileNotFoundError if the prediction TSV hasn't been generated.
    """
    import os, pandas as pd
    if not os.path.exists(TMBED_ALT_ISOFORM_TSV):
        raise FileNotFoundError(
            f'expected tmbed predictions at {TMBED_ALT_ISOFORM_TSV}. '
            'Run the tmbed pipeline first (see predictions/tmbed/).'
        )
    df = pd.read_csv(TMBED_ALT_ISOFORM_TSV, sep='\t')
    df['has_tm'] = df['has_tm'].astype(str).str.lower().isin(['true', '1'])
    df['parent_canonical_has_tm'] = df['parent_canonical_has_tm'].astype(str).str.lower().isin(['true', '1'])
    return df


def alt_isoforms_no_tm(refresh=False):
    """Alt isoforms predicted by tmbed to have no transmembrane segment.

    Tmbed runs ProtT5+CNN+CRF and labels per-residue topology
    (H/B/S/i/o). 'no TM' = no H or B in the predicted topology.
    Includes splicing-introduced soluble alts of canonical-TM proteins.
    """
    alts = reviewed_human_alt_isoforms(refresh=refresh)
    pred = _load_tmbed_alt_predictions()
    no_tm_entries = set(pred.loc[~pred['has_tm'], 'Entry'])
    return alts[alts['Entry'].isin(no_tm_entries)].reset_index(drop=True)


def alt_isoforms_lost_tm(refresh=False):
    """Alt isoforms whose parent canonical has UniProt TM but tmbed predicts
    the alt as TM-less. Splicing-induced loss-of-TM candidates — likely
    soluble forms of membrane proteins.
    """
    alts = reviewed_human_alt_isoforms(refresh=refresh)
    pred = _load_tmbed_alt_predictions()
    mask_entries = set(pred.loc[(~pred['has_tm']) & pred['parent_canonical_has_tm'], 'Entry'])
    return alts[alts['Entry'].isin(mask_entries)].reset_index(drop=True)


def signal_peptide(strict=False, refresh=False):
    """Reviewed human proteins with a signal peptide annotation.

    strict=True: at least one signal peptide annotation has experimental evidence.
    """
    df = reviewed_human(refresh=refresh)
    col = 'Signal peptide'
    if col not in df.columns:
        raise KeyError(f'expected column {col!r}; have {list(df.columns)[:10]}...')
    mask = _has_value(df, col)
    if strict:
        mask &= df[col].str.contains(STRICT_CODE, na=False)
    return df[mask].reset_index(drop=True)


def transmem(strict=False, refresh=False):
    """Reviewed human proteins with a transmembrane annotation."""
    df = reviewed_human(refresh=refresh)
    col = 'Transmembrane'
    if col not in df.columns:
        raise KeyError(f'expected column {col!r}; have {list(df.columns)[:10]}...')
    mask = _has_value(df, col)
    if strict:
        mask &= df[col].str.contains(STRICT_CODE, na=False)
    return df[mask].reset_index(drop=True)


def secretome_uniprot(refresh=False):
    """Signal peptide yes, transmembrane no. UniProt annotation-based secretome."""
    df = reviewed_human(refresh=refresh)
    sp_col, tm_col = 'Signal peptide', 'Transmembrane'
    mask = _has_value(df, sp_col) & ~_has_value(df, tm_col)
    return df[mask].reset_index(drop=True)


def sp_strict_no_tm(refresh=False):
    """Experimental signal peptide, no transmembrane. Strict-evidence secretome."""
    df = reviewed_human(refresh=refresh)
    sp_col, tm_col = 'Signal peptide', 'Transmembrane'
    mask = _has_value(df, sp_col) & df[sp_col].str.contains(STRICT_CODE, na=False) & ~_has_value(df, tm_col)
    return df[mask].reset_index(drop=True)


def no_sp_no_tm(refresh=False):
    """No signal peptide AND no transmembrane. Clean negatives for SP classification."""
    df = reviewed_human(refresh=refresh)
    sp_col, tm_col = 'Signal peptide', 'Transmembrane'
    mask = ~_has_value(df, sp_col) & ~_has_value(df, tm_col)
    return df[mask].reset_index(drop=True)


def secretome_hpa(refresh=False):
    """Human secretome as defined by Human Protein Atlas (Predicted secreted proteins).

    Returns reviewed UniProt human entries whose accession is in HPA's curated
    "Predicted secreted proteins" list. We intersect with reviewed_human so the
    columns (Sequence, Pfam, etc.) come from the UniProt cache and downstream
    set/embedding code can treat it like any other proteome filter.
    """
    df = reviewed_human(refresh=refresh)
    hpa = hpa_secretome_uniprots(refresh=refresh)
    return df[df['Entry'].isin(hpa)].reset_index(drop=True)


def secretome(refresh=False):
    """Union of HPA-predicted secretome and UniProt SP+/TM- annotated secretome.

    The default umbrella "secretome" — broader than either source alone:
    everything HPA flags as predicted-secreted, plus everything UniProt
    annotates with a signal peptide and no transmembrane segment.
    """
    hpa_df = secretome_hpa(refresh=refresh)
    sp_df = secretome_uniprot(refresh=refresh)
    entries = set(hpa_df['Entry']).union(sp_df['Entry'])
    df = reviewed_human(refresh=refresh)
    return df[df['Entry'].isin(entries)].reset_index(drop=True)


def secretome_with_sp(refresh=False):
    """Union of HPA secretome and all UniProt signal-peptide proteins (TM allowed)."""
    hpa_df = secretome_hpa(refresh=refresh)
    sp_df = signal_peptide(refresh=refresh)
    entries = set(hpa_df['Entry']).union(sp_df['Entry'])
    df = reviewed_human(refresh=refresh)
    return df[df['Entry'].isin(entries)].reset_index(drop=True)


def no_secretome_no_tm(refresh=False):
    """reviewed_human minus secretome (HPA ∪ UniProt SP+/TM-) minus transmem.

    Strictest clean negatives for SP/secretion classifiers: removes anything
    flagged secreted by either source AND anything with a transmembrane
    annotation. Equivalent to (no_sp_no_tm) ∩ (not in HPA secretome).
    """
    df = reviewed_human(refresh=refresh)
    sec = set(secretome(refresh=refresh)['Entry'])
    tm = set(transmem(refresh=refresh)['Entry'])
    drop = sec | tm
    return df[~df['Entry'].isin(drop)].reset_index(drop=True)


def surfaceome(confidence='trusted', refresh=False):
    """Human cell-surface proteome from CSPA (Bausch-Fluck et al. 2015, PLoS One).

    Intersected with reviewed_human so downstream code gets the full UniProt
    annotation (Sequence, Pfam, etc.) for each hit.

    confidence : 'high' → CSPA category 1 only
                 'trusted' → categories 1+2 (excludes '3 - unspecific'), default
                 'all' → all three categories
    """
    df = reviewed_human(refresh=refresh)
    cspa = cspa_uniprots(organism='Human', confidence=confidence, refresh=refresh)
    return df[df['Entry'].isin(cspa)].reset_index(drop=True)


def transcription_factors(candidates=False, refresh=False):
    """Curated human transcription factors (Lambert et al. 2018, humantfs.ccbr.utoronto.ca).

    Intersected with reviewed_human so downstream code gets full UniProt
    annotation for each TF. This is the field-standard human TF census,
    cleaner than any single UniProt keyword / GO term.

    candidates : if False (default) only the 1,639 curated bona fide TFs;
                 if True, all 2,765 assessed candidate genes (includes the
                 1,126 manually rejected as non-TFs).
    """
    df = reviewed_human(refresh=refresh)
    tf = humantfs_uniprots(only_tfs=not candidates, refresh=refresh)
    return df[df['Entry'].isin(tf)].reset_index(drop=True)


# --- Ortholog sets: every ortholog of every member of a base set ---

def _orthologs_for_set(base_fn, **kw):
    """Return DataFrame of orthologs for every member of `base_fn(**kw)` that has an OMA group.

    Columns: human_uniprot, oma_group, oma_id, uniprot (mapped, may be empty), sequence.
    """
    import pandas as pd
    base = base_fn(**kw)
    cov = _oma.load_coverage()
    have = cov[cov['uniprot'].isin(base['Entry']) & cov['oma_group'].notna()]
    mapping = _oma.load_uniprot_mapping().drop_duplicates('oma_id').set_index('oma_id')['uniprot']

    rows = []
    for _, r in have.iterrows():
        gid = int(r['oma_group'])
        try:
            msa = _oma.load_msa(gid, fetch=False)
        except FileNotFoundError:
            continue
        for oid, aligned in msa.items():
            rows.append({
                'human_uniprot': r['uniprot'],
                'oma_group':     gid,
                'oma_id':        oid,
                'uniprot':       mapping.get(oid, ''),
                'sequence':      aligned.replace('-', ''),
            })
    return pd.DataFrame(rows)


def transmem_orthologs(strict=False, refresh=False):
    """All orthologs of every reviewed human transmembrane protein with an OMA group."""
    return _orthologs_for_set(transmem, strict=strict, refresh=refresh)


# Registry: set name → (function, kwargs)
NAMED_SETS = {
    'reviewed_human':              (reviewed_human, {}),
    'reviewed_human_alt_isoforms': (reviewed_human_alt_isoforms, {}),
    'alt_isoforms_no_tm':          (alt_isoforms_no_tm, {}),
    'alt_isoforms_lost_tm':        (alt_isoforms_lost_tm, {}),
    'signal_peptide':         (signal_peptide, {'strict': False}),
    'signal_peptide_strict':  (signal_peptide, {'strict': True}),
    'transmem':               (transmem, {'strict': False}),
    'transmem_strict':        (transmem, {'strict': True}),
    'secretome':              (secretome, {}),
    'secretome_hpa':          (secretome_hpa, {}),
    'secretome_uniprot':      (secretome_uniprot, {}),
    'secretome_with_sp':      (secretome_with_sp, {}),
    'sp_strict_no_tm':        (sp_strict_no_tm, {}),
    'no_sp_no_tm':            (no_sp_no_tm, {}),
    'no_secretome_no_tm':     (no_secretome_no_tm, {}),
    'surfaceome':             (surfaceome, {'confidence': 'trusted'}),
    'surfaceome_strict':      (surfaceome, {'confidence': 'high'}),
    'surfaceome_all':         (surfaceome, {'confidence': 'all'}),
    'transcription_factors':  (transcription_factors, {}),
    'transmem_orthologs':     (transmem_orthologs, {'strict': False}),
}


def resolve_dynamic(name, refresh=False):
    """Resolve a named set to (ids, sequences_dict, df). Returns None if name not registered."""
    if name not in NAMED_SETS:
        return None
    fn, kw = NAMED_SETS[name]
    df = fn(**kw, refresh=refresh) if 'refresh' in fn.__code__.co_varnames else fn(**kw)
    ids = df['Entry'].tolist() if 'Entry' in df.columns else df['accession'].tolist()
    seqs = {}
    seq_col = 'Sequence' if 'Sequence' in df.columns else 'sequence'
    if seq_col in df.columns:
        seqs = dict(zip(ids, df[seq_col].tolist()))
    return ids, seqs, df


def list_dynamic():
    """Return the list of registered dynamic set names."""
    return sorted(NAMED_SETS.keys())


def dedupe(set_or_df, level=50, keep='first'):
    """Reduce a protein set to one representative per UniRef cluster.

    Parameters
    ----------
    set_or_df : str or DataFrame
        Named set (e.g. ``'signal_peptide_strict'``) or a DataFrame that
        already carries the relevant ``UniRef{level}`` column.
    level : int
        UniRef identity level: 50, 90, or 100.
    keep : 'first' | 'longest'
        Which row to keep per cluster. ``'longest'`` keeps the longest
        sequence (requires the ``Length`` column).

    Notes
    -----
    Rows with an empty cluster ID are kept as-is — UniProt's mapping does
    leave a small number of recently-added entries unclustered at the
    lower-identity levels; treating them as singletons is correct.
    """
    if level not in (50, 90, 100):
        raise ValueError(f'level must be 50, 90, or 100; got {level}')
    col = f'UniRef{level}'

    if isinstance(set_or_df, str):
        result = resolve_dynamic(set_or_df)
        if result is None:
            raise KeyError(f'unknown set {set_or_df!r}; see list_dynamic()')
        df = result[2]
    else:
        df = set_or_df
    if col not in df.columns:
        raise KeyError(
            f'expected {col!r}; load_proteome() merges it in by default. '
            f'Have: {[c for c in df.columns if "UniRef" in c]}'
        )

    has_id = df[col] != ''
    clustered = df[has_id]
    singletons = df[~has_id]

    if keep == 'first':
        clustered = clustered.drop_duplicates(subset=col, keep='first')
    elif keep == 'longest':
        if 'Length' not in clustered.columns:
            raise KeyError("keep='longest' requires a 'Length' column")
        order = clustered['Length'].astype(int).sort_values(ascending=False).index
        clustered = clustered.loc[order].drop_duplicates(subset=col, keep='first')
    else:
        raise ValueError(f"keep must be 'first' or 'longest'; got {keep!r}")

    import pandas as pd
    return pd.concat([clustered, singletons]).sort_index().reset_index(drop=True)
