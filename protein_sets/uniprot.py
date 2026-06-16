"""Pull and cache UniProt proteome metadata + sequences as a single TSV.

One download per (proteome, reviewed) → one cached file with as many columns
as possible. All set definitions filter against this cache locally.

UniRef cluster IDs (50/90/100) live in `uniref.py` — UniProtKB's TSV endpoint
doesn't expose them as fields, so they're fetched via the ID Mapping API and
cached as separate sidecars under `/ewsc/jlaffy/protein_sets/uniref/`.
`load_proteome(with_uniref=True)` merges them in at load time.
"""

import os
import re
import time
import urllib.request
import urllib.parse

import pandas as pd

from ._cache import cached_pull

PROTEOME_CACHE = '/ewsc/jlaffy/protein_sets/uniprot/proteomes'
DEFAULT_REFRESH_DAYS = 30

# As many useful fields as we can ask for in a single TSV. Add more here as
# new set definitions need them — single download covers everything.
FIELDS = [
    'accession', 'id', 'reviewed', 'protein_existence',
    'protein_name', 'gene_primary', 'gene_names',
    'organism_id', 'organism_name',
    'length', 'sequence',
    # subcellular / topology
    'cc_subcellular_location',
    'ft_signal', 'ft_transit', 'ft_transmem', 'ft_topo_dom', 'ft_intramem',
    'ft_chain', 'ft_propep', 'ft_init_met', 'ft_peptide',
    'ft_motif',
    # PTM / modifications
    'ft_disulfid', 'ft_carbohyd', 'ft_lipid', 'ft_mod_res',
    # function / disease
    'cc_function', 'cc_disease',
    # cross-refs
    'xref_pfam', 'xref_interpro', 'xref_alphafolddb',
    # ontology
    'go', 'go_c', 'go_f', 'go_p',
    'keyword',
]


def _cache_path(proteome, reviewed=True, root=None):
    root = root or PROTEOME_CACHE
    suffix = '_reviewed' if reviewed else '_all'
    return os.path.join(root, f'{proteome}{suffix}.tsv')


def pull_proteome(proteome='UP000005640', reviewed=True,
                  fields=None, root=None, refresh=False, max_days=DEFAULT_REFRESH_DAYS):
    """Download a UniProt proteome's metadata + sequences as TSV. Cached.

    Returns the path to the cached file. Auto-refreshes if older than `max_days`.
    On refresh failure with a prior cache present, returns the stale cache and
    emits a RuntimeWarning rather than crashing.
    """
    out = _cache_path(proteome, reviewed, root)
    fields = fields or FIELDS
    query = f'proteome:{proteome}'
    if reviewed:
        query += ' AND reviewed:true'
    params = {
        'query': query,
        'format': 'tsv',
        'fields': ','.join(fields),
        'compressed': 'false',
    }
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode(params)

    def _fetch(tmp):
        print(f'pulling {proteome} from UniProt ({len(fields)} fields)...')
        t0 = time.time()
        with urllib.request.urlopen(url, timeout=600) as r:
            with open(tmp, 'wb') as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
        print(f'  saved {out} ({os.path.getsize(tmp) / 1e6:.1f} MB, {time.time() - t0:.1f}s)')

    return cached_pull(out, _fetch, refresh=refresh, max_days=max_days)


def load_proteome(proteome='UP000005640', reviewed=True, refresh=False, root=None,
                  with_uniref=True):
    """Load cached proteome metadata as a DataFrame. Pulls if missing or stale.

    If `with_uniref` (default True), merges in `UniRef50/90/100` cluster IDs as
    additional columns from `uniref.load_uniref_mapping()`. Sidecar files share
    the 30-day TTL refresh logic.
    """
    p = pull_proteome(proteome, reviewed=reviewed, root=root, refresh=refresh)
    df = pd.read_csv(p, sep='\t', dtype=str).fillna('')
    df.columns = [c.strip() for c in df.columns]
    if with_uniref:
        from . import uniref
        for level in (50, 90, 100):
            m = uniref.load_uniref_mapping(proteome, level=level, reviewed=reviewed,
                                            refresh=refresh)
            df = df.merge(m, on='Entry', how='left').fillna('')
    return df


# --- Alternative isoforms (non-canonical) ---
#
# UniProt's `includeIsoform=true` flag returns one row per isoform. The
# canonical isoform appears under the bare accession (e.g. `P09038`); each
# alternative isoform appears under a `-N` suffix (e.g. `P09038-1`,
# `P09038-2`). The canonical's own `-N` form is suppressed in the output, so
# filtering to rows where Entry contains '-' cleanly isolates the alts.
#
# Note: the canonical's logical IsoId number is NOT always `-1`. About 4% of
# multi-isoform entries have canonical = `-2`/`-3`/`-N` for some other N. The
# bare accession always carries the displayed (canonical) sequence regardless.
# We parse the canonical IsoId from the `Alternative products` field and store
# it as a `canonical_isoid` column for downstream callers that need to know.

ALT_ISOFORM_FIELDS = list(FIELDS) + ['cc_alternative_products']

# Captures the whole IsoId list (may be comma-joined: "P61978-1, Q07244-1")
# for the entry whose Sequence is marked Displayed.
_DISPLAYED_ISOID_RE = re.compile(r'IsoId=([^;]+);\s*Sequence=Displayed')


def _parse_canonical_isoid(text, parent=None):
    """Extract the IsoId marked `Sequence=Displayed` from a UniProt
    `Alternative products (isoforms)` field.

    Some entries share an IsoId block across paralogous accessions
    (e.g. `IsoId=P61978-1, Q07244-1; Sequence=Displayed`). When `parent`
    is supplied, picks the IsoId belonging to that accession; otherwise
    returns the first one. Empty string if no match.
    """
    if not text:
        return ''
    m = _DISPLAYED_ISOID_RE.search(text)
    if not m:
        return ''
    isoids = [s.strip() for s in m.group(1).split(',') if s.strip()]
    if not isoids:
        return ''
    if parent is None:
        return isoids[0]
    for iso in isoids:
        if iso == parent or iso.startswith(parent + '-'):
            return iso
    return isoids[0]


def _alt_isoform_cache_path(proteome, root=None):
    root = root or PROTEOME_CACHE
    return os.path.join(root, f'{proteome}_reviewed_alt_isoforms.tsv')


def pull_alt_isoforms(proteome='UP000005640', root=None, refresh=False,
                     max_days=DEFAULT_REFRESH_DAYS):
    """Download all *non-canonical* isoforms for a reviewed proteome. Cached.

    Uses UniProt's `includeIsoform=true` and filters to rows whose `Entry` has
    a `-N` suffix. The canonical of each parent already lives under the bare
    accession in the regular proteome cache (`reviewed_human`); this set is
    strictly the additional alt isoforms.

    Adds a `canonical_isoid` column (parsed from `Alternative products`) so
    each alt row knows which `-N` was canonical for its parent.
    """
    out = _alt_isoform_cache_path(proteome, root)
    fields = ALT_ISOFORM_FIELDS
    query = f'proteome:{proteome} AND reviewed:true'
    params = {
        'query':          query,
        'format':         'tsv',
        'fields':         ','.join(fields),
        'compressed':     'false',
        'includeIsoform': 'true',
    }
    url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode(params)

    def _fetch(tmp):
        print(f'pulling {proteome} alt isoforms from UniProt ({len(fields)} fields)...')
        t0 = time.time()
        with urllib.request.urlopen(url, timeout=600) as r:
            with open(tmp, 'wb') as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
        # Filter to alt-isoform rows + parse canonical_isoid.
        df = pd.read_csv(tmp, sep='\t', dtype=str).fillna('')
        n_total = len(df)
        df = df[df['Entry'].str.contains('-', regex=False)].reset_index(drop=True)
        ap_col = 'Alternative products (isoforms)'
        if ap_col in df.columns:
            df['canonical_isoid'] = df.apply(
                lambda r: _parse_canonical_isoid(
                    r[ap_col], parent=r['Entry'].split('-', 1)[0]
                ),
                axis=1,
            )
        else:
            df['canonical_isoid'] = ''
        df.to_csv(tmp, sep='\t', index=False)
        size_mb = os.path.getsize(tmp) / 1e6
        print(f'  saved {out} ({n_total} rows → {len(df)} alt isoforms, '
              f'{size_mb:.1f} MB, {time.time() - t0:.1f}s)')

    return cached_pull(out, _fetch, refresh=refresh, max_days=max_days)


def load_alt_isoforms(proteome='UP000005640', root=None, refresh=False):
    """Load cached alt-isoform DataFrame. Pulls if missing/stale.

    One row per non-canonical isoform with full sequence. `Entry` carries
    the `-N` suffix; the parent accession can be recovered by stripping
    everything from `-` onwards. `canonical_isoid` records which `-N` was
    canonical for that parent (so callers don't have to assume `-1`).
    """
    p = pull_alt_isoforms(proteome, root=root, refresh=refresh)
    df = pd.read_csv(p, sep='\t', dtype=str).fillna('')
    df.columns = [c.strip() for c in df.columns]
    return df
