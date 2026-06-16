"""Pull and cache Human Protein Atlas tables (secretome, subcellular location).

All HPA TSVs come from the search-as-download endpoint, which returns one row
per protein with whatever columns you ask for via ?columns=...
"""

import os
import time
import urllib.request
import urllib.parse

import pandas as pd

from ._cache import cached_pull
from ._config import data_path

HPA_CACHE = data_path('hpa')
DEFAULT_REFRESH_DAYS = 30

# HPA "Predicted secreted proteins" search → TSV download.
# Their endpoint accepts ?format=tsv&download=yes; no compression flag needed.
HPA_SECRETED_QUERY = 'protein_class:Predicted secreted proteins'

# Subcellular: pull every entry, request only the columns we need.
# Empty query string = all HPA entries; we filter to rows with any subloc
# annotation at load time.
HPA_SUBLOC_QUERY = ''
HPA_SUBLOC_COLUMNS = 'g,gs,up,scml,scal,scel,relsc'


def _cache_path(name='predicted_secreted', root=None):
    root = root or HPA_CACHE
    return os.path.join(root, f'{name}.tsv')


def pull_hpa_secretome(name='predicted_secreted', root=None,
                      refresh=False, max_days=DEFAULT_REFRESH_DAYS):
    """Download HPA's predicted secreted proteins list as TSV. Cached.

    Returns the path to the cached file. Auto-refreshes if older than `max_days`.
    On refresh failure with a prior cache present, returns the stale cache and
    emits a RuntimeWarning rather than crashing.
    """
    out = _cache_path(name, root)
    params = {'format': 'tsv', 'download': 'yes'}
    url = ('https://www.proteinatlas.org/search/'
           + urllib.parse.quote(HPA_SECRETED_QUERY)
           + '?' + urllib.parse.urlencode(params))

    def _fetch(tmp):
        print(f'pulling HPA Predicted Secreted Proteins...')
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


def load_hpa_secretome(name='predicted_secreted', root=None, refresh=False):
    """Load cached HPA secretome TSV as a DataFrame. Pulls if missing or stale."""
    p = pull_hpa_secretome(name, root=root, refresh=refresh)
    df = pd.read_csv(p, sep='\t', dtype=str).fillna('')
    df.columns = [c.strip() for c in df.columns]
    return df


def pull_hpa_subloc(name='subcellular_location', root=None,
                    refresh=False, max_days=DEFAULT_REFRESH_DAYS):
    """Download HPA subcellular location TSV via the search-as-download endpoint.

    Returns the path to the cached file. Columns requested: Gene, Gene synonym,
    Uniprot, Subcellular main location, Subcellular additional location,
    Subcellular extracellular location, Subcellular location reliability.
    """
    out = _cache_path(name, root)
    params = {
        'format': 'tsv',
        'download': 'yes',
        'columns': HPA_SUBLOC_COLUMNS,
    }
    url = ('https://www.proteinatlas.org/search/'
           + urllib.parse.quote(HPA_SUBLOC_QUERY)
           + '?' + urllib.parse.urlencode(params))

    def _fetch(tmp):
        print(f'pulling HPA subcellular location...')
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


def load_hpa_subloc(name='subcellular_location', root=None, refresh=False,
                    only_annotated=True):
    """Load cached HPA subcellular location TSV as a DataFrame.

    Columns (after column-name normalization): Gene, Gene synonym, Uniprot,
    Subcellular main location, Subcellular additional location,
    Subcellular extracellular location, Subcellular location reliability.

    With only_annotated=True (default), drop rows where every subcellular
    column is empty.
    """
    p = pull_hpa_subloc(name, root=root, refresh=refresh)
    df = pd.read_csv(p, sep='\t', dtype=str).fillna('')
    df.columns = [c.strip() for c in df.columns]
    if only_annotated:
        loc_cols = [c for c in df.columns if c.startswith('Subcellular')
                    and 'reliability' not in c.lower()]
        if loc_cols:
            mask = df[loc_cols].apply(lambda s: s.str.len() > 0).any(axis=1)
            df = df[mask].reset_index(drop=True)
    return df


_HPA_SUBLOC_SLIM_COLS = [
    'Uniprot', 'Gene',
    'Subcellular main location',
    'Subcellular additional location',
    'Subcellular location',
    'Reliability (IF)',
]


def hpa_subloc(query=None, slim=True, refresh=False):
    """Convenience view over the HPA subcellular table.

    `query`: None (all annotated entries), a list of UniProt accessions, or a
    DataFrame carrying an `Entry`/`uniprot` column to intersect against.

    `slim=True` returns just the subloc-relevant columns. Set False for the
    full HPA proteinatlas table (108 columns).
    """
    df = load_hpa_subloc(refresh=refresh)
    if slim:
        cols = [c for c in _HPA_SUBLOC_SLIM_COLS if c in df.columns]
        df = df[cols].copy()
    if query is not None:
        if hasattr(query, 'columns'):
            ids = set(query['Entry'] if 'Entry' in query.columns else query['uniprot'])
        else:
            ids = set(query)
        df = df[df['Uniprot'].isin(ids)].reset_index(drop=True)
    return df


def hpa_secretome_uniprots(refresh=False):
    """Return the set of UniProt accessions in HPA's predicted secretome.

    HPA rows can have empty or comma-separated UniProt fields; this flattens them.
    """
    df = load_hpa_secretome(refresh=refresh)
    if 'Uniprot' not in df.columns:
        raise KeyError(f'expected column "Uniprot"; have {list(df.columns)[:10]}...')
    out = set()
    for u in df['Uniprot']:
        for x in u.replace(' ', '').split(','):
            if x:
                out.add(x)
    return out
