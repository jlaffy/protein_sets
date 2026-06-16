"""Pull and cache the Human Transcription Factors database (Lambert et al. 2018).

Source: Lambert, Jolma, Campitelli, et al. (2018) "The Human Transcription
Factors." Cell 172: 650-665.  http://humantfs.ccbr.utoronto.ca

The curated database extract (v1.01) is a single CSV of 2,765 candidate genes,
each adjudicated by manual curation into `Is TF?` = Yes (1,639 bona fide
sequence-specific DNA-binding TFs) or No. This is the field-standard human TF
census — substantially cleaner than any single UniProt keyword or GO term.

The file is keyed by Ensembl gene ID + HGNC symbol (no UniProt accession), so
`humantfs_uniprots()` maps gene symbols onto the reviewed-human proteome
(primary symbol first, then any gene-name synonym). ~1,632 of the 1,639 curated
TFs map; the handful that don't are novel/read-through loci with no reviewed
UniProt entry.
"""

import os
import time
import urllib.request

import pandas as pd

from ._cache import cached_pull
from ._config import data_path

HUMANTFS_CACHE = data_path('humantfs')
DEFAULT_REFRESH_DAYS = 30

HUMANTFS_URL = ('http://humantfs.ccbr.utoronto.ca/download/v_1.01/'
                'DatabaseExtract_v_1.01.csv')

IS_TF_COL = 'Is TF?'
SYMBOL_COL = 'HGNC symbol'
ENSEMBL_COL = 'Ensembl ID'


def _csv_path(root=None):
    return os.path.join(root or HUMANTFS_CACHE, 'DatabaseExtract_v_1.01.csv')


def pull_humantfs(root=None, refresh=False, max_days=DEFAULT_REFRESH_DAYS):
    """Download the humantfs v1.01 database extract CSV. Cached under HUMANTFS_CACHE.

    Returns the path to the cached CSV. Auto-refreshes if older than `max_days`;
    falls back to a stale cache on refresh failure.
    """
    out = _csv_path(root)

    def _fetch(tmp):
        print(f'pulling humantfs from {HUMANTFS_URL} ...')
        t0 = time.time()
        req = urllib.request.Request(HUMANTFS_URL, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=600) as r:
            with open(tmp, 'wb') as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
        print(f'  saved {out} ({os.path.getsize(tmp) / 1e6:.1f} MB, {time.time() - t0:.1f}s)')

    return cached_pull(out, _fetch, refresh=refresh, max_days=max_days)


def load_humantfs(root=None, refresh=False, only_tfs=False):
    """Load the humantfs database extract as a DataFrame (2,765 candidate rows).

    only_tfs : if True, return only the 1,639 curated bona fide TFs (`Is TF?` == 'Yes').
    The leading unnamed index column is dropped.
    """
    p = pull_humantfs(root=root, refresh=refresh)
    df = pd.read_csv(p)
    df = df.drop(columns=[c for c in df.columns if c.startswith('Unnamed')], errors='ignore')
    df.columns = [c.strip() for c in df.columns]
    if only_tfs:
        df = df[df[IS_TF_COL].astype(str).str.strip() == 'Yes'].reset_index(drop=True)
    return df


def _symbol_to_entry(refresh=False):
    """Map every gene symbol (primary + synonyms) in reviewed_human to a UniProt Entry."""
    from .uniprot import load_proteome
    prot = load_proteome(refresh=refresh)
    primary, synonym = {}, {}
    for entry, prim, names in zip(prot['Entry'],
                                  prot['Gene Names (primary)'].fillna(''),
                                  prot['Gene Names'].fillna('')):
        if prim:
            primary.setdefault(prim, entry)
        for n in names.split():
            synonym.setdefault(n, entry)
    return primary, synonym


def humantfs_uniprots(only_tfs=True, refresh=False):
    """Return the set of UniProt accessions for human TFs (Lambert et al. 2018).

    only_tfs : if True (default), only the 1,639 curated bona fide TFs;
               if False, all 2,765 assessed candidate genes.

    Gene symbols are mapped onto reviewed_human by primary symbol, falling back
    to any gene-name synonym. Symbols with no reviewed-human match are dropped.
    """
    tf = load_humantfs(refresh=refresh, only_tfs=only_tfs)
    primary, synonym = _symbol_to_entry(refresh=refresh)
    out = set()
    for sym in tf[SYMBOL_COL].dropna().astype(str):
        sym = sym.strip()
        e = primary.get(sym) or synonym.get(sym)
        if e:
            out.add(e)
    return out
