"""Pull and cache the Cell Surface Protein Atlas (CSPA).

Source: Bausch-Fluck et al. (2015) PLoS One 10: e0121314.
Single xlsx (S2_File.xlsx) hosted at wlab.ethz.ch with three sheets:
  - Table A:  1,492 human surface proteins
  - Table B:  1,296 mouse surface proteins
  - Table_C: ~13,942 N-glycopeptides (both organisms)

Each protein row carries a CSPA confidence category:
  1 - high confidence   (the trustworthy call)
  2 - putative          (identified but less certain)
  3 - unspecific        (detected but likely noise)
"""

import os
import time
import urllib.request

import pandas as pd

from ._cache import cached_pull

CSPA_CACHE = '/ewsc/jlaffy/protein_sets/cspa'
DEFAULT_REFRESH_DAYS = 30

CSPA_URL = 'https://wlab.ethz.ch/CSPA/data/S2_File.xlsx'

SHEET_HUMAN = 'Table A'
SHEET_MOUSE = 'Table B'
SHEET_GLYCOPEPTIDES = 'Table_C'

# Category labels as they appear in the xlsx.
CATEGORY_HIGH = '1 - high confidence'
CATEGORY_PUTATIVE = '2 - putative'
CATEGORY_UNSPECIFIC = '3 - unspecific'


def _xlsx_path(root=None):
    return os.path.join(root or CSPA_CACHE, 'S2_File.xlsx')


def pull_cspa(root=None, refresh=False, max_days=DEFAULT_REFRESH_DAYS):
    """Download the CSPA S2_File.xlsx. Cached under CSPA_CACHE.

    Returns the path to the cached xlsx. Auto-refreshes if older than `max_days`.
    On refresh failure with a prior cache present, returns the stale cache and
    emits a RuntimeWarning rather than crashing.
    """
    out = _xlsx_path(root)

    def _fetch(tmp):
        print(f'pulling CSPA from {CSPA_URL} ...')
        t0 = time.time()
        with urllib.request.urlopen(CSPA_URL, timeout=600) as r:
            with open(tmp, 'wb') as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
        print(f'  saved {out} ({os.path.getsize(tmp) / 1e6:.1f} MB, {time.time() - t0:.1f}s)')

    return cached_pull(out, _fetch, refresh=refresh, max_days=max_days)


def _load_sheet(sheet, root=None, refresh=False):
    p = pull_cspa(root=root, refresh=refresh)
    df = pd.read_excel(p, sheet_name=sheet)
    df.columns = [c.strip() for c in df.columns]
    return df


def load_cspa_human(root=None, refresh=False):
    """Human surfaceome table (1,492 proteins). UniProt accession is in `ID_link`."""
    return _load_sheet(SHEET_HUMAN, root=root, refresh=refresh)


def load_cspa_mouse(root=None, refresh=False):
    """Mouse surfaceome table (1,296 proteins). UniProt accession is in `ID_link`."""
    return _load_sheet(SHEET_MOUSE, root=root, refresh=refresh)


def load_cspa_glycopeptides(root=None, refresh=False):
    """N-glycopeptides table (both organisms, ~13,942 rows). UniProt accession is in `ID link`."""
    return _load_sheet(SHEET_GLYCOPEPTIDES, root=root, refresh=refresh)


def cspa_uniprots(organism='Human', confidence='high', refresh=False):
    """Return the set of UniProt accessions in CSPA for the given organism.

    organism : 'Human' | 'Mouse'
    confidence : 'high' → category 1 only
                 'trusted' → categories 1+2 (high + putative), excludes unspecific
                 'all' → all three categories including unspecific
    """
    if organism == 'Human':
        df = load_cspa_human(refresh=refresh)
    elif organism == 'Mouse':
        df = load_cspa_mouse(refresh=refresh)
    else:
        raise ValueError(f'organism must be Human or Mouse, got {organism!r}')

    if confidence == 'high':
        df = df[df['CSPA category'] == CATEGORY_HIGH]
    elif confidence == 'trusted':
        df = df[df['CSPA category'].isin([CATEGORY_HIGH, CATEGORY_PUTATIVE])]
    elif confidence == 'all':
        pass
    else:
        raise ValueError(f'confidence must be high|trusted|all, got {confidence!r}')

    return set(df['ID_link'].dropna().astype(str))
