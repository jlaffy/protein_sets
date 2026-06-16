"""Load and cache InterPro domain annotations for human proteins.

Two cached files:
  - human_protein2ipr.tsv  — per-protein domain boundaries from all member DBs
  - interpro_entry_types.tsv — InterPro entry classification (Domain, Family, etc.)

Currently loaded from a pre-filtered local copy. A future pull_interpro() will
fetch from the InterPro REST API for auto-refresh.
"""

import os
from datetime import datetime

import pandas as pd

INTERPRO_CACHE = '/ewsc/jlaffy/protein_sets/interpro'

_PROTEIN2IPR_COLS = ['uniprot', 'interpro_id', 'interpro_name',
                     'member_db_acc', 'start', 'end']

# Member DB prefix → human-readable name
DB_NAMES = {
    'PF': 'pfam', 'PS': 'prosite', 'SM': 'smart', 'SSF': 'superfamily',
    'PR': 'prints', 'G': 'gene3d', 'PTHR': 'panther', 'cd': 'cdd',
    'PIRSF': 'pirsf', 'TIGR': 'tigrfams', 'MF': 'hamap',
    'SFLDF': 'sfld_family', 'SFLDS': 'sfld_superfamily', 'SFLDG': 'sfld_group',
    'NF': 'ncbifam',
}

# Reverse: name → prefix(es)
_NAME_TO_PREFIX = {}
for _pfx, _name in DB_NAMES.items():
    _NAME_TO_PREFIX.setdefault(_name, []).append(_pfx)


def cache_age_days(root=None):
    """Days since the protein2ipr cache was last modified, or +inf if missing."""
    p = os.path.join(root or INTERPRO_CACHE, 'human_protein2ipr.tsv')
    if not os.path.exists(p):
        return float('inf')
    age = datetime.now() - datetime.fromtimestamp(os.path.getmtime(p))
    return age.total_seconds() / 86400


def load_protein2ipr(root=None):
    """Load human_protein2ipr.tsv as a DataFrame.

    Columns: uniprot, interpro_id, interpro_name, member_db_acc, start, end
    """
    p = os.path.join(root or INTERPRO_CACHE, 'human_protein2ipr.tsv')
    if not os.path.exists(p):
        raise FileNotFoundError(
            f'InterPro cache not found at {p}. '
            f'Copy human_protein2ipr.tsv into {root or INTERPRO_CACHE}/')
    df = pd.read_csv(p, sep='\t', header=None, names=_PROTEIN2IPR_COLS, dtype=str)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    return df


def load_entry_types(root=None):
    """Load interpro_entry_types.tsv → DataFrame with interpro_id, entry_type, entry_name."""
    p = os.path.join(root or INTERPRO_CACHE, 'interpro_entry_types.tsv')
    if not os.path.exists(p):
        raise FileNotFoundError(f'Entry types file not found at {p}')
    df = pd.read_csv(p, sep='\t', dtype=str)
    df.columns = ['interpro_id', 'entry_type', 'entry_name']
    return df


def load_interpro(root=None, with_entry_types=True):
    """Load protein2ipr with member_db and (optionally) entry_type columns added.

    Extra columns:
        member_db   — human-readable DB name (pfam, smart, gene3d, ...)
        entry_type  — InterPro classification (Domain, Family, Homologous_superfamily, ...)
    """
    df = load_protein2ipr(root)
    df['member_db'] = (df['member_db_acc']
                       .str.extract(r'^([A-Za-z]+)', expand=False)
                       .map(DB_NAMES))
    if with_entry_types:
        et = load_entry_types(root)
        df = df.merge(et[['interpro_id', 'entry_type']], on='interpro_id', how='left')
    return df


def db_prefixes(db_name):
    """Return the member_db_acc prefix(es) for a given DB name.

    E.g. db_prefixes('pfam') → ['PF'], db_prefixes('gene3d') → ['G']
    """
    name = db_name.lower()
    if name not in _NAME_TO_PREFIX:
        raise KeyError(f'unknown DB {db_name!r}; known: {sorted(_NAME_TO_PREFIX)}')
    return _NAME_TO_PREFIX[name]
