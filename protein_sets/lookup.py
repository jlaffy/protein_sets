"""Resolve ID inputs to a flat list of UniProt accessions.

Resolution order:
  1. list/tuple/set                                  → as-is
  2. file path that exists                           → read line-by-line
  3. dynamic set (e.g. 'transmem')                   → run filter on cached proteome
  4. registered study name (e.g. 'schaeffer22')      → parent UniProts from rows.parquet
  5. comma-separated literal                         → split
  6. single token                                    → treat as a UniProt ID
"""

import os

from .definitions import NAMED_SETS, resolve_dynamic


def resolve(query, return_sequences=False):
    """Return a list of UniProt IDs (and optionally a {id: sequence} dict)."""
    if isinstance(query, (list, tuple, set)):
        ids = list(query)
        return (ids, {}) if return_sequences else ids

    q = str(query).strip()

    if os.path.exists(q):
        ids = _read_id_file(q)
        return (ids, {}) if return_sequences else ids

    if q in NAMED_SETS:
        ids, seqs, _ = resolve_dynamic(q)
        return (ids, seqs) if return_sequences else ids

    from . import studies as _studies
    if q in _studies.available():
        df = _studies.load(q)
        ids = df['uniprot'].drop_duplicates().tolist()
        seqs = {}  # study rows carry region sequences, not parent — fetch on demand
        return (ids, seqs) if return_sequences else ids

    if ',' in q:
        ids = [s.strip() for s in q.split(',') if s.strip()]
        return (ids, {}) if return_sequences else ids

    return ([q], {}) if return_sequences else [q]


def resolve_with_seqs(query):
    """Like resolve() but always returns sequences. Falls back to load_proteome()
    for raw IDs / file inputs that don't carry sequences.

    Returns (ids, seqs_dict).
    """
    ids, seqs = resolve(query, return_sequences=True)
    missing = [i for i in ids if i not in seqs]
    if missing:
        from .uniprot import load_proteome
        df = load_proteome()
        m = dict(zip(df['Entry'], df['Sequence']))
        for uid in missing:
            if uid in m:
                seqs[uid] = m[uid]
    still_missing = [i for i in ids if i not in seqs]
    if still_missing:
        raise KeyError(f'no sequence for {len(still_missing)} ID(s): {still_missing[:5]}')
    return ids, seqs


def _read_id_file(path):
    """Read a txt/tsv file: one ID per line, optional header skipped, takes first column."""
    ids = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.lower().startswith('uniprot'):
                continue
            ids.append(line.split('\t')[0].split(',')[0])
    return ids
