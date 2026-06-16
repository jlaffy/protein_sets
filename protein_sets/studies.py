"""Bespoke studies — external datasets registered as named entities.

A study is a directory under ``STUDIES_ROOT`` with the standard layout:

    studies/<name>/
        raw/                # original downloads, untouched
        parsed/
            rows.parquet    # standardized row contract (required at load time)
        metadata.yaml       # required
        parse.py            # raw → parsed, re-runnable, not imported at runtime

Once a study has ``parsed/rows.parquet`` and ``metadata.yaml``, it's first-class:
``studies.load(name)`` returns rows on the same contract as ``regions``, so it
flows through ``esm_embed.embed_regions``, structure fetching, and clustering
with no special-casing. The resolver in ``lookup.resolve`` also accepts study
names, returning parent UniProt IDs.

Each study's parse logic stays bespoke (lives in its own ``parse.py``); the
package only standardizes the output shape.
"""

import os

import pandas as pd
import yaml

from .layout import attach_provenance
from ._config import data_path


STUDIES_ROOT = data_path('studies')


def _study_dir(name, root=None):
    return os.path.join(root or STUDIES_ROOT, name)


def available(root=None):
    """Return sorted study names that have both rows.parquet and metadata.yaml."""
    r = root or STUDIES_ROOT
    if not os.path.isdir(r):
        return []
    names = []
    for n in sorted(os.listdir(r)):
        d = os.path.join(r, n)
        if (os.path.isfile(os.path.join(d, 'metadata.yaml'))
                and os.path.isfile(os.path.join(d, 'parsed', 'rows.parquet'))):
            names.append(n)
    return names


def path(name, root=None):
    """Directory path for `name`."""
    return _study_dir(name, root)


def info(name, root=None):
    """Parsed metadata.yaml for `name` as a dict."""
    with open(os.path.join(_study_dir(name, root), 'metadata.yaml')) as f:
        return yaml.safe_load(f)


def load(name, root=None, **filters):
    """Return rows.parquet for `name`, optionally filtered by meta columns.

    Examples
    --------
    >>> df = load('schaeffer22')
    >>> df = load('schaeffer22', cell_type='HEK293')
    """
    meta = info(name, root)
    df = pd.read_parquet(os.path.join(_study_dir(name, root), 'parsed', 'rows.parquet'))
    for col, val in filters.items():
        if isinstance(val, (list, tuple, set)):
            df = df[df[col].isin(val)]
        else:
            df = df[df[col] == val]
    df = df.reset_index(drop=True)

    prov = meta.get('provenance') or {}
    attach_provenance(
        df,
        set=name,
        operation='studies.load',
        source_category=prov.get('source_category', 'studies'),
        subgroup=prov.get('subgroup', name),
        kind=prov.get('kind'),
        args={'filters': filters} if filters else {},
        universe=meta.get('parent_universe', 'human'),
        study=name,
    )
    return df
