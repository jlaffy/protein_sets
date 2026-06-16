"""Central data-root configuration.

All cached data for this package lives under a single root directory. The root
is resolved once, at import, in this priority order:

1. ``$PROTEIN_SETS_DATA`` — explicit override (used by off-cluster installs).
2. ``/ewsc/jlaffy/protein_sets`` — the original shared NFS cache, IF it exists.
   This keeps the canonical cluster setup behaving exactly as before: same
   paths, same data, no re-download.
3. ``~/.cache/protein_sets`` — XDG-style default for a fresh install on a
   machine with neither the env var nor the shared cache.

Every per-source cache (uniprot, uniref, hpa, cspa, humantfs, interpro, oma,
studies) derives from this root via ``data_path(...)`` so there is a single
place to point the whole package at a different location.

The predictions root (tmbed alt-isoform TSVs) is separate because it lives
outside the protein_sets cache tree on the cluster; it follows the same
override pattern via ``$PROTEIN_SETS_PREDICTIONS``.
"""

import os

_LEGACY_DATA_ROOT = '/ewsc/jlaffy/protein_sets'
_LEGACY_PREDICTIONS_ROOT = '/ewsc/jlaffy/predictions'


def _resolve(env_var, legacy_root, default_subdir):
    override = os.environ.get(env_var)
    if override:
        return os.path.expanduser(override)
    if os.path.isdir(legacy_root):
        return legacy_root
    return os.path.join(os.path.expanduser('~'), '.cache', default_subdir)


DATA_ROOT = _resolve('PROTEIN_SETS_DATA', _LEGACY_DATA_ROOT, 'protein_sets')
PREDICTIONS_ROOT = _resolve(
    'PROTEIN_SETS_PREDICTIONS', _LEGACY_PREDICTIONS_ROOT, 'protein_sets_predictions'
)


def data_path(*parts):
    """Join `parts` onto the resolved data root."""
    return os.path.join(DATA_ROOT, *parts)


def predictions_path(*parts):
    """Join `parts` onto the resolved predictions root."""
    return os.path.join(PREDICTIONS_ROOT, *parts)
