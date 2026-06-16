"""Common cache+refresh wrapper for pull_* functions.

Every external data source in this package follows the same pattern:
"if a cached file is missing or older than N days, refetch it from the
upstream API." This helper centralizes that pattern, plus two safety nets:

1. **Atomic writes.** Every refresh writes to ``<path>.tmp`` first, then
   ``os.replace``s onto the final path. A network failure mid-stream can
   no longer leave a half-written cache file in place.
2. **Stale-fallback on refresh failure.** If the fetch raises *and* a
   previous valid cache exists, log a RuntimeWarning and return the
   stale cache. Stale data is far better than a cascade where one
   upstream-API outage breaks every consumer of every pull function.
   First-time pulls (no prior cache) still raise, since there is
   nothing to fall back to.
"""

import os
import warnings
from datetime import datetime


def cache_age_days(path):
    """Days since `path` was last modified, or +inf if missing."""
    if not os.path.exists(path):
        return float('inf')
    age = datetime.now() - datetime.fromtimestamp(os.path.getmtime(path))
    return age.total_seconds() / 86400


def cached_pull(out_path, fetch_fn, refresh=False, max_days=30):
    """Resilient cache-or-fetch.

    Parameters
    ----------
    out_path : str
        Final path of the cached file.
    fetch_fn : callable(tmp_path)
        Writes the freshly-fetched content to `tmp_path`. Caller does not
        need to handle directories, atomic renames, or fallback logic.
    refresh : bool
        If True, ignore cache and refetch.
    max_days : int
        Cache TTL.

    Returns
    -------
    str
        Path to a usable cache file (either the freshly-pulled one, or
        the prior stale one if the refresh failed).
    """
    if not refresh and cache_age_days(out_path) < max_days:
        return out_path
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp = out_path + '.tmp'
    try:
        fetch_fn(tmp)
        os.replace(tmp, out_path)
    except Exception as e:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass
        if os.path.exists(out_path):
            warnings.warn(
                f'refresh failed ({type(e).__name__}: {e}); '
                f'using stale cache at {out_path}',
                RuntimeWarning, stacklevel=2,
            )
            return out_path
        raise
    return out_path
