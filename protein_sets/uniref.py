"""UniRef cluster IDs and (eventually) full UniRef cluster catalog.

Phase 1 (here today): per-accession cluster lookup via UniProt's ID Mapping
API. For each UniProtKB accession in a proteome, fetch its UniRef50/90/100
cluster IDs and cache as a small 2-column TSV (Entry, UniRef{level}). Used by
`load_proteome(with_uniref=True)` to add cluster-id columns and by
`dedupe(set, level=...)` to collapse paralogous duplicates.

Phase 2 (planned): cross-species cluster catalog from UniRef XML / idmapping
bulk dumps — full member lists, cluster representative metadata, lowest
common ancestor (LCA) taxonomy. Will live in this same module and share the
`UNIREF_CACHE` root.
"""

import json
import os
import re
import time
import urllib.parse
import urllib.request

import pandas as pd

from ._cache import cached_pull
from ._config import data_path

UNIREF_CACHE = data_path('uniref')
DEFAULT_REFRESH_DAYS = 30
IDMAPPING_BASE = 'https://rest.uniprot.org/idmapping'


def _cache_path(proteome, level, reviewed=True, root=None):
    root = root or UNIREF_CACHE
    suffix = '_reviewed' if reviewed else '_all'
    return os.path.join(root, f'{proteome}{suffix}_uniref{level}.tsv')


def _idmapping_run(ids, to_db, poll_timeout=600, poll_interval=2, page_size=500):
    """Submit a UniProt ID-mapping job and return [(from, to), ...]."""
    data = urllib.parse.urlencode({
        'ids': ','.join(ids),
        'from': 'UniProtKB_AC-ID',
        'to': to_db,
    }).encode('utf-8')
    req = urllib.request.Request(f'{IDMAPPING_BASE}/run', data=data, method='POST')
    with urllib.request.urlopen(req, timeout=120) as r:
        job = json.loads(r.read())
    if 'jobId' not in job:
        raise RuntimeError(f'idmapping submit failed: {job}')
    job_id = job['jobId']

    # Poll status until done. The /status endpoint 303-redirects to /results
    # when finished; urllib follows redirects, so a "done" response may arrive
    # either as {"jobStatus": "FINISHED"} or as the results payload itself.
    t0 = time.time()
    while time.time() - t0 < poll_timeout:
        with urllib.request.urlopen(f'{IDMAPPING_BASE}/status/{job_id}', timeout=60) as r:
            status = json.loads(r.read())
        if status.get('jobStatus') == 'FINISHED' or 'results' in status:
            break
        if status.get('jobStatus') in ('RUNNING', 'NEW', None):
            time.sleep(poll_interval)
            continue
        raise RuntimeError(f'idmapping job {job_id} failed: {status}')
    else:
        raise TimeoutError(f'idmapping job {job_id} did not finish in {poll_timeout}s')

    pairs = []
    url = f'{IDMAPPING_BASE}/results/{job_id}?format=tsv&size={page_size}'
    while url:
        with urllib.request.urlopen(url, timeout=poll_timeout) as r:
            text = r.read().decode('utf-8')
            link = r.headers.get('Link', '')
        lines = text.strip().split('\n')
        start = 1 if lines and lines[0].startswith('From\t') else 0
        for line in lines[start:]:
            if '\t' in line:
                pairs.append(tuple(line.split('\t', 1)))
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else None
    return pairs


def fetch_cluster_ids(accessions, level):
    """Look up UniRef{level} cluster IDs for an arbitrary accession list.

    Lower-level helper for callers that have their own list of accessions and
    don't want the proteome-cache-bound `pull_uniref_mapping` semantics.
    """
    if level not in (50, 90, 100):
        raise ValueError(f'level must be 50, 90, or 100; got {level}')
    return _idmapping_run(list(accessions), f'UniRef{level}')


def pull_uniref_mapping(proteome='UP000005640', level=50, reviewed=True,
                         root=None, refresh=False, max_days=DEFAULT_REFRESH_DAYS):
    """Pull UniRef{level} cluster IDs for all accessions in a proteome.

    Uses the UniProt ID Mapping API. Returns path to a cached 2-column TSV
    (`Entry`, `UniRef{level}`). Auto-refreshes if older than `max_days`. On
    refresh failure with a prior cache present, returns the stale cache and
    emits a RuntimeWarning rather than crashing.
    """
    if level not in (50, 90, 100):
        raise ValueError(f'level must be 50, 90, or 100; got {level}')
    out = _cache_path(proteome, level, reviewed, root)

    def _fetch(tmp):
        # Late import to avoid circular: uniprot.load_proteome uses us back.
        from .uniprot import pull_proteome
        proteome_path = pull_proteome(proteome, reviewed=reviewed)
        accessions = pd.read_csv(proteome_path, sep='\t', dtype=str, usecols=['Entry'])['Entry'].tolist()
        print(f'mapping {len(accessions):,} accessions -> UniRef{level}...')
        t0 = time.time()
        pairs = _idmapping_run(accessions, f'UniRef{level}')
        out_df = pd.DataFrame(pairs, columns=['Entry', f'UniRef{level}'])
        out_df.to_csv(tmp, sep='\t', index=False)
        print(f'  saved {out} ({len(out_df):,} rows, {time.time() - t0:.1f}s)')

    return cached_pull(out, _fetch, refresh=refresh, max_days=max_days)


def load_uniref_mapping(proteome='UP000005640', level=50, reviewed=True,
                         refresh=False, root=None):
    """Load cached UniRef{level} mapping as a DataFrame. Pulls if missing or stale."""
    p = pull_uniref_mapping(proteome, level=level, reviewed=reviewed,
                             root=root, refresh=refresh)
    return pd.read_csv(p, sep='\t', dtype=str).fillna('')
