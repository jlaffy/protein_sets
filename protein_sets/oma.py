"""OMA orthologs: bulk UniProt mapping, group coverage, group metadata, HOG sizes,
and per-group MSAs.

Design mirrors uniprot.py / hpa.py / interpro.py: pull_X downloads + caches,
load_X reads cache. Plus convenience helpers for the common ortholog queries:
given a human UniProt, what's its OMA group, who are the orthologs, where's the MSA.

Cache layout under /ewsc/jlaffy/protein_sets/oma/:
    mapping/    oma-uniprot.txt.gz, oma-human-uniprot.txt
    coverage/   oma_coverage.csv          (uniprot, omaid, oma_group, oma_hog_id)
    metadata/   oma_group_metadata.csv    (description, n_members, species_codes)
    hog/        hog_sizes.csv
    msa/        oma_group_<gid>_aln.fasta (one per OMA group)
"""

import gzip
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm

from ._config import data_path

OMA_CACHE        = Path(data_path('oma'))
MAPPING_DIR      = OMA_CACHE / 'mapping'
COVERAGE_PATH    = OMA_CACHE / 'coverage' / 'oma_coverage.csv'
METADATA_PATH    = OMA_CACHE / 'metadata' / 'oma_group_metadata.csv'
HOG_PATH         = OMA_CACHE / 'hog'      / 'hog_sizes.csv'
MSA_DIR          = OMA_CACHE / 'msa'

BULK_MAPPING_URL = 'https://omabrowser.org/All/oma-uniprot.txt.gz'
OMA_BASE         = 'https://omabrowser.org'
MSA_URL          = OMA_BASE + '/oma/omagroup/{}/msa/'
OMA_GROUP_API    = OMA_BASE + '/api/group/{}/'
OMA_HOG_API      = OMA_BASE + '/api/hog/{}/members/'

DEFAULT_WORKERS  = 20


# --- HTTP session with retries (shared pattern across pulls) ---

def _session():
    s = requests.Session()
    r = Retry(total=4, backoff_factor=2,
              status_forcelist=[429, 500, 502, 503, 504, 522])
    s.mount('https://', HTTPAdapter(max_retries=r))
    return s


# --- Bulk OMA→UniProt mapping ---

def pull_uniprot_mapping(refresh=False):
    """Download the bulk OMA→UniProt mapping file (~74MB). Cached."""
    out = MAPPING_DIR / 'oma-uniprot.txt.gz'
    if out.exists() and not refresh:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    print(f'[oma] downloading bulk mapping from {BULK_MAPPING_URL} ...')
    t0 = time.time()
    r = _session().get(BULK_MAPPING_URL, stream=True, timeout=600)
    r.raise_for_status()
    with open(out, 'wb') as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)
    print(f'[oma] saved {out} ({out.stat().st_size / 1e6:.1f} MB, {time.time()-t0:.1f}s)')
    return out


def load_uniprot_mapping(refresh=False):
    """Load the OMA→UniProt mapping as a DataFrame (oma_id, uniprot, is_accession).

    The OMA bulk file contains multiple rows per OMA ID — entry names
    (e.g. RHBL2_HUMAN), primary accessions (e.g. Q9NX52), and secondary
    accessions. AFDB lookup requires accessions (no underscore), so we expose
    is_accession as a column and downstream code can prefer accession rows.
    """
    p = pull_uniprot_mapping(refresh=refresh)
    rows = []
    with gzip.open(p, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            oma, uni = line.rstrip('\n').split('\t', 1)
            rows.append((oma, uni, '_' not in uni))
    return pd.DataFrame(rows, columns=['oma_id', 'uniprot', 'is_accession'])


def best_accession_per_oma(refresh=False):
    """Return Series {oma_id: uniprot_accession} preferring accession over entry name.

    For each OMA ID, picks the first row where is_accession=True. Falls back to
    the first row (entry name) if no accession row exists for that OMA ID.
    """
    m = load_uniprot_mapping(refresh=refresh)
    pref = m.sort_values('is_accession', ascending=False)  # accessions first
    return pref.drop_duplicates('oma_id').set_index('oma_id')['uniprot']


# --- OMA group coverage (per human UniProt → group/hog assignment) ---

def load_coverage(refresh=False):
    """Load the per-uniprot OMA coverage table: uniprot, omaid, oma_group, oma_hog_id.

    `pull_coverage()` (re)builds it via OMA + UniProt APIs — slow, so we ship
    a pre-computed cache. Pass refresh=True to rebuild from scratch.
    """
    if refresh or not COVERAGE_PATH.exists():
        raise NotImplementedError(
            'coverage rebuild not ported yet — see the original '
            'conservation_scores/claude/scripts/get_oma_coverage.py'
        )
    return pd.read_csv(COVERAGE_PATH)


# --- OMA group metadata (description, members, species) ---

def pull_group_metadata(group_ids, workers=DEFAULT_WORKERS, append=False):
    """Fetch OMA group metadata for a list of group IDs and merge into the coverage cache."""
    cov = load_coverage()
    rows = {}
    s = _session()

    def _one(gid):
        try:
            r = s.get(OMA_GROUP_API.format(int(gid)), timeout=20)
            if r.status_code != 200:
                return {'oma_group': gid, 'description': None, 'n_members': None, 'species_codes': None}
            d = r.json()
            members = d.get('members', [])
            return {
                'oma_group':     gid,
                'description':   d.get('description'),
                'n_members':     len(members),
                'species_codes': ','.join(m['species']['code'] for m in members),
            }
        except Exception:
            return {'oma_group': gid, 'description': None, 'n_members': None, 'species_codes': None}

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_one, gid): gid for gid in group_ids}
        for f in tqdm(as_completed(futs), total=len(futs), desc='oma group metadata'):
            res = f.result()
            rows[res['oma_group']] = res

    new = pd.DataFrame(rows.values())
    out = cov.merge(new, on='oma_group', how='left')
    METADATA_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(METADATA_PATH, index=False)
    return METADATA_PATH


def load_group_metadata():
    """Load the per-protein OMA group metadata table (joined to coverage)."""
    return pd.read_csv(METADATA_PATH)


# --- HOG sizes ---

def pull_hog_sizes(hog_ids, workers=DEFAULT_WORKERS):
    """Fetch member counts for a list of HOG IDs."""
    s = _session()

    def _one(hid):
        try:
            r = s.get(OMA_HOG_API.format(hid), timeout=15)
            if r.status_code == 200:
                return hid, len(r.json().get('members', []))
        except Exception:
            pass
        return hid, None

    rows = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_one, h): h for h in hog_ids}
        for f in tqdm(as_completed(futs), total=len(futs), desc='oma hog sizes'):
            hid, n = f.result()
            rows[hid] = n
    df = pd.DataFrame({'oma_hog_id': list(rows), 'hog_n_members': list(rows.values())})
    HOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(HOG_PATH, index=False)
    return HOG_PATH


def load_hog_sizes():
    return pd.read_csv(HOG_PATH)


# --- Per-group MSAs ---

def pull_msa(group_id, refresh=False):
    """Fetch one OMA group's pre-computed MSA from the OMA Browser, save to cache."""
    out = MSA_DIR / f'oma_group_{int(group_id)}_aln.fasta'
    if out.exists() and not refresh:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    s = _session()
    html = s.get(MSA_URL.format(int(group_id)), timeout=30).text
    m = re.search(r'var url = "(/media/msa/[^"]+)"', html)
    if not m:
        raise ValueError(f'could not find MSA media URL for group {group_id}')
    fasta = s.get(OMA_BASE + m.group(1), timeout=30).text
    out.write_text(fasta)
    return out


def pull_msas(group_ids, workers=DEFAULT_WORKERS, max_rounds=5, refresh=False):
    """Drive pull_msa across many groups in parallel, with retries.

    Returns (n_done, list_of_failed_group_ids).
    """
    MSA_DIR.mkdir(parents=True, exist_ok=True)

    def _missing(gids):
        return [g for g in gids
                if not (MSA_DIR / f'oma_group_{int(g)}_aln.fasta').exists()]

    todo = group_ids if refresh else _missing(group_ids)
    for r in range(1, max_rounds + 1):
        if not todo:
            break
        n_workers = workers if r == 1 else max(4, workers // 4)
        print(f'[oma:pull_msas] round {r}: {len(todo)} groups ({n_workers} workers)')
        with ThreadPoolExecutor(max_workers=n_workers) as ex:
            futs = {ex.submit(pull_msa, g, refresh=refresh): g for g in todo}
            for f in tqdm(as_completed(futs), total=len(futs), desc='oma msa'):
                try:
                    f.result()
                except Exception:
                    pass
        todo = _missing(group_ids)

    n_done = len(group_ids) - len(todo)
    return n_done, todo


def msa_path(group_id):
    """Path to the cached MSA fasta for an OMA group ID. Does not pull."""
    return MSA_DIR / f'oma_group_{int(group_id)}_aln.fasta'


def load_msa(group_id, fetch=True):
    """Load an OMA group MSA as {ortholog_oma_id: aligned_sequence}.

    If `fetch=True` and the MSA isn't cached, pull it first.
    """
    p = msa_path(group_id)
    if not p.exists():
        if not fetch:
            raise FileNotFoundError(f'{p} not cached (fetch=False)')
        pull_msa(group_id)
    seqs, sid, parts = {}, None, []
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            if sid is not None:
                seqs[sid] = ''.join(parts)
            sid = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if sid is not None:
        seqs[sid] = ''.join(parts)
    return seqs


# --- Convenience API: starting from a UniProt accession ---

def oma_group_id(uniprot):
    """Return the OMA group ID for a human UniProt accession, or None if not in coverage."""
    cov = load_coverage()
    hit = cov[cov['uniprot'] == uniprot]
    if hit.empty:
        return None
    gid = hit.iloc[0].get('oma_group')
    if pd.isna(gid):
        return None
    return int(gid)


def oma_msa(uniprot, fetch=True):
    """Return the OMA group MSA for a human UniProt accession as {oma_id: aligned_seq}.

    Raises KeyError if the protein isn't assigned to an OMA group.
    """
    gid = oma_group_id(uniprot)
    if gid is None:
        raise KeyError(f'no OMA group assigned to {uniprot!r}')
    return load_msa(gid, fetch=fetch)


def oma_msa_file(uniprot):
    """Return the Path to the cached MSA file for a human UniProt's OMA group."""
    gid = oma_group_id(uniprot)
    if gid is None:
        raise KeyError(f'no OMA group assigned to {uniprot!r}')
    return msa_path(gid)


def orthologs(uniprot, with_seqs=False, fetch=True):
    """Return the list of ortholog OMA IDs in the same OMA group as `uniprot`.

    If with_seqs=True, returns a DataFrame with columns
    [oma_id, uniprot (mapped, may be empty), aligned_seq, sequence].
    """
    msa = oma_msa(uniprot, fetch=fetch)
    if not with_seqs:
        return list(msa.keys())

    # Join through the bulk OMA→UniProt mapping
    mapping = load_uniprot_mapping()
    first_uni = mapping.drop_duplicates('oma_id').set_index('oma_id')['uniprot']

    rows = []
    for oid, aligned in msa.items():
        rows.append({
            'oma_id':       oid,
            'uniprot':      first_uni.get(oid, ''),
            'aligned_seq':  aligned,
            'sequence':     aligned.replace('-', ''),
        })
    return pd.DataFrame(rows)
