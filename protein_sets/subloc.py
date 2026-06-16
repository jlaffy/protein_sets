"""Parse the UniProt `Subcellular location [CC]` free-text column.

UniProt encodes subcellular locations as one text blob per protein, with
isoform tags, topology annotations, ECO evidence codes, and free-text notes.
This module turns it into a long-format DataFrame: one row per
(uniprot, location), with topology + evidence kept alongside.

No "main vs additional" claim is made — locations appear in the order UniProt
wrote them, without semantic interpretation. Sublocations stay attached
("Nucleus, PML body" is one row). Isoform blocks are preserved with an
`isoform` column (canonical = empty string).
"""

import re
import pandas as pd

from .uniprot import load_proteome
from .lookup import resolve

# ECO code → human-readable tier. Order = priority when a single location
# carries multiple ECO codes (best tier wins).
ECO_TIERS = [
    ('ECO:0000269', 'experimental'),
    ('ECO:0000250', 'by_similarity'),
    ('ECO:0000305', 'inferred'),
    ('ECO:0000255', 'predicted'),
    ('ECO:0000303', 'non_traceable'),
    ('ECO:0007744', 'combinatorial'),
]
TIER_BY_CODE = dict(ECO_TIERS)
TIER_ORDER = [t for _, t in ECO_TIERS] + ['none']

_PREFIX = re.compile(r'^\s*SUBCELLULAR LOCATION:\s*', re.IGNORECASE)
_ISOFORM_HEADER = re.compile(r'\[Isoform [^\]]+\]:')
_EVIDENCE_BLOCK = re.compile(r'\{ECO:[^}]+\}')
_NOTE_BLOCK = re.compile(r'Note=.*$', re.DOTALL)


def _best_tier(eco_codes):
    if not eco_codes:
        return 'none'
    for code, tier in ECO_TIERS:
        if code in eco_codes:
            return tier
    return 'none'


def _parse_evidence(evidence_blocks):
    """Extract ECO codes + PubMed IDs from a list of `{ECO:...}` strings."""
    eco_codes = []
    pmids = []
    for block in evidence_blocks:
        inner = block[1:-1] if block.startswith('{') else block
        for clause in inner.split(','):
            clause = clause.strip()
            parts = clause.split('|')
            eco = parts[0].strip()
            if eco.startswith('ECO:'):
                eco_codes.append(eco)
            if len(parts) > 1:
                src = parts[1].strip()
                if src.startswith('PubMed:'):
                    pmids.append(src.split(':', 1)[1])
    eco_codes = list(dict.fromkeys(eco_codes))
    pmids = list(dict.fromkeys(pmids))
    return eco_codes, pmids


def _split_on_periods(text):
    """Split on periods that aren't inside `{...}` braces."""
    out, buf, depth = [], [], 0
    for ch in text:
        if ch == '{':
            depth += 1
            buf.append(ch)
        elif ch == '}':
            depth -= 1
            buf.append(ch)
        elif ch == '.' and depth == 0:
            out.append(''.join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        out.append(''.join(buf))
    return out


def _parse_one_location(entry):
    """Parse one location entry into (location, topology, ecos, pmids, tier)."""
    blocks = _EVIDENCE_BLOCK.findall(entry)
    eco_codes, pmids = _parse_evidence(blocks)
    text = _EVIDENCE_BLOCK.sub('', entry).strip().rstrip('. \t')
    if ';' in text:
        loc, _, topo = text.partition(';')
        loc, topo = loc.strip(), topo.strip().rstrip('. \t')
    else:
        loc, topo = text.strip(), ''
    return loc, topo, eco_codes, pmids, _best_tier(eco_codes)


def _split_isoform_blocks(text):
    """Yield (isoform_tag, body) pairs. Canonical block has tag = ''."""
    headers = list(_ISOFORM_HEADER.finditer(text))
    if not headers:
        yield ('', text)
        return
    if headers[0].start() > 0:
        head = text[:headers[0].start()].strip()
        if head:
            yield ('', head)
    for i, m in enumerate(headers):
        tag = m.group(0)[1:-2].strip()  # "[Isoform 2]:" -> "Isoform 2"
        start = m.end()
        end = headers[i + 1].start() if i + 1 < len(headers) else len(text)
        yield (tag, text[start:end].strip())


def _parse_subloc_field(uniprot, gene, text):
    """One free-text field → list of row dicts."""
    if not text or not text.strip():
        return []
    text = _PREFIX.sub('', text)
    rows = []
    for isoform_tag, body in _split_isoform_blocks(text):
        body = _NOTE_BLOCK.sub('', body).strip()
        if not body:
            continue
        for entry in _split_on_periods(body):
            entry = entry.strip()
            if not entry:
                continue
            loc, topo, ecos, pmids, tier = _parse_one_location(entry)
            if not loc:
                continue
            rows.append({
                'gene': gene,
                'uniprot': uniprot,
                'location': loc,
                'eco_codes': ecos,
                'evidence_tier': tier,
                'pubmed_ids': pmids,
                'topology': topo,
                'isoform': isoform_tag,
            })
    return rows


_LONG_COLS = ['gene', 'uniprot', 'location', 'eco_codes', 'evidence_tier',
              'pubmed_ids', 'topology', 'isoform']


def uniprot_subloc(query=None, refresh=False):
    """Parse UniProt Subcellular location [CC] into a long DataFrame.

    One row per (uniprot, location). The evidence on each row applies to that
    specific location only. Columns (in order):
        gene, uniprot, location, eco_codes, evidence_tier, pubmed_ids,
        topology, isoform.

    `query`: None (all reviewed human), a list of accessions, a named-set
        string, or a DataFrame with an `Entry`/`uniprot` column.
    """
    df = load_proteome(refresh=refresh)
    if query is not None:
        ids = set(resolve(query))
        df = df[df['Entry'].isin(ids)]
    col = 'Subcellular location [CC]'
    gene_col = 'Gene Names (primary)'
    rows = []
    for entry, gene, txt in zip(df['Entry'].values,
                                df[gene_col].values,
                                df[col].values):
        rows.extend(_parse_subloc_field(entry, gene, txt))
    return pd.DataFrame(rows, columns=_LONG_COLS)


def uniprot_subloc_wide(query=None, refresh=False):
    """Pivot `uniprot_subloc` to one row per protein with list-valued cols.

    Columns: gene, uniprot, locations, eco_codes, evidence_tiers,
    pubmed_ids, topologies, isoforms.
    """
    long = uniprot_subloc(query=query, refresh=refresh)
    wide_cols = ['gene', 'uniprot', 'locations', 'eco_codes', 'evidence_tiers',
                 'pubmed_ids', 'topologies', 'isoforms']
    if long.empty:
        return pd.DataFrame(columns=wide_cols)
    grouped = long.groupby(['uniprot', 'gene'], sort=False).agg({
        'location': list,
        'eco_codes': list,
        'evidence_tier': list,
        'pubmed_ids': list,
        'topology': list,
        'isoform': list,
    }).reset_index()
    grouped = grouped[['gene', 'uniprot', 'location', 'eco_codes',
                       'evidence_tier', 'pubmed_ids', 'topology', 'isoform']]
    grouped.columns = wide_cols
    return grouped
