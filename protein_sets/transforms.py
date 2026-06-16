"""Synthetic-sequence transforms over real proteins.

Each function takes a query (UniProt ID, gene name, named set, list, or file
path), resolves it via the package's own resolve(), applies a sequence
transformation, and returns a DataFrame matching the region-table contract
used downstream by esm_embed:

    item_id, uniprot, source, source_id, start, end, length, sequence,
    + transform-specific columns

`df.attrs['provenance']` records the upstream universe + the operation, so
esm_embed appends the embedding step on top.

Current transforms:
    mutation_scan  — window-substitution scan
    to_fasta       — write a transform DataFrame to FASTA

Future transforms will follow the same contract: truncate, graft, chimera,
delete_region, alanine_scan, designed library injection, etc.
"""

import pandas as pd

from .lookup import resolve_with_seqs as _resolve_with_seqs


def mutation_scan(query, alphabet=('A', 'R', 'E', 'F'), window=5):
    """Window-substitution scan over each position of each WT protein.

    For every position i in [0, L), the residues in the symmetric window
    [max(0, i - half), min(L, i + half + 1)] (half = window // 2) are
    replaced with the target amino acid repeated. The window truncates at
    sequence ends.

    Per WT: 1 WT row + L * len(alphabet) mutant rows.

    item_id format:
        '{uid}_wt'                  (the unmodified WT)
        '{uid}_pos{pos:04d}_{aa}'   (mutants; pos zero-padded to 4 digits)
    """
    ids, seqs = _resolve_with_seqs(query)

    half = window // 2
    rows = []
    for uid in ids:
        seq = seqs[uid]
        L = len(seq)
        rows.append({
            'item_id':   f'{uid}_wt',
            'uniprot':   uid,
            'source':    'mutation_scan',
            'source_id': 'wt',
            'start':     1,
            'end':       L,
            'length':    L,
            'sequence':  seq,
            'kind':      'wt',
            'pos':       -1,
            'aa':        '',
            'win_start': 0,
            'win_end':   0,
            'win_len':   0,
        })
        for pos in range(L):
            ws = max(0, pos - half)
            we = min(L, pos + half + 1)
            wl = we - ws
            for aa in alphabet:
                mut = seq[:ws] + aa * wl + seq[we:]
                rows.append({
                    'item_id':   f'{uid}_pos{pos:04d}_{aa}',
                    'uniprot':   uid,
                    'source':    'mutation_scan',
                    'source_id': aa,
                    'start':     1,
                    'end':       L,
                    'length':    L,
                    'sequence':  mut,
                    'kind':      'mut',
                    'pos':       pos,
                    'aa':        aa,
                    'win_start': ws,
                    'win_end':   we,
                    'win_len':   wl,
                })

    out = pd.DataFrame(rows)
    from .layout import attach_provenance
    attach_provenance(
        out,
        set=query if isinstance(query, str) else 'custom',
        operation='mutation_scan',
        args={'alphabet': list(alphabet), 'window': window},
        source_category='transforms',
        subgroup='mutscan',
        kind=None,
    )
    return out


def to_fasta(df, path, line_width=0):
    """Write a transform DataFrame (item_id, sequence columns) to FASTA.

    line_width: 0 = single line per sequence; >0 = wrap to that width.
    Returns the number of records written.
    """
    n = 0
    with open(path, 'w') as f:
        for _, r in df.iterrows():
            f.write(f'>{r["item_id"]}\n')
            seq = r['sequence']
            if line_width > 0:
                for i in range(0, len(seq), line_width):
                    f.write(seq[i:i + line_width] + '\n')
            else:
                f.write(seq + '\n')
            n += 1
    return n
