"""Universe-path layout — single source of truth for the embedding cache tree.

Each protein_sets builder attaches a `provenance` dict to its output DataFrame
that includes three fields specifying where embeddings of that DataFrame should
live:

  source_category : top-level group (e.g. 'uniprot_features', 'interpro',
                    'transforms'). Determines the parent directory.
  subgroup        : feature/column name (e.g. 'signal_peptide', 'lipidation',
                    'motif', 'mutscan'). Determines the leaf directory or its
                    prefix.
  kind            : within-subgroup class (e.g. 'gpi', 'phospho', 'nls') or
                    None for single-class features. Appended to leaf when set.

The mode (fresh vs extract) is appended to the source category at embed time:

  embeddings/<model>/<source_category>_<mode>/<leaf>/

where leaf = subgroup if kind is None, else f"{subgroup}_{kind}".

Special case: synthetic-variant categories like 'transforms' don't take a mode
suffix because they are always produced fresh (no parent to slice from).

This module is the canonical translation between provenance and universe path.
Add new categories here, not in builder code.
"""

# Categories whose universes are always fresh (no extract path possible).
# Skip the mode suffix on these.
_FIXED_MODE_CATEGORIES = frozenset({'transforms'})


# Default subgroup names from a UniProt feature column. Used by uniprot_regions
# to auto-derive subgroup when a wrapper doesn't pass one explicitly.
SUBGROUP_FROM_UNIPROT_FEATURE = {
    'Signal peptide':       'signal_peptide',
    'Transit peptide':      'transit_peptide',
    'Transmembrane':        'transmembrane',
    'Intramembrane':        'intramembrane',
    'Topological domain':   'topological_domain',
    'Propeptide':           'propeptide',
    'Chain':                'chain',
    'Peptide':              'peptide_active',
    'Initiator methionine': 'initiator_methionine',
    'Lipidation':           'lipidation',
    'Disulfide bond':       'disulfide',
    'Glycosylation':        'glycosylation',
    'Modified residue':     'ptm',
    'Motif':                'motif',
}


# Default subgroup names from an InterPro member-db filter.
SUBGROUP_FROM_INTERPRO_DB = {
    'pfam':         'pfam',
    'superfamily':  'superfamily',
    'gene3d':       'gene3d',
    'smart':        'smart',
    'cdd':          'cdd',
    'prosite':      'prosite',
    'panther':      'panther',
    'pirsf':        'pirsf',
    'tigrfams':     'tigrfams',
    'hamap':        'hamap',
    'ncbifam':      'ncbifam',
}


def universe_path(provenance, mode='extract'):
    """Return the universe path corresponding to a provenance dict.

    >>> universe_path({'source_category': 'uniprot_features',
    ...                'subgroup': 'lipidation', 'kind': 'gpi'}, mode='extract')
    'uniprot_features_extract/lipidation_gpi'
    >>> universe_path({'source_category': 'uniprot_features',
    ...                'subgroup': 'signal_peptide', 'kind': None}, mode='fresh')
    'uniprot_features_fresh/signal_peptide'
    >>> universe_path({'source_category': 'transforms',
    ...                'subgroup': 'mutscan', 'kind': None}, mode='fresh')
    'transforms/mutscan'
    """
    cat = provenance['source_category']
    base = cat if cat in _FIXED_MODE_CATEGORIES else f'{cat}_{mode}'
    sub = provenance['subgroup']
    kind = provenance.get('kind')
    leaf = f'{sub}_{kind}' if kind else sub
    return f'{base}/{leaf}'


def attach_provenance(df, *, set, operation, source_category, subgroup,
                     kind=None, args=None, universe='human', **extra):
    """Attach extended provenance dict to `df.attrs['provenance']`.

    All builders should call this rather than constructing the dict by hand.
    The 3 layout-relevant fields (source_category, subgroup, kind) are required
    keyword args to make sure no builder forgets them.
    """
    df.attrs['provenance'] = {
        'universe':        universe,
        'set':             set,
        'operation':       operation,
        'args':            args or {},
        'source_category': source_category,
        'subgroup':        subgroup,
        'kind':            kind,
        **extra,
    }
    return df
