---
name: protein-sets
description: "Canonical lens for any question about subsets, annotations, or categories of proteins (default universe: reviewed_human). Five layers: (1) sequence sources — proteomes, future ORF universes; (2) named sets and set algebra — secretome, transmem, intersections, dedupe; (3) regions within sequences — domains, motifs, sites, recurring annotations from InterPro and UniProt features; (4) transforms — synthetic-sequence experiments built on top: mutation scans, future graft / chimera / truncate / delete-region; (5) studies — bespoke per-project datasets registered under /ewsc/jlaffy/protein_sets/studies/ that emit region or protein rows (e.g. schaeffer22 cleavage-site 15mers tagged by cell type). All five share the same `(item_id, uniprot, sequence, ...)` row contract and feed esm_embed / protein-structures / protein-clusters downstream. Study names are also resolver-addressable like named sets. When the user references a group of proteins or an annotation without specifying source, route the conversation through this skill first — check NAMED_SETS, region builders, studies.available(), and load_proteome() columns. If a category isn't formalized yet, surface the underlying data and propose a definition before answering ad-hoc."
---

You are helping the user work with protein sequences — selecting, filtering, subsetting, and transforming them. The code lives in the `protein_sets` package at `/home/unix/jlaffy/protein_sets` (`pip install -e ~/protein_sets`). All edits — new sources, set definitions, region types, transforms — go there.

## Universe-path layout (single source of truth)

Every builder attaches an extended `provenance` dict to its DataFrame. The
three layout-relevant fields determine where embeddings of that DataFrame land:

```
df.attrs['provenance'] = {
    ...,
    'source_category': 'uniprot_features' | 'interpro' | 'transforms' | ...,
    'subgroup':        'lipidation' | 'motif' | 'signal_peptide' | ...,
    'kind':            'gpi' | 'nls' | None,    # None for single-class features
}
```

The helper `protein_sets.layout.universe_path(prov, mode)` (also exported as
`protein_sets.universe_path`) translates this into a path:

```
embeddings/<model>/<source_category>_<mode>/<leaf>/
```
where `leaf = subgroup if kind is None else f'{subgroup}_{kind}'`.

Special case: `source_category='transforms'` skips the `_<mode>` suffix
(synthetic variants are always fresh — there's no parent to slice from).

### When you write a new builder

1. Pick `source_category`: `uniprot_features` (UniProt TSV columns), `interpro`
   (InterPro xrefs), `transforms` (synthetic variants), or a new top-level if
   the data origin is genuinely new.
2. Pick `subgroup`: usually the column name (`signal_peptide`, `lipidation`,
   `motif`) for UniProt features, or the member db (`pfam`) for InterPro.
   Defaults are in `SUBGROUP_FROM_UNIPROT_FEATURE` and
   `SUBGROUP_FROM_INTERPRO_DB`.
3. Pick `kind`:
   - `None` if the feature is single-class (e.g. `signal_peptide`,
     `disulfide` — leaf becomes just the subgroup).
   - A class string if the feature has multiple distinct sub-classes (e.g.
     `lipidation` → kind ∈ {`gpi`, `palmitoyl`, ...}; `motif` → kind ∈ {`nls`,
     `nes`, ...}). Leaf becomes `subgroup_kind`.
4. Call `attach_provenance(df, set=..., operation=..., args=...,
   source_category=..., subgroup=..., kind=...)` instead of building the dict
   by hand. The required-keyword signature stops you forgetting fields.

### When the kind decision is ambiguous

Use the discover-and-propose protocol below: surface the column's distinct
`/note=` values, propose a sensible split, confirm with the user before
writing the classifier + wrappers.

## Default behavior on ambiguous "subset" / "annotation" questions

When the user references a group of proteins or an annotation without specifying
source (e.g. "GPI proteins", "the kinases", "things with a propeptide"), follow
this protocol *before* writing analysis code:

1. **Default universe**: assume `reviewed_human`. There is no other universe
   currently in use (no unreviewed, no cross-species).
2. **Check what already exists**, in this order:
   a. `NAMED_SETS.keys()` and `list_dynamic()` — registered named sets.
   b. Region builders in `regions.py` (`pfam_regions`, `signal_peptide_regions`,
      `propeptide_regions`, `uniprot_regions`, etc.).
   c. Transforms in `transforms.py`.
3. **If nothing fits, discover before answering**:
   a. Surface candidate columns in `load_proteome()` (e.g. `Keywords`,
      `Lipidation`, `Subcellular location [CC]`).
   b. List the distinct categories/values present in those columns.
   c. Propose a concrete definition (named set or region builder) and confirm
      with the user before either (i) registering it in `NAMED_SETS` /
      `regions.py`, or (ii) answering ad-hoc with a one-shot filter.
4. Prefer formalizing over ad-hoc filtering when the category looks reusable.
   Ad-hoc is fine only for genuinely one-off exploratory questions.

**Anti-pattern**: do not write `/tmp` regex scripts on `load_proteome()` to
answer recurring annotation questions. If the analysis would benefit from being
reusable across sessions, propose adding it as a named set or region builder
first, then use that.

## Sites vs bonds vs spans — when to prompt

Some features have multiple table shapes available (e.g. disulfides have both
`disulfide_sites` and `disulfide_bonds`). Pick by intent:

| user wants | shape | function family |
|---|---|---|
| embeddings / per-position analysis | site rows (1 residue each) | `*_sites` |
| counts, loop lengths, pair lookups, "which X is paired with Y" | bond/pair rows (2 positions per row) | `*_bonds` |
| sliced subsequences | span rows | `*_regions` |

**Prompt when ambiguous.** If the user asks for "disulfide bonds" but the
surrounding context implies embeddings (e.g. "generate embeddings for the
disulfides", "extract from full-length"), surface the choice — they may have
meant sites. Likewise if they say "disulfide sites" but ask about loop lengths
or pair structure (which sites tables don't carry), redirect to bonds. The
function names are precise; the user's word may not be.

## Vocabulary

Four kinds of things this skill produces, all sharing the same row contract `(item_id, uniprot, sequence, [start, end, ...])` plus `df.attrs['provenance']`. Use the precise term when discussing — they don't mean the same thing.

- **Universe** — a persisted `(id, sequence)` population from an authoritative source. Universes come *only* from UniProt populations (reviewed_human, future unreviewed_human) or future ORF libraries. Each universe has its own ID space. Default for almost everything: `reviewed_human`.
- **Named set** — a logically-defined filter inside a universe. Examples: `secretome`, `transmem`, `signal_peptide_strict`, `surfaceome`. Same sequences as the parent universe — just a subset of items. Composes via intersection / union / `dedupe()`.
- **Region** — annotation-derived boundaries *within* sequences (sub-strings). Examples: `pfam_regions(...)`, `signal_peptide_regions(...)`. When materialized, regions form a *new universe* because the sequences differ from the parent (they're sub-strings).
- **Transform** — a sequence-level operation that produces synthetic sequences from a query (any of the above). Examples: `mutation_scan(...)`. When materialized, transforms also form a new universe — sequences are modified.
- **Study** — a registered bespoke dataset under `/ewsc/jlaffy/protein_sets/studies/<name>/`. Each study has its own `raw/`, `parse.py` (re-runnable, not imported), `parsed/rows.parquet`, and `metadata.yaml`. Output may be region rows (sites, windowed peptides — e.g. `schaeffer22`) or, in future, protein rows (full-sequence external libraries). The package ships the contract + loader; per-study `parse.py` is bespoke. Resolver-addressable (`resolve("schaeffer22")` → parent UniProts).

**Rule that unifies everything**: same sequences as parent → it's a *named set* within the parent universe. Different sequences → it's a *new universe* (derived via slicing, transform, or external study). esm-embed records each universe's lineage in its `universe.yaml` (parent + operations chain).

Naming convention disambiguates layer at a glance:
```
signal_peptide()         → named set (whole proteins)
signal_peptide_regions() → region builder (sub-sequences)
mutation_scan()          → transform (synthetic sequences)
```

Three phrasings, three meanings:
- "the secretome **set**" → named set (registered, in `NAMED_SETS`)
- "a **set** that has X and Y" → ad-hoc filter described inline (not registered)
- "the **set**" / "this **set**" → generic, whatever's currently in hand

## Layered model

Five conceptual layers in the package, all sharing the row contract:

1. **Sequence sources** — proteomes (UniProt reviewed human; future: unreviewed, ORF universes). The "what's real" layer.
2. **Sets** — named subsets and set algebra (`secretome`, `transmem`, `dedupe(...)`, intersections, unions).
3. **Regions** — annotation-derived boundaries within sequences (Pfam domains, signal peptides, transmembrane spans, sites). Recurring annotation types reusable across many proteins.
4. **Transforms** — synthetic-sequence experiments built on top of the above: mutation scans (and in the future: graft, chimera, truncate, delete-region). Take a query + a transformation regime, return a DataFrame of synthetic `(item_id, sequence)` rows that feed esm_embed exactly like real sequences do.
5. **Studies** — bespoke per-project datasets (external experiments, paper supplementary tables) registered under `/ewsc/jlaffy/protein_sets/studies/<name>/`. Output is region rows (e.g. `schaeffer22` cut-site 15mers per cell type) or, for future "new-sequence" studies, protein rows. Per-study `parse.py` is bespoke; once parsed, rows flow through the same downstream pipes as any other layer.

The system is designed to grow — new sources, named sets, region types, transforms, and studies slot in following the existing patterns.

## Outputs

Three concrete output forms, regardless of which layer you're querying:

- **DataFrame** — the canonical output of every query. Carries `df.attrs['provenance']`.
- **IDs** — just the `Entry` column when accessions are all you need. `resolve(query)` returns a list directly.
- **FASTA** — `to_fasta(df, path)` writes records for handing off to non-Python tools.

## Downstream tools

This skill defines *inputs*. Downstream tools consume them and produce their own kinds of output:

```
protein-sets         (defines: universes, sets, regions, transforms)
    │
    ▼  DataFrame or query (resolved via protein_sets.resolve)
esm-embed            (output Y: embeddings)
protein-structures   (output Y: structures, 3di, contacts)
<future tool>        (output Y: predictions, alignments, ...)
```

Every downstream tool should accept either a DataFrame from this package or a query string `resolve(...)` understands. Same row contract everywhere.

## Import

```python
from protein_sets import (
    # proteome
    load_proteome, pull_proteome, PROTEOME_CACHE,
    # UniRef cluster IDs
    load_uniref_mapping, pull_uniref_mapping, fetch_cluster_ids, UNIREF_CACHE,
    # HPA
    load_hpa_secretome, hpa_secretome_uniprots, HPA_CACHE,
    # CSPA
    load_cspa_human, load_cspa_glycopeptides, cspa_uniprots, CSPA_CACHE,
    # humantfs (Lambert 2018 human transcription factors)
    load_humantfs, humantfs_uniprots, HUMANTFS_CACHE,
    # InterPro
    load_interpro, load_entry_types, INTERPRO_CACHE,
    # named sets
    NAMED_SETS, resolve_dynamic, list_dynamic,
    reviewed_human, signal_peptide, transmem, secretome, surfaceome,
    dedupe,
    # regions
    interpro_domain_regions, pfam_regions,
    uniprot_regions, uniprot_regions_remove, signal_peptide_regions,
    # generic external-table → region builder (used by studies)
    from_external_table,
    # transforms (synthetic sequences derived from real proteins)
    mutation_scan, to_fasta,
)
# bespoke studies (registered under /ewsc/jlaffy/protein_sets/studies/)
from protein_sets import studies
studies.available()                   # → ['schaeffer22', ...]
studies.info('schaeffer22')           # → parsed metadata.yaml
studies.load('schaeffer22')           # → rows on the standard contract
studies.load('schaeffer22', cell_type='HEK293')  # filter by meta col
```

## Data sources

| source | cached at | refresh | what it provides |
|---|---|---|---|
| UniProt reviewed human proteome | `/ewsc/jlaffy/protein_sets/uniprot/proteomes/UP000005640_reviewed.tsv` | auto, 30-day TTL | 20,416 proteins: accession, gene name, sequence, signal peptide, transmembrane, intramembrane, lipidation, PTMs, topology, subcellular location, GO, keywords, Pfam xrefs, and more |
| UniRef cluster IDs (50 / 90 / 100) | `/ewsc/jlaffy/protein_sets/uniref/UP000005640_reviewed_uniref{50,90,100}.tsv` | auto, 30-day TTL | per-accession UniRef cluster IDs for every reviewed human protein. Pulled via the UniProt ID Mapping API as 2-column sidecars and merged into `load_proteome()` as the `UniRef50/90/100` columns. Powers `dedupe()`. **Phase 2 (planned, not built)**: full cross-species cluster catalog (member lists, LCA taxonomy, cluster representatives) under the same `uniref/` cache root, no sequences. |
| HPA Predicted Secreted Proteins | `/ewsc/jlaffy/protein_sets/hpa/predicted_secreted.tsv` | auto, 30-day TTL | curated secreted protein list from Human Protein Atlas |
| CSPA Cell Surface Protein Atlas | `/ewsc/jlaffy/protein_sets/cspa/S2_File.xlsx` | auto, 30-day TTL | Bausch-Fluck et al. (2015) PLoS One e0121314 — MS-derived human surfaceome (1,492 entries, fed into the `surfaceome*` named sets after intersecting with reviewed_human). The file also ships ~13,942 N-glycopeptides accessible via `load_cspa_glycopeptides()` if needed. |
| humantfs (Human Transcription Factors) | `/ewsc/jlaffy/protein_sets/humantfs/DatabaseExtract_v_1.01.csv` | auto, 30-day TTL | Lambert et al. (2018) Cell 172:650 — the field-standard curated human TF census (humantfs.ccbr.utoronto.ca): 2,765 assessed candidates, 1,639 curated bona fide TFs (`Is TF?` == Yes), with DBD family, binding mode, and motif status. Keyed by Ensembl/HGNC symbol (no UniProt), so `humantfs_uniprots()` maps symbols onto reviewed_human (~1,632/1,639 map); feeds the `transcription_factors` named set. Cleaner than any single UniProt keyword or GO term. |
| InterPro domain boundaries | `/ewsc/jlaffy/protein_sets/interpro/human_protein2ipr.tsv` + `interpro_entry_types.tsv` | manual | per-protein domain start/end positions from all member DBs (Pfam, SMART, SUPERFAMILY, Gene3D, CDD, etc.) with InterPro entry type classification |

Additional sources can be added over time following the same pattern: a pull function, a cache directory under `/ewsc/jlaffy/protein_sets/`, and a load function.

### The UniProt TSV

This is the single source of truth for protein-level annotation. Load it as a DataFrame:

```python
df = load_proteome()  # 20,416 rows
```

Actual columns (from `df.columns`):

| column | content |
|---|---|
| `Entry` | UniProt accession (e.g. P01308) |
| `Gene Names (primary)` | primary gene symbol |
| `Sequence` | full amino acid sequence |
| `Length` | sequence length |
| `Signal peptide` | signal peptide feature annotations (structured text with positions + evidence) |
| `Transmembrane` | transmembrane region annotations |
| `Intramembrane` | intramembrane region annotations |
| `Lipidation` | lipid modifications including GPI anchors (e.g. `LIPID 228; /note="GPI-anchor amidated serine"`) |
| `Topological domain` | inside/outside/membrane topology |
| `Subcellular location [CC]` | free-text subcellular location |
| `Chain` | mature chain after processing |
| `Propeptide` | propeptide regions |
| `Initiator methionine` | initiator methionine |
| `Peptide` | active peptides |
| `Disulfide bond` | disulfide bond positions |
| `Glycosylation` | glycosylation sites |
| `Modified residue` | other PTM sites |
| `Function [CC]` | function description |
| `Involvement in disease` | disease associations |
| `Gene Ontology (GO)` | all GO terms |
| `Gene Ontology (cellular component)` | GO CC terms |
| `Gene Ontology (molecular function)` | GO MF terms |
| `Gene Ontology (biological process)` | GO BP terms |
| `Keywords` | UniProt keywords (semicolon-separated, e.g. `GPI-anchor;Membrane;Signal`) |
| `Pfam` | Pfam cross-references |
| `InterPro` | InterPro cross-references |
| `AlphaFoldDB` | AlphaFold cross-references |
| `UniRef50` / `UniRef90` / `UniRef100` | UniRef cluster IDs (e.g. `UniRef50_P01308`) — merged in by `load_proteome(with_uniref=True)` (default) from sidecar TSVs in `/ewsc/jlaffy/protein_sets/uniref/`. Each protein gets one cluster ID per level; cluster IDs are unique grouping keys for `dedupe()` |

All columns are strings. Empty string means no annotation. Feature columns (Signal peptide, Transmembrane, Lipidation, etc.) contain structured text with positions and evidence codes — these would need parsing to extract start/end coordinates.

## Named protein sets

| name | definition | ~count |
|---|---|---|
| `reviewed_human` | all reviewed human UniProt entries | 20,416 |
| `signal_peptide` | has signal peptide annotation | 3,614 |
| `signal_peptide_strict` | signal peptide with experimental evidence (ECO:0000269) | 738 |
| `transmem` | has transmembrane annotation | 5,228 |
| `transmem_strict` | transmembrane with experimental evidence | 285 |
| `sp_no_tm` | signal peptide yes, transmembrane no | 2,123 |
| `secretome` | HPA Predicted Secreted Proteins intersected with reviewed human | 1,870 |
| `surfaceome` | CSPA human surface proteins (categories 1+2: high + putative), intersected with reviewed human | ~1,241 |
| `surfaceome_strict` | CSPA category 1 (high confidence) only, intersected with reviewed human | ~996 |
| `surfaceome_all` | all CSPA human entries incl. category 3 (unspecific), intersected with reviewed human | ~1,444 |
| `transcription_factors` | Lambert et al. 2018 curated human TFs (`Is TF?` == Yes), mapped to reviewed human | ~1,631 |

### Resolving a named set

```python
# Get the full DataFrame (all UniProt columns preserved)
df = transmem()              # 5,228 rows, all columns
df = signal_peptide()        # 3,614 rows
df = secretome()             # 1,870 rows
df = surfaceome()            # ~1,241 rows (CSPA categories 1+2 ∩ reviewed_human)

# Or via the registry (useful when the set name is a variable)
ids, seqs, df = resolve_dynamic('secretome')
# ids:  list of UniProt accessions
# seqs: dict {accession: sequence}
# df:   full DataFrame

# List available set names
list_dynamic()  # ['reviewed_human', 'secretome', 'signal_peptide', ...]
```

### Redundancy-reduced sets (UniRef dedup)

Use `dedupe()` to collapse a set down to one representative per UniRef cluster — useful for ML train/test splits, MSA construction, or any analysis that should not be biased by paralogous duplicates.

```python
from protein_sets import dedupe

# Drop near-duplicates within a set at 50% identity
d = dedupe('signal_peptide_strict', level=50)   # 738 → 711
d = dedupe('secretome', level=50)               # 1,870 → 1,744
d = dedupe('transmem', level=50)                # 5,228 → 4,863

# Other identity levels
d = dedupe('secretome', level=90)               # 1,870 → 1,812 (90% identity)
d = dedupe('secretome', level=100)              # 1,870 → 1,859 (exact dedup)

# Pick the longest sequence per cluster instead of first
d = dedupe('transmem', level=50, keep='longest')

# Works on a DataFrame too (must already carry the UniRef{level} column)
df = secretome()
d = dedupe(df, level=50)
```

A small number of entries (~3 in the reviewed human proteome) have no UniRef50/90 assignment in the current UniProt release — `dedupe()` keeps these as singletons rather than dropping them.

The cluster IDs are pulled via the UniProt ID Mapping API by `pull_uniref_mapping()`, cached as 2-column sidecars under `/ewsc/jlaffy/protein_sets/uniref/`, and merged into `load_proteome()` automatically. Currently human-proteome-anchored: a forward map per accession, no full cluster membership. Phase 2 (planned) will add the cross-species cluster catalog (member lists, LCA taxonomy) under the same cache root.

### Adding a new named set

Add a function to `~/protein_sets/protein_sets/definitions.py` and register it in `NAMED_SETS`:

```python
def gpi_anchored(refresh=False):
    """Reviewed human proteins with a GPI-anchor annotation."""
    df = reviewed_human(refresh=refresh)
    mask = df['Keywords'].str.contains('GPI-anchor', case=False, na=False)
    return df[mask].reset_index(drop=True)

NAMED_SETS['gpi_anchored'] = (gpi_anchored, {})
```

The pattern is always: load the proteome (or another source), filter on annotation columns, return the filtered DataFrame. Downstream consumers (esm-embed, foldseek, analysis scripts) can then use the set by name.

New sets can combine multiple annotations (e.g. an `all_membrane` set that unions transmembrane, intramembrane, and GPI-anchored proteins) or bring in new sources.

## Regions

Regions are subsequences or sites within proteins, defined by `(uniprot, start, end)` coordinates. Currently sourced from InterPro via `interpro_domain_regions()`.

### InterPro regions

`interpro_domain_regions()` extracts regions from InterPro annotations. It supports filtering by `db` (source database) and `entry_type` (InterPro classification). The two filters are independent — you can use either, both, or neither.

**Entry types** — Domain and Homologous_superfamily are both domain-level regions in practice; the others are sites, repeats, or family-level assignments:

| entry_type | rows | what it is |
|---|---|---|
| `Domain` | 284,129 | structural/functional domains (Pfam, SMART, CDD, etc.) |
| `Homologous_superfamily` | 137,004 | structural domain classifications (SUPERFAMILY, Gene3D) — domain-level in practice |
| `Family` | 129,850 | whole-protein family assignments |
| `Repeat` | 53,287 | repeated structural units |
| `Conserved_site` | 14,251 | short conserved sequence motifs |
| `Binding_site` | 3,381 | ligand/substrate binding positions |
| `Active_site` | 2,864 | catalytic residue positions |
| `PTM` | 1,318 | post-translational modification sites |

**Database filters**: `pfam`, `smart`, `cdd`, `gene3d`, `superfamily`, `prosite`, `panther`, `pirsf`, `tigrfams`, `hamap`, `ncbifam`

```python
# Pfam domains for the secretome
regions = pfam_regions('secretome')

# All domain-level regions (Domain + Homologous_superfamily entry types use different DBs)
regions = interpro_domain_regions('secretome', entry_type='Domain')

# Binding sites across transmem proteins
regions = interpro_domain_regions('transmem', entry_type='Binding_site')
```

Region table columns: `item_id`, `uniprot`, `source`, `source_id`, `interpro_id`, `entry_type`, `start`, `end`, `length`, `sequence`.

Coordinates are 1-indexed inclusive (UniProt convention). `item_id` format: `{uniprot}_{member_db_acc}_{start}_{end}`.

### UniProt feature regions

The UniProt TSV feature columns contain positional annotations in a uniform `TYPE start..end; /qualifiers...` format (multiple features semicolon-separated). Two functions extract them, sharing the same output contract as `interpro_domain_regions()`:

```python
uniprot_regions(set, feature='Signal peptide')          # the annotated region
uniprot_regions_remove(set, feature='Signal peptide')   # complement (one row per parent, item_id = parent uniprot)
```

`feature` accepts the column name (`'Signal peptide'`), the keyword (`'SIGNAL'`), or the short form (`'sp'`). `_remove` requires the complement to be a **single contiguous span per protein** — unambiguous by design. For N-terminal features (SP, PROPEP) this is `(end+1, length)`. For interior features (TM helices, etc.) the complement splits and the function **raises** with examples; you must add an explicit `remove_policy` to that feature's `UNIPROT_FEATURE_KEYWORDS` entry or use a different operation. Item ID is the parent UniProt accession — meant for populating `<embeddings>/<model>/human_<short>_removed/` universes.

Per-feature config in `UNIPROT_FEATURE_KEYWORDS` carries: `keyword` (UniProt's literal token), optional `short` (colloquial name like `sp`/`tm`), and boundary policy `fill_unknown_start` / `fill_unknown_end` (how to handle `?` coordinates). Currently only `Signal peptide` has full policy filled in; others grow organically.

Columns with positional region data: `Signal peptide`, `Transmembrane`, `Intramembrane`, `Topological domain`, `Chain`, `Propeptide`, `Peptide`, `Lipidation`, `Disulfide bond`, `Glycosylation`, `Modified residue`.

### Provenance

Every region function attaches a provenance dict to `df.attrs['provenance']`:

```python
{universe: 'human', set: '<named_set>', operation: '<function_name>', args: {...}}
```

Downstream caches (e.g. `esm_embed.embed_regions`) inherit this and append their own operation, persisting the chain to `<universe>/universe.yaml`. Any artifact derived from the same protein-set call has identical upstream provenance — fully traceable.

### Choosing the right region source

Some region types clearly come from one source. Others overlap between InterPro and UniProt — when ambiguous, ask the user.

| region type | source | why |
|---|---|---|
| structural domains (Pfam, SMART, etc.) | InterPro | this is what InterPro is for |
| homologous superfamilies | InterPro | SUPERFAMILY/Gene3D structural classifications |
| repeats | InterPro | repeat unit boundaries |
| signal peptide / transit peptide spans | UniProt features | processing features, not in InterPro |
| transmembrane helix positions | UniProt features | per-helix spans (InterPro has the domain *containing* TM helices, different thing) |
| chain / propeptide / peptide boundaries | UniProt features | maturation/processing |
| topological domains (extra/intracellular) | UniProt features | topology segments |
| active sites | **ask** | both InterPro (ProSite) and UniProt annotate these, different coverage |
| binding sites | **ask** | both sources, different granularity |
| PTMs | **ask** | InterPro has some; UniProt has Modified residue, Glycosylation, Lipidation, Disulfide bond — very different annotation types |
| conserved sites | **ask** | InterPro has them, UniProt may have evidence codes |

## Transforms (synthetic-sequence experiments)

Transforms take real proteins (resolved via the same query language used everywhere else in this skill) and apply a sequence-level transformation to produce a DataFrame of synthetic `(item_id, uniprot, sequence, ...)` rows. Output matches the region-table contract, so `esm_embed.embed_regions(...)` consumes them identically to real-protein regions.

Today: `mutation_scan` only. Future transforms (truncate, graft, chimera, delete-region, alanine-scan, designed-library injection) follow the same contract.

### `mutation_scan(query, alphabet=('A','R','E','F'), window=5)`

Window-substitution scan. For each position `i ∈ [0, L)` in each WT, replaces the symmetric window `[max(0, i - half), min(L, i + half + 1)]` (half = `window // 2`) with the target amino acid repeated. Window truncates at sequence ends.

Per WT, returns `1` WT row + `L * len(alphabet)` mutant rows.

`item_id` format:
- `'{uid}_wt'` — the unmodified WT
- `'{uid}_pos{pos:04d}_{aa}'` — mutants. Position is zero-padded to 4 digits so lexicographic sort = numeric sort.

```python
from protein_sets import mutation_scan, to_fasta

df = mutation_scan('P01308')                              # insulin, 4 AAs × 110 positions = 441 rows
df = mutation_scan('secretome', alphabet=('A','R'))       # whole named set
to_fasta(df, '/tmp/mutscan.fa')
```

Plus all 14 columns of metadata: `kind` (`wt`|`mut`), `pos`, `aa`, `win_start`, `win_end`, `win_len`.

### Wiring into esm_embed

Two paths, both supported:

**A. DataFrame → `embed_regions`** (recommended — preserves provenance via `df.attrs['provenance']`):

```python
from esm_embed.api import embed_regions, summarize
embed_regions(df, model='esmc_600m_1152',
              universe='mutscan', set_name='P01308',
              source='fresh', yes=True)
summarize('esmc_600m_1152', 'mutscan', 'P01308', form='cls')
summarize('esmc_600m_1152', 'mutscan', 'P01308', form='mean')
```

**B. FASTA → `embed_fasta`** (when you only have a file):

```python
from esm_embed.api import embed_fasta
embed_fasta('/tmp/mutscan.fa', universe='mutscan', set_name='P01308',
            model='esmc_600m_1152', yes=True)
```

### Where the embeddings land (esm-embed config)

The actual embedding step happens in `esm_embed`, which writes its output to `$ESM_EMBED_CACHE_ROOT/<model>/<universe>/`. Convention:
- `/ewsc/jlaffy/embeddings/` — canonical real-proteome universes (`human`, `human_sp_removed`, …). The default.
- `/ewsc/jlaffy/embeddings_variants/` — synthetic-variant universes produced by transforms (`mutscan`, future `mutscan_chimera`, etc.).

Keep them separate so the canonical embeddings directory isn't diluted with derived/synthetic data:

```bash
ESM_EMBED_CACHE_ROOT=/ewsc/jlaffy/embeddings_variants python my_scan.py
```

### Adding a new transform

1. Add a function to `~/protein_sets/protein_sets/transforms.py`. Same input contract as the existing functions (a query), same row contract output, set `df.attrs['provenance']`.
2. Re-export from `protein_sets/__init__.py`.
3. Document the transform in this section — its inputs, the columns it adds, and the `item_id` naming scheme it uses.
4. (Optional) Add a smoke test that exercises one named-set query end-to-end.

The naming convention for `item_id` should encode the transform's parameters so the FASTA / manifest / per_residue files are self-describing.

## Studies (bespoke external datasets)

Studies are how external/per-project data joins the row contract. The package ships the contract + the loader + a generic windowing helper; each study lives outside the package, in a known on-disk layout.

### On-disk layout

```
/ewsc/jlaffy/protein_sets/studies/<name>/
  raw/                     # original downloads, untouched
  parsed/
    rows.parquet           # REQUIRED — standardized row contract
  metadata.yaml            # REQUIRED — see below
  parse.py                 # raw → parsed. Re-runnable. NOT imported at load time.
  docs/  scripts/          # optional, for completeness
```

Only `rows.parquet` + `metadata.yaml` are loader-required; everything else exists for reproducibility / independence.

### Loader API

```python
from protein_sets import studies
studies.available()                          # ['schaeffer22', ...]
studies.info('schaeffer22')                  # parsed metadata.yaml as dict
studies.load('schaeffer22')                  # rows.parquet + provenance attached
studies.load('schaeffer22', cell_type='HEK293')  # filter on any meta column
studies.path('schaeffer22')                  # on-disk directory
```

`studies.load(...)` attaches `df.attrs['provenance']` using the `provenance` block in `metadata.yaml` (defaults: `source_category='studies'`, `subgroup=<name>`).

### Row contract for studies

Region-style studies (e.g. `schaeffer22` — 15mers around cut sites):

| col | meaning |
|---|---|
| `item_id` | `{study}__{uniprot}__{position}` plus a `__{meta}` slug when needed for uniqueness |
| `uniprot` | parent accession (keyed against `parent_universe`, usually `reviewed_human`) |
| `sequence` | windowed peptide |
| `start`, `end` | 1-indexed inclusive parent coords |
| `length` | `end - start + 1` |
| `study` | study name |
| `position` | anchor (e.g. cut site), parent coords |
| meta cols | carried through (e.g. `cell_type`, `schaeffer_cut_topo`, ...) |

Future protein-style studies (full-sequence external libraries) would emit one row per protein — same `item_id` / `uniprot` / `sequence`, no `start/end/position`.

### Generic builder: `from_external_table`

```python
from protein_sets import from_external_table

rows = from_external_table(
    positions,                  # df with cols: uniprot, position, **meta
    study='schaeffer22',
    window=15,
    anchor='center',            # or 'n_term' | 'c_term'
    kind='cleavage_w15',        # provenance kind
)
```

Handles parent-sequence lookup (via `resolve_with_seqs`), windowing, item_id minting, OOB drops, and provenance. All `meta_cols` flow through and are included in the item_id for uniqueness — pre-trim `positions` to the columns that matter for uniqueness, then `parse.py` can re-merge richer per-event metadata afterwards.

### Resolver integration

Study names are first-class in `resolve()`:

```python
resolve('schaeffer22')           # → list of parent UniProts (deduped)
```

So `/protein-structures` and any other UniProt-level pipe works with a study name. For region-level downstream work (embedding the actual 15mers), pass the loaded DataFrame: `embed_regions(studies.load('schaeffer22'))`.

### Adding a new study

1. `mkdir /ewsc/jlaffy/protein_sets/studies/<name>/{raw,parsed}`
2. Drop raw downloads under `raw/`
3. Write `parse.py`: load raw → build positions df → `from_external_table(...)` → optionally merge richer metadata → write `parsed/rows.parquet`
4. Write `metadata.yaml` (see `schaeffer22` for template): `name`, `kind` (`region` | `protein`), `description`, `source` (citation/DOI/URL), `parent_universe`, `window` (if region), `columns`, `provenance` (source_category/subgroup/kind)
5. Run `python parse.py`
6. Verify: `studies.available()` shows it; `studies.load(name)` returns expected rows

### What's bespoke vs. shared

Bespoke (per-study): `raw/` contents, `parse.py` logic, `metadata.yaml` contents.
Shared (in package, don't duplicate): the row contract, `from_external_table`, `studies.load/info/available`, resolver integration, provenance/layout conventions.

## How to interpret requests

| user says | do |
|---|---|
| "how many proteins have signal peptides?" | `len(signal_peptide())` |
| "give me secretome IDs" | `ids, _, _ = resolve_dynamic('secretome')` |
| "give me the surfaceome" (cell surface proteins) | `surfaceome()` — CSPA high+putative by default; use `surfaceome_strict` for high-confidence only, or `load_cspa_human()` for the raw table with CD markers, confidence categories, detection counts |
| "CSPA N-glycopeptides" | `load_cspa_glycopeptides()` — peptide-level table with modified-form notation (N[115] for N-glycosylation site) |
| "transcription factors" / "human TF list" / "is X a TF?" | `transcription_factors()` — Lambert et al. 2018 curated census (the field standard); `load_humantfs()` for the raw table incl. DBD family / binding mode / `Is TF?` flag. Prefer this over UniProt keywords (`Transcription regulation` is broad/noisy) or GO:0003700 |
| "which transmembrane proteins are also in the secretome?" | intersect `transmem()` and `secretome()` DataFrames on `Entry` |
| "proteins with GPI anchors but no transmembrane domain" | filter `load_proteome()`: `Keywords` contains `GPI-anchor` AND `Transmembrane` is empty |
| "what Pfam domains does P01308 have?" | `load_interpro()[lambda d: d.uniprot == 'P01308']` or check the `Pfam` column in `load_proteome()` |
| "extract Pfam domain sequences for transmem proteins" | `pfam_regions('transmem')` — DataFrame with sequences |
| "count proteins by subcellular location" | `load_proteome()` then parse/group the `Subcellular location [CC]` column |
| "write a FASTA for the secretome" | resolve IDs + sequences, write `>accession\nsequence\n` format |
| "explore membrane proteins beyond just transmembrane" | combine transmem, intramembrane, GPI-anchor (Keywords + Lipidation columns) from the proteome |
| "what InterPro entry types are there?" | `load_entry_types()` or `load_interpro()['entry_type'].value_counts()` |
| "give me UniRef50-deduplicated `<set>`" / "non-redundant `<set>` at 50% identity" | `dedupe('<set>', level=50)` — group by the `UniRef50` column, keep one row per cluster |
| "redundancy-reduced training set" / "dedup before train/test split" | `dedupe(set, level=50)` — UniRef50 is the standard for ML splits; UniRef90 for milder dedup |
| "what UniRef50 cluster is `P01308` in?" | `load_proteome()` then `df.loc[df.Entry == 'P01308', 'UniRef50']`; or `fetch_cluster_ids(['P01308'], level=50)` for an arbitrary accession not in the human proteome |
| "expand `<set>` to cross-species cluster members" | not yet supported — that's UniRef Phase 2; today only the forward map (human accession → cluster_id) is cached, not the reverse (cluster_id → all members) |
| "5×AA window mutation scan on insulin" | `mutation_scan('P01308', window=5)` then `embed_regions(df, ...)` |
| "scan the secretome with mutations to A and R only" | `mutation_scan('secretome', alphabet=('A','R'))` |
| "write the mutants to FASTA" | `to_fasta(df, path)` |
| "embed the variants on /ewsc" | set `ESM_EMBED_CACHE_ROOT=/ewsc/jlaffy/embeddings_variants` then `embed_regions(...)` |
| "alanine scan", "truncation series", "graft this region from A onto B" | not yet implemented — would be added as new transforms in `transforms.py` following the same contract |
| "what studies do I have?" / "list bespoke datasets" | `studies.available()` |
| "load the schaeffer cut sites" / "give me schaeffer22 15mers" | `studies.load('schaeffer22')` |
| "just the HEK293 cut sites" | `studies.load('schaeffer22', cell_type='HEK293')` |
| "structures for schaeffer22 proteins" | `resolve('schaeffer22')` → list of parent UniProts → pass to /protein-structures |
| "register a new study X" / "add this paper's cut sites as a universe" | follow the layout in the Studies section: `studies/<name>/{raw,parsed,metadata.yaml,parse.py}`; use `from_external_table` for site-style studies |

## What NOT to do

- Don't generate embeddings from this skill — that's the esm-embed skill.
- Don't hardcode protein counts — they change when UniProt refreshes. Always compute from the data.
- Don't assume column names — verify with `df.columns`. The actual names are listed above but UniProt can change them across releases.
- Don't re-download data unnecessarily — `load_proteome()` auto-pulls if stale, no need to call `pull_proteome()` explicitly unless forcing a refresh.
