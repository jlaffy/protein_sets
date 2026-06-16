# protein-sets

Resolve, filter, and subset human protein universes by annotation — a single
row contract `(item_id, uniprot, sequence, ...)` shared across **named sets**
(secretome, transmembrane, signal-peptide, surfaceome, TFs…), **regions**
(domains, signal peptides, sites), **transforms** (mutation scans), and
**studies** (your own bespoke datasets).

All reference data **auto-downloads on first use** and is cached locally — no
manual file wrangling.

## Install

```bash
pip install git+https://github.com/jlaffy/protein_sets
```

## Quickstart

```python
from protein_sets import load_proteome, signal_peptide, secretome, resolve

load_proteome()          # full reviewed-human UniProt table (auto-downloads once)
signal_peptide()         # proteins with a signal-peptide annotation
secretome()              # HPA predicted secreted ∩ reviewed human
resolve("secretome")     # → list of UniProt accessions
```

The first call fetches what it needs (UniProt, HPA, CSPA, humantfs, UniRef) into
a local cache; every later call is instant. The 30-day TTL refreshes stale data
automatically.

## Where data is cached

Resolved once at import, in priority order:

1. `$PROTEIN_SETS_DATA` — set this to put the cache wherever you like.
2. `/ewsc/jlaffy/protein_sets` — the shared NFS cache, **if it exists** (cluster
   users transparently reuse it — nothing to configure).
3. `~/.cache/protein_sets` — default for a fresh install.

```bash
export PROTEIN_SETS_DATA=~/my_protein_cache   # optional
```

## What auto-downloads vs. what needs data

| capability | data |
|---|---|
| proteome, named sets, UniProt-feature regions, transforms, UniRef dedup | **auto-downloads** ✅ |
| InterPro domain regions | needs a local InterPro file (auto-download planned for 2.0) |
| OMA orthologs | bulk OMA download on first use |
| studies | bring your own — drop a dataset under `<cache>/studies/<name>/` |

## Adding your own study

```
<cache>/studies/<name>/
    raw/                 # original downloads
    parsed/rows.parquet  # standardized rows (required)
    metadata.yaml        # required
    parse.py             # raw → parsed, re-runnable
```

Then `studies.load("<name>")` returns rows on the standard contract, usable
everywhere a named set is.

## Claude Code skill

`SKILL.md` ships in this repo — copy it to `~/.claude/skills/protein-sets/` (or
a project's `.claude/skills/`) to get `/protein-sets` routing in Claude Code.
