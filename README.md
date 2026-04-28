# getKEGGpathways

A Python package for batch-fetching KEGG pathway gene lists as pandas DataFrames.

## Installation

```bash
pip install /path/to/getkeggpathways
```

Dependencies: `requests`, `pandas`, `tqdm`

## Quick Start

```python
from getKEGGpathways import KEGGpathways

# Fetch all pathway genes for mouse (mmu)
df = KEGGpathways.get("mmu")
print(df.head())
```

Output format:

```
  organism_id pathway_id              pathway_name gene_id gene_name
0         mmu   mmu00010      Glycolysis / Gluconeogenesis   103988   Gck
1         mmu   mmu00010      Glycolysis / Gluconeogenesis   106557   Ldhal6b
...
```

The returned `DataFrame` has 5 columns:

| Column | Description |
|---|---|
| `organism_id` | KEGG organism code |
| `pathway_id` | Pathway identifier |
| `pathway_name` | Pathway description |
| `gene_id` | Entrez gene ID |
| `gene_name` | Gene symbol |

## API

### `KEGGpathways.get()`

```python
df = KEGGpathways.get(
    organism="mmu",       # organism code, e.g. hsa, mmu, dre
    *,
    latency=0.3,          # seconds between batch requests (avoids rate-limiting)
    batch_size=10,        # pathway IDs per GET request
    timeout=30,           # HTTP request timeout in seconds
    max_retries=3,        # max retries on network / 5xx errors
    n_pathways=None,      # limit to first N pathways (useful for testing)
    force_refresh=False,  # if True, ignore cached pickle and re-fetch
    no_org_names=True     #Remove the comment suffix ( - organism_name ...) for pathway names
    source_col=None,      # copy the specified column as a new column `source` (useful for decoupler)
    target_col=None       # copy the specified column as a new column `target` (useful for decoupler)
)
```

## Caching

- **Organism list**: cached at `~/.cache/getkeggpathways/organisms.txt` after first request, valid for 30 days.
- **Query results**: each successful `get()` call caches the full result at `~/.cache/getkeggpathways/{organism}.pkl`, valid for 30 days. Subsequent calls read from cache with no network requests.

Use `force_refresh=True` to bypass the cache and re-fetch from the KEGG API.

## Supported Organisms

To browse all available organism codes:

```python
from getKEGGpathways.kegg import _ORGANISM_CACHE_FILE
import pandas as pd

df = pd.read_csv(_ORGANISM_CACHE_FILE, sep="\t", names=["id", "name"])
print(df.head())
```

Common organisms:

| Code | Species |
|---|---|
| `hsa` | Homo sapiens (human) |
| `mmu` | Mus musculus (mouse) |
| `rno` | Rattus norvegicus (rat) |
| `dre` | Danio rerio (zebrafish) |
