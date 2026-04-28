import re
import time
from pathlib import Path
from warnings import warn

import pandas as pd
import requests
from tqdm import tqdm

_BASE_URL = "https://rest.kegg.jp"
_CACHE_DIR = Path.home() / ".cache" / "getkeggpathways"
_ORGANISM_CACHE_FILE = _CACHE_DIR / "organisms.txt"
_CACHE_TTL_DAYS = 30

_DEFAULT_BATCH_SIZE = 10
_DEFAULT_LATENCY = 0.3
_DEFAULT_TIMEOUT = 30
_DEFAULT_MAX_RETRIES = 3

_COLUMNS = ["organism_id", "pathway_id", "pathway_name", "gene_id", "gene_name"]


def _pickle_path(organism):
    return _CACHE_DIR / f"{organism}.pkl"


class KEGGpathways:

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @staticmethod
    def get(
        organism="hsa",
        *, # Keyword-only
        latency=None,
        batch_size=None,
        timeout=None,
        max_retries=None,
        n_pathways=None,
        force_refresh=False,
        no_org_names=True,
        source_col=None,
        target_col=None,
    ):
        """Return a DataFrame of KEGG pathway genes for *organism*.

        Parameters
        ----------
        organism : str
            KEGG organism code (e.g. ``"mmu"``, ``"hsa"``).
        latency : float, optional
            Seconds to sleep between batch requests.
        batch_size : int, optional
            Number of pathway IDs sent per GET request.
        timeout : int, optional
            Seconds before each HTTP request times out.
        max_retries : int, optional
            Maximum retries on network / 5xx errors.
        n_pathways : int, optional
            Limit fetching to the first *n_pathways* (useful for testing).
        force_refresh : bool, optional
            If True, ignore cached pickle and re-fetch from KEGG API.
        no_org_names: bool, default: True
            Remove the organism suffix ( - organism_name ...) for pathway names
        source_col: str, optional
            Add copy the specified column as a new column `source` (useful for decoupler)
        target_col: str, optional
            Add copy the specified column as a new column `target` (useful for decoupler)
        """
        if latency is None:
            latency = _DEFAULT_LATENCY
        if batch_size is None:
            batch_size = _DEFAULT_BATCH_SIZE
        if timeout is None:
            timeout = _DEFAULT_TIMEOUT
        if max_retries is None:
            max_retries = _DEFAULT_MAX_RETRIES

        KEGGpathways._ensure_cache_dir()

        # Return cached result if available and fresh
        ppath = _pickle_path(organism)
        if not force_refresh and ppath.exists():
            age_days = (time.time() - ppath.stat().st_mtime) / 86400
            if age_days < _CACHE_TTL_DAYS:
                df = pd.read_pickle(ppath)
                if n_pathways and n_pathways > 0:
                    # Filter cached result to first n_pathways' worth
                    first_ids = df["pathway_id"].unique()[:n_pathways]
                    df = df[df["pathway_id"].isin(first_ids)]
                return df

        organisms = KEGGpathways._load_organism_cache(timeout, max_retries)
        if organism not in organisms:
            raise ValueError(
                f"Unknown organism code: {organism!r}. "
                f"Available codes: {', '.join(sorted(organisms)[:20])}..."
            )

        pathways = KEGGpathways._fetch_pathways(organism, timeout, max_retries)
        if not pathways:
            return pd.DataFrame(columns=_COLUMNS)
        if n_pathways and n_pathways > 0:
            pathways = pathways[0:n_pathways]

        pathway_ids = [p[0] for p in pathways]
        all_records = []
        batches = list(KEGGpathways._batch(pathway_ids, batch_size))

        for batch in tqdm(batches, desc=f"Fetching {organism} pathways", unit="batch"):
            text = KEGGpathways._fetch_pathway_details(batch, timeout, max_retries)
            all_records.extend(KEGGpathways._parse_get_response(text, organism))
            if len(batch) == batch_size:
                time.sleep(latency)

        if not all_records:
            return pd.DataFrame(columns=_COLUMNS)
        df = pd.DataFrame(all_records, columns=_COLUMNS)
        if not n_pathways:
            df.to_pickle(ppath)
        
        if no_org_names:
            df["pathway_name"] = df['pathway_name']\
                .str.replace(r' - .*', '', regex=True)
        if target_col:
            if target_col in df.columns:
                df.insert(0, "target", df[target_col])
            else:
                warn(f"column {target_col} not found in data.")
        if source_col:
            if source_col in df.columns:
                df.insert(0, "source", df[source_col])
            else:
                warn(f"column {source_col} not found in data.")
        return df

    # ------------------------------------------------------------------
    # Cache helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _ensure_cache_dir():
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _load_organism_cache(timeout, max_retries):
        if _ORGANISM_CACHE_FILE.exists():
            age_days = (time.time() - _ORGANISM_CACHE_FILE.stat().st_mtime) / 86400
            if age_days < _CACHE_TTL_DAYS:
                content = _ORGANISM_CACHE_FILE.read_text(encoding="utf-8")
                organisms = KEGGpathways._parse_organism_tsv(content)
                if organisms:
                    return organisms

        resp = KEGGpathways._request_with_retry(
            f"{_BASE_URL}/list/organism", timeout, max_retries
        )
        _ORGANISM_CACHE_FILE.write_text(resp.text, encoding="utf-8")
        return KEGGpathways._parse_organism_tsv(resp.text)

    # ------------------------------------------------------------------
    # HTTP layer
    # ------------------------------------------------------------------

    @staticmethod
    def _request_with_retry(url, timeout, max_retries):
        last_exc = None
        for attempt in range(max_retries + 1):
            try:
                resp = requests.get(url, timeout=timeout)
                resp.raise_for_status()
                return resp
            except requests.HTTPError:
                if resp.status_code < 500:
                    raise
                last_exc = RuntimeError(
                    f"KEGG API returned {resp.status_code} for {url}"
                )
            except requests.RequestException as exc:
                last_exc = exc

            if attempt < max_retries:
                time.sleep(2 ** attempt)
        raise ConnectionError(
            f"Failed to fetch {url} after {max_retries + 1} attempts"
        ) from last_exc

    @staticmethod
    def _fetch_organisms(timeout, max_retries):
        url = f"{_BASE_URL}/list/organism"
        resp = KEGGpathways._request_with_retry(url, timeout, max_retries)
        return KEGGpathways._parse_organism_tsv(resp.text)

    @staticmethod
    def _fetch_pathways(organism, timeout, max_retries):
        url = f"{_BASE_URL}/list/pathway/{organism}"
        resp = KEGGpathways._request_with_retry(url, timeout, max_retries)
        pathways = []
        for line in resp.text.strip().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                pathways.append((parts[0].strip(), parts[1].strip()))
        return pathways

    @staticmethod
    def _fetch_pathway_details(pathway_ids, timeout, max_retries):
        url = f"{_BASE_URL}/get/{'+'.join(pathway_ids)}"
        resp = KEGGpathways._request_with_retry(url, timeout, max_retries)
        return resp.text

    # ------------------------------------------------------------------
    # Parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_organism_tsv(text):
        organisms = {}
        for line in text.strip().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                org_id = parts[1].strip()
                org_name = parts[2].strip()
                organisms[org_id] = org_name
        return organisms

    @staticmethod
    def _parse_get_response(text, organism):
        records = []
        entries = text.split("///")

        for entry in entries:
            if not entry.strip():
                continue

            pathway_id = None
            pathway_name = None
            in_gene_section = False
            gene_lines = []

            for line in entry.splitlines():
                if not line.strip():
                    continue

                # Match section headers
                if not line[0].isspace():
                    in_gene_section = False
                    if line.startswith("ENTRY"):
                        m = re.match(r"ENTRY\s+(\w+)", line)
                        if m:
                            pathway_id = m.group(1)
                    elif line.startswith("NAME"):
                        m = re.match(r"NAME\s+(.+)", line)
                        if m:
                            pathway_name = m.group(1).strip()
                    elif line.startswith("GENE"):
                        in_gene_section = True
                        gene_lines.append(line)
                    continue

                # Continuation lines
                if in_gene_section:
                    gene_lines.append(line)

            if not pathway_id:
                continue

            for gene_line in gene_lines:
                # Remove section header "GENE" prefix if present
                clean = re.sub(r"^GENE\s+", "", gene_line).strip()
                if not clean:
                    continue
                parts = clean.split()
                if len(parts) >= 2:
                    gene_id = parts[0]
                    gene_name = parts[1].rstrip(";")
                    records.append({
                        "organism_id": organism,
                        "pathway_id": pathway_id,
                        "pathway_name": pathway_name or "",
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                    })

        return records

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------

    @staticmethod
    def _batch(iterable, size):
        batch = []
        for item in iterable:
            batch.append(item)
            if len(batch) == size:
                yield batch
                batch = []
        if batch:
            yield batch
