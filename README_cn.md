# getKEGGpathways

批量从 [KEGG](https://www.genome.jp/kegg/) 数据库获取通路基因列表的 Python 包。

## 安装

```bash
pip install -e /path/to/getkeggpathways
```

依赖：`requests`, `pandas`, `tqdm`

## 快速开始

```python
from getKEGGpathways import KEGGpathways

# 获取小鼠 (mmu) 所有通路的基因列表
df = KEGGpathways.get("mmu")
print(df.head())
```

输出格式：

```
  organism_id pathway_id              pathway_name gene_id gene_name
0         mmu   mmu00010      Glycolysis / Gluconeogenesis   103988   Gck
1         mmu   mmu00010      Glycolysis / Gluconeogenesis   106557   Ldhal6b
...
```

返回的 `DataFrame` 包含 5 列：

| 列名 | 说明 |
|---|---|
| `organism_id` | 物种 KEGG 缩写 |
| `pathway_id` | 通路 ID |
| `pathway_name` | 通路名称 |
| `gene_id` | 基因 Entrez ID |
| `gene_name` | 基因简称 |

## API

### `KEGGpathways.get()`

```python
df = KEGGpathways.get(
    organism="mmu",       # 物种缩写，如 hsa, mmu, dre 等
    *,
    latency=0.3,          # 请求间隔（秒），避免 IP 被封
    batch_size=10,        # 每次请求的通路数量
    timeout=30,           # HTTP 超时（秒）
    max_retries=3,        # 网络错误最大重试次数
    n_pathways=None,      # 限制只获取前 N 个通路（测试用）
    force_refresh=False,  # 强制重新获取，忽略缓存
    no_org_names=True     # 去掉pathway_name内容后面的物种名称
    source_col=None,      # 把指定的列作为“source”列添加到返回数据的左边
    target_col=None       # 把指定的列作为“target”列添加到返回数据的左边
)
```

## 数据缓存

- **物种列表**：首次请求后缓存到 `~/.cache/getkeggpathways/organisms.txt`，有效期 30 天
- **查询结果**：每次成功的 `get()` 调用会将完整结果缓存到 `~/.cache/getkeggpathways/{物种}.pkl`，有效期 30 天。下次调用时直接读取缓存，无需联网

使用 `force_refresh=True` 可强制忽略缓存重新获取。

## 支持物种

查看所有支持的物种缩写：

```python
from getKEGGpathways.kegg import _ORGANISM_CACHE_FILE
import pandas as pd

df = pd.read_csv(_ORGANISM_CACHE_FILE, sep="\t", names=["id", "name"])
print(df.head())
```

常用物种：

| 缩写 | 物种 |
|---|---|
| `hsa` | Homo sapiens (human) |
| `mmu` | Mus musculus (mouse) |
| `rno` | Rattus norvegicus (rat) |
| `dre` | Danio rerio (zebrafish) |
| `dme` | Drosophila melanogaster (fruit fly) |
| `cel` | Caenorhabditis elegans (nematode) |
| `sce` | Saccharomyces cerevisiae (yeast) |
| `eco` | Escherichia coli |
