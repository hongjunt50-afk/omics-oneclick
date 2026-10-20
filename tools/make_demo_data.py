"""生成合成表达谱用于自测：含真实差异基因与共表达模块结构。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HALLMARK_LIKE = [
    "Cell cycle", "Apoptosis", "Inflammatory response", "Oxidative phosphorylation",
    "Epithelial mesenchymal transition", "Glycolysis", "DNA repair", "Hypoxia",
    "TNF signaling", "Interferon response", "p53 pathway", "Wnt signaling",
]


def _symbol_pool(seed: int) -> list[str]:
    """优先用内置 KEGG 基因集里的真实基因符号，让富集/GSEA 能真正命中。"""
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    try:
        from app.core.genesets import asset_dir

        gmt = asset_dir() / "KEGG_2021_Human.gmt"
        if gmt.exists():
            symbols: list[str] = []
            seen = set()
            with open(gmt, encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    for gene in line.rstrip("\n").split("\t")[2:]:
                        gene = gene.strip().upper()
                        if gene and gene not in seen:
                            seen.add(gene)
                            symbols.append(gene)
            if len(symbols) >= 1500:
                return symbols
    except Exception:
        pass
    return [f"GENE{i:05d}" for i in range(20000)]


def make_dataset(
    out_dir: Path,
    *,
    n_genes: int = 2000,
    per_group: int = 5,
    groups: tuple[str, ...] = ("Control", "Treat"),
    n_de: int = 80,
    seed: int = 7,
) -> tuple[Path, Path]:
    rng = np.random.default_rng(seed)
    total = per_group * len(groups)
    samples = [f"{g}_{i + 1}" for g in groups for i in range(per_group)]
    labels = [g for g in groups for _ in range(per_group)]

    pool = _symbol_pool(seed)
    chosen = list(rng.choice(pool, size=min(n_genes, len(pool)), replace=False))

    # 12 个共表达模块，模块内基因共享潜在因子
    n_modules = 12
    per_module = n_genes // n_modules
    latent = rng.normal(size=(n_modules, total))
    loadings = rng.normal(0.8, 0.35, size=(n_modules, per_module))

    blocks = []
    for m in range(n_modules):
        block = loadings[m][:, None] * latent[m][None, :] + rng.normal(0, 0.55, size=(per_module, total))
        blocks.append(block)

    extra = n_genes - per_module * n_modules
    if extra > 0:
        blocks.append(rng.normal(0, 0.6, size=(extra, total)))

    values = np.vstack(blocks)

    # 注入差异基因
    treat_idx = np.array([i for i, g in enumerate(labels) if g == groups[-1]])
    ctrl_idx = np.array([i for i, g in enumerate(labels) if g == groups[0]])
    de_rows = rng.choice(values.shape[0], size=n_de, replace=False)
    for k, row in enumerate(de_rows):
        shift = 1.6 if k % 2 == 0 else -1.6
        values[row, treat_idx] += shift
        values[row, ctrl_idx] -= shift * 0.35

    # 转成对数正态表达值
    expr = np.exp(values * 0.9 + 5.0)
    frame = pd.DataFrame(expr, index=chosen[: values.shape[0]], columns=samples).round(3)
    frame = frame[~frame.index.duplicated()]
    frame.index.name = "gene_symbol"

    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = out_dir / "demo_expression_matrix.csv"
    frame.to_csv(matrix_path)

    phenotype = pd.DataFrame({"sample": samples, "condition": labels})
    pheno_path = out_dir / "demo_phenotype.csv"
    phenotype.to_csv(pheno_path, index=False)
    return matrix_path, pheno_path


if __name__ == "__main__":
    import sys

    target = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    matrix, pheno = make_dataset(target)
    print(f"写入 {matrix}")
    print(f"写入 {pheno}")
