"""命令行入口：便于自动化测试与批量跑批，GUI 之外的第二通道。"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.core.data_loader import load_dataset
from app.core.pipeline import AnalysisConfig, Runner
from app.core.report import build_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="omicsone", description="OmicsOne 一键组学分析（命令行）")
    parser.add_argument("data", help="表达矩阵文件路径")
    parser.add_argument("-o", "--output", default="OmicsOne_Results", help="结果输出目录")
    parser.add_argument("--groups", help="分组表（样本名,组名）")
    parser.add_argument("--group-a", help="比较组 A（分子）")
    parser.add_argument("--group-b", help="比较组 B（分母）")
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--lfc", type=float, default=1.0)
    parser.add_argument("--normalization", default="quantile", choices=["quantile", "median", "none"])
    parser.add_argument("--ora", nargs="*", default=["KEGG_2021_Human", "GO_Biological_Process_2023"])
    parser.add_argument("--gsea", default="KEGG_2021_Human")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--no-wgcna", action="store_true")
    parser.add_argument("--steps", nargs="*", default=None, help="只跑指定步骤：qc deg pca enrichment gsea wgcna")
    args = parser.parse_args(argv)

    def progress(fraction: float, message: str) -> None:
        bar = "█" * int(fraction * 28)
        print(f"\r[{bar:<28}] {fraction * 100:5.1f}% {message[:40]:<40}", end="", flush=True)

    def log(level: str, message: str) -> None:
        print(f"\n[{level}] {message}")

    dataset = load_dataset(args.data, group_file=args.groups)
    config = AnalysisConfig(
        output_dir=Path(args.output),
        group_a=args.group_a,
        group_b=args.group_b,
        padj_cut=args.padj,
        lfc_cut=args.lfc,
        normalization=args.normalization,
        ora_libraries=list(args.ora),
        gsea_library=args.gsea,
        gsea_permutations=args.permutations,
        run_wgcna=not args.no_wgcna,
    )
    runner = Runner(dataset, config, progress=progress, log=log)
    result = runner.run_all(args.steps)
    report = build_report(result, dataset, config)
    progress(1.0, "完成")
    print(f"\n结果目录：{result.project_dir}")
    print(f"HTML 报告：{report}")
    if result.errors:
        print("部分步骤失败：")
        for err in result.errors:
            print(f"  - {err}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
