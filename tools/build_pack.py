"""对已有的结果目录重新生成提交材料包（不重跑分析）。

用法：
    python tools/build_pack.py <结果目录> [-o 输出目录] [--screenshot 界面截图.png]

为什么不重跑：材料包只读 tables/ 与 summary.json，重跑既慢又可能改变随机结果。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.core.submission import build_submission_pack  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="从结果目录生成 28 项提交材料包")
    parser.add_argument("result_dir", help="OmicsOne_<项目名>_<时间戳> 目录")
    parser.add_argument("-o", "--output", default=None, help="材料包输出目录（默认 <结果目录>/提交材料）")
    parser.add_argument("--screenshot", default=None, help="软件导入/下载界面截图（对应材料 1-1）")
    parser.add_argument("--heatmap-top", type=int, default=50, help="Top N 差异基因热图的 N（默认 50）")
    args = parser.parse_args(argv)

    pack = build_submission_pack(
        args.result_dir,
        args.output,
        screenshot=args.screenshot,
        heatmap_top=args.heatmap_top,
    )
    print(f"材料目录：{pack.out_dir}")
    print(pack.summary_line())
    for row in pack.checks:
        print(f"  {'OK ' if row.ok else '!! '} {row.code:<4} {row.name:<30} {row.measured:>12} / {row.target}")
    for warning in pack.warnings:
        print(f"  [提示] {warning}")
    return 0 if pack.ok_count == len(pack.checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
