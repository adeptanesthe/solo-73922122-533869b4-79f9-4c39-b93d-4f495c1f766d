"""从 samples/demo.fa 生成单样本冲突核查的演示 VCF（坐标由参考自动计算，保证 REF 合法）。

运行：python3 -I scripts/make_conflict_demo.py
生成三份文件：
  conflict_feasible.vcf    可行：含同碱基两个剂量 1 SNV、纯合 SNV、缺失与锚点 SNV 共存
  conflict_odd_cycle.vcf   无解：同一参考碱基上三个剂量 1 SNV（任意两条可安排，三条无解）
  conflict_hom.vcf         无解：剂量 2 变异与另一变异共享被替换/删除的参考碱基
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import parse_fasta  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def find(seq: str, motif: str, start: int = 1) -> int:
    """返回 motif 在参考中的 1-based 起点。"""
    idx = seq.find(motif, start - 1)
    if idx < 0:
        raise SystemExit(f"motif {motif!r} not found in reference")
    return idx + 1


def row(pos, ref, alt, gt):
    return f"chrDemo\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}"


def main() -> None:
    ref = parse_fasta((ROOT / "samples" / "demo.fa").read_bytes())
    seq = ref.sequence
    out_dir = ROOT / "samples"

    # 找一个参考碱基 C（SNV 演示位）与其相邻的 G（另一个独立位点）
    p_cg = find(seq, "CG")
    c_pos = p_cg          # 参考 C
    g_pos = p_cg + 1      # 参考 G
    # 再找一个 A 串，构造纯缺失（锚点 + 被删 A）与锚点 SNV
    a_run = find(seq, "AAAA")
    anchor_pos = a_run - 1
    anchor_base = seq[anchor_pos - 1]

    # 1) 可行场景 ---------------------------------------------------------------
    feasible = [
        # 同一碱基 c_pos 上两个不同 ALT，各剂量 1：分占两条拷贝即可
        row(c_pos, "C", "G", "0/1"),
        row(c_pos, "C", "T", "0/1"),
        # 另一位点的纯合 SNV（剂量 2），与前者碱基不同，不冲突
        row(g_pos, "G", "A", "1/1"),
        # 纯缺失（剂量 1，删除 A 串第一个 A）与锚点 SNV：锚点不计入被删碱基，可共存
        row(anchor_pos, anchor_base + "A", anchor_base, "0/1"),
        row(anchor_pos, anchor_base, "C" if anchor_base != "C" else "G", "0/1"),
    ]
    (out_dir / "conflict_feasible.vcf").write_text(
        HEADER + "\n".join(feasible) + "\n", encoding="utf-8"
    )

    # 2) 奇环无解场景：同一参考碱基上三个剂量 1 SNV ------------------------------
    odd_cycle = [
        row(c_pos, "C", "G", "0/1"),
        row(c_pos, "C", "T", "0/1"),
        row(c_pos, "C", "A", "0/1"),
        # 附带一个完全无关的可行变异，证明证据只取奇环节点
        row(g_pos, "G", "T", "0/1"),
    ]
    (out_dir / "conflict_odd_cycle.vcf").write_text(
        HEADER + "\n".join(odd_cycle) + "\n", encoding="utf-8"
    )

    # 3) 纯合占位无解场景：剂量 2 SNV 与同碱基另一剂量 1 SNV ----------------------
    hom = [
        row(c_pos, "C", "G", "1/1"),
        row(c_pos, "C", "T", "0/1"),
    ]
    (out_dir / "conflict_hom.vcf").write_text(
        HEADER + "\n".join(hom) + "\n", encoding="utf-8"
    )

    print(f"wrote {len(feasible)} records -> samples/conflict_feasible.vcf")
    print(f"wrote {len(odd_cycle)} records -> samples/conflict_odd_cycle.vcf")
    print(f"wrote {len(hom)} records -> samples/conflict_hom.vcf")


if __name__ == "__main__":
    main()
