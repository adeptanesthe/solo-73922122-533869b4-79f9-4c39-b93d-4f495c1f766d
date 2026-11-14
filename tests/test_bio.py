"""核心逻辑测试，仅依赖标准库：python -I -m unittest -v"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import (  # noqa: E402
    CONFLICT_ALLOWED_KINDS,
    CheckError,
    compare,
    normalize,
    parse_fasta,
    parse_vcf,
    check_conflicts,
    affected_bases,
    build_items,
)

# 1-based: 1A 2C 3G 4T | 5A 6A 7A 8A | 9C 10G 11T 12T
SEQ = "ACGTAAAACGTT"
FASTA = f">chrT\n{SEQ}\n".encode()


def vcf(records, sample="S1", extra_header_cols=""):
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample
        + extra_header_cols,
    ]
    lines.extend(records)
    return "\n".join(lines) + "\n"


def rec(pos, ref, alt, gt, chrom="chrT", fmt="GT"):
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t{fmt}\t{gt}"


def parse(data):
    refobj = parse_fasta(FASTA)
    calls, _ = parse_vcf(data.encode(), "test.vcf", refobj)
    return refobj, calls


class NormalizeTests(unittest.TestCase):
    def test_insertion_in_homopolymer_different_writing(self):
        # 在 5-8 的 A 串中插入一个 A：pos4 T>TA 与 pos5 A>AA 应归一为同一变异
        a = normalize(4, "T", "TA", SEQ)
        b = normalize(5, "A", "AA", SEQ)
        self.assertEqual(a, b)
        self.assertEqual(a, (4, "T", "TA"))

    def test_deletion_in_homopolymer_different_writing(self):
        a = normalize(4, "TA", "T", SEQ)   # 删 pos5 的 A
        b = normalize(5, "AA", "A", SEQ)   # 删 pos6 的 A，同属 A 串
        self.assertEqual(a, b)
        self.assertEqual(a, (4, "TA", "T"))

    def test_snv_not_moved(self):
        self.assertEqual(normalize(2, "C", "G", SEQ), (2, "C", "G"))

    def test_different_indel_lengths_stay_distinct(self):
        one = normalize(4, "TA", "T", SEQ)
        two = normalize(4, "TAA", "T", SEQ)
        self.assertNotEqual(one, two)

    def test_repeat_run_tt_at_end_deletion(self):
        # 末端 TT 串（pos11-12）中删一个 T：两种写法归一到最左锚点
        a = normalize(10, "GT", "G", SEQ)    # 删 pos11 的 T
        b = normalize(11, "TT", "T", SEQ)    # 删 pos12 的 T
        self.assertEqual(a, b)
        self.assertEqual(a[0], 10)

    def test_insertion_shift_stops_when_it_would_become_mnv(self):
        # pos12 后插入 G：再左移就不再是纯插入，规范结果保持在原位
        self.assertEqual(normalize(12, "T", "TG", SEQ), (12, "T", "TG"))


class ParseAndCompareTests(unittest.TestCase):
    def test_snv_consistent_with_line_numbers(self):
        _, ca = parse(vcf([rec(2, "C", "G", "0/1")]))
        _, cb = parse(vcf([rec(2, "C", "G", "0/1")]))
        refobj = parse_fasta(FASTA)
        groups = compare(refobj, ca, cb)
        self.assertEqual(len(groups["consistent"]), 1)
        entry = groups["consistent"][0]
        self.assertEqual(entry["left"].dosage, 1)
        self.assertEqual(entry["left"].line_numbers, [3])
        self.assertEqual(entry["right"].line_numbers, [3])

    def test_indel_written_differently_treated_as_same_variant(self):
        _, ca = parse(vcf([rec(4, "T", "TA", "0/1")]))
        _, cb = parse(vcf([rec(5, "A", "AA", "0/1")]))
        refobj = parse_fasta(FASTA)
        groups = compare(refobj, ca, cb)
        self.assertEqual(len(groups["consistent"]), 1)
        self.assertEqual(groups["consistent"][0]["variant"].label, "4:TinsTA")
        # 原始行号与原始写法都要保留
        self.assertEqual(groups["consistent"][0]["right"].raw_descriptions, ["5:A>AA"])

    def test_multi_alt_dosages_per_allele(self):
        _, calls = parse(vcf([rec(2, "C", "G,T", "1/2")]))
        keys = sorted(calls, key=lambda v: v.alt)
        self.assertEqual([k.alt for k in keys], ["G", "T"])
        self.assertTrue(all(c.dosage == 1 for c in calls.values()))

    def test_multi_alt_homozygous_alt_and_ignore_zero_dose(self):
        _, calls = parse(vcf([rec(2, "C", "G,T", "1/1")]))
        self.assertEqual(len(calls), 1)  # T 剂量 0 被忽略
        (only,) = calls.values()
        self.assertEqual(only.dosage, 2)

    def test_dosage_mismatch_group(self):
        _, ca = parse(vcf([rec(2, "C", "G", "1/1")]))
        _, cb = parse(vcf([rec(2, "C", "G", "0/1")]))
        groups = compare(parse_fasta(FASTA), ca, cb)
        self.assertEqual(len(groups["dosage_mismatch"]), 1)
        self.assertEqual(groups["consistent"], [])

    def test_single_side_groups(self):
        _, ca = parse(vcf([rec(2, "C", "G", "0/1")]))
        _, cb = parse(vcf([rec(3, "G", "T", "0/1")]))
        groups = compare(parse_fasta(FASTA), ca, cb)
        self.assertEqual(len(groups["left_only"]), 1)
        self.assertEqual(len(groups["right_only"]), 1)

    def test_two_insertion_alts_of_different_length_are_distinct(self):
        # 同一记录的两个 ALT：插 1 个 A 与插 2 个 A，归一后是两个不同变异，
        # 各剂量 1，行号各自只记一次
        data = vcf(["chrT\t5\t.\tA\tAA,AAA\t.\t.\t.\tGT\t1/2"])
        _, calls = parse(data)
        labels = {v.label: c for v, c in calls.items()}
        self.assertEqual(labels["4:TinsTA"].dosage, 1)
        self.assertEqual(labels["4:TinsTA"].line_numbers, [3])
        ins2 = [v for v in calls if v.kind == "INS" and len(v.alt) - 1 == 2]
        self.assertEqual(len(ins2), 1)
        self.assertEqual(calls[ins2[0]].dosage, 1)

    def test_same_variant_two_heterozygous_records_merge_to_dosage_two(self):
        data = vcf([rec(4, "T", "TA", "0/1"), rec(5, "A", "AA", "0/1")])
        _, calls = parse(data)
        self.assertEqual(len(calls), 1)
        (only,) = calls.values()
        self.assertEqual(only.dosage, 2)
        self.assertEqual(only.line_numbers, [3, 4])

    def test_ref_mismatch_fails_whole_check(self):
        with self.assertRaises(CheckError) as ctx:
            parse(vcf([rec(2, "G", "T", "0/1")]))  # pos2 参考实为 C
        self.assertIn("REF", str(ctx.exception))

    def test_contradictory_dosage_in_same_file_fails(self):
        data = vcf([rec(2, "C", "G", "0/1"), rec(2, "C", "G", "1/1")])
        with self.assertRaises(CheckError) as ctx:
            parse(data)
        self.assertIn("矛盾", str(ctx.exception))

    def test_homozygous_plus_extra_record_contradiction(self):
        data = vcf([rec(2, "C", "G", "1/1"), rec(2, "C", "G", "0/1")])
        with self.assertRaises(CheckError):
            parse(data)

    def test_phased_gt_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G", "0|1")]))

    def test_missing_gt_rejected(self):
        for bad in ("./1", "0/.", "./."):
            with self.subTest(gt=bad):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", "G", bad)]))

    def test_wrong_ploidy_rejected(self):
        for bad in ("0", "0/1/1"):
            with self.subTest(gt=bad):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", "G", bad)]))

    def test_gt_index_out_of_range(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G,T", "0/3")]))

    def test_insertion_too_long_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(1, "A", "A" + "C" * 21, "0/1")]))  # 插入 21 bp

    def test_deletion_too_long_rejected(self):
        # 构造 21 bp 缺失：需要足够参考长度，直接用 classify 路径验证
        long_ref = "A" + "C" * 21
        seq2 = long_ref + "GG"
        fasta2 = b">chrL\n" + seq2.encode()
        refobj = parse_fasta(fasta2)
        data = vcf([rec(1, long_ref, "A", "0/1")])
        with self.assertRaises(CheckError):
            parse_vcf(data.encode(), "x.vcf", refobj)

    def test_mnv_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(1, "AC", "GT", "0/1")]))

    def test_symbolic_and_non_acgt_alleles_rejected(self):
        for alt in (".", "G*", "<DEL>", "g"):
            with self.subTest(alt=alt):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", alt, "0/1")]))

    def test_chrom_mismatch_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G", "0/1", chrom="chrOther")]))

    def test_pos_out_of_range_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(99, "C", "G", "0/1")]))

    def test_multi_sample_vcf_rejected(self):
        data = vcf([rec(2, "C", "G", "0/1")], extra_header_cols="\tS2")
        with self.assertRaises(CheckError):
            parse(data)

    def test_missing_header_rejected(self):
        with self.assertRaises(CheckError):
            parse(rec(2, "C", "G", "0/1") + "\n")

    def test_format_not_leading_gt_rejected(self):
        data = vcf([rec(2, "C", "G", "30", fmt="DP:GT")])
        with self.assertRaises(CheckError):
            parse(data)

    def test_format_with_gt_plus_other_fields_ok(self):
        data = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT:DP\t0/1:30"])
        _, calls = parse(data)
        self.assertEqual(len(calls), 1)

    def test_ref_only_gt_produces_no_variant(self):
        _, calls = parse(vcf([rec(2, "C", "G", "0/0")]))
        self.assertEqual(calls, {})

    def test_comment_and_blank_lines_ignored(self):
        data = vcf(["##source=unit-test", "", rec(2, "C", "G", "0/1")])
        _, calls = parse(data)
        self.assertEqual(len(calls), 1)
        (only,) = calls.values()
        self.assertEqual(only.line_numbers, [5])  # 含两条 ## 头与空行的物理行号

    def test_crlf_line_endings_keep_physical_line_numbers(self):
        data = ("##fileformat=VCFv4.2\r\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\r\n"
                + rec(2, "C", "G", "0/1") + "\r\n")
        _, calls = parse(data)
        (only,) = calls.values()
        self.assertEqual(only.line_numbers, [3])

    def test_missing_gt_in_multi_alt_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G,T", "./2")]))

    def test_phased_multi_alt_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G,T", "1|2")]))


def parse_conflict(data):
    """按冲突核查规则解析（仅 SNV/DEL），返回 calls。"""
    refobj = parse_fasta(FASTA)
    calls, _ = parse_vcf(
        data.encode(), "conflict.vcf", refobj, allowed_kinds=CONFLICT_ALLOWED_KINDS
    )
    return calls


def verdict_of(data):
    return check_conflicts(parse_conflict(data))


class ConflictParseTests(unittest.TestCase):
    def test_insertion_rejected_in_conflict_mode(self):
        with self.assertRaises(CheckError) as ctx:
            parse_conflict(vcf([rec(4, "T", "TA", "0/1")]))
        self.assertIn("插入", str(ctx.exception))

    def test_mnv_rejected_in_conflict_mode(self):
        with self.assertRaises(CheckError):
            parse_conflict(vcf([rec(1, "AC", "GT", "0/1")]))

    def test_snv_and_del_accepted(self):
        calls = parse_conflict(
            vcf([
                rec(2, "C", "G", "0/1"),
                rec(4, "TA", "T", "0/1"),
            ])
        )
        kinds = sorted(v.kind for v in calls)
        self.assertEqual(kinds, ["DEL", "SNV"])

    def test_conflict_mode_keeps_size_and_gt_rules(self):
        # 分相 GT 仍被拒绝
        with self.assertRaises(CheckError):
            parse_conflict(vcf([rec(2, "C", "G", "0|1")]))
        # REF 不符仍被拒绝
        with self.assertRaises(CheckError):
            parse_conflict(vcf([rec(2, "G", "T", "0/1")]))

    def test_default_parse_still_accepts_insertions(self):
        # 原双 VCF 路径不受 allowed_kinds 影响
        _, calls = parse(vcf([rec(4, "T", "TA", "0/1")]))
        self.assertEqual(len(calls), 1)


class ConflictCheckTests(unittest.TestCase):
    def test_affected_bases_snv_and_deletion_excludes_anchor(self):
        calls = parse_conflict(vcf([rec(4, "TA", "T", "0/1")]))
        (variant,) = calls
        self.assertEqual(affected_bases(variant), (5,))  # 锚点 pos4 不计

        calls_snv = parse_conflict(vcf([rec(2, "C", "G", "0/1")]))
        (snv,) = calls_snv
        self.assertEqual(affected_bases(snv), (2,))

    def test_two_snv_different_bases_feasible(self):
        verdict = verdict_of(
            vcf([
                rec(2, "C", "G", "0/1"),
                rec(3, "G", "T", "0/1"),
            ])
        )
        self.assertTrue(verdict.feasible)
        self.assertEqual(verdict.reason, "ok")
        self.assertEqual(verdict.witness, ())

    def test_two_snv_same_base_feasible_on_two_copies(self):
        verdict = verdict_of(
            vcf([
                rec(2, "C", "G", "0/1"),
                rec(2, "C", "T", "0/1"),
            ])
        )
        self.assertTrue(verdict.feasible)

    def test_three_snv_same_base_odd_cycle_infeasible(self):
        # 任意两条可安排，三条合起来无解（三个不同 ALT，均为剂量 1）
        verdict = verdict_of(
            vcf([
                rec(2, "C", "G", "0/1"),
                rec(2, "C", "T", "0/1"),
                rec(2, "C", "A", "0/1"),
            ])
        )
        self.assertFalse(verdict.feasible)
        self.assertEqual(verdict.reason, "odd_cycle")
        labels = sorted(item.variant.label for item in verdict.witness)
        self.assertEqual(labels, ["2:C>A", "2:C>G", "2:C>T"])
        # 证据条目带行号与原始写法
        g = next(i for i in verdict.witness if i.variant.label == "2:C>G")
        self.assertEqual(g.line_numbers, [3])
        self.assertEqual(g.dosage, 1)

    def test_multi_alt_plus_separate_record_odd_cycle(self):
        # 多 ALT 记录 1/2 给出两个剂量 1，加上第三条记录的另一个 ALT 构成三角
        verdict = verdict_of(
            vcf([
                rec(2, "C", "G,T", "1/2"),   # G、T 各剂量 1
                rec(2, "C", "A", "0/1"),     # A 剂量 1
            ])
        )
        self.assertFalse(verdict.feasible)
        self.assertEqual(verdict.reason, "odd_cycle")
        self.assertEqual(len(verdict.witness), 3)

    def test_larger_graph_witness_is_minimal_cycle(self):
        # pos2 上三个剂量 1 SNV 构成三角奇环；另加 pos3 上的无关 SNV，
        # 证据只应包含三角上的 3 个变异，不含 pos3
        verdict = verdict_of(
            vcf([
                rec(2, "C", "G", "0/1"),
                rec(2, "C", "T", "0/1"),
                rec(2, "C", "A", "0/1"),
                rec(3, "G", "T", "0/1"),
            ])
        )
        self.assertFalse(verdict.feasible)
        self.assertEqual(verdict.reason, "odd_cycle")
        self.assertEqual(len(verdict.witness), 3)
        self.assertTrue(all(item.variant.pos == 2 for item in verdict.witness))

    def test_dosage_two_with_neighbor_infeasible_hom_reason(self):
        verdict = verdict_of(
            vcf([
                rec(2, "C", "G", "1/1"),
                rec(2, "C", "T", "0/1"),
            ])
        )
        self.assertFalse(verdict.feasible)
        self.assertEqual(verdict.reason, "hom_occupied")
        self.assertEqual(len(verdict.witness), 2)
        labels = {item.variant.label for item in verdict.witness}
        self.assertEqual(labels, {"2:C>G", "2:C>T"})

    def test_dosage_two_isolated_is_feasible(self):
        verdict = verdict_of(
            vcf([
                rec(2, "C", "G", "1/1"),
                rec(3, "G", "T", "1/1"),
            ])
        )
        self.assertTrue(verdict.feasible)

    def test_deletion_anchor_snv_not_conflicting(self):
        verdict = verdict_of(
            vcf([
                rec(4, "TA", "T", "0/1"),   # 删 pos5
                rec(4, "T", "G", "0/1"),    # 锚点 pos4 替换
            ])
        )
        self.assertTrue(verdict.feasible)

    def test_deletion_conflicts_chain_and_homozygous_collision(self):
        # 缺失 4:TAA>T 删 pos5,6；归一到 5:AA>A 的单碱基缺失删 pos6：
        # 两者共享被删碱基 pos6，但剂量各 1，可分占两条拷贝
        v1 = verdict_of(vcf([rec(4, "TAA", "T", "0/1"), rec(5, "AA", "A", "0/1")]))
        self.assertTrue(v1.feasible)

        # 同一对缺失中任一变为剂量 2，即占满两条拷贝，与另一个无解
        v2 = verdict_of(vcf([rec(4, "TAA", "T", "1/1"), rec(5, "AA", "A", "0/1")]))
        self.assertFalse(v2.feasible)
        self.assertEqual(v2.reason, "hom_occupied")

    def test_deletion_triangle_odd_cycle(self):
        # D1=4:TAA>T 删 {5,6}，D2=3:GTA>G 删 {4,5}（共享被删碱基 5），
        # S=5:A>C 替换碱基 5：D1-D2、D1-S、D2-S 三边构成三角奇环
        verdict = verdict_of(
            vcf([
                rec(4, "TAA", "T", "0/1"),   # 删 5,6
                rec(3, "GTA", "G", "0/1"),   # 删 4,5
                rec(5, "A", "C", "0/1"),     # 替换 5
            ])
        )
        self.assertFalse(verdict.feasible)
        self.assertEqual(verdict.reason, "odd_cycle")
        self.assertEqual(len(verdict.witness), 3)

    def test_homozygous_deletion_collision(self):
        verdict = verdict_of(
            vcf([
                rec(4, "TA", "T", "1/1"),   # 剂量 2 缺失，占两条拷贝
                rec(5, "A", "C", "0/1"),    # 同一被删碱基上的 SNV
            ])
        )
        self.assertFalse(verdict.feasible)
        self.assertEqual(verdict.reason, "hom_occupied")

    def test_square_cycle_even_is_feasible(self):
        # 4 个变异构成偶环（K2,2 式）：两个碱基各两个 SNV，跨碱基不冲突，
        # 同碱基两两冲突——实际是两条独立的边，必然可二染色
        verdict = verdict_of(
            vcf([
                rec(2, "C", "G", "0/1"),
                rec(2, "C", "T", "0/1"),
                rec(3, "G", "A", "0/1"),
                rec(3, "G", "C", "0/1"),
            ])
        )
        self.assertTrue(verdict.feasible)

    def test_normalized_deletion_writings_merge_then_conflict(self):
        # 同一缺失的两种写法先归一合并为剂量 2，再与被删碱基上的 SNV 冲突
        verdict = verdict_of(
            vcf([
                rec(4, "TA", "T", "0/1"),
                rec(5, "AA", "A", "0/1"),   # 与上者归一为同一缺失 -> 剂量 2
                rec(5, "A", "C", "0/1"),    # pos5 SNV
            ])
        )
        self.assertFalse(verdict.feasible)
        self.assertEqual(verdict.reason, "hom_occupied")

    def test_build_items_sorted_with_bases_and_line_numbers(self):
        calls = parse_conflict(
            vcf([
                rec(3, "G", "T", "0/1"),
                rec(2, "C", "G", "1/1"),
            ])
        )
        items = build_items(calls)
        self.assertEqual([item.variant.pos for item in items], [2, 3])
        self.assertEqual(items[0].dosage, 2)
        self.assertEqual(items[0].bases, (2,))
        self.assertEqual(items[0].line_numbers, [4])

    def test_empty_calls_feasible(self):
        verdict = check_conflicts({})
        self.assertTrue(verdict.feasible)
        self.assertEqual(build_items({}), [])


class FastaTests(unittest.TestCase):
    def test_bad_base_fails(self):
        with self.assertRaises(CheckError):
            parse_fasta(b">c\nACGN\n")

    def test_missing_header_fails(self):
        with self.assertRaises(CheckError):
            parse_fasta(b"ACGT\n")

    def test_two_records_fails(self):
        with self.assertRaises(CheckError):
            parse_fasta(b">a\nAC\n>b\nGT\n")

    def test_wrapped_sequence_ok(self):
        refobj = parse_fasta(b">chrT\nACGT\nAAAA\nCGTT\n")
        self.assertEqual(refobj.sequence, SEQ)


if __name__ == "__main__":
    unittest.main()
