import os
import sys
import argparse
import polars as pl
import intervaltree
from collections import defaultdict

BED_TEST_COLS = ("chrom", "st", "end", "name")
BED_TRUTH_COLS = (
    "chrom",
    "st",
    "end",
    "name",
    "score",
    "strand",
    "tst",
    "tend",
    "item_rgb",
)
DEFAULT_READ_CSV_PARAMS = dict(
    separator="\t",
    new_columns=BED_TEST_COLS,
    has_header=False,
    truncate_ragged_lines=True,
)
# Different than simulated because Q100 also wants to fix homopolymer errors.
GOOD_MTYPES = {"correct", "good", "Hap"}


def calculate_precision_recall(
    df_bed_misassemblies: pl.DataFrame,
    df_bed_ctrl: pl.DataFrame,
    df_bed_truth: pl.DataFrame,
) -> tuple[float, float, pl.DataFrame]:
    # Create range where we liftover the control intervals to match new misassembled intervals
    itree_update_ranges: dict[str, intervaltree.IntervalTree] = {}
    for grp, df_grp in df_bed_truth.group_by("chrom"):
        grp = grp[0]
        df_grp = (
            df_grp.with_columns(length=(pl.col("end") - pl.col("st")).cum_sum())
            .select(
                st=pl.col("end"),
                end=pl.col("st").shift(-1),
                length=pl.col("length"),
                name=pl.col("name"),
            )
            .fill_null(value=sys.maxsize)
        )
        itree_update_ranges[grp] = intervaltree.IntervalTree(
            intervaltree.Interval(
                itv["st"], itv["end"] + 1, (itv["length"], itv["name"])
            )
            for itv in df_grp.iter_rows(named=True)
        )

    # To calculate the number of FPs
    # * iterate through the control bed
    # * For each misassembly, shift coordinates (misjoin +, dupe -, etc.) to match misassembly bed.
    # * Store coordinates
    # * Iterate thru misassemblies in simulated case and count number of intervals where no overlap
    itree_ctrl_misassemblies_misasim_coords = defaultdict(intervaltree.IntervalTree)
    for row in df_bed_ctrl.iter_rows(named=True):
        try:
            ovl = itree_update_ranges[row["chrom"]].overlap(row["st"], row["end"])
        except KeyError:
            ovl = None
        bp_adj = 0
        if ovl:
            final_itv: intervaltree.Interval = max(ovl, key=lambda x: x.data[0])
            bp_adj, mtype = final_itv.data
            if mtype == "misjoin":
                bp_adj = -bp_adj
            elif mtype == "false_duplication":
                bp_adj = bp_adj
            elif mtype == "inversion":
                bp_adj = 0
        new_st = max(0, row["st"] + bp_adj)
        new_end = row["end"] + bp_adj
        itree_ctrl_misassemblies_misasim_coords[row["chrom"]].addi(
            new_st, new_end, row["name"]
        )

    itree_misassemblies = defaultdict(intervaltree.IntervalTree)
    itree_new_misassemblies = defaultdict(intervaltree.IntervalTree)
    for row in df_bed_misassemblies.iter_rows(named=True):
        ovl = itree_ctrl_misassemblies_misasim_coords[row["chrom"]].overlap(
            row["st"], row["end"]
        )
        if not ovl:
            itree_new_misassemblies[row["chrom"]].addi(
                row["st"], row["end"], row["name"]
            )
        itree_misassemblies[row["chrom"]].addi(row["st"], row["end"], row["name"])

    true_positive = 0
    false_negative = 0
    itree_simulated_misassemblies = defaultdict(intervaltree.IntervalTree)
    rows = []
    for row in df_bed_truth.iter_rows(named=True):
        ovl_st, ovl_end = row["tst"] - 5, row["tend"] + 5
        ovl = itree_misassemblies[row["chrom"]].overlap(ovl_st, ovl_end)
        if ovl:
            true_positive += 1
            rows.append((row["chrom"], ovl_st, ovl_end, "true_positive"))
        else:
            false_negative += 1
            rows.append((row["chrom"], ovl_st, ovl_end, "false_negative"))
        itree_simulated_misassemblies[row["chrom"]].addi(
            row["tst"], row["tend"] + 1, row["name"]
        )

    recall = true_positive / (true_positive + false_negative)

    false_positive = 0
    for chrom, itree in itree_simulated_misassemblies.items():
        itree_chrom_new_misassemblies = itree_new_misassemblies[chrom]
        for itv in sorted(itree_chrom_new_misassemblies.iter()):
            # Is a simulated misassembly
            if itree.overlaps(itv):
                continue
            # Otherwise, is false positive
            rows.append((chrom, itv.begin, itv.end, "false_positive"))
            false_positive += 1
    precision = true_positive / (true_positive + false_positive)

    df_missing = pl.DataFrame(
        rows, orient="row", schema=["#chrom", "st", "end", "type"]
    )
    return precision, recall, df_missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-a",
        "--input_test_bed",
        required=True,
        help="Input bed with misassembly calls.",
    )
    ap.add_argument(
        "-b",
        "--input_truth_bed",
        required=True,
        help="Input bed with truth misassembly positions.",
    )
    ap.add_argument(
        "--columns_test",
        nargs="+",
        default=BED_TEST_COLS,
        help="Expected columns in misassembly calls file.",
    )
    ap.add_argument(
        "--columns_truth",
        nargs="+",
        default=BED_TRUTH_COLS,
        help="Expected columns in misassembly truth file.",
    )
    ap.add_argument(
        "--output_dir_missed_calls",
        default=None,
        help="Output directory for missing calls.",
    )
    args = ap.parse_args()

    df_test = pl.read_csv(
        args.input_test_bed,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=args.columns_test,
    ).filter(~pl.col("name").is_in(GOOD_MTYPES))
    dfs_truth = pl.read_csv(
        args.input_truth_bed,
        separator="\t",
        has_header=False,
        new_columns=args.columns_truth,
    )
    header = ["precision", "recall"]
    print("\t".join(header))

    output_dir_missed_calls = args.output_dir_missed_calls
    precision, recall, missing_calls = calculate_precision_recall(
        df_test, pl.DataFrame(schema=df_test.schema), dfs_truth
    )

    row = [precision, recall]
    if output_dir_missed_calls:
        os.makedirs(output_dir_missed_calls, exist_ok=True)
        output_path = os.path.join(output_dir_missed_calls, "missed_calls.bed")
        missing_calls.write_csv(output_path, separator="\t", include_header=True)
    print("\t".join(str(elem) for elem in row))


if __name__ == "__main__":
    raise SystemExit(main())
