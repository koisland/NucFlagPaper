import os
import sys
import csv
import glob
import argparse
from typing import Any
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
GOOD_MTYPES = {"correct", "good", "Hap"}

RGX_DOWNSAMPLE = r"(?<downsample>0\.33|0\.5)"
RGX_SM_DTYPE = r"(?<sample>[^/]*?)_(?<dtype>hifi|ont_r10|ont_r9)"
RGX_MTYPE_NUM_LEN = (
    r"(?<mtype>false_duplication|misjoin|inversion)-(?<num>\d+)-(?<len>\d+)"
)
DEFAULT_READ_CSV_PARAMS = dict(
    separator="\t",
    new_columns=BED_TEST_COLS,
    has_header=False,
    truncate_ragged_lines=True,
)


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
    for row in df_bed_ctrl.filter(~pl.col("name").is_in(GOOD_MTYPES)).iter_rows(
        named=True
    ):
        try:
            ovl = itree_update_ranges[row["chrom"]].overlap(row["st"], row["end"])
        except KeyError:
            ovl = None
        bp_adj = 0
        if ovl:
            final_itv: intervaltree.Interval = max(ovl, key=lambda x: x.data[0])
            bp_adj, mtype = final_itv.data
            # Is simulated intervals.
            if mtype == "misjoin":
                bp_adj = -bp_adj
            elif mtype == "false_duplication":
                bp_adj = bp_adj
        new_st = max(0, row["st"] + bp_adj)
        new_end = row["end"] + bp_adj
        itree_ctrl_misassemblies_misasim_coords[row["chrom"]].addi(
            new_st, new_end, row["name"]
        )

    itree_misassemblies = defaultdict(intervaltree.IntervalTree)
    itree_new_misassemblies = defaultdict(intervaltree.IntervalTree)
    for row in df_bed_misassemblies.filter(
        ~pl.col("name").is_in(GOOD_MTYPES)
    ).iter_rows(named=True):
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
    missing_rows = []
    for row in df_bed_truth.iter_rows(named=True):
        ovl_st, ovl_end = row["tst"] - 5, row["tend"] + 5
        ovl = itree_misassemblies[row["chrom"]].overlap(ovl_st, ovl_end)
        if any(ovl_itv.data not in GOOD_MTYPES for ovl_itv in ovl):
            # print(f"{row['chrom']}:{row['tst']}-{row['tend']}_{row['name']}")
            # print(ovl)
            true_positive += 1
        else:
            false_negative += 1
            missing_rows.append((row["chrom"], ovl_st, ovl_end, "false_negative"))
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
            # print(chrom, itv)
            missing_rows.append((chrom, itv.begin, itv.end, "false_positive"))
            false_positive += 1
    precision = true_positive / (true_positive + false_positive)

    df_missing = pl.DataFrame(
        missing_rows, orient="row", schema=["#chrom", "st", "end", "type"]
    )
    return precision, recall, df_missing


def read_files(fglob: str, expected_columns: tuple[str]) -> pl.DataFrame:
    dfs: list[pl.DataFrame] = []
    sniffer = csv.Sniffer()

    for file in glob.glob(fglob):
        with open(file, "rt") as fh:
            try:
                header = next(fh)
                has_header = sniffer.has_header(header)
            except StopIteration:
                has_header = False

        if has_header:
            skip_rows = 1
        else:
            skip_rows = 0

        try:
            df = (
                pl.read_csv(
                    file,
                    separator="\t",
                    skip_rows=skip_rows,
                    has_header=False,
                    new_columns=expected_columns,
                    comment_prefix="#",
                    columns=range(len(expected_columns)),
                    schema_overrides={"st": pl.Int64, "end": pl.Int64},
                    ignore_errors=True,
                )
                .with_columns(fname=pl.lit(file))
                .filter(~pl.col("st").is_null() & ~pl.col("end").is_null())
            )
            dfs.append(df)
        except pl.exceptions.NoDataError:
            print(f"No calls in {file}.", file=sys.stderr)

    df = (
        pl.concat(dfs)
        .with_columns(
            downsample=pl.col("fname").str.extract_groups(RGX_DOWNSAMPLE),
            sm_dtype=pl.col("fname").str.extract_groups(RGX_SM_DTYPE),
            mtype_num_len=pl.col("fname").str.extract_groups(RGX_MTYPE_NUM_LEN),
        )
        .unnest("downsample", "sm_dtype", "mtype_num_len")
    )
    df = df.with_columns(is_control=pl.col("mtype").is_null()).with_columns(
        # Replace misspelling
        mtype=pl.when(pl.col("mtype") == "false_dupe")
        .then(pl.lit("false_duplication"))
        .otherwise(pl.col("mtype"))
    )
    return df


def read_truth_files(
    fglob_truth: str,
    expected_columns: list[str],
) -> dict[tuple[str | Any, str | Any, str | Any], pl.DataFrame]:
    dfs = {}
    for file in glob.glob(fglob_truth):
        # seed, _ = os.path.splitext(os.path.basename(file))
        df = (
            pl.read_csv(
                file,
                separator="\t",
                has_header=False,
                new_columns=expected_columns,
                comment_prefix="#",
            )
            .with_columns(
                # Replace misspelling
                name=pl.when(pl.col("name") == "false_dupe")
                .then(pl.lit("false_duplication"))
                .otherwise(pl.col("name"))
            )
            .filter(~pl.col("name").is_in(GOOD_MTYPES))
        )
        unique_mtypes = df["name"].unique()
        # Per chrom
        numbers = df.group_by(["chrom"]).agg(pl.col("name").len())["name"].mode()
        lengths = (
            df.with_columns(length=pl.col("end") - pl.col("st"))
            .group_by(["name"])
            .agg(pl.col("length").unique())
            .explode("length")["length"]
        )

        mtype = ",".join(unique_mtypes)
        num = ",".join(str(n) for n in numbers)
        length = ",".join(str(lt) for lt in lengths)

        df = df.with_columns(
            mtype=pl.lit(mtype),
            num=pl.lit(num),
            len=pl.lit(length),
        )
        dfs[(mtype, num, length)] = df
    return dfs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-a",
        "--input_test_dir",
        required=True,
        help="Input dir with misassembly and control calls.",
    )
    ap.add_argument("--glob_test", required=True, help="File glob to find above calls.")
    ap.add_argument(
        "-b",
        "--input_truth_dir",
        required=True,
        help="Input dir with truth misassembly positions.",
    )
    ap.add_argument(
        "--glob_truth",
        required=True,
        help="File glob to find with truth misassemblies.",
    )
    ap.add_argument("-d", "--dtype", required=True, help="Data type to filter for.")
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
    dtype = args.dtype

    df_res = read_files(
        fglob=os.path.join(args.input_test_dir, args.glob_test),
        expected_columns=args.columns_test,
    ).filter(pl.col("dtype") == args.dtype)
    dfs_truth = read_truth_files(
        fglob_truth=os.path.join(args.input_truth_dir, args.glob_truth),
        expected_columns=args.columns_truth,
    )

    header = [
        "sample",
        "downsample",
        "dtype",
        "mtype",
        "num",
        "len",
        "precision",
        "recall",
    ]
    print("\t".join(header))

    output_dir_missed_calls = args.output_dir_missed_calls
    os.makedirs(output_dir_missed_calls, exist_ok=True)

    all_groups = df_res["sample", "downsample", "dtype"].unique()
    for group in all_groups.iter_rows():
        sample, downsample, dtype = group
        if not downsample:
            is_downsampled = pl.col("downsample").is_null()
        else:
            is_downsampled = pl.col("downsample") == downsample

        dfs_group_test = df_res.filter(
            (pl.col("sample") == sample)
            & is_downsampled
            & (pl.col("dtype") == dtype)
            & ~pl.col("is_control")
        ).partition_by(["mtype", "num", "len"], as_dict=True)

        df_group_control = df_res.filter(
            (pl.col("sample") == sample)
            & is_downsampled
            & (pl.col("dtype") == dtype)
            & pl.col("is_control")
        )
        for mtype_group, df_groups_test_mtypes in dfs_group_test.items():
            df_mtype_truth = dfs_truth.get(mtype_group)
            if not isinstance(df_mtype_truth, pl.DataFrame):
                print(f"Skipping {mtype_group}...", file=sys.stderr)
                continue
            precision, recall, missing_calls = calculate_precision_recall(
                df_groups_test_mtypes, df_group_control, df_mtype_truth
            )

            row = [*group, *mtype_group, precision, recall]
            if output_dir_missed_calls:
                output_path = os.path.join(
                    output_dir_missed_calls,
                    "_".join(e if e else "None" for e in [*group, *mtype_group])
                    + ".bed",
                )
                missing_calls.write_csv(
                    output_path, separator="\t", include_header=True
                )
            print("\t".join(str(elem) for elem in row))


if __name__ == "__main__":
    raise SystemExit(main())
