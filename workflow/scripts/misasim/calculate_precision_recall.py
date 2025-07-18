import os
import sys
import argparse
import polars as pl
import intervaltree
from collections import defaultdict

BED9_COLS = ["chrom", "st", "end", "name", "score", "strand", "tst", "tend", "item_rgb"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--bed_misassemblies", required=True)
    ap.add_argument("-c", "--bed_ctrl", required=True)
    ap.add_argument("-s", "--bed_simulated", required=True, help="Regions where misassembly simulated.")

    args = ap.parse_args()

    df_bed_misassemblies = pl.read_csv(args.bed_misassemblies, separator="\t", new_columns=BED9_COLS, has_header=False)
    df_bed_ctrl = pl.read_csv(args.bed_ctrl, separator="\t", new_columns=BED9_COLS, has_header=False)
    df_bed_simulated = (
        pl.read_csv(args.bed_simulated, separator="\t", new_columns=BED9_COLS, has_header=False)
        .filter(pl.col("name") != "good")
    )

    # Create range where we liftover the control intervals to match new misassembled intervals
    itree_update_ranges: dict[str, intervaltree.IntervalTree] = {}
    for grp, df_grp in df_bed_simulated.group_by("chrom"):
        grp = grp[0]
        df_grp = (
            df_grp
            .with_columns(length=(pl.col("end")-pl.col("st")).cum_sum())
            .select(st=pl.col("end"), end=pl.col("st").shift(-1), length=pl.col("length"), name=pl.col("name"))
            .fill_null(value=sys.maxsize)
        )
        itree_update_ranges[grp] = intervaltree.IntervalTree(
            intervaltree.Interval(itv["st"], itv["end"] + 1, (itv["length"], itv["name"]))
            for itv in df_grp.iter_rows(named=True)
        )

    # To calculate the number of FPs
    # * iterate through the control bed
    # * For each misassembly, shift coordinates (misjoin +, dupe -, etc.) to match misassembly bed.
    # * Store coordinates
    # * Iterate thru misassemblies in simulated case and count number of intervals where no overlap
    itree_ctrl_misassemblies_misasim_coords = defaultdict(intervaltree.IntervalTree)
    for row in df_bed_ctrl.filter(pl.col("name") != "good").iter_rows(named=True):
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
            elif mtype == "false_dupe":
                bp_adj = bp_adj
        new_st = max(0, row["st"] + bp_adj)
        new_end = row["end"] + bp_adj
        itree_ctrl_misassemblies_misasim_coords[row["chrom"]].addi(new_st, new_end, row["name"])

    itree_misassemblies = defaultdict(intervaltree.IntervalTree)
    itree_new_misassemblies = defaultdict(intervaltree.IntervalTree)
    for row in df_bed_misassemblies.filter(pl.col("name") != "good").iter_rows(named=True):
        ovl = itree_ctrl_misassemblies_misasim_coords[row["chrom"]].overlap(row["st"], row["end"])
        if not ovl:
            itree_new_misassemblies[row["chrom"]].addi(row["st"], row["end"], row["name"])
        itree_misassemblies[row["chrom"]].addi(row["st"], row["end"], row["name"])

    true_positive = 0
    false_negative = 0
    itree_simulated_misassemblies = defaultdict(intervaltree.IntervalTree)
    for row in df_bed_simulated.iter_rows(named=True):
        ovl_st, ovl_end = row["tst"] - 5, row["tend"] + 5
        ovl = itree_misassemblies[row["chrom"]].overlap(ovl_st, ovl_end)
        if any(ovl_itv.data != "good" for ovl_itv in ovl):
            # print(f"{row['chrom']}:{row['tst']}-{row['tend']}_{row['name']}")
            # print(ovl)
            true_positive += 1
        else:
            false_negative += 1
        itree_simulated_misassemblies[row["chrom"]].addi(row["tst"], row["tend"] + 1, row["name"])

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
            false_positive += 1
    precision = true_positive / (true_positive + false_positive)
    print(f"{precision}\t{recall}")

 
if __name__ == "__main__":
    raise SystemExit(main())
