import sys
import polars as pl
import intervaltree as it

from typing import Any
from collections import defaultdict


def main():
    infile = sys.argv[1]
    bp_merge = int(sys.argv[2])
    if infile == "-":
        infile = sys.stdin

    df_infile = pl.read_csv(
        infile,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "st",
            "end",
            "rchrom",
            "rst",
            "rend",
            "itv",
            "ignore",
            "strand",
            "rstrand",
        ],
    )

    itrees = defaultdict(it.IntervalTree)

    for row in df_infile.iter_rows(named=True):
        itrees[row["chrom"]].add(
            it.Interval(row["st"] - bp_merge, row["end"] + bp_merge, row)
        )

    def reduce_data(i1: dict[str, Any], i2: dict[str, Any]):
        length_i1 = i1["end"] - i1["st"]
        length_i2 = i2["end"] - i2["st"]

        return {
            "chrom": i1["chrom"],
            "st": min(i1["st"], i2["st"]),
            "end": max(i1["end"], i2["end"]),
            "rchrom": i1["rchrom"],
            "rst": min(i1["rst"], i2["rst"]),
            "rend": max(i1["rend"], i2["rend"]),
            "itv": i1["chrom"],
            "ignore": i1["ignore"],
            "strand": i1["strand"] if length_i1 > length_i2 else i2["strand"],
            "rstrand": i1["rstrand"],
        }

    for chrom, itree in itrees.items():
        itree.merge_overlaps(data_reducer=reduce_data)
        for st, end, data in sorted(itree.iter()):
            print(
                chrom,
                st + bp_merge,
                end - bp_merge,
                data["rchrom"],
                data["rst"],
                data["rend"],
                data["itv"],
                data["ignore"],
                data["strand"],
                data["rstrand"],
                sep="\t",
            )


if __name__ == "__main__":
    raise SystemExit(main())
