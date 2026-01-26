import json
import argparse
import sys
import matplotlib.colors

import polars as pl
import matplotlib.pyplot as plt

from typing import TextIO
from intervaltree import Interval, IntervalTree
from collections import defaultdict, Counter
from supervenn import supervenn, make_sets_from_chunk_sizes


RENAME_TOOLS = {
    "nucflag_no_homopolymers": "NucFlag v1.0 (No homopolymers)",
    "nucflag": "NucFlag v1.0",
    "deepvariant": "DeepVariant v1.9 (FILTER=='PASS')",
    "flagger": "HMM-Flagger v1.1.0",
    "inspector": "Inspector v1.3",
}
COLORS = dict(
    zip(
        RENAME_TOOLS.values(),
        ["blue", "purple", "maroon", "magenta", "teal"],
        strict=True,
    )
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", type=str, help="JSON string of call maps.")
    ap.add_argument(
        "-o", "--output_prefix", type=str, default="out", help="Output prefix."
    )

    args = ap.parse_args()
    output_prefix = args.output_prefix

    calls: dict[str, str] = json.loads(args.input)

    # {chrom: IntervalTree[Interval[start, end, name]]}
    itree_calls: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
    for tool, file in calls.items():
        df_call = pl.read_csv(
            file,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            columns=list(range(0, 4)),
            schema={
                "#chrom": pl.String,
                "chromStart": pl.UInt64,
                "chromEnd": pl.UInt64,
                "name": pl.String,
            },
            truncate_ragged_lines=True,
        )
        df_call = df_call.filter(~pl.col("name").is_in(["hap", "Hap", "correct"]))
        for call in df_call.iter_rows(named=True):
            itree_calls[call["#chrom"]].add(
                Interval(
                    call["chromStart"],
                    call["chromEnd"],
                    (call["name"], tool),
                )
            )
        if tool == "nucflag":
            df_homopolymers_only = df_call.filter(
                ~pl.col("name").eq("homopolymer")
            ).with_columns(len=pl.col("chromEnd") - pl.col("chromStart"))
            print(
                f"Homopolymer length total: {df_homopolymers_only['len'].sum()}",
                file=sys.stderr,
            )
            for call in df_call.filter(~pl.col("name").eq("homopolymer")).iter_rows(
                named=True
            ):
                itree_calls[call["#chrom"]].add(
                    Interval(
                        call["chromStart"],
                        call["chromEnd"],
                        (call["name"], f"{tool}_no_homopolymers"),
                    )
                )

    ovl_counts: Counter[str] = Counter()
    all_ovls = set()
    file_handles: dict[str, TextIO] = {}
    for chrom, itrees in itree_calls.items():
        itv: Interval
        ovl: set[Interval]
        for itv in itrees.iter():
            ovl = itrees.overlap(itv)
            # Sort all intervals into immutable tuple
            # Check if used combination.
            ovl_sorted: tuple[Interval, ...] = tuple(sorted(set(itv for itv in ovl)))
            if ovl_sorted in all_ovls:
                continue
            # Then create overlap name.
            ovl_name: str = "-".join(sorted(set(itv.data[1] for itv in ovl_sorted)))
            if file_handles.get(ovl_name):
                fh = file_handles[ovl_name]
                print(
                    chrom,
                    *[
                        f"{oitv.begin}\t{oitv.end}\t{oitv.data[0]}\t{oitv.data[1]}"
                        for oitv in ovl
                    ],
                    sep="\t",
                    file=fh,
                )
            else:
                file_handles[ovl_name] = open(f"{output_prefix}_{ovl_name}.bed", "wt")
            ovl_counts[ovl_name] += 1
            all_ovls.add(ovl_sorted)

    for fh in file_handles.values():
        fh.close()

    rows_ovl_counts = []
    for ovl_name, cnt in ovl_counts.items():
        row = {tool: False for tool in calls.keys()}
        for tool in ovl_name.split("-"):
            row[tool] = True
        # Rename
        row = {RENAME_TOOLS[tool]: value for tool, value in row.items()}
        row["size"] = cnt
        rows_ovl_counts.append(row)

    df_ovl_counts = pl.DataFrame(rows_ovl_counts, orient="row").to_pandas()
    # Reorder.
    df_ovl_counts = df_ovl_counts[[*COLORS.keys(), "size"]]
    fig, ax = plt.subplots(figsize=(16, 8), dpi=600, layout="constrained")
    sets, labels = make_sets_from_chunk_sizes(df_ovl_counts)
    colors = [matplotlib.colors.to_rgba(COLORS[label], alpha=0.5) for label in labels]
    supervenn(
        sets,
        labels,
        ax=ax,
        min_width_for_annotation=1200,
        fontsize=12,
        color_cycle=colors,
    )
    fig.savefig(f"{output_prefix}_venn.png", bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
