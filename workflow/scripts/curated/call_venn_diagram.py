import json
import argparse

import polars as pl
import matplotlib.pyplot as plt

from intervaltree import Interval, IntervalTree
from collections import defaultdict, Counter
from supervenn import supervenn, make_sets_from_chunk_sizes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", type=str, help="JSON string of call maps.")
    ap.add_argument("-o", "--output", type=str, help="Output plot.")

    args = ap.parse_args()

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

    ovl_counts: Counter[str] = Counter()
    all_ovls = set()
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
            ovl_counts[ovl_name] += 1
            all_ovls.add(ovl_sorted)

    rows_ovl_counts = []
    for ovl_name, cnt in ovl_counts.items():
        row = {tool: False for tool in calls.keys()}
        for tool in ovl_name.split("-"):
            row[tool] = True
        row["size"] = cnt
        rows_ovl_counts.append(row)

    df_ovl_counts = pl.DataFrame(rows_ovl_counts, orient="row").to_pandas()

    fig, ax = plt.subplots(figsize=(16, 8), dpi=600, layout="constrained")
    sets, labels = make_sets_from_chunk_sizes(df_ovl_counts)
    supervenn(sets, labels, ax=ax)
    fig.savefig(args.output, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
