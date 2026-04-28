import sys
import json
import argparse
import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from typing import TextIO
from intervaltree import Interval, IntervalTree
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.patches import Rectangle
from collections import defaultdict, Counter
from upsetplot import UpSet, from_indicators

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

plt.rcParams["font.family"] = "Arial"


def minimalize_ax(ax: Axes, *, remove_ticks: bool = False) -> None:
    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        ax.tick_params(
            axis="both",
            left=False,
            top=False,
            right=False,
            bottom=False,
            labelleft=False,
            labeltop=False,
            labelright=False,
            labelbottom=False,
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

    rows_ovl_counts: list[dict[str, bool | int]] = []
    for ovl_name, cnt in ovl_counts.items():
        row = {tool: False for tool in calls.keys()}
        for tool in ovl_name.split("-"):
            row[tool] = True
        # Rename
        row = {RENAME_TOOLS[tool]: value for tool, value in row.items()}
        row["value"] = cnt
        rows_ovl_counts.append(row)

    df_pl_ovl_counts = pl.DataFrame(rows_ovl_counts, orient="row").sort(
        descending=True, by="value"
    )
    # Reorder.
    df_ovl_counts = df_pl_ovl_counts.to_pandas().fillna(value=np.nan)[
        [*COLORS.keys(), "value"]
    ]

    fig, ax = plt.subplots(figsize=(16, 4), dpi=600, layout="tight")
    minimalize_ax(ax, remove_ticks=True)

    df_ovl_counts = from_indicators(list(COLORS.keys()), data=df_ovl_counts)
    upset = UpSet(
        data=df_ovl_counts,
        show_counts=True,
        sum_over="value",
        sort_by="cardinality",
        sort_categories_by="input",
    )
    for lbl, color in COLORS.items():
        upset.style_categories(
            categories=[lbl],
            bar_facecolor=color,
            shading_facecolor=to_rgba(color, alpha=0.2),
        )

    subplots = upset.plot(fig=fig)

    # Add stacked bar. We cannot use existing implementation because we don't have a grouping var.
    ax_intersections: Axes = subplots["intersections"]
    intersections = sorted(rows_ovl_counts, key=lambda x: x["value"], reverse=True)
    x_coords = []
    widths = []
    heights = defaultdict(list)
    for container in ax_intersections.containers:
        for intersection, ptch in zip(intersections, container.patches):
            ptch: Rectangle
            x_coords.append(ptch.get_x())
            widths.append(ptch.get_width())

            intersection_size = sum(
                1 for tool in COLORS.keys() if intersection.get(tool)
            )
            # Divide height by intersection size so even height prop
            for tool in COLORS.keys():
                heights[tool].append(
                    ptch.get_height() / intersection_size
                    if intersection.get(tool)
                    else 0
                )
            ptch.remove()

    bottom = np.zeros(len(intersections))
    for tool, color in COLORS.items():
        heights_tool = heights[tool]
        ax_intersections.bar(
            x=x_coords,
            height=heights_tool,
            width=widths,
            bottom=bottom,
            color=color,
            align="edge",
        )
        bottom += heights_tool

    fig.savefig(f"{output_prefix}_venn.png", bbox_inches="tight")
    fig.savefig(f"{output_prefix}_venn.pdf", bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
