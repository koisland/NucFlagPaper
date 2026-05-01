from matplotlib.patches import Patch
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
    "nucflag": "NucFlag v1.0",
    "deepvariant": "DeepVariant v1.9 (FILTER=='PASS')",
    "flagger": "HMM-Flagger v1.1.0",
    "inspector": "Inspector v1.3",
}
COLORS = dict(
    zip(
        RENAME_TOOLS.values(),
        ["purple", "maroon", "magenta", "teal"],
        strict=True,
    )
)
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
    alignment="left",
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

    ovl_counts: Counter[str] = Counter()
    all_ovls = set()
    file_handles: dict[str, TextIO] = {}
    nucflag_non_homopolymer_ovl_name_counts: Counter[tuple[str, ...]] = Counter()
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
            ovl_name = sorted(set(itv.data[1] for itv in ovl_sorted))
            ovl_name_joined: str = "-".join(ovl_name)
            if file_handles.get(ovl_name_joined):
                fh = file_handles[ovl_name_joined]
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
                file_handles[ovl_name_joined] = open(
                    f"{output_prefix}_{ovl_name_joined}.bed", "wt"
                )
            # Count when non-homopolymer
            if itv.data[1] == "nucflag" and itv.data[0] != "homopolymer":
                nucflag_non_homopolymer_ovl_name_counts[
                    tuple(RENAME_TOOLS[tool] for tool in ovl_name)
                ] += 1

            ovl_counts[ovl_name_joined] += 1
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
        sort_by="-degree",
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
    # Jank but use height of bar to determine which intersection
    height_to_intersection_map = {it["value"]: it for it in rows_ovl_counts}
    x_coords = []
    widths = []
    heights = defaultdict(list)
    all_intersections = []
    non_homopolymer_heights = []
    for container in ax_intersections.containers:
        for ptch in container.patches:
            ptch: Rectangle
            height = ptch.get_height()
            x_coords.append(ptch.get_x())
            widths.append(ptch.get_width())
            intersection = height_to_intersection_map[height]
            non_homopolymer_height = nucflag_non_homopolymer_ovl_name_counts.get(
                tuple(sorted(k for k, v in intersection.items() if v and k != "value")),
                0,
            )
            non_homopolymer_heights.append(non_homopolymer_height)
            intersection_size = sum(
                1 for tool in COLORS.keys() if intersection.get(tool)
            )
            # Divide height by intersection size so even height prop
            for tool in COLORS.keys():
                heights[tool].append(
                    height / intersection_size if intersection.get(tool) else 0
                )
            ptch.remove()
            all_intersections.append(intersection)

    bottom = np.zeros(len(rows_ovl_counts))
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

    # Add final hatched bar
    non_homopolymer_bars = ax_intersections.bar(
        x=x_coords,
        height=non_homopolymer_heights,
        width=widths,
        color=COLORS["NucFlag v1.0"],
        hatch="///",
        edgecolor="black",
        align="edge",
    )
    non_homopolymer_bar_labels = [
        cnt.get_height() if cnt.get_height() else "" for cnt in non_homopolymer_bars
    ]
    ax_intersections.bar_label(non_homopolymer_bars, non_homopolymer_bar_labels)
    ax_intersections.set_ylabel("# of error calls")
    ax_intersections.legend(
        loc="upper left",
        handles=[Patch(edgecolor="black", fill=False, hatch="////")],
        labels=["Non-homopolymers"],
        **LEGEND_KWARGS,
    )
    fig.savefig(f"{output_prefix}_venn.png", bbox_inches="tight")
    fig.savefig(f"{output_prefix}_venn.pdf", bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
