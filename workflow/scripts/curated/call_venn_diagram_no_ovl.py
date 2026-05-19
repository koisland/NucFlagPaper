import json
import argparse
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from typing import TextIO
from intervaltree import Interval, IntervalTree
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
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
        df_call = df_call.filter(
            ~pl.col("name").is_in(["het_or_mismap", "Hap", "correct"])
        )
        for call in df_call.iter_rows(named=True):
            itree_calls[call["#chrom"]].add(
                Interval(
                    call["chromStart"],
                    call["chromEnd"],
                    (call["name"], tool),
                )
            )

    ovl_counts: Counter[str] = Counter()
    file_handles: dict[str, TextIO] = {}
    nucflag_homopolymer_ovl_name_counts: Counter[tuple[str, ...]] = Counter()
    for chrom, itrees in itree_calls.items():
        itv: Interval
        ovl: set[Interval]
        ignore_itvs = set()
        for itv in itrees.iter():
            if itv in ignore_itvs:
                continue
            ovl = itrees.overlap(itv)
            # Check if used itv.
            filt_ovl: set[Interval] = set(itv for itv in ovl if itv not in ignore_itvs)
            if not filt_ovl:
                continue
            ignore_itvs.update(filt_ovl)

            # Then create overlap name.
            ovl_name = sorted(set(itv.data[1] for itv in filt_ovl))
            ovl_name_joined: str = "-".join(ovl_name)
            if not file_handles.get(ovl_name_joined):
                file_handles[ovl_name_joined] = open(
                    f"{output_prefix}_{ovl_name_joined}.bed", "wt"
                )
            fh = file_handles[ovl_name_joined]
            print(
                chrom,
                *[
                    f"{oitv.begin}\t{oitv.end}\t{oitv.data[0]}\t{oitv.data[1]}"
                    for oitv in filt_ovl
                ],
                sep="\t",
                file=fh,
            )
            # Count when homopolymer
            if itv.data[1] == "nucflag" and itv.data[0] == "homopolymer":
                nucflag_homopolymer_ovl_name_counts[
                    tuple(RENAME_TOOLS[tool] for tool in ovl_name)
                ] += 1

            ovl_counts[ovl_name_joined] += 1

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

    fig, ax = plt.subplots(figsize=(30, 4), dpi=600, layout="tight")
    minimalize_ax(ax, remove_ticks=True)

    df_ovl_counts = from_indicators(list(COLORS.keys()), data=df_ovl_counts)
    upset = UpSet(
        data=df_ovl_counts,
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

    # Redraw labels on totals since we need to draw intersections manually
    ax_totals: Axes = subplots["totals"]
    total_labels = []
    total_bar_cont = ax_totals.containers[0]
    for ptch in total_bar_cont.patches:
        width = ptch.get_width()
        total_labels.append(int(width))
    ax_totals.bar_label(
        total_bar_cont,
        total_labels,
        fontsize=8,
    )

    # Add stacked bar. We cannot use existing implementation because we don't have a grouping var.
    ax_intersections: Axes = subplots["intersections"]
    # Jank but use height of bar to determine which intersection
    height_to_intersection_map = {it["value"]: it for it in rows_ovl_counts}
    x_coords = []
    widths = []
    heights = defaultdict(list)
    all_intersections = []
    homopolymer_heights = []
    for container in ax_intersections.containers:
        for ptch in container.patches:
            ptch: Rectangle
            width = ptch.get_width()
            height = ptch.get_height()
            # Offset by half width to fit both bars
            x_coords.append(ptch.get_x() - (width / 2))
            widths.append(width)
            intersection = height_to_intersection_map[height]
            homopolymer_height = nucflag_homopolymer_ovl_name_counts.get(
                tuple(sorted(k for k, v in intersection.items() if v and k != "value")),
                0,
            )
            homopolymer_heights.append(homopolymer_height)
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

    x_coords = np.array(x_coords)
    bottom = np.zeros(len(rows_ovl_counts))
    for tool, color in COLORS.items():
        heights_tool = heights[tool]
        bars = ax_intersections.bar(
            x=x_coords,
            height=heights_tool,
            width=widths,
            bottom=bottom,
            color=color,
            align="edge",
        )
        bottom += heights_tool
        # Only show label if top of bar
        bar_labels = []
        for cnt in bars.patches:
            height = cnt.get_height()
            if height == 0:
                label = ""
            else:
                height = int(cnt.get_y() + height)
                label = height if height in height_to_intersection_map else ""
            bar_labels.append(label)

        ax_intersections.bar_label(
            bars,
            bar_labels,
            fontsize=7.5,
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
        )

    # Add final hatched bar
    homopolymer_bars = ax_intersections.bar(
        # Offset by half width
        x=x_coords + widths[0],
        height=homopolymer_heights,
        width=widths,
        color=COLORS["NucFlag v1.0"],
        hatch="///",
        edgecolor="black",
        align="edge",
    )
    homopolymer_bar_labels = [
        int(cnt.get_height()) if cnt.get_height() else "" for cnt in homopolymer_bars
    ]
    ax_intersections.bar_label(
        homopolymer_bars,
        homopolymer_bar_labels,
        fontsize=7.5,
        path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
    )
    ax_intersections.set_ylabel("# of error calls")
    ax_intersections.legend(
        loc="upper left",
        handles=[
            Patch(edgecolor="black", fill=False),
            Patch(edgecolor="black", fill=False, hatch="////"),
        ],
        labels=["Total calls", "Homopolymer calls"],
        **LEGEND_KWARGS,
    )
    fig.savefig(f"{output_prefix}_venn.png", bbox_inches="tight")
    fig.savefig(f"{output_prefix}_venn.pdf", bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
