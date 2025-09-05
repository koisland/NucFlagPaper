import argparse
from matplotlib.colors import rgb2hex
import polars as pl
import pyideogram as pyid
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Patch

COLOR_KEY = {
    "truth": "black",
    "false_positive": "red",
    "true_positive": "blue",
    "false_negative": "orange",
}


def minimize_ax(ax: Axes, *, remove_ticks: bool = False):
    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        ax.set_xticks([], [])
        ax.set_yticks([], [])


def rgb_to_hex(srs: pl.Series) -> pl.Series:
    color_hex = []
    for elem in srs:
        if elem.startswith("#"):
            color_hex.append(elem)
        else:
            rgb = tuple(int(e) / 255 for e in elem.split(","))
            assert len(rgb) == 3, f"Invalid item_rgb format for {rgb}"
            color_hex.append(rgb2hex(rgb))
    return pl.Series(name="color", values=color_hex)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--outfile", default="out.png", help="Output file.")
    args = ap.parse_args()

    cytobands = pyid.dataloader.load_cytobands(
        "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/cytoBand.hg002v1.0.sorted.bed"
    )
    dfs_calls = [
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/data/v1.0.1_truth.bed",
            separator="\t",
            columns=[0, 1, 2, 3],
            new_columns=["chrom", "st", "end", "patch"],
        )
        .with_columns(type=pl.lit("truth"))
        .filter(pl.col("chrom") == "chrY_PATERNAL"),
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/nucflag_v1.0.1_missed/missed_calls.tsv",
            separator="\t",
        ).filter(pl.col("chrom") == "chrY_PATERNAL"),
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/inspector_v1.0.1_missed/missed_calls.tsv",
            separator="\t",
        ).filter(pl.col("chrom") == "chrY_PATERNAL"),
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/flagger_v1.0.1_missed/missed_calls.tsv",
            separator="\t",
        ).filter(pl.col("chrom") == "chrY_PATERNAL"),
    ]

    df_segdup = (
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/HG002.SDs.010624.45col.bed",
            separator="\t",
            has_header=False,
            columns=range(0, 10),
            new_columns=[
                "chrom",
                "st",
                "end",
                "name",
                "score",
                "strand",
                "tst",
                "tend",
                "item_rgb",
            ],
        )
        .with_columns(
            pl.col("item_rgb").map_batches(rgb_to_hex, return_dtype=pl.String)
        )
        .filter(pl.col("chrom") == "chrY_PATERNAL")
    )
    df_censat = (
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/hg002v1.0.1.cenSatv2.0.noheader.bed",
            separator="\t",
            has_header=False,
            new_columns=[
                "chrom",
                "st",
                "end",
                "name",
                "score",
                "strand",
                "tst",
                "tend",
                "item_rgb",
            ],
        )
        .with_columns(
            pl.col("item_rgb").map_batches(rgb_to_hex, return_dtype=pl.String)
        )
        .filter(pl.col("chrom") == "chrY_PATERNAL")
    )

    color_key = COLOR_KEY
    tool_calls = ["truth", "nucflag", "inspector", "flagger"]
    height_ratios = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.5]
    fig, axes = plt.subplots(
        ncols=1,
        nrows=len(height_ratios),
        figsize=(20, 4),
        height_ratios=height_ratios,
    )

    indices_calls = range(0, 4)
    idx_censat = 4
    idx_segdup = 5
    idx_chrom = 6
    ax: Axes = axes[idx_chrom]

    ax.xaxis.set_tick_params(which="both", length=0, labelleft=False)
    ax.yaxis.set_tick_params(which="both", length=0)
    # pyideogram removes xyticks
    minimize_ax(ax)

    # Add sequence context.
    ax_censat: Axes = axes[idx_censat]
    ax_segdup: Axes = axes[idx_segdup]
    minimize_ax(ax_censat, remove_ticks=True)
    minimize_ax(ax_segdup, remove_ticks=True)
    ax_censat.set_ylabel(
        "censat", rotation=0, ha="right", va="center", fontsize="x-small"
    )
    ax_segdup.set_ylabel(
        "segdup", rotation=0, ha="right", va="center", fontsize="x-small"
    )

    for row in df_censat.iter_rows(named=True):
        ax_censat.axvspan(xmin=row["st"], xmax=row["end"], color=row["item_rgb"])
    for row in df_segdup.iter_rows(named=True):
        ax_segdup.axvspan(xmin=row["st"], xmax=row["end"], color=row["item_rgb"])

    # Add calls
    for i, idx in enumerate(indices_calls):
        ax_other: Axes = axes[idx]
        df_calls = dfs_calls[i]
        ax_other.set_ylim(0, 1)
        minimize_ax(ax_other, remove_ticks=True)

        ax_other.set_ylabel(
            tool_calls[i], rotation=0, ha="right", va="center", fontsize="x-small"
        )

        for row in df_calls.iter_rows(named=True):
            color = color_key[row["type"]]
            ax_other.axvspan(xmin=row["st"], xmax=row["end"], color=color)

    pyid.ideogramh("chrY_PATERNAL", cytobands, ax)

    ax_legend: Axes = axes[-1]
    minimize_ax(ax_legend, remove_ticks=True)
    ax_legend.legend(
        handles=[Patch(facecolor=color, label=lbl) for lbl, color in color_key.items()],
        bbox_to_anchor=(0.5, 0.5),
        loc="center",
        ncol=len(color_key),
    )
    fig.savefig(args.outfile, bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
