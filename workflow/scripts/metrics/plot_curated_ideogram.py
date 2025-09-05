import argparse
import polars as pl
import pyideogram as pyid
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Patch

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"] + ["chrY"]
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-f", "--fp", action="store_true", help='Include "false-positives".'
    )
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
        ).with_columns(type=pl.lit("truth")),
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/nucflag_v1.0.1_missed/missed_calls.tsv",
            separator="\t",
        ),
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/inspector_v1.0.1_missed/missed_calls.tsv",
            separator="\t",
        ),
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlagPaper/results/curated/summary/flagger_v1.0.1_missed/missed_calls.tsv",
            separator="\t",
        ),
    ]
    if not args.fp:
        dfs_calls = [df.filter(pl.col("type") != "false_positive") for df in dfs_calls]
        color_key = COLOR_KEY
        color_key.pop("false_positive")
    else:
        color_key = COLOR_KEY

    tool_calls = ["truth", "nucflag", "inspector", "flagger"]

    df_fai = (
        pl.read_csv(
            "results/curated/data/asm/v1.0.1.fa.gz.fai",
            separator="\t",
            columns=[0, 1],
            new_columns=["chrom", "length"],
            has_header=False,
        )
        .with_columns(
            chrom_split=pl.col("chrom")
            .str.split_exact(by="_", n=1)
            .struct.rename_fields(["chrom_name", "hap"])
        )
        .unnest("chrom_split")
        .filter(~pl.col("chrom").is_in(["chrM", "chrEBV"]))
        .with_columns(pl.col("chrom_name").cast(pl.Enum(CHROMS)))
        .sort(by="chrom_name")
    )
    max_length = df_fai["length"].max()

    chroms = df_fai["chrom"]
    base_height_ratios = [0.1, 0.1, 0.1, 0.1, 0.5]
    num_tracks = len(base_height_ratios)
    height_ratios = base_height_ratios * len(chroms)
    height_ratios.append(0.5)
    fig, axes = plt.subplots(
        ncols=1,
        # 1 additional track for legend
        nrows=len(chroms) * num_tracks + 1,
        figsize=(40, len(chroms)),
        height_ratios=height_ratios,
    )

    for i, chrom in enumerate(chroms):
        idx_chrom = i + ((num_tracks - 1) * (i + 1))
        idx_st = idx_chrom - (num_tracks - 1)
        indices_other = range(idx_st, idx_chrom)
        ax: Axes = axes[idx_chrom]

        ax.xaxis.set_tick_params(which="both", length=0, labelleft=False)
        ax.yaxis.set_tick_params(which="both", length=0)
        # pyideogram removes xyticks
        minimize_ax(ax)

        for i, idx in enumerate(indices_other):
            ax_other: Axes = axes[idx]
            df_calls = dfs_calls[i].filter(pl.col("chrom") == chrom)
            ax_other.set_xlim(0, max_length)
            ax_other.set_ylim(0, 1)
            minimize_ax(ax_other, remove_ticks=True)

            ax_other.set_ylabel(
                tool_calls[i], rotation=0, ha="right", va="center", fontsize="x-small"
            )

            for row in df_calls.iter_rows(named=True):
                color = color_key[row["type"]]
                ax_other.axvspan(xmin=row["st"], xmax=row["end"], color=color)

        ax.set_xlim(0, max_length)

        pyid.ideogramh(chrom, cytobands, ax)

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
