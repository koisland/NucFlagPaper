import argparse
import polars as pl
import pyideogram as pyid
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.patches import Patch

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrXY"]
CALL_COLOR_KEY = {
    "truth": "black",
    "false_positive": "gray",
    "true_positive": "blue",
    "false_negative": "red",
}
CALL_NAMES = {
    "truth": "Truth",
    "false_positive": "False Positive",
    "true_positive": "True Positive",
    "false_negative": "False Negative",
}
BAND_COLORS = pyid.BANDCOL | {"none": (1.0, 1.0, 1.0)}
LBL_KWARGS = dict(rotation=0, ha="right", va="center")
TOOL_NAMES = ["Truth", "NucFlag v1.0", "Inspector v1.3", "HMM-Flagger v1.1.0"]
TOOL_COLORS = ["black", "purple", "teal", "magenta"]


def minimize_ax(ax: Axes, *, remove_ticks: bool = False):
    for spine in ["left", "right", "bottom", "top"]:
        ax.spines[spine].set_visible(False)
    if remove_ticks:
        ax.set_xticks([], [])
        ax.set_yticks([], [])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cytobands")
    ap.add_argument("--fai")
    ap.add_argument("--truth")
    ap.add_argument("--nucflag")
    ap.add_argument("--flagger")
    ap.add_argument("--inspector")
    ap.add_argument("--segdups")
    ap.add_argument("--censat")
    ap.add_argument(
        "-f", "--fp", action="store_true", help='Include "false-positives".'
    )
    ap.add_argument("-o", "--output_prefix", default="ideogram", help="Output file.")
    args = ap.parse_args()

    cytobands = pyid.dataloader.load_cytobands(args.cytobands)
    missed_calls_kwargs = dict(
        comment_prefix="#",
        separator="\t",
        has_header=False,
        new_columns=["chrom", "st", "end", "type"],
    )
    dfs_calls = [
        pl.read_csv(
            args.truth,
            separator="\t",
            columns=[0, 1, 2, 3],
            new_columns=["chrom", "st", "end", "patch"],
        ).with_columns(type=pl.lit("truth")),
        pl.read_csv(args.nucflag, **missed_calls_kwargs),
        pl.read_csv(args.inspector, **missed_calls_kwargs),
        pl.read_csv(args.flagger, **missed_calls_kwargs),
    ]

    if not args.fp:
        dfs_calls = [
            df.filter(pl.col("type") != "false_positive").with_columns(
                pl.col("type").cast(pl.Enum(CALL_COLOR_KEY.keys()))
            )
            for df in dfs_calls
        ]
        color_key = CALL_COLOR_KEY
        color_key.pop("false_positive")
    else:
        dfs_calls = [
            df.with_columns(pl.col("type").cast(pl.Enum(CALL_COLOR_KEY.keys())))
            for df in dfs_calls
        ]
        color_key = CALL_COLOR_KEY

    df_fai = (
        pl.read_csv(
            args.fai,
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
        .with_columns(
            chrom_name=pl.when(pl.col("chrom_name").is_in(["chrX", "chrY"]))
            .then(pl.lit("chrXY"))
            .otherwise(pl.col("chrom_name"))
        )
        .with_columns(pl.col("chrom_name").cast(pl.Enum(CHROMS)))
        .sort(by="chrom_name")
    )
    max_length = df_fai["length"].max()

    # Draw overview
    outfile_overview = f"{args.output_prefix}_overview.png"
    fig_overview, axes_overview = plt.subplots(
        nrows=len(TOOL_NAMES),
        figsize=(5, 10),
        layout="constrained",
    )

    def func(pct, total):
        num = round((pct / 100) * total)
        return f"{pct:.1f}%\n({num})"

    for i, (tool, color) in enumerate(zip(TOOL_NAMES, TOOL_COLORS)):
        ax_overview: Axes = axes_overview[i]
        df_calls = dfs_calls[i]
        type_counts = dict(df_calls["type"].value_counts().iter_rows())
        ax_overview.pie(
            type_counts.values(),
            colors=[CALL_COLOR_KEY[typ] for typ in type_counts.keys()],
            radius=1.0,
            autopct=lambda pct: func(pct, dfs_calls[0].shape[0]),
            textprops=dict(color="w", fontsize=12),
        )
        ax_overview.set_ylabel(tool, color=color, rotation=0, ha="right", fontsize=20)

    fig_overview.legend(
        handles=[
            Patch(facecolor=color, label=CALL_NAMES[lbl])
            for lbl, color in color_key.items()
        ],
        bbox_to_anchor=(0.5, -0.05),
        loc="center",
        ncol=len(color_key),
        frameon=False,
        edgecolor="black",
        fontsize=20,
        handlelength=0.7,
        handleheight=0.7,
    )
    fig_overview.savefig(outfile_overview, bbox_inches="tight", dpi=600)

    chrom_names = df_fai["chrom_name"].unique(maintain_order=True)
    base_height_ratios = [0.5, 0.25, 0.25, 0.25, 0.25, 0.25]
    width_ratios = [0.0025, 0.5, 0.075, 0.0025, 0.5]
    num_tracks = len(base_height_ratios)
    height_ratios = base_height_ratios * len(chrom_names)
    height_ratios.append(0.5)
    fig, axes = plt.subplots(
        ncols=len(width_ratios),
        # 1 additional track for legend
        nrows=len(chrom_names) * num_tracks + 1,
        figsize=(20, len(chrom_names)),
        height_ratios=height_ratios,
        width_ratios=width_ratios,
    )
    fig.subplots_adjust(wspace=0.02)

    ax_col_indices_stats = [0, 3]
    ax_col_idx_spacer = 2
    ax_col_indices_chrom = [1, 4]

    for chrom_name_idx, chrom_name in enumerate(chrom_names):
        ax_row_idx_chrom = chrom_name_idx * num_tracks
        ax_row_indices_tracks = range(
            ax_row_idx_chrom + 1, ax_row_idx_chrom + num_tracks
        )
        # print(ax_row_idx_chrom)
        # print(tuple(ax_row_indices_tracks))
        for ax_col_idx_chrom, ax_col_idx_stats, hap in zip(
            ax_col_indices_chrom, ax_col_indices_stats, ("MATERNAL", "PATERNAL")
        ):
            if hap == "MATERNAL" and chrom_name == "chrXY":
                new_chrom_name = "chrX"
            elif hap == "PATERNAL" and chrom_name == "chrXY":
                new_chrom_name = "chrY"
            else:
                new_chrom_name = chrom_name
            chrom_name_hap = f"{new_chrom_name}_{hap}"
            print(chrom_name_hap)

            ax_chrom: Axes = axes[ax_row_idx_chrom, ax_col_idx_chrom]
            ax_chrom_stats: Axes = axes[ax_row_idx_chrom, ax_col_idx_stats]
            ax_chrom_spacer: Axes = axes[ax_row_idx_chrom, ax_col_idx_spacer]

            ax_chrom.xaxis.set_tick_params(which="both", length=0, labelleft=False)
            ax_chrom.yaxis.set_tick_params(which="both", length=0)
            # pyideogram removes xyticks
            minimize_ax(ax_chrom)
            minimize_ax(ax_chrom_stats, remove_ticks=True)
            minimize_ax(ax_chrom_spacer, remove_ticks=True)

            for idx_data, ax_row_idx in enumerate(ax_row_indices_tracks):
                ax_track: Axes = axes[ax_row_idx, ax_col_idx_chrom]
                ax_stats: Axes = axes[ax_row_idx, ax_col_idx_stats]
                ax_spacer: Axes = axes[ax_row_idx, ax_col_idx_spacer]

                ax_track.set_xlim(0, max_length)
                ax_track.set_ylim(0, 1)
                minimize_ax(ax_track, remove_ticks=True)
                minimize_ax(ax_stats, remove_ticks=True)
                minimize_ax(ax_spacer, remove_ticks=True)

                # Last track is spacer
                if ax_row_idx == ax_row_indices_tracks[-1]:
                    continue

                df_calls = dfs_calls[idx_data].filter(pl.col("chrom") == chrom_name_hap)

                color = TOOL_COLORS[idx_data]
                label = TOOL_NAMES[idx_data]
                # Write stats
                type_counts = dict(df_calls["type"].value_counts().iter_rows())
                ax_stats.pie(
                    type_counts.values(),
                    colors=[CALL_COLOR_KEY[typ] for typ in type_counts.keys()],
                    radius=3,
                )

                # Write regions
                ax_stats.set_ylabel(
                    label, **LBL_KWARGS, color=color, fontsize="x-small"
                )
                for row in df_calls.iter_rows(named=True):
                    color = color_key[row["type"]]
                    ax_track.axvspan(xmin=row["st"], xmax=row["end"], color=color)

            ax_chrom.set_xlim(0, max_length)
            pyid.ideogramh(
                chrom=chrom_name_hap,
                bands=cytobands,
                ax=ax_chrom,
                color=BAND_COLORS,
                label="",
            )
            ax_chrom.set_yticks([], [])
            ax_chrom.set_ylabel(chrom_name_hap, **LBL_KWARGS)

    ax_legend_spacer: Axes = axes[len(chrom_names) * num_tracks, ax_col_idx_spacer]
    minimize_ax(ax_legend_spacer, remove_ticks=True)

    # Add legend for each tool.
    for ax_col_idx_stats, ax_col_idx_chrom in zip(
        ax_col_indices_stats, ax_col_indices_chrom
    ):
        ax_legend: Axes = axes[len(chrom_names) * num_tracks, ax_col_idx_chrom]
        minimize_ax(ax_legend, remove_ticks=True)
        minimize_ax(
            axes[len(chrom_names) * num_tracks, ax_col_idx_stats], remove_ticks=True
        )
        ax_legend.legend(
            handles=[
                Patch(facecolor=color, label=CALL_NAMES[lbl])
                for lbl, color in color_key.items()
                if lbl != "truth"
            ],
            bbox_to_anchor=(0.5, 0.5),
            loc="center",
            ncol=len(color_key),
            frameon=False,
            edgecolor="black",
            handlelength=0.7,
            handleheight=0.7,
        )
    # Reduce white space between haps
    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=600)
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
