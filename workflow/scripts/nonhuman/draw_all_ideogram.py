import sys
import argparse
import polars as pl
import pyideogram as pyid
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex
from matplotlib.patches import Patch
from collections import defaultdict
from matplotlib.lines import Line2D

BED9P_COLS = [
    ("#chrom", pl.String),
    ("chromStart", pl.UInt64),
    ("chromEnd", pl.UInt64),
    ("name", pl.String),
    ("score", pl.UInt64),
    ("strand", pl.String),
    ("thickStart", pl.UInt64),
    ("thickEnd", pl.UInt64),
    ("itemRgb", pl.String),
    ("zscore", pl.Float32),
    ("af", pl.Float32),
    ("asm", pl.String),
]
IGNORE_CALLS = {
    "correct",
    "homopolymer",
    "dinucleotide",
    "simple_repeat",
    "het_or_mismap",
}
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
    alignment="left",
    edgecolor="black",
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


BAND_COLORS = pyid.BANDCOL | {"none": (1.0, 1.0, 1.0), "acen": (0.95, 0.95, 0.95)}
LBL_KWARGS = dict(rotation=0, ha="right", va="center")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--infile", help="Input misassemblies.")
    ap.add_argument("-y", "--cytobands", help="Input cytoband file.")
    ap.add_argument("-r", "--rename_key", help="Rename chromosomes key-value map.")
    ap.add_argument("-o", "--output_prefix", help="Output prefix.")
    args = ap.parse_args()

    df_calls = (
        pl.read_csv(
            args.infile,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            schema=dict(BED9P_COLS),
            truncate_ragged_lines=True,
        )
        .with_columns(length=pl.col("chromEnd") - pl.col("chromStart"))
        .filter(~pl.col("name").is_in(IGNORE_CALLS))
    )
    color_key = {
        name: rgb2hex([int(e) / 255.0 for e in itemRgb.split(",")])
        if not itemRgb.startswith("#")
        else itemRgb
        for name, itemRgb in df_calls.select("name", "itemRgb").unique().iter_rows()
    }

    df_fai = df_calls.group_by(["#chrom"]).agg(
        length=pl.col("chromEnd").max() - pl.col("chromStart").min()
    )
    fai_map = dict(df_fai.iter_rows())
    cytobands = pyid.dataloader.load_cytobands(args.cytobands)

    cytoband_cens: defaultdict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, itvs in cytobands.items():
        # Format of pyideogram.dataloader.load_cytobands
        # chrom: ([end], [pos], [band])
        row_itvs: list[tuple[int, str, str]] = list(zip(*itvs))
        for i, (end, pos, band) in enumerate(row_itvs):
            if band != "acen":
                continue
            try:
                prev_band = row_itvs[max(i - 1, 0)]
                st = prev_band[0]
            except IndexError:
                st = 0
            cytoband_cens[chrom].append((st, end))

    rename_key = dict(
        pl.read_csv(
            args.rename_key, separator="\t", has_header=False, columns=[0, 1]
        ).iter_rows()
    )
    chrom_names = (
        df_fai.filter(pl.col("length") > 1_000_000)
        .sort(by="length", descending=True)["#chrom"]
        .unique(maintain_order=True)
    )
    if chrom_names.is_empty():
        return 1

    # chrom, calls, spacer
    base_height_ratios = [0.1, 1.0]
    num_tracks = len(base_height_ratios)
    height_ratios = base_height_ratios * len(chrom_names)

    fig, axes = plt.subplots(
        ncols=1,
        # 1 additional track for legend
        nrows=len(chrom_names) * num_tracks,
        figsize=(8, len(chrom_names) * 1),
        height_ratios=height_ratios,
        sharex=True,
    )
    fig.subplots_adjust(wspace=0.02)

    for chrom_name_idx, chrom_name in enumerate(sorted(chrom_names)):
        df_chrom_calls = df_calls.filter(pl.col("#chrom") == chrom_name)
        chrom_length = fai_map[chrom_name]

        print(
            f"On contig #{chrom_name_idx + 1} {chrom_name} ({chrom_length // 1_000_000} Mbp)...",
            file=sys.stderr,
        )
        ax_row_idx_tracks: int = chrom_name_idx * num_tracks
        ax_row_idx_chrom: int = chrom_name_idx * num_tracks + 1
        ax_chrom: Axes = axes[ax_row_idx_chrom]
        ax_track: Axes = axes[ax_row_idx_tracks]

        ax_chrom.xaxis.set_tick_params(which="both", length=0, labelleft=False)
        ax_chrom.yaxis.set_tick_params(which="both", length=0)

        # pyideogram removes xyticks
        minimalize_ax(ax_chrom)
        minimalize_ax(ax_track, remove_ticks=True)

        # Draw bar above cen
        for st, end in cytoband_cens[chrom_name]:
            ax_track.axvspan(xmin=st, xmax=end, color="red")
        # Write regions
        for row in df_chrom_calls.iter_rows(named=True):
            color = color_key[row["name"]]
            ax_chrom.axvspan(
                xmin=row["chromStart"],
                xmax=row["chromEnd"],
                ymin=0.05,
                ymax=0.95,
                color=color,
                zorder=2,
            )

        pyid.ideogramh(
            chrom=chrom_name,
            bands=cytobands,
            ax=ax_chrom,
            color=BAND_COLORS,
            label="",
        )

        ax_chrom.set_yticks([], [])
        new_chrom_name = rename_key[chrom_name]
        ax_chrom.set_ylabel(new_chrom_name, **LBL_KWARGS)

    # Add legend.
    fig_legend, axes_legend = plt.subplots(layout="constrained", nrows=2)
    ax_legend: Axes = axes_legend[0]
    ax_legend_cen: Axes = axes_legend[1]
    minimalize_ax(ax_legend, remove_ticks=True)
    minimalize_ax(ax_legend_cen, remove_ticks=True)
    ax_legend.legend(
        handles=[
            Patch(edgecolor="black", facecolor=color, label=lbl)
            for lbl, color in color_key.items()
        ],
        ncol=6,
        title="Calls",
        loc="center left",
        **LEGEND_KWARGS,
    )
    ax_legend_cen.legend(
        handles=[Line2D([0], [0], linewidth=2, color="red")],
        labels=["Centromere"],
        borderaxespad=0,
        fancybox=False,
        frameon=False,
        alignment="left",
        loc="center left",
    )

    # Reduce white space between haps
    fig.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight", dpi=600)
    fig.savefig(f"{args.output_prefix}.png", bbox_inches="tight", dpi=600)
    fig_legend.savefig(f"{args.output_prefix}_legend.pdf", bbox_inches="tight", dpi=600)
    fig_legend.savefig(f"{args.output_prefix}_legend.png", bbox_inches="tight", dpi=600)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
