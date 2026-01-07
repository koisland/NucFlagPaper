import argparse
from matplotlib.patches import Patch
import polars as pl
import pyideogram as pyid
import matplotlib.pyplot as plt

from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex

CALL_COLOR_KEY = {
    "truth": "black",
    "false_positive": "gray",
    "true_positive": "blue",
    "false_negative": "red",
}
CALL_NAMES = {
    "truth": "Truth",
    "false_positive": "Unclear",
    "true_positive": "True Positive",
    "false_negative": "False Negative",
}
LBL_KWARGS = dict(rotation=0, ha="right", va="center", fontsize="medium")
TOOL_NAMES = [
    "Truth",
    "NucFlag v1.0",
    "Inspector v1.3",
    "HMM-Flagger v1.1.0",
    "DeepVariant v1.9.0",
]
TOOL_COLORS = ["black", "purple", "teal", "magenta", "maroon"]
SEGDUP_ORDERING = (
    r"Less than 90% similarity",
    r"90 - 98% similarity",
    r"98 - 99% similarity",
    r"Greater than 99% similarity",
)
DF_CENSAT_COLORS = pl.DataFrame(
    {
        "item_rgb": [
            "#00994c",
            "#e0e0e0",
            "#a9a9a9",
            "#ac33c7",
            "#fa99ff",
            "#333366",
            "#000000",
            "#fa0000",
            "#ffcc99",
            "#00cccc",
            "#78a1bb",
            "#990000",
            "#91ff00",
            "#ff9200",
        ],
        "name": [
            "hsat1B",
            "ct",
            "GAP",
            "gSat",
            "bSat",
            "HSat2",
            "rDNA",
            "active_hor",
            "mon",
            "cenSat(Other)",
            "HSat3",
            "dhor",
            "hsat1A",
            "hor",
        ],
    }
)
SEGDUP_COLORS = dict(
    zip(
        [
            r"Less than 90% similarity",
            r"90 - 98% similarity",
            r"98 - 99% similarity",
            r"Greater than 99% similarity",
        ],
        ["#800080", "#808080", "#ffff00ff", "#ffa500"],
    )
)


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
    ap.add_argument("--cytobands")
    ap.add_argument("--truth")
    ap.add_argument("--nucflag")
    ap.add_argument("--flagger")
    ap.add_argument("--inspector")
    ap.add_argument("--deepvariant")
    ap.add_argument("--segdups")
    ap.add_argument("--censat")
    ap.add_argument("-o", "--outfile", default="out.png", help="Output file.")
    ap.add_argument(
        "-f", "--fp", action="store_true", help='Include "false-positives".'
    )
    ap.add_argument("-c", "--chrom", help="Chromosome name.", default="chrY_PATERNAL")
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
        )
        .with_columns(type=pl.lit("truth"))
        .filter(pl.col("chrom") == args.chrom),
        pl.read_csv(args.nucflag, **missed_calls_kwargs).filter(
            pl.col("chrom") == args.chrom
        ),
        pl.read_csv(args.inspector, **missed_calls_kwargs).filter(
            pl.col("chrom") == args.chrom
        ),
        pl.read_csv(args.flagger, **missed_calls_kwargs).filter(
            pl.col("chrom") == args.chrom
        ),
        pl.read_csv(args.deepvariant, **missed_calls_kwargs).filter(
            pl.col("chrom") == args.chrom
        ),
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

    df_segdup = (
        pl.read_csv(
            args.segdups,
            separator="\t",
            has_header=False,
            columns=[*range(0, 9), 35],
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
                "filter_score",
            ],
        )
        .with_columns(
            pl.col("item_rgb").map_batches(rgb_to_hex, return_dtype=pl.String),
        )
        .filter(pl.col("chrom") == args.chrom)
        .with_columns(
            name=pl.when(pl.col("filter_score") < 0.9)
            .then(pl.lit(r"Less than 90% similarity"))
            .when(pl.col("filter_score").is_between(0.9, 0.98))
            .then(pl.lit(r"90 - 98% similarity"))
            .when(pl.col("filter_score").is_between(0.98, 0.99))
            .then(pl.lit(r"98 - 99% similarity"))
            .otherwise(pl.lit(r"Greater than 99% similarity"))
            .cast(pl.Enum(SEGDUP_ORDERING))
        )
    )
    df_censat = (
        pl.read_csv(
            args.censat,
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
        .drop("name")
        .filter(pl.col("chrom") == args.chrom)
        .join(DF_CENSAT_COLORS, on="item_rgb")
    )

    color_key = CALL_COLOR_KEY
    height_ratios = [0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5]  #
    fig, axes = plt.subplots(
        ncols=1,
        nrows=len(height_ratios),
        figsize=(20, 4),
        height_ratios=height_ratios,
        sharex=True,
    )

    indices_calls = range(3, 8)  #
    idx_chrom = 0
    idx_censat = 1
    idx_segdup = 2
    idx_censat_legend = 8  #
    idx_segdup_legend = 9  #
    idx_legend = 10  #
    ax: Axes = axes[idx_chrom]

    ax.xaxis.set_tick_params(which="both", length=0, labelleft=False)
    ax.yaxis.set_tick_params(which="both", length=0)
    # pyideogram removes xyticks
    minimize_ax(ax)

    pyid.ideogramh(args.chrom, cytobands, ax)

    # Add sequence context.
    ax_censat: Axes = axes[idx_censat]
    ax_segdup: Axes = axes[idx_segdup]
    minimize_ax(ax_censat, remove_ticks=True)
    minimize_ax(ax_segdup, remove_ticks=True)
    ax_censat.set_ylabel("Centromere/Satellites", **LBL_KWARGS)
    ax_segdup.set_ylabel("Segmental duplications", **LBL_KWARGS)

    for row in df_censat.iter_rows(named=True):
        ax_censat.axvspan(
            xmin=row["st"], xmax=row["end"], color=row["item_rgb"], label=row["name"]
        )
    for row in df_segdup.iter_rows(named=True):
        ax_segdup.axvspan(
            xmin=row["st"], xmax=row["end"], color=row["item_rgb"], label=row["name"]
        )

    # Add calls
    for i, idx in enumerate(indices_calls):
        ax_other: Axes = axes[idx]
        df_calls = dfs_calls[i].sort(by="st")
        ax_other.set_ylim(0, 1)
        minimize_ax(ax_other, remove_ticks=True)
        tool_name = TOOL_NAMES[i]
        ax_other.set_ylabel(tool_name, **LBL_KWARGS, color=TOOL_COLORS[i])
        print(tool_name)
        for row in df_calls.iter_rows(named=True):
            color = color_key[row["type"]]
            # Plot on top.
            if row["type"] == "true_positive":
                zorder = 1
            else:
                zorder = 0
            ax_other.axvspan(
                xmin=row["st"], xmax=row["end"], color=color, zorder=zorder
            )

    for ax_ref, i in (
        (ax_censat, idx_censat_legend),
        (ax_segdup, idx_segdup_legend),
        (None, idx_legend),
    ):
        if ax_ref:
            handles, labels = ax_ref.get_legend_handles_labels()
            labels_handles = dict(zip(labels, handles))
        else:
            labels_handles = {
                CALL_NAMES[call]: Patch(color=color)
                for call, color in CALL_COLOR_KEY.items()
                if call != "truth"
            }

        ax_legend: Axes = axes[i]
        minimize_ax(ax_legend, remove_ticks=True)

        if i == idx_segdup_legend:
            labels_handles = {
                categ: Patch(color=SEGDUP_COLORS[categ]) for categ in SEGDUP_ORDERING
            }

        ax_legend.legend(
            handles=labels_handles.values(),
            labels=labels_handles.keys(),
            bbox_to_anchor=(0.5, 0.5),
            loc="center",
            ncol=len(labels_handles),
            frameon=False,
            edgecolor="black",
            handlelength=0.7,
            handleheight=0.7,
        )

    fig.savefig(args.outfile, bbox_inches="tight", dpi=600)


if __name__ == "__main__":
    raise SystemExit(main())
