import argparse

import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

from matplotlib.text import Text
from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex
from matplotlib.patches import Patch


COLS = [
    "chrom",
    "st",
    "end",
    "name",
    "ochrom",
    "ost",
    "oend",
    "oname",
    "oitem_rgb",
]
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
            "#808080",
        ],
        "oname": [
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
            "No Overlap",
        ],
    }
)
CENSAT_COLORS = dict(DF_CENSAT_COLORS.select("oname", "item_rgb").iter_rows())
DF_SEGDUP_COLORS = pl.DataFrame(
    {
        "item_rgb": ["#800080", "#808080", "#ffff00ff", "#ffa500", "#000000"],
        "oname": [
            r"Less than 90% similarity",
            r"90 - 98% similarity",
            r"98 - 99% similarity",
            r"Greater than 99% similarity",
            "No Overlap",
        ],
    }
)
SEGDUP_COLORS = dict(DF_SEGDUP_COLORS.select("oname", "item_rgb").iter_rows())
CORRECT_CALLS = {"good", "correct", "Hap"}
TOOL_COLORS = {
    "NucFlag v1.0.0": "purple",
    "Inspector v1.3": "teal",
    "HMM-Flagger v1.1.0": "magenta",
    "DeepVariant v1.9.0": "maroon",
}
LEGEND_KWARGS = dict(
    handlelength=1.0,
    handleheight=1.0,
    borderaxespad=0,
    fancybox=False,
    frameon=False,
)


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


# https://stackoverflow.com/a/67594395
def update_bars(ax: Axes, round_to: int):
    for c in ax.containers:
        labels = [f"{round(v.get_height() / 1_000_000, round_to)}" for v in c]
        ax.bar_label(c, labels=labels, label_type="edge")

    ylim = ax.get_ylim()
    yticks = range(int(ylim[0]), int(ylim[1]), 10_000_000)
    ax.set_yticks(yticks, [f"{round(tick / 1_000_000, round_to)}" for tick in yticks])
    ax.set_xlabel(None)


def update_deepvariant_label(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(pl.col("name").str.split("-")).with_columns(
        name=pl.when(
            pl.col("name").list[0].str.len_chars()
            < pl.col("name").list[1].str.len_chars()
        )
        .then(pl.lit("insertion"))
        .when(
            pl.col("name").list[0].str.len_chars()
            > pl.col("name").list[1].str.len_chars()
        )
        .then(pl.lit("deletion"))
        .when(pl.col("name").list[0] != pl.col("name").list[1])
        .then(pl.lit("snv"))
        .otherwise(pl.lit("other"))
    )


def draw_bar_all(df_all: pl.DataFrame, outfile: str):
    g = sns.catplot(
        data=df_all,
        x="name",
        y="length",
        col="lbl",
        hue="lbl",
        sharex=False,
        kind="bar",
        legend=None,
        aspect=1.5,
        palette=TOOL_COLORS,
    )
    for ax in g.axes.ravel():
        update_bars(ax, round_to=2)

    g.set_titles(row_template="{row_name}", col_template="{col_name}")

    g.tick_params(axis="x", rotation=45)
    g.set_ylabels("Length (Mbp)")
    g.savefig(outfile, dpi=600, bbox_inches="tight")


def draw_bar_annot(df_annot: pl.DataFrame, outfile: str):
    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        layout="constrained",
        figsize=(20, 5),
        sharex=False,
    )

    for ax, typ, colors in (
        (axes[0], "Centromere Satellites (cenSat)", CENSAT_COLORS),
        (axes[1], "Segmental Duplications", SEGDUP_COLORS),
    ):
        ax: Axes
        sns.barplot(
            df_annot.filter(pl.col("typ") == typ),
            x="oname",
            y="length",
            hue="lbl",
            order=colors.keys(),
            palette=TOOL_COLORS,
            legend=None,
            ax=ax,
        )
        for lbl in ax.get_xticklabels():
            lbl: Text
            # Color and stroke around text
            lbl.set_color(colors[lbl.get_text()])
            lbl.set_path_effects([pe.Stroke(linewidth=0.2, foreground="black")])

        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        ax.set_ylabel("Length (Mbp)")
        update_bars(ax, round_to=1)

    fig.legend(
        labels=TOOL_COLORS.keys(),
        handles=[Patch(color=color) for color in TOOL_COLORS.values()],
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        title=None,
        **LEGEND_KWARGS,
    )
    fig.savefig(outfile, dpi=600, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_calls", nargs="+", help="Input calls")
    ap.add_argument("-l", "--labels", nargs="+", help="Labels")
    ap.add_argument("-c", "--chrom", required=True, help="Chromosome.")
    ap.add_argument(
        "-s",
        "--input_censat",
        nargs="+",
        help="Input calls intersected with censat annotations.",
    )
    ap.add_argument(
        "-d",
        "--input_segdup",
        nargs="+",
        help="Input calls intersected with segdup annotations.",
    )
    ap.add_argument("-o", "--output_prefix", help="Output prefix.")
    args = ap.parse_args()
    chrom = args.chrom

    dfs_all_annot: list[pl.DataFrame] = []
    dfs_all: list[pl.DataFrame] = []
    for lbl, call, censat, segdup in zip(
        args.labels, args.input_calls, args.input_censat, args.input_segdup, strict=True
    ):
        df_censat = (
            pl.read_csv(
                censat,
                separator="\t",
                new_columns=COLS,
                has_header=False,
            )
            .filter(~pl.col("name").is_in(CORRECT_CALLS))
            .unique(subset=["chrom", "st", "end", "name", "oname", "oitem_rgb"])
            .with_columns(
                oname=pl.when(pl.col("oname") == "no_overlap")
                .then(pl.lit("No Overlap"))
                .otherwise(pl.col("oname")),
                item_rgb=pl.col("oitem_rgb").map_batches(
                    rgb_to_hex, return_dtype=pl.String
                ),
            )
        )

        df_segdup = (
            pl.read_csv(
                segdup,
                separator="\t",
                new_columns=COLS,
                has_header=False,
            )
            # Oname is percent ident
            .with_columns(
                oname=pl.when(pl.col("oname") == "no_overlap")
                .then(pl.lit(0.0))
                .otherwise(pl.col("oname"))
                .cast(pl.Float64)
            )
            .filter(~pl.col("name").is_in(CORRECT_CALLS))
            .with_columns(
                oname=pl.when(pl.col("oname") == 0.0)
                .then(pl.lit("No Overlap"))
                .when(pl.col("oname") < 0.9)
                .then(pl.lit(r"Less than 90% similarity"))
                .when(pl.col("oname").is_between(0.9, 0.98))
                .then(pl.lit(r"90 - 98% similarity"))
                .when(pl.col("oname").is_between(0.98, 0.99))
                .then(pl.lit(r"98 - 99% similarity"))
                .otherwise(pl.lit(r"Greater than 99% similarity")),
            )
            # Just take any if multiple. Probably not the best solution.
            .unique(subset=["chrom", "st", "end", "name", "oname"])
            .join(DF_SEGDUP_COLORS, on="oname")
        )
        kwargs = dict(
            separator="\t",
            new_columns=COLS[0:4],
            has_header=False,
        )
        try:
            df_misassemblies = pl.read_csv(
                call, comment_prefix="#", **kwargs, columns=range(0, 4)
            ).filter(~pl.col("name").is_in(CORRECT_CALLS))
        except (pl.exceptions.ShapeError, pl.exceptions.OutOfBoundsError):
            # flagger
            df_misassemblies = pl.read_csv(
                call, comment_prefix="track", **kwargs, columns=range(0, 4)
            ).filter(~pl.col("name").is_in(CORRECT_CALLS))

        if lbl == "DeepVariant v1.9.0":
            df_censat = update_deepvariant_label(df_censat)
            df_segdup = update_deepvariant_label(df_segdup)
            df_misassemblies = update_deepvariant_label(df_misassemblies)

        if chrom != "all":
            df_censat = df_censat.filter(pl.col("chrom") == chrom)
            df_segdup = df_segdup.filter(pl.col("chrom") == chrom)
            df_misassemblies = df_misassemblies.filter(pl.col("chrom") == chrom)

        df_censat_total = (
            df_censat.drop("oname")
            .join(DF_CENSAT_COLORS, on="item_rgb")
            .group_by(["oname"])
            .agg(
                length=(pl.col("end") - pl.col("st")).sum(),
                item_rgb=pl.col("item_rgb").first(),
            )
            .sort(by=["oname", "length"])
            .with_columns(typ=pl.lit("Centromere Satellites (cenSat)"), lbl=pl.lit(lbl))
        )

        df_segdup_total = (
            df_segdup.group_by(["oname"])
            .agg(
                length=(pl.col("end") - pl.col("st")).sum(),
                item_rgb=pl.col("item_rgb").first(),
            )
            .sort(by=["oname", "length"])
            .with_columns(typ=pl.lit("Segmental Duplications"), lbl=pl.lit(lbl))
        )
        df_annot = pl.concat([df_segdup_total, df_censat_total])
        dfs_all_annot.append(df_annot)

        df_misassemblies_total = (
            df_misassemblies.group_by(["name"])
            .agg(length=(pl.col("end") - pl.col("st")).sum())
            .sort(by=["name", "length"])
            .with_columns(lbl=pl.lit(lbl))
        )
        dfs_all.append(df_misassemblies_total)

    df_all = pl.concat(dfs_all)
    draw_bar_all(df_all, f"{args.output_prefix}.png")

    df_all_annot = pl.concat(dfs_all_annot)
    draw_bar_annot(df_all_annot, f"{args.output_prefix}_ctx.png")


if __name__ == "__main__":
    raise SystemExit(main())
