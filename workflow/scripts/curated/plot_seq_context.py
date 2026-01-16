import sys
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import rgb2hex

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
            "no_overlap",
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
            "Other",
        ],
    }
)
SEGDUP_COLORS = dict(DF_SEGDUP_COLORS.select("oname", "item_rgb").iter_rows())
CORRECT_CALLS = {"good", "correct", "Hap"}


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
    chrom = sys.argv[1]
    calls = sys.argv[2]
    censat = sys.argv[3]
    segdup = sys.argv[4]
    output_prefix = sys.argv[5]

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
        .with_columns(
            oname=pl.when(pl.col("oname") == "no_overlap")
            .then(pl.lit(0.0))
            .otherwise(pl.col("oname"))
            .cast(pl.Float64)
        )
        .filter(~pl.col("name").is_in(CORRECT_CALLS))
        .with_columns(
            oname=pl.when(pl.col("oname") == 0.0)
            .then(pl.lit("Other"))
            .when(pl.col("oname") < 0.9)
            .then(pl.lit(r"Less than 90% similarity"))
            .when(pl.col("oname").is_between(0.9, 0.98))
            .then(pl.lit(r"90 - 98% similarity"))
            .when(pl.col("oname").is_between(0.98, 0.99))
            .then(pl.lit(r"98 - 99% similarity"))
            .otherwise(pl.lit(r"Greater than 99% similarity")),
        )
        .unique(subset=["chrom", "st", "end", "name", "oname", "oitem_rgb"])
        .join(DF_SEGDUP_COLORS, on="oname")
    )
    kwargs = dict(
        separator="\t",
        new_columns=COLS[0:4],
        has_header=False,
    )
    try:
        df_misassemblies = pl.read_csv(
            calls, comment_prefix="#", **kwargs, columns=range(0, 4)
        ).filter(~pl.col("name").is_in(CORRECT_CALLS))
    except (pl.exceptions.ShapeError, pl.exceptions.OutOfBoundsError):
        # flagger
        df_misassemblies = pl.read_csv(
            calls, comment_prefix="track", **kwargs, columns=range(0, 4)
        ).filter(~pl.col("name").is_in(CORRECT_CALLS))

    if chrom != "all":
        df_censat = df_censat.filter(pl.col("chrom") == chrom)
        df_segdup = df_segdup.filter(pl.col("chrom") == chrom)
        df_misassemblies = df_misassemblies.filter(pl.col("chrom") == chrom)

    df_censat_total = (
        df_censat.drop("oname")
        .join(DF_CENSAT_COLORS, on="item_rgb")
        .group_by(["oname", "name"])
        .agg(
            length=(pl.col("end") - pl.col("st")).sum(),
            item_rgb=pl.col("item_rgb").first(),
        )
        .sort(by="length")
        .pivot(index="oname", on="name", values="length")
        .fill_null(0)
        .to_pandas()
        .set_index("oname")
    )
    df_censat_total = df_censat_total.reindex(index=DF_CENSAT_COLORS["oname"]).fillna(0)

    df_segdup_total = (
        df_segdup.group_by(["oname", "name"])
        .agg(
            length=(pl.col("end") - pl.col("st")).sum(),
            item_rgb=pl.col("item_rgb").first(),
        )
        .sort(by="length")
        .pivot(index="oname", on="name", values="length")
        .fill_null(0)
        .to_pandas()
        .set_index("oname")
    )
    df_segdup_total = df_segdup_total.reindex(index=DF_SEGDUP_COLORS["oname"]).fillna(0)

    fig_censat, ax_censat = plt.subplots(figsize=(16, 8))
    sns.heatmap(
        df_censat_total,
        ax=ax_censat,
        annot=True,
        linewidth=0.5,
        fmt=",.0f",
        cmap="mako",
    )
    for lbl in ax_censat.get_yticklabels():
        lbl.set_rotation(0)
        lbl.set_color(CENSAT_COLORS[lbl.get_text()])

    ax_censat.set_xticklabels(
        ax_censat.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor"
    )

    ax_censat.set_ylabel("Sequence type")
    ax_censat.set_xlabel("Misassembly")

    fig_segdup, ax_segdup = plt.subplots(figsize=(16, 8))
    sns.heatmap(
        df_segdup_total,
        ax=ax_segdup,
        annot=True,
        linewidth=0.5,
        fmt=",.0f",
        cmap="mako",
    )
    for lbl in ax_segdup.get_yticklabels():
        lbl.set_rotation(0)
        lbl.set_color(SEGDUP_COLORS[lbl.get_text()])

    ax_segdup.set_xticklabels(
        ax_segdup.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor"
    )

    ax_segdup.set_ylabel("Segmental duplication similarity")
    ax_segdup.set_xlabel("Misassembly")

    fig, ax = plt.subplots(figsize=(4, 4))
    df_misassemblies_total = (
        df_misassemblies.group_by(["name"])
        .agg(length=(pl.col("end") - pl.col("st")).sum())
        .sort(by="length")
    )

    plot_bar = sns.barplot(
        df_misassemblies_total,
        x="name",
        y="length",
        ax=ax,
    )
    for container in plot_bar.containers:
        ax.bar_label(container, fmt=lambda v: f"{v:,.0f}")

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.set_ylabel("Length (bp)")
    ax.set_xlabel(None)
    ax.set_xticklabels(
        ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor"
    )

    df_misassemblies_total.write_csv(
        f"{output_prefix}_total_{chrom}.tsv", separator="\t"
    )
    df_segdup_total.to_csv(f"{output_prefix}_segdup_{chrom}.tsv", sep="\t", index=True)
    df_censat_total.to_csv(f"{output_prefix}_censat_{chrom}.tsv", sep="\t", index=True)

    fig.savefig(f"{output_prefix}_total_{chrom}.png", bbox_inches="tight")
    fig_segdup.savefig(f"{output_prefix}_segdup_{chrom}.png", bbox_inches="tight")
    fig_censat.savefig(f"{output_prefix}_censat_{chrom}.png", bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
